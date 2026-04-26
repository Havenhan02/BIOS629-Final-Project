# Load libraries
library(tidyverse)
library(survival)
library(survminer)
library(randomForestSRC)
library(survcomp)
library(caret)

# Set seed
set.seed(629)

# Read data
cat("=== LOADING DATA ===\n")
data <- read.csv("/Users/morgan/Desktop/Biostat629_Interim_project/dataset/msk_chord_2024_NSCLC_data.csv", skip = 1)
cat("Dataset dimensions:", dim(data)[1], "rows x", dim(data)[2], "columns\n")

# Select and clean variables
analysis_data <- data %>%
  select(
    age = `Current.Age`,
    sex = Sex,
    stage = `Stage..Highest.Recorded.`,
    fga = `Fraction.Genome.Altered`,
    smoking = `Smoking.History..NLP.`,
    mutation_count = `Mutation.Count`,
    os_time = `Overall.Survival..Months.`,
    os_status = `Overall.Survival.Status`,
    tmb = `TMB..nonsynonymous.`,
    tumor_purity = `Tumor.Purity`
  )

# Clean and transform
cat("\n=== DATA CLEANING ===\n")
analysis_data <- analysis_data %>%
  mutate(
    event = ifelse(grepl("DECEASED|1:DECEASED", os_status), 1, 0),
    stage = as.factor(stage),
    sex = as.factor(sex),
    age = as.numeric(age),
    fga = as.numeric(fga),
    mutation_count = as.numeric(mutation_count),
    os_time = as.numeric(os_time),
    tmb = as.numeric(tmb),
    tumor_purity = as.numeric(tumor_purity)
  ) %>%
  filter(!is.na(os_time) & os_time > 0)

analysis_data$smoking <- factor(
  analysis_data$smoking,
  levels = c("Never", "Former/Current Smoker", "Unknown")
)

cat("Initial sample size:", nrow(analysis_data), "\n")
cat("Events:", sum(analysis_data$event), "\n")

# Complete cases
analysis_complete <- analysis_data %>%
  select(-tumor_purity) %>%
  filter(
    !is.na(age), !is.na(sex), !is.na(stage), !is.na(fga), 
    !is.na(smoking), !is.na(mutation_count), !is.na(tmb)
  )

cat("Final complete case sample size:", nrow(analysis_complete), "\n")
cat("Events in complete data:", sum(analysis_complete$event), "\n")

# Create output directory
if(!dir.exists("final_output")) dir.create("final_output")

# Save descriptive stats
cat("\n=== DESCRIPTIVE STATISTICS ===\n")
desc_stats <- summary(analysis_complete)
print(desc_stats)

# 5-fold CV
cat("\n=== 5-FOLD CROSS-VALIDATION ===\n")
folds <- createFolds(analysis_complete$event, k = 5, list = TRUE)

cv_results <- data.frame(
  Fold = integer(),
  Model = character(),
  C_index = numeric(),
  stringsAsFactors = FALSE
)

all_predictions <- data.frame(
  row_index = integer(),
  fold = integer(),
  cox_clinical_lp = numeric(),
  cox_full_lp = numeric(),
  rsf_risk = numeric(),
  os_time = numeric(),
  event = numeric()
)

for(i in 1:5) {
  cat("Processing Fold", i, "...\n")
  
  test_idx <- folds[[i]]
  train_data <- analysis_complete[-test_idx, ]
  test_data <- analysis_complete[test_idx, ]
  
  # Clinical Cox
  cox_clinical <- coxph(Surv(os_time, event) ~ age + sex + stage + smoking, data = train_data)
  lp_clinical <- predict(cox_clinical, newdata = test_data, type = "lp")
  cox_clinical_cindex <- survConcordance(Surv(test_data$os_time, test_data$event) ~ lp_clinical)$concordance
  
  # Full Cox
  cox_full <- coxph(Surv(os_time, event) ~ age + sex + stage + smoking + fga + mutation_count + tmb, data = train_data)
  lp_full <- predict(cox_full, newdata = test_data, type = "lp")
  cox_full_cindex <- survConcordance(Surv(test_data$os_time, test_data$event) ~ lp_full)$concordance
  
  # RSF
  rsf_model <- rfsrc(Surv(os_time, event) ~ age + sex + stage + smoking + fga + mutation_count + tmb, 
                     data = train_data, ntree = 500, forest = TRUE)
  rsf_pred <- predict(rsf_model, newdata = test_data)
  rsf_cindex <- survConcordance(Surv(test_data$os_time, test_data$event) ~ rsf_pred$predicted)$concordance
  
  cv_results <- rbind(cv_results, data.frame(Fold=i, Model="Clinical Cox", C_index=cox_clinical_cindex))
  cv_results <- rbind(cv_results, data.frame(Fold=i, Model="Clinical + Genomic Cox", C_index=cox_full_cindex))
  cv_results <- rbind(cv_results, data.frame(Fold=i, Model="Random Survival Forest", C_index=rsf_cindex))
  
  all_predictions <- rbind(all_predictions, data.frame(
    row_index = test_idx, fold = i,
    cox_clinical_lp = lp_clinical, cox_full_lp = lp_full, rsf_risk = rsf_pred$predicted,
    os_time = test_data$os_time, event = test_data$event
  ))
}

cat("\n=== CV RESULTS ===\n")
cv_summary <- cv_results %>%
  group_by(Model) %>%
  summarise(Mean_C_index = mean(C_index), SD_C_index = sd(C_index), .groups='drop') %>%
  arrange(desc(Mean_C_index))
print(cv_summary)

print(cv_results %>% pivot_wider(names_from=Model, values_from=C_index))

# Train final models
cat("\n=== FINAL MODELS ===\n")
final_cox_clinical <- coxph(Surv(os_time, event) ~ age + sex + stage + smoking, data = analysis_complete)
final_cox_full <- coxph(Surv(os_time, event) ~ age + sex + stage + smoking + fga + mutation_count + tmb, data = analysis_complete)
final_rsf <- rfsrc(Surv(os_time, event) ~ age + sex + stage + smoking + fga + mutation_count + tmb, 
                   data = analysis_complete, ntree = 1000, importance = TRUE)

cat("\nClinical Cox Model:\n")
print(summary(final_cox_clinical)$concordance)
cat("\nFull Cox Model:\n")
print(summary(final_cox_full)$concordance)
cat("\nRSF OOB Error:", final_rsf$err.rate[length(final_rsf$err.rate)], "\n")
cat("RSF C-index:", 1 - final_rsf$err.rate[length(final_rsf$err.rate)], "\n")

# Variable importance
cat("\n=== VARIABLE IMPORTANCE ===\n")
vimp_result <- vimp(final_rsf)
importance_df <- data.frame(Variable = names(vimp_result$importance), Importance = as.numeric(vimp_result$importance)) %>%
  arrange(desc(Importance))
print(importance_df)

# Risk stratification
cat("\n=== RISK STRATIFICATION ===\n")
analysis_complete$row_index <- 1:nrow(analysis_complete)
analysis_risk <- analysis_complete %>%
  left_join(all_predictions %>% select(row_index, cox_clinical_lp, cox_full_lp, rsf_risk), by = "row_index")

missing_idx <- which(is.na(analysis_risk$rsf_risk))
if(length(missing_idx) > 0) {
  analysis_risk$cox_clinical_lp[missing_idx] <- predict(final_cox_clinical, type = "lp")[missing_idx]
  analysis_risk$cox_full_lp[missing_idx] <- predict(final_cox_full, type = "lp")[missing_idx]
  analysis_risk$rsf_risk[missing_idx] <- final_rsf$predicted.oob[missing_idx]
}

analysis_risk <- analysis_risk %>%
  mutate(
    rsf_risk_group = cut(rsf_risk, breaks = quantile(rsf_risk, probs = c(0, 1/3, 2/3, 1)),
                         labels = c("Low Risk", "Intermediate Risk", "High Risk"), include.lowest = TRUE)
  )

cat("Risk group distribution:\n")
print(table(analysis_risk$rsf_risk_group))

# Survival by risk group
risk_survival <- analysis_risk %>%
  group_by(rsf_risk_group) %>%
  summarise(n = n(), median_os = median(os_time), events = sum(event), 
            event_rate = round(100 * sum(event) / n(), 1), .groups = 'drop')
cat("\nSurvival by risk group:\n")
print(risk_survival)

# Save results
saveRDS(cv_results, "final_output/cv_results.rds")
saveRDS(cv_summary, "final_output/cv_summary.rds")
saveRDS(importance_df, "final_output/importance_df.rds")
saveRDS(analysis_risk, "final_output/analysis_risk.rds")
write.csv(cv_results, "final_output/cv_results.csv", row.names = FALSE)
write.csv(cv_summary, "final_output/cv_summary.csv", row.names = FALSE)
write.csv(importance_df, "final_output/variable_importance.csv", row.names = FALSE)
write.csv(risk_survival, "final_output/risk_group_summary.csv", row.names = FALSE)

cat("\n=== RESULTS SAVED TO final_output/ ===\n")
