library(tidyverse)
library(survival)
library(survminer)
library(randomForestSRC)
library(caret)

set.seed(629)

data <- read.csv("/Users/morgan/Desktop/Biostat629_Interim_project/dataset/msk_chord_2024_NSCLC_data.csv", skip = 1)
analysis_data <- data %>%
  select(
    age = `Current.Age`, sex = Sex, stage = `Stage..Highest.Recorded.`,
    fga = `Fraction.Genome.Altered`, smoking = `Smoking.History..NLP.`,
    mutation_count = `Mutation.Count`, os_time = `Overall.Survival..Months.`,
    os_status = `Overall.Survival.Status`, tmb = `TMB..nonsynonymous.`
  ) %>%
  mutate(
    event = ifelse(grepl("DECEASED|1:DECEASED", os_status), 1, 0),
    stage = as.factor(stage), sex = as.factor(sex),
    age = as.numeric(age), fga = as.numeric(fga),
    mutation_count = as.numeric(mutation_count), os_time = as.numeric(os_time),
    tmb = as.numeric(tmb)
  ) %>% filter(!is.na(os_time) & os_time > 0)

analysis_data$smoking <- factor(analysis_data$smoking, levels = c("Never", "Former/Current Smoker", "Unknown"))
analysis_complete <- analysis_data %>% select(-os_status) %>%
  filter(!is.na(age) & !is.na(sex) & !is.na(stage) & !is.na(fga) & !is.na(smoking) & !is.na(mutation_count) & !is.na(tmb))

cat("Sample size:", nrow(analysis_complete), "\n")

folds <- createFolds(analysis_complete$event, k = 5, list = TRUE)
all_preds <- data.frame(row_index = integer(), rsf_risk = numeric())

for(i in 1:5) {
  cat("Fold", i, "...\n")
  test_idx <- folds[[i]]
  train_data <- analysis_complete[-test_idx, ]
  test_data <- analysis_complete[test_idx, ]
  rsf_model <- rfsrc(Surv(os_time, event) ~ age + sex + stage + smoking + fga + mutation_count + tmb,
                     data = train_data, ntree = 500, forest = TRUE)
  rsf_pred <- predict(rsf_model, newdata = test_data)
  all_preds <- rbind(all_preds, data.frame(row_index = test_idx, rsf_risk = rsf_pred$predicted))
}

analysis_complete$row_index <- seq_len(nrow(analysis_complete))
analysis_risk <- analysis_complete %>% left_join(all_preds, by = "row_index")

missing_idx <- which(is.na(analysis_risk$rsf_risk))
if(length(missing_idx) > 0) {
  final_rsf <- rfsrc(Surv(os_time, event) ~ age + sex + stage + smoking + fga + mutation_count + tmb,
                     data = analysis_complete, ntree = 500, importance = FALSE)
  analysis_risk$rsf_risk[missing_idx] <- final_rsf$predicted.oob[missing_idx]
}

analysis_risk <- analysis_risk %>%
  mutate(rsf_risk_group = cut(rsf_risk, breaks = quantile(rsf_risk, probs = c(0, 1/3, 2/3, 1)),
                              labels = c("Low Risk", "Intermediate Risk", "High Risk"), include.lowest = TRUE))

cat("\nRisk group distribution:\n")
print(table(analysis_risk$rsf_risk_group))

risk_survival <- analysis_risk %>%
  group_by(rsf_risk_group) %>%
  summarise(n = n(), median_os = median(os_time), events = sum(event),
            event_rate = round(100 * sum(event)/n(), 1), .groups='drop')
cat("\nRisk group survival:\n")
print(risk_survival)

km_risk <- survfit(Surv(os_time, event) ~ rsf_risk_group, data = analysis_risk)
p3 <- ggsurvplot(km_risk, data = analysis_risk, pval = TRUE, pval.method = TRUE,
  risk.table = TRUE, title = "Survival by RSF Risk Groups",
  subtitle = "Patients stratified by Random Survival Forest predicted risk",
  xlab = "Time (Months)", ylab = "Survival Probability",
  palette = c("#00A087", "#E69F00", "#E64B35"),
  risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,
  ggtheme = theme_minimal(base_size = 12),
  legend.title = "Risk Group", legend.labs = c("Low Risk", "Intermediate Risk", "High Risk"))

png("final_output/fig6_km_risk_groups.png", width=1200, height=800, res=150)
print(p3)
dev.off()
cat("Saved fig6_km_risk_groups.png\n")

saveRDS(analysis_risk, "final_output/analysis_risk.rds")
write.csv(risk_survival, "final_output/risk_group_summary.csv", row.names = FALSE)
cat("Done!\n")
