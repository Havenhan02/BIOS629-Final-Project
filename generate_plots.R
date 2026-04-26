library(tidyverse)
library(survival)
library(survminer)
library(randomForestSRC)

# Load results
cv_results <- readRDS("final_output/cv_results.rds")
cv_summary <- readRDS("final_output/cv_summary.rds")
importance_df <- readRDS("final_output/importance_df.rds")

cat("Generating plots...\n")

# 1. CV Box Plot
cv_plot <- ggplot(cv_results, aes(x = Model, y = C_index, fill = Model)) +
  geom_boxplot(alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.7) +
  scale_fill_manual(values = c("#4DBBD5", "#00A087", "#E64B35")) +
  coord_cartesian(ylim = c(0.60, 0.72)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "5-Fold Cross-Validation: Model Performance Comparison",
    subtitle = "C-index distribution across folds",
    y = "C-index",
    x = NULL
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )
ggsave("final_output/fig1_cv_boxplot.png", cv_plot, width = 10, height = 7, dpi = 300)
cat("Saved fig1_cv_boxplot.png\n")

# 2. CV Bar Plot with error bars
cv_bar_plot <- ggplot(cv_summary, aes(x = reorder(Model, Mean_C_index), y = Mean_C_index, fill = Model)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_C_index - SD_C_index, ymax = Mean_C_index + SD_C_index), 
                width = 0.2, linewidth = 1) +
  scale_fill_manual(values = c("#4DBBD5", "#00A087", "#E64B35")) +
  coord_cartesian(ylim = c(0.60, 0.72)) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean C-index from 5-Fold Cross-Validation",
    subtitle = "Error bars represent standard deviation across folds",
    y = "C-index",
    x = NULL
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  geom_text(aes(label = sprintf("%.3f", Mean_C_index)), vjust = -1, size = 5, fontface = "bold")
ggsave("final_output/fig2_cv_barplot.png", cv_bar_plot, width = 10, height = 6, dpi = 300)
cat("Saved fig2_cv_barplot.png\n")

# 3. Variable Importance Plot
importance_plot <- ggplot(importance_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Variable Importance from Random Survival Forest",
    subtitle = "Higher values indicate greater importance for prediction",
    x = "Variable",
    y = "Importance Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
ggsave("final_output/fig3_variable_importance.png", importance_plot, width = 10, height = 6, dpi = 300)
cat("Saved fig3_variable_importance.png\n")

cat("All plots saved to final_output/\n")
