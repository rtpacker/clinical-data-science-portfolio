# =============================================================================
# Residential Segregation & Thrombo-Inflammatory Biomarkers
# Forest Plot Figures
# =============================================================================
# Description: Reads regression result Excel files and produces three forest
#              plots faceted by segregation index (Dissimilarity, Isolation,
#              Interaction):
#
#   Plot 1 — Reverse-T % change (log-transformed biomarkers):
#             CRP, D-Dimer, IFN, TNF, IL-6
#   Plot 2 — Beta coefficient (untransformed biomarkers):
#             E-Selectin, FIX
#   Plot 3 — Odds Ratio (IL-1b tertile comparisons):
#             IL-1b T1 vs T2, IL-1b T1 vs T3
#
# 
# =============================================================================


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(openxlsx)


# -----------------------------------------------------------------------------
# Paths — update before running
# -----------------------------------------------------------------------------
input.path  <- "path/to/regression_results/"
output.path <- "path/to/output/figures/"


# =============================================================================
# Helper Functions
# =============================================================================

# -----------------------------------------------------------------------------
# load_biomarker()
# Reads all three index sheets from one regression result workbook,
# filters to the target model row, and tags Group and Biomarker columns.
# -----------------------------------------------------------------------------
load_biomarker <- function(filename, biomarker_label, model_row = "Age and gender adjusted") {
  sheets <- c(Dissimilarity = "Dissimilarity_Results",
              Isolation     = "Isolation_Results",
              Interaction   = "Interaction_Results")

  lapply(names(sheets), function(group) {
    read.xlsx(paste0(input.path, filename), sheet = sheets[[group]]) %>%
      filter(Model == model_row) %>%
      mutate(Biomarker = biomarker_label, Group = group)
  }) %>% bind_rows()
}

# -----------------------------------------------------------------------------
# parse_ci_reverse_t()
# Parses the CI string for log-transformed biomarkers and converts
# Reverse_T from a "%" string to a scaled numeric (divide by 100).
# -----------------------------------------------------------------------------
parse_ci_reverse_t <- function(df) {
  df %>%
    mutate(Reverse_T_t = as.numeric(sub("%", "", Reverse_T)) / 100,
           CI          = gsub("[()]", "", CI)) %>%
    separate(CI, into = c("CI_lower", "CI_upper"), sep = ",") %>%
    mutate(CI_lower = as.numeric(gsub("%", "", CI_lower)) / 100,
           CI_upper = as.numeric(gsub("%", "", CI_upper)) / 100)
}

# -----------------------------------------------------------------------------
# parse_ci_beta()
# Parses the 95CI string for untransformed biomarkers (Beta scale).
# -----------------------------------------------------------------------------
parse_ci_beta <- function(df) {
  df %>%
    mutate(Beta = as.numeric(Beta),
           CI   = gsub("[()]", "", `95CI`)) %>%
    separate(CI, into = c("CI_lower", "CI_upper"), sep = ",") %>%
    mutate(CI_lower = as.numeric(gsub("%", "", CI_lower)),
           CI_upper = as.numeric(gsub("%", "", CI_upper)))
}

# -----------------------------------------------------------------------------
# parse_ci_or()
# Parses the OR_CI_95 string for logistic (IL-1b tertile) models.
# -----------------------------------------------------------------------------
parse_ci_or <- function(df) {
  df %>%
    mutate(Odds_Ratio = as.numeric(Odds_Ratio),
           CI         = gsub("[()]", "", OR_CI_95)) %>%
    separate(CI, into = c("CI_lower", "CI_upper"), sep = ",") %>%
    mutate(CI_lower = as.numeric(gsub("%", "", CI_lower)),
           CI_upper = as.numeric(gsub("%", "", CI_upper)))
}

# Shared color scale and Group factor order for all plots
index_colors <- c("Dissimilarity" = "red", "Isolation" = "blue", "Interaction" = "purple")
index_levels <- c("Dissimilarity", "Isolation", "Interaction")


# =============================================================================
# PLOT 1: Reverse-T Forest Plot — Log-Transformed Biomarkers
# =============================================================================
combined_results <- bind_rows(
  load_biomarker("crpresultsregressindexweighted.xlsx",    "CRP"),
  load_biomarker("ddimerresultsregressindexweighted.xlsx", "D-Dimer"),
  load_biomarker("ifnresultsregressindexweighted.xlsx",    "IFN"),
  load_biomarker("tnfresultsregressindexweighted.xlsx",    "TNF"),
  load_biomarker("il6resultsregressindexweighted.xlsx",    "IL-6")
) %>%
  parse_ci_reverse_t() %>%
  mutate(Group = factor(Group, levels = index_levels))

# Main faceted forest plot
plot1 <- ggplot(combined_results, aes(x = Reverse_T_t, y = Biomarker, color = Group)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper),
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_text(aes(label = round(Reverse_T_t, 3)),
            position = position_dodge(width = 0.5),
            vjust = -1, size = 3, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ Group, ncol = 1) +
  scale_color_manual(values = index_colors) +
  labs(title = "Forest Plot: Percent Difference by Biomarker and Segregation Index",
       x = "Percent Difference", y = "Biomarker", color = "Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Companion text label panels (one per index group)
make_label_panel <- function(data, group_name) {
  data %>%
    filter(Group == group_name) %>%
    mutate(label_text = paste0(Biomarker, ": ", round(Reverse_T_t, 3),
                               " (", round(CI_lower, 3), ", ", round(CI_upper, 3), ")")) %>%
    ggplot(aes(x = 1, y = Biomarker, label = label_text)) +
    geom_text(size = 4, hjust = 0.5, vjust = 0.5, color = "black") +
    labs(title = paste0(group_name, ": Reverse-T (lower CI, upper CI)")) +
    theme_void() +
    theme(
      plot.title   = element_text(size = 12, face = "bold", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin  = margin(10, 10, 10, 10)
    )
}

final_plot_reverseT <- plot1 | (make_label_panel(combined_results, "Dissimilarity") /
                                make_label_panel(combined_results, "Isolation")     /
                                make_label_panel(combined_results, "Interaction"))
final_plot_reverseT

ggsave(paste0(output.path, "forest_plot_reverseT.png"), final_plot_reverseT,
       width = 14, height = 10, dpi = 300)


# =============================================================================
# PLOT 2: Beta Forest Plot — Untransformed Biomarkers (E-Selectin, FIX)
# =============================================================================
combined_beta_results <- bind_rows(
  load_biomarker("eselectinresultsregressindexweighted.xlsx", "E-Selectin"),
  load_biomarker("fixresultsregressindexweighted.xlsx",       "FIX")
) %>%
  parse_ci_beta() %>%
  mutate(Group = factor(Group, levels = index_levels))

plot_beta <- ggplot(combined_beta_results, aes(x = Beta, y = Biomarker, color = Group)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper),
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_text(aes(label = round(Beta, 3)),
            position = position_dodge(width = 0.5),
            vjust = -1, size = 3, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ Group, ncol = 1) +
  scale_color_manual(values = index_colors) +
  labs(title = "Forest Plot: Beta Coefficient by Biomarker and Segregation Index",
       x = "Beta", y = "Biomarker", color = "Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_beta

ggsave(paste0(output.path, "forest_plot_beta.png"), plot_beta,
       width = 10, height = 8, dpi = 300)


# =============================================================================
# PLOT 3: Odds Ratio Forest Plot — IL-1b Tertile Comparisons
# =============================================================================
combined_tertile_results <- bind_rows(
  load_biomarker("il1bt1vt2resultsregressindexweighted.xlsx", "IL-1b T1 vs T2"),
  load_biomarker("il1bt1vt3resultsregressindexweighted.xlsx", "IL-1b T1 vs T3")
) %>%
  parse_ci_or() %>%
  mutate(Group = factor(Group, levels = index_levels))

plot_or <- ggplot(combined_tertile_results, aes(x = Odds_Ratio, y = Biomarker, color = Group)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper),
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_text(aes(label = round(Odds_Ratio, 3)),
            position = position_dodge(width = 0.5),
            vjust = -1, size = 3, show.legend = FALSE) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  facet_wrap(~ Group, ncol = 1) +
  scale_color_manual(values = index_colors) +
  labs(title = "Forest Plot: Odds Ratio by IL-1b Tertile and Segregation Index",
       x = "Odds Ratio", y = "Biomarker", color = "Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_or

ggsave(paste0(output.path, "forest_plot_odds_ratio.png"), plot_or,
       width = 10, height = 8, dpi = 300)
