# Supplementary analysis:
# Diagnostic performance of family history as a marker of PGV carrier status

library(haven)
library(dplyr)
library(tidyr)
library(purrr)
library(writexl)

input_file <- "ACC_PV_FH_controls_20260803.dta"
output_file <- "Supplementary_Table_S6_FH_diagnostic_performance.xlsx"

data <- read_dta(input_file) %>%
  mutate(
    case_control = case_when(
      as.character(CaCo_status) %in% c("0", "control", "Control") ~ "control",
      as.character(CaCo_status) %in% c("1", "case", "Case") ~ "case",
      TRUE ~ NA_character_
    ),
    pgv_any = case_when(
      gene_clinlof == "No_carrier" ~ 0L,
      gene_clinlof == "QC out" ~ NA_integer_,
      is.na(gene_clinlof) ~ NA_integer_,
      TRUE ~ 1L
    )
  ) %>%
  filter(case_control == "control", !is.na(pgv_any))

if (nrow(data) != 11811) {
  warning(
    "The control count is ", nrow(data),
    "; the expected manuscript count is 11,811."
  )
}

brca_genes <- c("BRCA1", "BRCA2")

hrr_genes <- c(
  "ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1",
  "CHEK2", "PALB2", "RAD51C", "RAD51D"
)

hrr_without_brca_genes <- setdiff(hrr_genes, brca_genes)
mmr_genes <- c("MLH1", "MSH2", "MSH6", "PMS2")

data <- data %>%
  mutate(
    pgv_brca1_2 = as.integer(pgv_any == 1L & gene_clinlof %in% brca_genes),
    pgv_brca1 = as.integer(pgv_any == 1L & gene_clinlof == "BRCA1"),
    pgv_brca2 = as.integer(pgv_any == 1L & gene_clinlof == "BRCA2"),
    pgv_hrr = as.integer(pgv_any == 1L & gene_clinlof %in% hrr_genes),
    pgv_hrr_without_brca = as.integer(
      pgv_any == 1L & gene_clinlof %in% hrr_without_brca_genes
    ),
    pgv_mmr = as.integer(pgv_any == 1L & gene_clinlof %in% mmr_genes)
  )

pgv_groups <- tibble(
  pgv_var = c(
    "pgv_any", "pgv_brca1_2", "pgv_brca1", "pgv_brca2",
    "pgv_hrr", "pgv_hrr_without_brca", "pgv_mmr"
  ),
  PGV_group = c(
    "All PGVs", "BRCA1/2", "BRCA1", "BRCA2",
    "HRR", "HRR excluding BRCA1/2", "MMR"
  )
)

fh_categories <- tibble(
  fh_var = c(
    "fh_any_fdr", "fh_breast_fdr", "fh_ovarian_fdr",
    "fh_pancreas_fdr", "fh_prostate_fdr", "fh_colon_fdr",
    "fh_gastric_fdr", "fh_lung_fdr", "fh_other_fdr"
  ),
  Family_history = c(
    "Any cancer", "Breast cancer", "Ovarian cancer",
    "Pancreatic cancer", "Prostate cancer", "Colorectal cancer",
    "Gastric cancer", "Lung cancer", "Other cancers"
  )
)

required_variables <- c(pgv_groups$pgv_var, fh_categories$fh_var)
missing_variables <- setdiff(required_variables, names(data))

if (length(missing_variables) > 0) {
  stop(
    "The following required variables are missing: ",
    paste(missing_variables, collapse = ", ")
  )
}

safe_percent <- function(numerator, denominator) {
  ifelse(denominator > 0, 100 * numerator / denominator, NA_real_)
}

calculate_performance <- function(df, pgv_var, pgv_label, fh_var, fh_label) {
  analysis_df <- df %>%
    filter(!is.na(.data[[pgv_var]]), !is.na(.data[[fh_var]]))
  
  tab <- table(
    FH = factor(analysis_df[[fh_var]], levels = c(0, 1)),
    PGV = factor(analysis_df[[pgv_var]], levels = c(0, 1))
  )
  
  tn <- as.numeric(tab["0", "0"])
  fn <- as.numeric(tab["0", "1"])
  fp <- as.numeric(tab["1", "0"])
  tp <- as.numeric(tab["1", "1"])
  
  tibble(
    PGV_group = pgv_label,
    Family_history = fh_label,
    N_analyzed = nrow(analysis_df),
    TP = tp,
    FP = fp,
    FN = fn,
    TN = tn,
    Sensitivity_percent = safe_percent(tp, tp + fn),
    Specificity_percent = safe_percent(tn, tn + fp),
    PPV_percent = safe_percent(tp, tp + fp),
    NPV_percent = safe_percent(tn, tn + fn),
    False_negative_percent = safe_percent(fn, tp + fn)
  )
}

analysis_grid <- crossing(pgv_groups, fh_categories)

full_results <- analysis_grid %>%
  pmap_dfr(
    function(pgv_var, PGV_group, fh_var, Family_history) {
      calculate_performance(
        df = data,
        pgv_var = pgv_var,
        pgv_label = PGV_group,
        fh_var = fh_var,
        fh_label = Family_history
      )
    }
  )

display_table <- full_results %>%
  mutate(
    Sensitivity = ifelse(
      is.na(Sensitivity_percent), "NE", sprintf("%.1f", Sensitivity_percent)
    ),
    Specificity = ifelse(
      is.na(Specificity_percent), "NE", sprintf("%.1f", Specificity_percent)
    ),
    PPV = ifelse(is.na(PPV_percent), "NE", sprintf("%.2f", PPV_percent)),
    NPV = ifelse(is.na(NPV_percent), "NE", sprintf("%.2f", NPV_percent)),
    `False-negative proportion` = ifelse(
      is.na(False_negative_percent),
      "NE",
      sprintf("%.1f", False_negative_percent)
    )
  ) %>%
  select(
    `PGV group` = PGV_group,
    `Family history` = Family_history,
    TP, FP, FN, TN,
    `Sensitivity (%)` = Sensitivity,
    `Specificity (%)` = Specificity,
    `PPV (%)` = PPV,
    `NPV (%)` = NPV,
    `False-negative proportion (%)` = `False-negative proportion`
  )

qc_summary <- tibble(
  Check = c(
    "Controls analyzed", "Any PGV carriers", "BRCA1/2 carriers",
    "BRCA1 carriers", "BRCA2 carriers", "HRR carriers",
    "HRR carriers excluding BRCA1/2", "MMR carriers"
  ),
  N = c(
    nrow(data),
    sum(data$pgv_any == 1L, na.rm = TRUE),
    sum(data$pgv_brca1_2 == 1L, na.rm = TRUE),
    sum(data$pgv_brca1 == 1L, na.rm = TRUE),
    sum(data$pgv_brca2 == 1L, na.rm = TRUE),
    sum(data$pgv_hrr == 1L, na.rm = TRUE),
    sum(data$pgv_hrr_without_brca == 1L, na.rm = TRUE),
    sum(data$pgv_mmr == 1L, na.rm = TRUE)
  )
)

notes <- tibble(
  Item = c(
    "Population", "Test definition", "Outcome definition",
    "True positive", "False positive", "False negative", "True negative",
    "Sensitivity", "Specificity", "Positive predictive value",
    "Negative predictive value", "False-negative proportion",
    "HRR excluding BRCA1/2"
  ),
  Definition = c(
    "Cancer-free controls included in the final analysis dataset.",
    "Positive family history in the specified category.",
    "PGV carrier status in the specified gene group.",
    "FH-positive PGV carrier.",
    "FH-positive non-carrier.",
    "FH-negative PGV carrier.",
    "FH-negative non-carrier.",
    "TP / (TP + FN).",
    "TN / (TN + FP).",
    "TP / (TP + FP).",
    "TN / (TN + FN).",
    "FN / (TP + FN), equivalent to 100% minus sensitivity.",
    "ATM, BARD1, BRIP1, CHEK2, PALB2, RAD51C, and RAD51D."
  )
)

write_xlsx(
  list(
    Supplementary_Table_S6 = display_table,
    Full_numeric_results = full_results,
    Quality_control = qc_summary,
    Notes = notes
  ),
  output_file
)

cat("\nDiagnostic-performance analysis completed.\n")
cat("Input file: ", input_file, "\n", sep = "")
cat("Output file: ", output_file, "\n", sep = "")
cat("Controls analyzed: ", nrow(data), "\n\n", sep = "")

cat("Key BRCA1/2 examples:\n")

full_results %>%
  filter(
    PGV_group == "BRCA1/2",
    Family_history %in% c("Breast cancer", "Ovarian cancer")
  ) %>%
  select(
    PGV_group, Family_history, TP, FP, FN, TN,
    Sensitivity_percent, Specificity_percent,
    PPV_percent, NPV_percent, False_negative_percent
  ) %>%
  print(n = Inf)
