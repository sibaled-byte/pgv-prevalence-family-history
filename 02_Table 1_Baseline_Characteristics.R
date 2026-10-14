# Table 1: Baseline characteristics of cancer-free controls
#
# This script uses the final saved analysis dataset and reproduces
# the structure of the manuscript Table 1 with columns for:
#   - All controls
#   - Female controls
#   - Male controls
#
# Values are presented as n (%) unless otherwise specified.

library(haven)
library(dplyr)
library(tidyr)
library(writexl)

input_file <- "ACC_PV_FH_controls_20260803.dta"
output_file <- "Table1_baseline_characteristics.xlsx"

data <- read_dta(input_file) %>%
  mutate(
    case_control = case_when(
      as.character(CaCo_status) %in% c("0", "control", "Control") ~ "control",
      as.character(CaCo_status) %in% c("1", "case", "Case") ~ "case",
      TRUE ~ NA_character_
    ),
    sex_label = case_when(
      as.character(sex) %in% c("1", "Male", "male") ~ "Male",
      as.character(sex) %in% c("2", "Female", "female") ~ "Female",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(case_control == "control")

age_variable <- if ("fvage" %in% names(data)) {
  "fvage"
} else if ("age" %in% names(data)) {
  "age"
} else {
  stop("Neither fvage nor age was found in the analysis dataset.")
}

data <- data %>%
  mutate(
    age_value = as.numeric(.data[[age_variable]]),
    age_group = case_when(
      age_value < 40 ~ "<40",
      age_value >= 40 & age_value < 50 ~ "40–49",
      age_value >= 50 & age_value < 60 ~ "50–59",
      age_value >= 60 & age_value < 70 ~ "60–69",
      age_value >= 70 ~ "≥70",
      TRUE ~ NA_character_
    ),
    age_group = factor(
      age_group,
      levels = c("<40", "40–49", "50–59", "60–69", "≥70")
    ),
    pv_binary = case_when(
      gene_clinlof == "No_carrier" ~ 0L,
      gene_clinlof == "QC out" ~ NA_integer_,
      is.na(gene_clinlof) ~ NA_integer_,
      TRUE ~ 1L
    )
  ) %>%
  filter(!is.na(pv_binary))

brca_genes <- c("BRCA1", "BRCA2")

hrr_genes <- c(
  "ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1",
  "CHEK2", "PALB2", "RAD51C", "RAD51D"
)

mmr_genes <- c("MLH1", "MSH2", "MSH6", "PMS2")

data <- data %>%
  mutate(
    pgv_all = pv_binary,
    pgv_brca = as.integer(
      pv_binary == 1L & gene_clinlof %in% brca_genes
    ),
    pgv_hrr = as.integer(
      pv_binary == 1L & gene_clinlof %in% hrr_genes
    ),
    pgv_mmr = as.integer(
      pv_binary == 1L & gene_clinlof %in% mmr_genes
    )
  )

format_n_percent <- function(n, denominator) {
  if (is.na(denominator) || denominator == 0) {
    return(NA_character_)
  }
  
  paste0(
    format(n, big.mark = ",", scientific = FALSE, trim = TRUE),
    " (",
    sprintf("%.1f", 100 * n / denominator),
    ")"
  )
}

format_mean_sd <- function(x) {
  if (all(is.na(x))) {
    return(NA_character_)
  }
  
  paste0(
    sprintf("%.1f", mean(x, na.rm = TRUE)),
    " (",
    sprintf("%.1f", sd(x, na.rm = TRUE)),
    ")"
  )
}

summarise_binary <- function(df, variable) {
  denominator <- sum(!is.na(df[[variable]]))
  n_positive <- sum(df[[variable]] == 1L, na.rm = TRUE)
  
  format_n_percent(n_positive, denominator)
}

summarise_age_group <- function(df, level) {
  denominator <- sum(!is.na(df$age_group))
  n_group <- sum(df$age_group == level, na.rm = TRUE)
  
  format_n_percent(n_group, denominator)
}

summarise_any_fh_no <- function(df) {
  denominator <- sum(!is.na(df$fh_any_fdr))
  n_no <- sum(df$fh_any_fdr == 0L, na.rm = TRUE)
  
  format_n_percent(n_no, denominator)
}

summarise_any_fh_yes <- function(df) {
  denominator <- sum(!is.na(df$fh_any_fdr))
  n_yes <- sum(df$fh_any_fdr == 1L, na.rm = TRUE)
  
  format_n_percent(n_yes, denominator)
}

make_table1_column <- function(df, column_name) {
  tibble(
    Characteristic = c(
      "Age (years)",
      "<40",
      "40–49",
      "50–59",
      "60–69",
      "≥70",
      "Mean age, years (SD)",
      "Any family history of cancer",
      "No",
      "Yes",
      "Family history by site",
      "Breast cancer",
      "Ovarian cancer",
      "Pancreatic cancer",
      "Prostate cancer",
      "Colorectal cancer",
      "Gastric cancer",
      "Lung cancer",
      "Other cancers",
      "PGV carrier status",
      "Non-carriers",
      "All PGVs",
      "BRCA1/2",
      "HRR",
      "MMR"
    ),
    value = c(
      "",
      summarise_age_group(df, "<40"),
      summarise_age_group(df, "40–49"),
      summarise_age_group(df, "50–59"),
      summarise_age_group(df, "60–69"),
      summarise_age_group(df, "≥70"),
      format_mean_sd(df$age_value),
      "",
      summarise_any_fh_no(df),
      summarise_any_fh_yes(df),
      "",
      summarise_binary(df, "fh_breast_fdr"),
      summarise_binary(df, "fh_ovarian_fdr"),
      summarise_binary(df, "fh_pancreas_fdr"),
      summarise_binary(df, "fh_prostate_fdr"),
      summarise_binary(df, "fh_colon_fdr"),
      summarise_binary(df, "fh_gastric_fdr"),
      summarise_binary(df, "fh_lung_fdr"),
      summarise_binary(df, "fh_other_fdr"),
      "",
      format_n_percent(
        sum(df$pgv_all == 0L, na.rm = TRUE),
        sum(!is.na(df$pgv_all))
      ),
      summarise_binary(df, "pgv_all"),
      summarise_binary(df, "pgv_brca"),
      summarise_binary(df, "pgv_hrr"),
      summarise_binary(df, "pgv_mmr")
    )
  ) %>%
    rename(!!column_name := value)
}

all_controls <- data
female_controls <- data %>% filter(sex_label == "Female")
male_controls <- data %>% filter(sex_label == "Male")

table1 <- make_table1_column(
  all_controls,
  paste0("All (N = ", format(nrow(all_controls), big.mark = ","), ")")
) %>%
  left_join(
    make_table1_column(
      female_controls,
      paste0(
        "Female (N = ",
        format(nrow(female_controls), big.mark = ","),
        ")"
      )
    ),
    by = "Characteristic"
  ) %>%
  left_join(
    make_table1_column(
      male_controls,
      paste0(
        "Male (N = ",
        format(nrow(male_controls), big.mark = ","),
        ")"
      )
    ),
    by = "Characteristic"
  )


# P values comparing female and male participants
# Pearson's chi-squared or Fisher's exact test

get_p <- function(var) {
  tab <- table(data$sex_label, data[[var]])
  
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  
  if (any(chi$expected < 5)) {
    p <- fisher.test(tab)$p.value
  } else {
    p <- chi$p.value
  }
  
  return(p)
}

# Calculate P values
p_values <- c(
  "Age group" = get_p("age_group"),
  "Any family history" = get_p("fh_any_fdr"),
  "Breast cancer" = get_p("fh_breast_fdr"),
  "Ovarian cancer" = get_p("fh_ovarian_fdr"),
  "Pancreatic cancer" = get_p("fh_pancreas_fdr"),
  "Prostate cancer" = get_p("fh_prostate_fdr"),
  "Colorectal cancer" = get_p("fh_colon_fdr"),
  "Gastric cancer" = get_p("fh_gastric_fdr"),
  "Lung cancer" = get_p("fh_lung_fdr"),
  "Other cancers" = get_p("fh_other_fdr"),
  "All PGVs" = get_p("pgv_all"),
  "BRCA1/2" = get_p("pgv_brca"),
  "HRR" = get_p("pgv_hrr"),
  "MMR" = get_p("pgv_mmr")
)

# Format P values
format_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# Add P values to Table 1
table1$P_value <- ""

row_names <- c(
  "Age (years)", "Any family history of cancer",
  "Breast cancer", "Ovarian cancer",
  "Pancreatic cancer", "Prostate cancer",
  "Colorectal cancer", "Gastric cancer",
  "Lung cancer", "Other cancers",
  "All PGVs", "BRCA1/2", "HRR", "MMR"
)

stopifnot(length(row_names) == length(p_values))
stopifnot(all(row_names %in% table1$Characteristic))

table1$P_value[match(row_names, table1$Characteristic)] <-
  vapply(p_values, format_p, character(1))

# Mean age (SD) intentionally has no P value

missing_summary <- tibble(
  Variable = c(
    "Age",
    "Sex",
    "Any family history",
    "Breast family history",
    "Ovarian family history",
    "Pancreatic family history",
    "Prostate family history",
    "Colorectal family history",
    "Gastric family history",
    "Lung family history",
    "Other family history",
    "PGV carrier status"
  ),
  Missing_N = c(
    sum(is.na(data$age_value)),
    sum(is.na(data$sex_label)),
    sum(is.na(data$fh_any_fdr)),
    sum(is.na(data$fh_breast_fdr)),
    sum(is.na(data$fh_ovarian_fdr)),
    sum(is.na(data$fh_pancreas_fdr)),
    sum(is.na(data$fh_prostate_fdr)),
    sum(is.na(data$fh_colon_fdr)),
    sum(is.na(data$fh_gastric_fdr)),
    sum(is.na(data$fh_lung_fdr)),
    sum(is.na(data$fh_other_fdr)),
    sum(is.na(data$pgv_all))
  )
)

table_notes <- tibble(
  Item = c(
    "Population",
    "Any family history",
    "Site-specific family history",
    "Values",
    "PGV",
    "HRR",
    "MMR"
  ),
  Definition = c(
    "Cancer-free controls included in the final analysis dataset.",
    "Cancer reported in at least one parent or sibling on the baseline questionnaire.",
    "Cancer site reported for a parent or sibling and harmonised to three-character ICD-10 categories.",
    "n (%) unless otherwise specified.",
    "Pathogenic germline variant.",
    "Homologous recombination repair.",
    "Mismatch repair."
  )
)

write_xlsx(
  list(
    Table1 = table1,
    Missing_data = missing_summary,
    Notes = table_notes
  ),
  output_file
)

cat(
  "\nTable 1 completed.\n",
  "Input: ", input_file, "\n",
  "Output: ", output_file, "\n",
  "All controls: ", nrow(all_controls), "\n",
  "Female controls: ", nrow(female_controls), "\n",
  "Male controls: ", nrow(male_controls), "\n",
  sep = ""
)