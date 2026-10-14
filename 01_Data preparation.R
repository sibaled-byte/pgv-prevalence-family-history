# Pathogenic germline variant prevalence and association with
# site-specific family history
# Data preparation script

library(haven)
library(dplyr)
library(writexl)

# Input files
herpacc2_file   <- "herpacc2_familyhistory.dta"
herpacc3_file   <- "herpacc3_familyhistory.dta"
pv_file         <- "familyhistory_pv_hpc.dta"
environment_file <- "herpacc_envi.dta"
sex_check_file  <- "sex_check.dta"

# Output files
output_data <- "ACC_PV_FH_controls_20260803.dta"
output_qc      <- "ACC_PGV_FH_prevalence_QC_questionnaire_20260803.xlsx"
output_session <- "ACC_PGV_FH_prevalence_session_info_20260803.txt"

# Check that all input files are available
input_files <- c(
  herpacc2_file,
  herpacc3_file,
  pv_file,
  environment_file,
  sex_check_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing input file(s):",
      paste(missing_files, collapse = ", ")
    )
  )
}

# Clean ICD-10 values and keep the three-character cancer-site code
clean_icd10 <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "0")] <- NA_character_
  
  numeric_code <- grepl("^\\d{1,2}$", x)
  x[numeric_code] <- paste0("C", x[numeric_code])
  
  one_digit_code <- grepl("^C\\d$", x)
  x[one_digit_code] <- paste0(
    "C0",
    sub("^C", "", x[one_digit_code])
  )
  
  valid_code <- grepl("^C\\d{2}", x)
  
  result <- rep(NA_character_, length(x))
  result[valid_code] <- substr(x[valid_code], 1, 3)
  
  result
}

# Check whether any first-degree relative had a specified cancer
has_site <- function(x, codes) {
  any(
    unlist(x, use.names = FALSE) %in% codes,
    na.rm = TRUE
  )
}

# Check that participant IDs are unique
check_unique_acc <- function(data, name) {
  if (anyDuplicated(data$acc) > 0) {
    stop(paste("Duplicate acc values found in", name))
  }
}

# Check that a merge did not change the number of rows
check_same_rows <- function(before, after, name) {
  if (nrow(before) != nrow(after)) {
    stop(
      paste(
        name,
        "changed the row count from",
        nrow(before),
        "to",
        nrow(after)
      )
    )
  }
}

# ICD-10 cancer-site definitions
icd_breast   <- "C50"
icd_gastric  <- "C16"
icd_colon    <- paste0("C", 18:20)
icd_lung     <- "C34"
icd_ovarian  <- "C56"
icd_pancreas <- "C25"
icd_prostate <- "C61"

# Read HERPACC family-history data
herpacc2 <- read_dta(herpacc2_file) %>%
  mutate(hpc = 2L)

herpacc3 <- read_dta(herpacc3_file) %>%
  mutate(hpc = 3L)

check_unique_acc(herpacc2, "HERPACC-2 family-history data")
check_unique_acc(herpacc3, "HERPACC-3 family-history data")

source_counts <- tibble(
  dataset = c("HERPACC-2", "HERPACC-3"),
  participants = c(nrow(herpacc2), nrow(herpacc3))
)

# Create HERPACC-2 family-history variables
herpacc2 <- herpacc2 %>%
  mutate(
    father_cancer  = `父親がん`,
    mother_cancer  = `母親がん`,
    sibling_cancer = `兄弟姉妹のがん`,
    
    father_icd = clean_icd10(`父親がんCD`),
    mother_icd = clean_icd10(`母親がんCD`),
    sib1_icd   = clean_icd10(`兄弟姉妹_診断CD1`),
    sib2_icd   = clean_icd10(`兄弟姉妹_診断CD2`),
    sib3_icd   = clean_icd10(`兄弟姉妹_診断CD3`),
    sib4_icd   = clean_icd10(`兄弟姉妹_診断CD4`)
  ) %>%
  rowwise() %>%
  mutate(
    family_icd = list(c(
      father_icd,
      mother_icd,
      sib1_icd,
      sib2_icd,
      sib3_icd,
      sib4_icd
    )),
    fh_any_fdr = as.integer(
      father_cancer == 3 |
        mother_cancer == 3 |
        sibling_cancer == 2
    ),
    fh_breast_fdr    = as.integer(has_site(family_icd, icd_breast)),
    fh_gastric_fdr   = as.integer(has_site(family_icd, icd_gastric)),
    fh_colon_fdr     = as.integer(has_site(family_icd, icd_colon)),
    fh_lung_fdr      = as.integer(has_site(family_icd, icd_lung)),
    fh_ovarian_fdr   = as.integer(has_site(family_icd, icd_ovarian)),
    fh_pancreas_fdr  = as.integer(has_site(family_icd, icd_pancreas)),
    fh_prostate_fdr  = as.integer(has_site(family_icd, icd_prostate))
  ) %>%
  ungroup() %>%
  select(-family_icd)

# Create HERPACC-3 family-history variables
herpacc3 <- herpacc3 %>%
  mutate(
    father_cancer  = `父親_がん`,
    mother_cancer  = `母親_がん`,
    sibling_cancer = `同胞のがん`,
    
    father1_icd = clean_icd10(`父親_がん1_診断部位CD`),
    father2_icd = clean_icd10(`父親_がん2_診断部位CD`),
    mother1_icd = clean_icd10(`母親_がん1_診断部位CD`),
    mother2_icd = clean_icd10(`母親_がん2_診断部位CD`),
    sib1_1_icd  = clean_icd10(`同胞1_診断部位CD1`),
    sib1_2_icd  = clean_icd10(`同胞1_診断部位CD2`),
    sib2_1_icd  = clean_icd10(`同胞2_診断部位CD1`),
    sib2_2_icd  = clean_icd10(`同胞2_診断部位CD2`),
    sib3_1_icd  = clean_icd10(`同胞3_診断部位CD1`),
    sib3_2_icd  = clean_icd10(`同胞3_診断部位CD2`),
    sib4_1_icd  = clean_icd10(`同胞4_診断部位CD1`),
    sib4_2_icd  = clean_icd10(`同胞4_診断部位CD2`)
  ) %>%
  rowwise() %>%
  mutate(
    family_icd = list(c(
      father1_icd,
      father2_icd,
      mother1_icd,
      mother2_icd,
      sib1_1_icd,
      sib1_2_icd,
      sib2_1_icd,
      sib2_2_icd,
      sib3_1_icd,
      sib3_2_icd,
      sib4_1_icd,
      sib4_2_icd
    )),
    fh_any_fdr = as.integer(
      father_cancer == 3 |
        mother_cancer == 3 |
        sibling_cancer == 2
    ),
    fh_breast_fdr    = as.integer(has_site(family_icd, icd_breast)),
    fh_gastric_fdr   = as.integer(has_site(family_icd, icd_gastric)),
    fh_colon_fdr     = as.integer(has_site(family_icd, icd_colon)),
    fh_lung_fdr      = as.integer(has_site(family_icd, icd_lung)),
    fh_ovarian_fdr   = as.integer(has_site(family_icd, icd_ovarian)),
    fh_pancreas_fdr  = as.integer(has_site(family_icd, icd_pancreas)),
    fh_prostate_fdr  = as.integer(has_site(family_icd, icd_prostate))
  ) %>%
  ungroup() %>%
  select(-family_icd)

# FH coding follows the original analysis. Participants with completely
# missing parental and sibling FH data were not included in the analytic cohort.
# Family-size missingness is handled separately.
fh_variables <- c(
  "fh_any_fdr",
  "fh_breast_fdr",
  "fh_gastric_fdr",
  "fh_colon_fdr",
  "fh_lung_fdr",
  "fh_ovarian_fdr",
  "fh_pancreas_fdr",
  "fh_prostate_fdr"
)

# Apply the prespecified FH coding used for the manuscript.
herpacc2 <- herpacc2 %>%
  mutate(
    across(
      all_of(fh_variables),
      ~ ifelse(is.na(.x), 0L, .x)
    )
  )

herpacc3 <- herpacc3 %>%
  mutate(
    across(
      all_of(fh_variables),
      ~ ifelse(is.na(.x), 0L, .x)
    )
  )

# Combine HERPACC-2 and HERPACC-3
family_data <- bind_rows(herpacc2, herpacc3)

check_unique_acc(
  family_data,
  "combined HERPACC family-history data"
)

# Create the other-cancer family-history category
family_data <- family_data %>%
  mutate(
    fh_other_fdr = as.integer(
      fh_any_fdr == 1L &
        !(
          fh_breast_fdr == 1L |
            fh_gastric_fdr == 1L |
            fh_colon_fdr == 1L |
            fh_lung_fdr == 1L |
            fh_ovarian_fdr == 1L |
            fh_pancreas_fdr == 1L |
            fh_prostate_fdr == 1L
        )
    )
  )

# Check cleaned ICD-10 values
icd_columns <- names(family_data)[grepl("icd$", names(family_data))]

all_icd <- na.omit(
  unlist(
    family_data[icd_columns],
    use.names = FALSE
  )
)

invalid_icd <- unique(
  all_icd[!grepl("^C\\d{2}$", all_icd)]
)

if (length(invalid_icd) > 0) {
  stop(
    paste(
      "Unexpected ICD-10 values remain:",
      paste(invalid_icd, collapse = ", ")
    )
  )
}

icd_frequency <- as.data.frame(
  table(all_icd),
  stringsAsFactors = FALSE
) %>%
  rename(
    ICD10_site = all_icd,
    N = Freq
  ) %>%
  arrange(ICD10_site)

family_history_summary <- family_data %>%
  summarise(
    N_total = n(),
    Any_cancer_FH = sum(fh_any_fdr == 1L, na.rm = TRUE),
    Breast_FH = sum(fh_breast_fdr == 1L, na.rm = TRUE),
    Gastric_FH = sum(fh_gastric_fdr == 1L, na.rm = TRUE),
    Colorectal_FH = sum(fh_colon_fdr == 1L, na.rm = TRUE),
    Lung_FH = sum(fh_lung_fdr == 1L, na.rm = TRUE),
    Ovarian_FH = sum(fh_ovarian_fdr == 1L, na.rm = TRUE),
    Pancreatic_FH = sum(fh_pancreas_fdr == 1L, na.rm = TRUE),
    Prostate_FH = sum(fh_prostate_fdr == 1L, na.rm = TRUE),
    Other_FH = sum(fh_other_fdr == 1L, na.rm = TRUE)
  )

# Keep only the family-history variables needed for analysis
family_data <- family_data %>%
  select(
    acc,
    hpc,
    fh_any_fdr,
    fh_breast_fdr,
    fh_gastric_fdr,
    fh_colon_fdr,
    fh_lung_fdr,
    fh_ovarian_fdr,
    fh_pancreas_fdr,
    fh_prostate_fdr,
    fh_other_fdr
  )

# Read PV data and remove older FH variables
pv_data <- read_dta(pv_file)

check_unique_acc(pv_data, "PV source data")

old_fh_variables <- c(
  "any_family",
  "breast_family",
  "gastric_family",
  "colon_family",
  "lung_family",
  "ovarian_family",
  "pancreatic_family",
  "prostate_family",
  "other_family"
)

pv_data <- pv_data %>%
  select(-any_of(old_fh_variables))

study_data <- pv_data %>%
  left_join(
    family_data,
    by = "acc"
  )

check_same_rows(
  pv_data,
  study_data,
  "Family-history merge"
)

# Merge environmental variables
environment_data <- read_dta(environment_file) %>%
  select(
    acc,
    smoker,
    py,
    menopause,
    drinker,
    quant_drink_gr,
    bmi
  )

check_unique_acc(environment_data, "environmental data")

before_environment <- study_data

study_data <- study_data %>%
  left_join(
    environment_data,
    by = "acc"
  )

check_same_rows(
  before_environment,
  study_data,
  "Environmental-data merge"
)

# Merge sex-QC flag
sex_check <- read_dta(sex_check_file) %>%
  transmute(
    acc,
    sex_flag = trimws(as.character(sex_check))
  )

check_unique_acc(sex_check, "sex-check data")

before_sex_check <- study_data

study_data <- study_data %>%
  left_join(
    sex_check,
    by = "acc"
  )

check_same_rows(
  before_sex_check,
  study_data,
  "Sex-check merge"
)

# Apply sample-QC and sex-mismatch exclusions
study_data <- study_data %>%
  mutate(
    sampleqc_value = trimws(as.character(sampleqc)),
    exclude_sex_mismatch = coalesce(sex_flag == "NG", FALSE),
    exclude_sample_qc =
      is.na(sampleqc_value) |
      sampleqc_value != "QC clear",
    exclude = exclude_sex_mismatch | exclude_sample_qc
  )

qc_summary <- study_data %>%
  summarise(
    N_before_exclusions = n(),
    Sex_mismatch = sum(exclude_sex_mismatch),
    Missing_sex_check = sum(is.na(sex_flag)),
    Sample_QC_failure = sum(exclude_sample_qc),
    Total_excluded = sum(exclude),
    N_after_exclusions = sum(!exclude)
  )

study_data <- study_data %>%
  filter(!exclude) %>%
  select(
    -sampleqc_value,
    -exclude_sex_mismatch,
    -exclude_sample_qc,
    -exclude
  )

# Keep cancer-free participants only
study_data <- study_data %>%
  mutate(
    case_control = trimws(tolower(as.character(CaCo_status)))
  ) %>%
  filter(case_control %in% c("0", "control"))

# Define PGV carrier status
study_data <- study_data %>%
  mutate(
    gene_name = trimws(as.character(gene_clinlof)),
    pv_binary = case_when(
      gene_name == "No_carrier" ~ 0L,
      is.na(gene_name) |
        gene_name == "" |
        gene_name == "QC out" ~ NA_integer_,
      TRUE ~ 1L
    )
  )

# Gene-group definitions used in the manuscript
breast_cancer_genes <- c(
  "ATM", "BRCA1", "BRCA2", "CHEK2",
  "NF1", "PALB2", "PTEN", "TP53"
)

hrr_genes <- c(
  "ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1",
  "CHEK2", "PALB2", "RAD51C", "RAD51D"
)

mmr_genes <- c(
  "MLH1", "MSH2", "MSH6", "PMS2"
)

study_data <- study_data %>%
  mutate(
    pv_bc = case_when(
      is.na(pv_binary) ~ NA_integer_,
      gene_name %in% breast_cancer_genes ~ 1L,
      TRUE ~ 0L
    ),
    pv_hrr = case_when(
      is.na(pv_binary) ~ NA_integer_,
      gene_name %in% hrr_genes ~ 1L,
      TRUE ~ 0L
    ),
    pv_brca = case_when(
      is.na(pv_binary) ~ NA_integer_,
      gene_name %in% c("BRCA1", "BRCA2") ~ 1L,
      TRUE ~ 0L
    ),
    pv_mmr = case_when(
      is.na(pv_binary) ~ NA_integer_,
      gene_name %in% mmr_genes ~ 1L,
      TRUE ~ 0L
    ),
    fh_breast_or_ovarian_fdr = case_when(
      fh_breast_fdr == 1L | fh_ovarian_fdr == 1L ~ 1L,
      fh_breast_fdr == 0L & fh_ovarian_fdr == 0L ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# Prepare covariates
study_data <- study_data %>%
  mutate(
    smoker_value = trimws(tolower(as.character(smoker))),
    smoking_status = case_when(
      smoker_value %in% c("0", "never") ~ "Never",
      smoker_value %in% c("1", "former", "current") ~ "Ever",
      TRUE ~ NA_character_
    ),
    smoker_binary = case_when(
      smoking_status == "Never" ~ 0L,
      smoking_status == "Ever" ~ 1L,
      TRUE ~ NA_integer_
    ),
    pack_years = case_when(
      smoker_binary == 0L ~ 0,
      smoker_binary == 1L ~ suppressWarnings(as.numeric(py)),
      TRUE ~ NA_real_
    ),
    drinker_value = suppressWarnings(as.numeric(drinker)),
    drinking_status = case_when(
      drinker_value == 1 ~ "Ever",
      drinker_value %in% c(2, 3) ~ "Never",
      drinker_value == 9 ~ NA_character_,
      TRUE ~ NA_character_
    ),
    drinker_binary = case_when(
      drinking_status == "Never" ~ 0L,
      drinking_status == "Ever" ~ 1L,
      TRUE ~ NA_integer_
    ),
    quant_drink_gr = case_when(
      drinker_binary == 0L ~ 0,
      drinker_binary == 1L ~ suppressWarnings(
        as.numeric(quant_drink_gr)
      ),
      TRUE ~ NA_real_
    ),
    bmi = suppressWarnings(as.numeric(bmi))
  )

# Final checks
check_unique_acc(study_data, "final prevalence-analysis dataset")

if (!all(as.character(study_data$sampleqc) == "QC clear")) {
  stop("The final dataset contains samples that did not pass QC")
}

final_summary <- study_data %>%
  summarise(
    N_total = n(),
    Women = sum(sex %in% c(2, "Female"), na.rm = TRUE),
    Men = sum(sex %in% c(1, "Male"), na.rm = TRUE),
    HERPACC2 = sum(hpc == 2, na.rm = TRUE),
    HERPACC3 = sum(hpc == 3, na.rm = TRUE),
    PGV_carriers = sum(pv_binary == 1L, na.rm = TRUE),
    Missing_PGV = sum(is.na(pv_binary)),
    Any_FH_positive = sum(fh_any_fdr == 1L, na.rm = TRUE),
    Breast_FH = sum(fh_breast_fdr == 1L, na.rm = TRUE),
    Ovarian_FH = sum(fh_ovarian_fdr == 1L, na.rm = TRUE),
    Pancreatic_FH = sum(fh_pancreas_fdr == 1L, na.rm = TRUE),
    Prostate_FH = sum(fh_prostate_fdr == 1L, na.rm = TRUE),
    Colorectal_FH = sum(fh_colon_fdr == 1L, na.rm = TRUE),
    Gastric_FH = sum(fh_gastric_fdr == 1L, na.rm = TRUE),
    Lung_FH = sum(fh_lung_fdr == 1L, na.rm = TRUE),
    Other_FH = sum(fh_other_fdr == 1L, na.rm = TRUE),
    Missing_all_FH = sum(
      if_all(starts_with("fh_"), is.na)
    )
  )

fh_by_phase <- study_data %>%
  group_by(hpc) %>%
  summarise(
    N_total = n(),
    Any_cancer_FH = sum(fh_any_fdr == 1L, na.rm = TRUE),
    Breast_FH = sum(fh_breast_fdr == 1L, na.rm = TRUE),
    Gastric_FH = sum(fh_gastric_fdr == 1L, na.rm = TRUE),
    Colorectal_FH = sum(fh_colon_fdr == 1L, na.rm = TRUE),
    Lung_FH = sum(fh_lung_fdr == 1L, na.rm = TRUE),
    Ovarian_FH = sum(fh_ovarian_fdr == 1L, na.rm = TRUE),
    Pancreatic_FH = sum(fh_pancreas_fdr == 1L, na.rm = TRUE),
    Prostate_FH = sum(fh_prostate_fdr == 1L, na.rm = TRUE),
    Other_FH = sum(fh_other_fdr == 1L, na.rm = TRUE),
    .groups = "drop"
  )

# Save the final analysis dataset and QC information
write_dta(
  study_data,
  output_data
)

write_xlsx(
  list(
    Source_counts = source_counts,
    Family_history_summary = family_history_summary,
    ICD10_frequency = icd_frequency,
    QC_exclusions = qc_summary,
    Final_sample = final_summary,
    FH_summary_by_phase = fh_by_phase
  ),
  output_qc
)

writeLines(
  capture.output(sessionInfo()),
  output_session
)

print(source_counts)
print(qc_summary)
print(final_summary)

message(
  paste0(
    "Data preparation completed.\n",
    "Analysis dataset: ", output_data, "\n",
    "QC workbook: ", output_qc, "\n",
    "Session information: ", output_session
  )
)
