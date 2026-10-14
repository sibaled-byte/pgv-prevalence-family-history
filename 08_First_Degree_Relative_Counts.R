# The questionnaire records the total sibship size, including the participant.
# Therefore, the number of siblings excluding the participant is:
#
# n_siblings = reported sibship size - 1
#
# The number of first-degree relatives covered by the family history
# assessment equals two parents plus siblings excluding the participant.

library(haven)
library(dplyr)
library(writexl)

herpacc2_file <- "herpacc2_familyhistory.dta"
herpacc3_file <- "herpacc3_familyhistory.dta"
analysis_file <- "ACC_PV_FH_controls_20260803.dta"

output_file <- "First_Degree_Relative_Counts.xlsx"

herpacc2_raw <- read_dta(herpacc2_file)
herpacc3_raw <- read_dta(herpacc3_file)
analysis_data <- read_dta(analysis_file)

# Total sibship-size variables.
# These variables include the participant in the reported count.
h2_sibship_var <- "兄弟姉妹人数"
h3_sibship_var <- "兄弟姉妹"

required_h2 <- c("acc", h2_sibship_var)
required_h3 <- c("acc", h3_sibship_var)

missing_h2 <- setdiff(required_h2, names(herpacc2_raw))
missing_h3 <- setdiff(required_h3, names(herpacc3_raw))

if (length(missing_h2) > 0) {
  stop(
    "Missing HERPACC-2 variable(s): ",
    paste(missing_h2, collapse = ", ")
  )
}

if (length(missing_h3) > 0) {
  stop(
    "Missing HERPACC-3 variable(s): ",
    paste(missing_h3, collapse = ", ")
  )
}

# Extract sibship size and derive the number of siblings excluding
# the participant.

herpacc2_relatives <- herpacc2_raw %>%
  transmute(
    acc,
    questionnaire_phase = 2L,
    sibship_size_reported = suppressWarnings(
      as.numeric(.data[[h2_sibship_var]])
    )
  )

herpacc3_relatives <- herpacc3_raw %>%
  transmute(
    acc,
    questionnaire_phase = 3L,
    sibship_size_reported = suppressWarnings(
      as.numeric(.data[[h3_sibship_var]])
    )
  )

relative_data <- bind_rows(
  herpacc2_relatives,
  herpacc3_relatives
) %>%
  mutate(
    # Treat special codes and impossible values as missing.
    sibship_size_reported = case_when(
      is.na(sibship_size_reported) ~ NA_real_,
      sibship_size_reported %in% c(88, 98, 99, 888, 999) ~ NA_real_,
      sibship_size_reported < 1 ~ NA_real_,
      TRUE ~ sibship_size_reported
    ),
    
    # The reported sibship size includes the participant.
    n_siblings = case_when(
      is.na(sibship_size_reported) ~ NA_real_,
      TRUE ~ sibship_size_reported - 1
    ),
    
    # First-degree relatives included in the FH assessment:
    # two parents plus siblings excluding the participant.
    n_first_degree_relatives = case_when(
      is.na(n_siblings) ~ NA_real_,
      TRUE ~ 2 + n_siblings
    )
  )

if (anyDuplicated(relative_data$acc) > 0) {
  stop(
    "Duplicate participant IDs were found in the combined relative data."
  )
}

# Restrict to the cancer-free controls used in the manuscript.

controls <- analysis_data %>%
  mutate(
    case_control = case_when(
      as.character(CaCo_status) %in%
        c("0", "control", "Control") ~ "control",
      
      as.character(CaCo_status) %in%
        c("1", "case", "Case") ~ "case",
      
      TRUE ~ NA_character_
    ),
    
    sex_label = case_when(
      as.character(sex) %in%
        c("1", "Male", "male") ~ "Male",
      
      as.character(sex) %in%
        c("2", "Female", "female") ~ "Female",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(case_control == "control") %>%
  left_join(
    relative_data,
    by = "acc"
  )

if (nrow(controls) != 11811) {
  warning(
    "The control count is ",
    nrow(controls),
    "; the expected manuscript count is 11,811."
  )
}

# Confirm that questionnaire phase agrees with the analysis dataset.

phase_check <- controls %>%
  summarise(
    N_with_both_phase_values = sum(
      !is.na(hpc) &
        !is.na(questionnaire_phase)
    ),
    
    N_phase_disagreement = sum(
      hpc != questionnaire_phase,
      na.rm = TRUE
    )
  )

summarise_relatives <- function(data) {
  
  nonmissing <- data$n_first_degree_relatives[
    !is.na(data$n_first_degree_relatives)
  ]
  
  tibble(
    N_total = nrow(data),
    
    N_nonmissing = length(nonmissing),
    
    N_missing = sum(
      is.na(data$n_first_degree_relatives)
    ),
    
    Missing_percent = round(
      100 *
        sum(is.na(data$n_first_degree_relatives)) /
        nrow(data),
      1
    ),
    
    Mean = ifelse(
      length(nonmissing) > 0,
      mean(nonmissing),
      NA_real_
    ),
    
    SD = ifelse(
      length(nonmissing) > 1,
      sd(nonmissing),
      NA_real_
    ),
    
    Median = ifelse(
      length(nonmissing) > 0,
      median(nonmissing),
      NA_real_
    ),
    
    Q1 = ifelse(
      length(nonmissing) > 0,
      quantile(
        nonmissing,
        0.25,
        names = FALSE
      ),
      NA_real_
    ),
    
    Q3 = ifelse(
      length(nonmissing) > 0,
      quantile(
        nonmissing,
        0.75,
        names = FALSE
      ),
      NA_real_
    ),
    
    Minimum = ifelse(
      length(nonmissing) > 0,
      min(nonmissing),
      NA_real_
    ),
    
    Maximum = ifelse(
      length(nonmissing) > 0,
      max(nonmissing),
      NA_real_
    )
  )
}

overall_summary <- summarise_relatives(controls) %>%
  mutate(
    Group = "Overall",
    .before = 1
  )

sex_summary <- controls %>%
  filter(!is.na(sex_label)) %>%
  group_split(sex_label) %>%
  lapply(
    function(df) {
      summarise_relatives(df) %>%
        mutate(
          Group = unique(df$sex_label),
          .before = 1
        )
    }
  ) %>%
  bind_rows()

phase_summary <- controls %>%
  filter(!is.na(hpc)) %>%
  group_split(hpc) %>%
  lapply(
    function(df) {
      summarise_relatives(df) %>%
        mutate(
          Group = paste0(
            "HERPACC-",
            unique(df$hpc)
          ),
          .before = 1
        )
    }
  ) %>%
  bind_rows()

summary_table <- bind_rows(
  overall_summary,
  sex_summary,
  phase_summary
) %>%
  mutate(
    `Mean (SD)` = ifelse(
      is.na(Mean),
      NA_character_,
      paste0(
        sprintf("%.2f", Mean),
        " (",
        sprintf("%.2f", SD),
        ")"
      )
    ),
    
    `Median (IQR)` = ifelse(
      is.na(Median),
      NA_character_,
      paste0(
        sprintf("%.2f", Median),
        " (",
        sprintf("%.2f", Q1),
        "-",
        sprintf("%.2f", Q3),
        ")"
      )
    ),
    
    Range = ifelse(
      is.na(Minimum),
      NA_character_,
      paste0(
        Minimum,
        "-",
        Maximum
      )
    )
  ) %>%
  select(
    Group,
    N_total,
    N_nonmissing,
    N_missing,
    Missing_percent,
    `Mean (SD)`,
    `Median (IQR)`,
    Range
  )

missingness_by_sex_phase <- controls %>%
  group_by(
    sex_label,
    hpc
  ) %>%
  summarise(
    N_total = n(),
    
    N_missing = sum(
      is.na(n_first_degree_relatives)
    ),
    
    Missing_percent = round(
      100 * N_missing / N_total,
      1
    ),
    
    .groups = "drop"
  )

sibling_frequency <- controls %>%
  count(
    n_siblings,
    name = "N"
  ) %>%
  mutate(
    Percent = round(
      100 * N / nrow(controls),
      2
    )
  ) %>%
  arrange(n_siblings)

relative_frequency <- controls %>%
  count(
    n_first_degree_relatives,
    name = "N"
  ) %>%
  mutate(
    Percent = round(
      100 * N / nrow(controls),
      2
    )
  ) %>%
  arrange(n_first_degree_relatives)

quality_control <- tibble(
  Check = c(
    "Final control sample",
    "Matched to questionnaire relative data",
    "Missing reported sibship size",
    "Missing sibling count",
    "Missing first-degree-relative count",
    "Sibling count below zero",
    "First-degree-relative count below two",
    "HERPACC phase disagreement"
  ),
  
  N = c(
    nrow(controls),
    
    sum(
      !is.na(controls$questionnaire_phase)
    ),
    
    sum(
      is.na(controls$sibship_size_reported)
    ),
    
    sum(
      is.na(controls$n_siblings)
    ),
    
    sum(
      is.na(controls$n_first_degree_relatives)
    ),
    
    sum(
      controls$n_siblings < 0,
      na.rm = TRUE
    ),
    
    sum(
      controls$n_first_degree_relatives < 2,
      na.rm = TRUE
    ),
    
    phase_check$N_phase_disagreement
  )
)

write_xlsx(
  list(
    Summary = summary_table,
    Missingness_by_sex_phase = missingness_by_sex_phase,
    Sibling_distribution = sibling_frequency,
    Relative_distribution = relative_frequency,
    Quality_control = quality_control
  ),
  output_file
)

cat("\nSummary of first-degree relatives:\n")
print(summary_table, n = Inf)

cat("\nMissingness by sex and HERPACC phase:\n")
print(missingness_by_sex_phase, n = Inf)

cat("\nSibling distribution excluding the participant:\n")
print(sibling_frequency, n = Inf)

cat("\nQuality-control checks:\n")
print(quality_control, n = Inf)

cat(
  "\nResults saved to:\n",
  output_file,
  "\n"
)
