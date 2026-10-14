
# Supplementary Table S4
# Sensitivity analyses for BRCA1, BRCA2, and
# HRR excluding BRCA1/2 according to FH and sex


library(haven)
library(dplyr)
library(tidyr)
library(purrr)
library(writexl)


# 1. Files


input_file <- "ACC_PV_FH_controls_20260803.dta"
output_file <- "Supplementary_Table_S4.xlsx"


# 2. Read final controls-only dataset


data <- read_dta(input_file) %>%
  mutate(
    sex_label = case_when(
      as.character(sex) %in% c("1", "Male", "male") ~ "Male",
      as.character(sex) %in% c("2", "Female", "female") ~ "Female",
      TRUE ~ NA_character_
    ),
    
    pv_binary = case_when(
      gene_clinlof == "No_carrier" ~ 0L,
      gene_clinlof == "QC out" ~ NA_integer_,
      is.na(gene_clinlof) ~ NA_integer_,
      TRUE ~ 1L
    )
  ) %>%
  filter(
    !is.na(pv_binary),
    !is.na(sex_label)
  )

# Age variable
age_variable <- if ("fvage" %in% names(data)) {
  "fvage"
} else if ("age" %in% names(data)) {
  "age"
} else {
  stop("Neither fvage nor age was found.")
}


# 3. Define sensitivity-analysis PGV groups


hrr_no_brca_genes <- c(
  "ATM",
  "BARD1",
  "BRIP1",
  "CHEK2",
  "PALB2",
  "RAD51C",
  "RAD51D"
)

data <- data %>%
  mutate(
    BRCA1 = as.integer(
      pv_binary == 1L & gene_clinlof == "BRCA1"
    ),
    
    BRCA2 = as.integer(
      pv_binary == 1L & gene_clinlof == "BRCA2"
    ),
    
    HRR_no_BRCA = as.integer(
      pv_binary == 1L &
        gene_clinlof %in% hrr_no_brca_genes
    )
  )

# Carrier-count quality-control check
cat("\nCarrier counts:\n")
cat("BRCA1:", sum(data$BRCA1), "\n")
cat("BRCA2:", sum(data$BRCA2), "\n")
cat("HRR excluding BRCA1/2:", sum(data$HRR_no_BRCA), "\n")


# 4. FH variables

family_history_map <- tibble(
  fh_var = c(
    "fh_any_fdr",
    "fh_breast_fdr",
    "fh_ovarian_fdr",
    "fh_pancreas_fdr",
    "fh_prostate_fdr",
    "fh_colon_fdr",
    "fh_gastric_fdr",
    "fh_lung_fdr",
    "fh_other_fdr"
  ),
  
  Family_history = c(
    "Any cancer",
    "Breast cancer",
    "Ovarian cancer",
    "Pancreatic cancer",
    "Prostate cancer",
    "Colorectal cancer",
    "Gastric cancer",
    "Lung cancer",
    "Other cancers"
  )
)

pgv_map <- tibble(
  carrier_var = c(
    "BRCA1",
    "BRCA2",
    "HRR_no_BRCA"
  ),
  
  PGV_group = c(
    "BRCA1",
    "BRCA2",
    "HRR†"
  )
)


# 5. FH prevalence function


make_prevalence_result <- function(df, carrier_var, fh_var) {
  
  analysis_df <- df %>%
    filter(
      !is.na(.data[[carrier_var]]),
      !is.na(.data[[fh_var]])
    )
  
  tab <- table(
    factor(
      analysis_df[[carrier_var]],
      levels = c(0, 1)
    ),
    factor(
      analysis_df[[fh_var]],
      levels = c(0, 1)
    )
  )
  
  n_fh_negative <- as.numeric(tab["1", "0"])
  n_fh_positive <- as.numeric(tab["1", "1"])
  
  N_fh_negative <- sum(tab[, "0"])
  N_fh_positive <- sum(tab[, "1"])
  
  tibble(
    n_FH_negative = n_fh_negative,
    N_FH_negative = N_fh_negative,
    
    prevalence_FH_negative = ifelse(
      N_fh_negative > 0,
      100 * n_fh_negative / N_fh_negative,
      NA_real_
    ),
    
    n_FH_positive = n_fh_positive,
    N_FH_positive = N_fh_positive,
    
    prevalence_FH_positive = ifelse(
      N_fh_positive > 0,
      100 * n_fh_positive / N_fh_positive,
      NA_real_
    )
  )
}


# 6. Adjusted logistic model


fit_adjusted_model <- function(
    df,
    carrier_var,
    fh_var,
    include_sex
) {
  
  covariates <- c(
    fh_var,
    age_variable,
    "hpc"
  )
  
  if (include_sex) {
    covariates <- c(
      covariates,
      "sex"
    )
  }
  
  analysis_df <- df %>%
    filter(
      !is.na(.data[[carrier_var]]),
      !is.na(.data[[fh_var]]),
      !is.na(.data[[age_variable]]),
      !is.na(hpc)
    )
  
  if (include_sex) {
    analysis_df <- analysis_df %>%
      filter(!is.na(sex))
  }
  
  # Need both carrier and non-carrier,
  # and both FH+ and FH−
  if (
    nrow(analysis_df) == 0 ||
    n_distinct(analysis_df[[carrier_var]]) < 2 ||
    n_distinct(analysis_df[[fh_var]]) < 2
  ) {
    return(
      tibble(
        OR = NA_real_,
        LCL = NA_real_,
        UCL = NA_real_,
        P = NA_real_,
        Model_N = nrow(analysis_df)
      )
    )
  }
  
  model <- suppressWarnings(
    glm(
      reformulate(
        covariates,
        response = carrier_var
      ),
      data = analysis_df,
      family = binomial()
    )
  )
  
  coefficient_table <- coef(summary(model))
  
  if (!fh_var %in% rownames(coefficient_table)) {
    return(
      tibble(
        OR = NA_real_,
        LCL = NA_real_,
        UCL = NA_real_,
        P = NA_real_,
        Model_N = nrow(analysis_df)
      )
    )
  }
  
  estimate <- coefficient_table[
    fh_var,
    "Estimate"
  ]
  
  standard_error <- coefficient_table[
    fh_var,
    "Std. Error"
  ]
  
  p_value <- coefficient_table[
    fh_var,
    "Pr(>|z|)"
  ]
  
  # Mark extreme / unstable estimates as NE
  OR <- exp(estimate)
  LCL <- exp(estimate - 1.96 * standard_error)
  UCL <- exp(estimate + 1.96 * standard_error)
  
  unstable <- (
    is.na(OR) ||
      is.infinite(OR) ||
      is.infinite(LCL) ||
      is.infinite(UCL) ||
      standard_error > 10
  )
  
  if (unstable) {
    return(
      tibble(
        OR = NA_real_,
        LCL = NA_real_,
        UCL = NA_real_,
        P = NA_real_,
        Model_N = nrow(analysis_df)
      )
    )
  }
  
  tibble(
    OR = OR,
    LCL = LCL,
    UCL = UCL,
    P = p_value,
    Model_N = nrow(analysis_df)
  )
}


# 7. Interaction model
#    OR = female-to-male ratio of ORs


fit_interaction_model <- function(
    df,
    carrier_var,
    fh_var
) {
  
  analysis_df <- df %>%
    filter(
      !is.na(.data[[carrier_var]]),
      !is.na(.data[[fh_var]]),
      !is.na(.data[[age_variable]]),
      !is.na(sex),
      !is.na(hpc)
    ) %>%
    mutate(
      fh_factor = factor(
        .data[[fh_var]],
        levels = c(0, 1)
      ),
      
      # Male = reference, Female = comparison
      sex_factor = factor(
        sex_label,
        levels = c("Male", "Female")
      )
    )
  
  if (
    nrow(analysis_df) == 0 ||
    n_distinct(analysis_df[[carrier_var]]) < 2 ||
    n_distinct(analysis_df$fh_factor) < 2 ||
    n_distinct(analysis_df$sex_factor) < 2
  ) {
    return(
      tibble(
        Interaction_OR = NA_real_,
        Interaction_LCL = NA_real_,
        Interaction_UCL = NA_real_,
        P_interaction = NA_real_
      )
    )
  }
  
  model <- suppressWarnings(
    glm(
      as.formula(
        paste(
          carrier_var,
          "~ fh_factor * sex_factor +",
          age_variable,
          "+ hpc"
        )
      ),
      data = analysis_df,
      family = binomial()
    )
  )
  
  coefficient_table <- coef(summary(model))
  
  interaction_term <- grep(
    "fh_factor.*:sex_factor|sex_factor.*:fh_factor",
    rownames(coefficient_table),
    value = TRUE
  )
  
  if (length(interaction_term) == 0) {
    return(
      tibble(
        Interaction_OR = NA_real_,
        Interaction_LCL = NA_real_,
        Interaction_UCL = NA_real_,
        P_interaction = NA_real_
      )
    )
  }
  
  estimate <- coefficient_table[
    interaction_term[1],
    "Estimate"
  ]
  
  standard_error <- coefficient_table[
    interaction_term[1],
    "Std. Error"
  ]
  
  p_value <- coefficient_table[
    interaction_term[1],
    "Pr(>|z|)"
  ]
  
  OR <- exp(estimate)
  LCL <- exp(estimate - 1.96 * standard_error)
  UCL <- exp(estimate + 1.96 * standard_error)
  
  unstable <- (
    is.na(OR) ||
      is.infinite(OR) ||
      is.infinite(LCL) ||
      is.infinite(UCL) ||
      OR > 1000 ||
      UCL > 1000 ||
      standard_error > 10
  )
  
  if (unstable) {
    return(
      tibble(
        Interaction_OR = NA_real_,
        Interaction_LCL = NA_real_,
        Interaction_UCL = NA_real_,
        P_interaction = NA_real_
      )
    )
  }
  
  tibble(
    Interaction_OR = OR,
    Interaction_LCL = LCL,
    Interaction_UCL = UCL,
    P_interaction = p_value
  )
}


# 8. Run overall analysis


run_overall_analysis <- function(
    carrier_var,
    pgv_label,
    fh_var,
    fh_label
) {
  
  prevalence <- make_prevalence_result(
    data,
    carrier_var,
    fh_var
  )
  
  model <- fit_adjusted_model(
    data,
    carrier_var,
    fh_var,
    include_sex = TRUE
  )
  
  interaction <- fit_interaction_model(
    data,
    carrier_var,
    fh_var
  )
  
  bind_cols(
    prevalence,
    model,
    interaction
  ) %>%
    mutate(
      PGV_group = pgv_label,
      Family_history = fh_label
    )
}


# 9. Run sex-stratified analysis


run_sex_analysis <- function(
    sex_name,
    carrier_var,
    pgv_label,
    fh_var,
    fh_label
) {
  
  sex_data <- data %>%
    filter(
      sex_label == sex_name
    )
  
  prevalence <- make_prevalence_result(
    sex_data,
    carrier_var,
    fh_var
  )
  
  model <- fit_adjusted_model(
    sex_data,
    carrier_var,
    fh_var,
    include_sex = FALSE
  )
  
  bind_cols(
    prevalence,
    model
  ) %>%
    mutate(
      Sex = sex_name,
      PGV_group = pgv_label,
      Family_history = fh_label
    )
}


# 10. Run all combinations


analysis_grid <- crossing(
  pgv_map,
  family_history_map
)

overall_results <- analysis_grid %>%
  pmap_dfr(
    function(
    carrier_var,
    PGV_group,
    fh_var,
    Family_history
    ) {
      
      run_overall_analysis(
        carrier_var,
        PGV_group,
        fh_var,
        Family_history
      )
    }
  )

sex_stratified_results <- crossing(
  Sex = c("Female", "Male"),
  analysis_grid
) %>%
  pmap_dfr(
    function(
    Sex,
    carrier_var,
    PGV_group,
    fh_var,
    Family_history
    ) {
      
      run_sex_analysis(
        Sex,
        carrier_var,
        PGV_group,
        fh_var,
        Family_history
      )
    }
  )


# 11. Multiple-testing correction


# IMPORTANT:
# S4 contains 3 PGV groups × 9 FH categories = 27 tests overall.
overall_results <- overall_results %>%
  mutate(
    FDR_q = p.adjust(P, method = "BH",
                     n = sum(!is.na(P))),
    P_bonf = p.adjust(P, method = "bonferroni",
                      n = sum(!is.na(P)))
  )

sex_stratified_results <- sex_stratified_results %>%
  mutate(
    FDR_q = p.adjust(P, method = "BH",
                     n = sum(!is.na(P))),
    P_bonf = p.adjust(P, method = "bonferroni",
                      n = sum(!is.na(P)))
  )


# 12. Formatting


format_prevalence <- function(
    n,
    N,
    prevalence
) {
  
  paste0(
    format(
      n,
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    ),
    "/",
    format(
      N,
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    ),
    " (",
    sprintf("%.2f", prevalence),
    ")"
  )
}

format_or_ci <- function(
    OR,
    LCL,
    UCL,
    fdr_q = NA_real_
) {
  
  if (
    is.na(OR) ||
    is.na(LCL) ||
    is.na(UCL)
  ) {
    return("NE")
  }
  
  marker <- ifelse(
    !is.na(fdr_q) &&
      fdr_q < 0.05,
    "*",
    ""
  )
  
  paste0(
    sprintf("%.2f", OR),
    marker,
    " (",
    sprintf("%.2f", LCL),
    "–",
    sprintf("%.2f", UCL),
    ")"
  )
}

format_interaction <- function(
    OR,
    LCL,
    UCL
) {
  
  if (
    is.na(OR) ||
    is.na(LCL) ||
    is.na(UCL)
  ) {
    return("NE")
  }
  
  paste0(
    sprintf("%.2f", OR),
    " (",
    sprintf("%.2f", LCL),
    "–",
    sprintf("%.2f", UCL),
    ")"
  )
}

format_p <- function(p) {
  
  if (is.na(p)) {
    return("NE")
  }
  
  if (p < 0.001) {
    return("<0.001")
  }
  
  sprintf("%.3f", p)
}


# 13. Overall display


overall_display <- overall_results %>%
  transmute(
    PGV_group,
    Family_history,
    
    Overall_FH_positive = pmap_chr(
      list(
        n_FH_positive,
        N_FH_positive,
        prevalence_FH_positive
      ),
      format_prevalence
    ),
    
    Overall_FH_negative = pmap_chr(
      list(
        n_FH_negative,
        N_FH_negative,
        prevalence_FH_negative
      ),
      format_prevalence
    ),
    
    Overall_OR = pmap_chr(
      list(
        OR,
        LCL,
        UCL,
        FDR_q
      ),
      format_or_ci
    ),
    
    Overall_P = map_chr(
      P,
      format_p
    ),
    
    Overall_FDR_q = map_chr(
      FDR_q,
      format_p
    ),
    
    Overall_Pbonf = map_chr(
      P_bonf,
      format_p
    ),
    
    Interaction_OR = pmap_chr(
      list(
        Interaction_OR,
        Interaction_LCL,
        Interaction_UCL
      ),
      format_interaction
    ),
    
    P_interaction = map_chr(
      P_interaction,
      format_p
    )
  )


# 14. Sex-specific display


sex_display <- sex_stratified_results %>%
  transmute(
    Sex,
    PGV_group,
    Family_history,
    
    FH_positive = pmap_chr(
      list(
        n_FH_positive,
        N_FH_positive,
        prevalence_FH_positive
      ),
      format_prevalence
    ),
    
    FH_negative = pmap_chr(
      list(
        n_FH_negative,
        N_FH_negative,
        prevalence_FH_negative
      ),
      format_prevalence
    ),
    
    Adjusted_OR = pmap_chr(
      list(
        OR,
        LCL,
        UCL,
        FDR_q
      ),
      format_or_ci
    ),
    
    P = map_chr(
      P,
      format_p
    ),
    
    FDR_q = map_chr(
      FDR_q,
      format_p
    ),
    
    Pbonf = map_chr(
      P_bonf,
      format_p
    )
  ) %>%
  pivot_wider(
    names_from = Sex,
    values_from = c(
      FH_positive,
      FH_negative,
      Adjusted_OR,
      P,
      FDR_q,
      Pbonf
    ),
    names_glue = "{Sex}_{.value}"
  )


# 15. Final Supplementary Table S4


supplementary_table_s4 <- overall_display %>%
  left_join(
    sex_display,
    by = c(
      "PGV_group",
      "Family_history"
    )
  ) %>%
  mutate(
    PGV_group = factor(
      PGV_group,
      levels = c(
        "BRCA1",
        "BRCA2",
        "HRR†"
      )
    ),
    
    Family_history = factor(
      Family_history,
      levels = family_history_map$Family_history
    )
  ) %>%
  arrange(
    PGV_group,
    Family_history
  ) %>%
  transmute(
    `PGV group` = as.character(PGV_group),
    `Family history` = as.character(Family_history),
    
    `Overall FH+ n/N (%)` = Overall_FH_positive,
    `Overall FH− n/N (%)` = Overall_FH_negative,
    `Overall adjusted OR (95% CI)` = Overall_OR,
    `Overall P` = Overall_P,
    `Overall FDR q` = Overall_FDR_q,
    `Overall Pbonf` = Overall_Pbonf,
    
    `Female FH+ n/N (%)` = Female_FH_positive,
    `Female FH− n/N (%)` = Female_FH_negative,
    `Female adjusted OR (95% CI)` = Female_Adjusted_OR,
    `Female P` = Female_P,
    `Female FDR q` = Female_FDR_q,
    `Female Pbonf` = Female_Pbonf,
    
    `Male FH+ n/N (%)` = Male_FH_positive,
    `Male FH− n/N (%)` = Male_FH_negative,
    `Male adjusted OR (95% CI)` = Male_Adjusted_OR,
    `Male P` = Male_P,
    `Male FDR q` = Male_FDR_q,
    `Male Pbonf` = Male_Pbonf,
    
    `Interaction OR (95% CI)` = Interaction_OR,
    `Pint` = P_interaction
  )


# 16. Notes


notes <- tibble(
  Item = c(
    "Population",
    "Overall model",
    "Sex-stratified models",
    "Interaction model",
    "Overall multiple-testing correction",
    "Sex-stratified multiple-testing correction",
    "Asterisk",
    "HRR†",
    "Abbreviations"
  ),
  
  Description = c(
    "Cancer-free controls in the final analysis dataset.",
    "Adjusted for age, sex, and HERPACC recruitment phase.",
    "Adjusted for age and HERPACC recruitment phase.",
    "Family history × sex interaction model adjusted for age and HERPACC recruitment phase; interaction OR represents the female-to-male ratio of ORs.",
"Benjamini–Hochberg FDR and Bonferroni corrections applied across 23 estimable comparisons among 27 planned overall tests; non-estimable comparisons were excluded.",
"Benjamini–Hochberg FDR and Bonferroni corrections applied jointly across 40 estimable comparisons among 54 planned sex-stratified tests; non-estimable comparisons were excluded.",
    "* indicates FDR q < 0.05.",
    "HRR† denotes homologous recombination repair genes excluding BRCA1 and BRCA2: ATM, BARD1, BRIP1, CHEK2, PALB2, RAD51C, and RAD51D.",
    "FH, family history; PGV, pathogenic germline variant; OR, odds ratio; CI, confidence interval; HRR, homologous recombination repair; NE, not estimable."
  )
)

# 17. Export


write_xlsx(
  list(
    Supplementary_Table_S4 = supplementary_table_s4,
    Overall_full_results = overall_results,
    Sex_stratified_full_results = sex_stratified_results,
    Notes = notes
  ),
  output_file
)

cat("\nSupplementary Table S4 completed.\n")
cat("Input: ", input_file, "\n", sep = "")
cat("Output: ", output_file, "\n", sep = "")
cat("Controls analysed: ", nrow(data), "\n", sep = "")

cat("\nAny cancer FH distribution:\n")
print(
  table(
    data$fh_any_fdr,
    useNA = "ifany"
  )
)

cat("\nOther cancer FH distribution:\n")
print(
  table(
    data$fh_other_fdr,
    useNA = "ifany"
  )
)

sum(is.na(overall_results$OR))
sum(is.na(sex_stratified_results$OR))

