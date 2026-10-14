library(haven)
library(dplyr)
library(purrr)
library(tibble)
library(writexl)

input_file <- "ACC_PV_FH_controls_20260803.dta"
output_file <- "Supplementary_Table_S5_incremental_pseudo_R2.xlsx"

minimum_gene_carriers <- 10L

data <- read_dta(input_file)

if (nrow(data) != 11811) {
  stop(
    "Analysis dataset contains ",
    nrow(data),
    " participants; expected 11,811 cancer-free controls."
  )
}

age_variable <- if ("age" %in% names(data)) {
  "age"
} else if ("fvage" %in% names(data)) {
  "fvage"
} else {
  stop("Neither age nor fvage was found in the analysis dataset.")
}

required_variables <- c("sex", "hpc", "gene_clinlof")
missing_required <- setdiff(required_variables, names(data))

if (length(missing_required) > 0) {
  stop(
    "Missing required variable(s): ",
    paste(missing_required, collapse = ", ")
  )
}

family_history_map <- tribble(
  ~fh_var, ~Family_history,
  "fh_any_fdr", "Any cancer",
  "fh_breast_fdr", "Breast cancer",
  "fh_ovarian_fdr", "Ovarian cancer",
  "fh_pancreas_fdr", "Pancreatic cancer",
  "fh_prostate_fdr", "Prostate cancer",
  "fh_colon_fdr", "Colorectal cancer",
  "fh_gastric_fdr", "Gastric cancer",
  "fh_lung_fdr", "Lung cancer",
  "fh_other_fdr", "Other cancers"
)

missing_fh <- setdiff(family_history_map$fh_var, names(data))

if (length(missing_fh) > 0) {
  stop(
    "Missing family-history variable(s): ",
    paste(missing_fh, collapse = ", ")
  )
}

data <- data %>%
  mutate(
    age_analysis = as.numeric(.data[[age_variable]]),
    pgv_all = case_when(
      gene_clinlof == "No_carrier" ~ 0L,
      gene_clinlof == "QC out" ~ NA_integer_,
      is.na(gene_clinlof) ~ NA_integer_,
      TRUE ~ 1L
    )
  )

brca_genes <- c("BRCA1", "BRCA2")

hrr_genes <- c(
  "ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1",
  "CHEK2", "PALB2", "RAD51C", "RAD51D"
)

mmr_genes <- c("MLH1", "MSH2", "MSH6", "PMS2")

data <- data %>%
  mutate(
    pgv_brca = case_when(
      is.na(pgv_all) ~ NA_integer_,
      pgv_all == 1L & gene_clinlof %in% brca_genes ~ 1L,
      TRUE ~ 0L
    ),
    pgv_hrr = case_when(
      is.na(pgv_all) ~ NA_integer_,
      pgv_all == 1L & gene_clinlof %in% hrr_genes ~ 1L,
      TRUE ~ 0L
    ),
    pgv_mmr = case_when(
      is.na(pgv_all) ~ NA_integer_,
      pgv_all == 1L & gene_clinlof %in% mmr_genes ~ 1L,
      TRUE ~ 0L
    )
  )

outcome_map <- tribble(
  ~outcome_var, ~Outcome_level, ~Outcome,
  "pgv_all", "Overall", "All PGVs",
  "pgv_brca", "Pathway/gene group", "BRCA1/2",
  "pgv_hrr", "Pathway/gene group", "HRR",
  "pgv_mmr", "Pathway/gene group", "MMR"
)

recode_gene_indicator <- function(x) {
  x_numeric <- suppressWarnings(as.numeric(x))

  case_when(
    x_numeric == 1 ~ 1L,
    x_numeric == 0 ~ 0L,
    x_numeric == 9 ~ NA_integer_,
    TRUE ~ NA_integer_
  )
}

lof_variables <- names(data)[grepl("_lof$", names(data))]
clin_variables <- names(data)[grepl("_clin$", names(data))]
genes <- sort(
  union(
    sub("_lof$", "", lof_variables),
    sub("_clin$", "", clin_variables)
  )
)

for (gene in genes) {
  lof_var <- paste0(gene, "_lof")
  clin_var <- paste0(gene, "_clin")

  lof_result <- if (lof_var %in% names(data)) {
    recode_gene_indicator(data[[lof_var]])
  } else {
    rep(0L, nrow(data))
  }

  clin_result <- if (clin_var %in% names(data)) {
    recode_gene_indicator(data[[clin_var]])
  } else {
    rep(0L, nrow(data))
  }

  data[[paste0("pgv_gene_", gene)]] <- case_when(
    lof_result == 1L | clin_result == 1L ~ 1L,
    lof_result == 0L & clin_result == 0L ~ 0L,
    TRUE ~ NA_integer_
  )
}

gene_outcome_map <- tibble(
  outcome_var = paste0("pgv_gene_", genes),
  Outcome_level = "Gene",
  Outcome = genes
)

outcome_map_all <- bind_rows(outcome_map, gene_outcome_map)

nagelkerke_r2 <- function(model, null_model) {
  n <- nobs(model)
  ll_model <- as.numeric(logLik(model))
  ll_null <- as.numeric(logLik(null_model))
  numerator <- 1 - exp((2 / n) * (ll_null - ll_model))
  denominator <- 1 - exp((2 / n) * ll_null)
  numerator / denominator
}

mcfadden_r2 <- function(model, null_model) {
  1 - as.numeric(logLik(model)) / as.numeric(logLik(null_model))
}

fit_incremental_models <- function(
    df,
    outcome_var,
    added_fh_vars,
    analysis_label,
    family_history_label
) {
  model_variables <- c(
    outcome_var,
    "age_analysis",
    "sex",
    "hpc",
    added_fh_vars
  )

  analysis_df <- df %>%
    select(all_of(model_variables)) %>%
    filter(if_all(everything(), ~ !is.na(.x)))

  n_total <- nrow(analysis_df)
  n_carriers <- sum(analysis_df[[outcome_var]] == 1L)
  n_noncarriers <- sum(analysis_df[[outcome_var]] == 0L)

  if (
    n_total == 0L ||
    n_carriers == 0L ||
    n_noncarriers == 0L ||
    any(map_lgl(added_fh_vars, ~ n_distinct(analysis_df[[.x]]) < 2L))
  ) {
    return(
      tibble(
        Analysis = analysis_label,
        Family_history = family_history_label,
        Model_N = n_total,
        Carrier_count = n_carriers,
        Base_Nagelkerke_R2 = NA_real_,
        Expanded_Nagelkerke_R2 = NA_real_,
        Incremental_Nagelkerke_R2 = NA_real_,
        Base_McFadden_R2 = NA_real_,
        Expanded_McFadden_R2 = NA_real_,
        Incremental_McFadden_R2 = NA_real_,
        Likelihood_ratio_P = NA_real_,
        Converged = FALSE,
        Note = "Not estimable because the outcome or an FH variable did not vary."
      )
    )
  }

  base_formula <- as.formula(
    paste(outcome_var, "~ age_analysis + sex + hpc")
  )

  expanded_formula <- reformulate(
    c("age_analysis", "sex", "hpc", added_fh_vars),
    response = outcome_var
  )

  null_model <- suppressWarnings(
    glm(
      reformulate("1", response = outcome_var),
      data = analysis_df,
      family = binomial()
    )
  )

  base_model <- suppressWarnings(
    glm(base_formula, data = analysis_df, family = binomial())
  )

  expanded_model <- suppressWarnings(
    glm(expanded_formula, data = analysis_df, family = binomial())
  )

  lrt <- anova(base_model, expanded_model, test = "LRT")
  lrt_p <- lrt$`Pr(>Chi)`[nrow(lrt)]
  converged <- isTRUE(base_model$converged) && isTRUE(expanded_model$converged)

  base_nagelkerke <- nagelkerke_r2(base_model, null_model)
  expanded_nagelkerke <- nagelkerke_r2(expanded_model, null_model)
  base_mcfadden <- mcfadden_r2(base_model, null_model)
  expanded_mcfadden <- mcfadden_r2(expanded_model, null_model)

  note <- if (!converged) {
    "At least one model did not converge."
  } else if (n_carriers < minimum_gene_carriers) {
    paste0(
      "Sparse outcome (fewer than ",
      minimum_gene_carriers,
      " carriers); estimate should not be interpreted."
    )
  } else {
    ""
  }

  tibble(
    Analysis = analysis_label,
    Family_history = family_history_label,
    Model_N = n_total,
    Carrier_count = n_carriers,
    Base_Nagelkerke_R2 = base_nagelkerke,
    Expanded_Nagelkerke_R2 = expanded_nagelkerke,
    Incremental_Nagelkerke_R2 = expanded_nagelkerke - base_nagelkerke,
    Base_McFadden_R2 = base_mcfadden,
    Expanded_McFadden_R2 = expanded_mcfadden,
    Incremental_McFadden_R2 = expanded_mcfadden - base_mcfadden,
    Likelihood_ratio_P = lrt_p,
    Converged = converged,
    Note = note
  )
}

individual_results <- pmap_dfr(
  outcome_map_all,
  function(outcome_var, Outcome_level, Outcome) {
    map2_dfr(
      family_history_map$fh_var,
      family_history_map$Family_history,
      ~ fit_incremental_models(
        df = data,
        outcome_var = outcome_var,
        added_fh_vars = .x,
        analysis_label = "Individual FH category",
        family_history_label = .y
      )
    ) %>%
      mutate(
        Outcome_level = Outcome_level,
        Outcome = Outcome,
        .before = 1
      )
  }
)

site_specific_fh_vars <- family_history_map %>%
  filter(fh_var != "fh_any_fdr") %>%
  pull(fh_var)

omnibus_results <- pmap_dfr(
  outcome_map_all,
  function(outcome_var, Outcome_level, Outcome) {
    fit_incremental_models(
      df = data,
      outcome_var = outcome_var,
      added_fh_vars = site_specific_fh_vars,
      analysis_label = "Omnibus FH block",
      family_history_label = "Eight site-specific FH categories"
    ) %>%
      mutate(
        Outcome_level = Outcome_level,
        Outcome = Outcome,
        .before = 1
      )
  }
)

all_results <- bind_rows(individual_results, omnibus_results)

primary_results <- all_results %>%
  filter(Outcome_level != "Gene")

gene_results <- all_results %>%
  filter(Outcome_level == "Gene")

gene_reportable <- gene_results %>%
  filter(Carrier_count >= minimum_gene_carriers, Converged)

gene_sparse <- gene_results %>%
  filter(Carrier_count < minimum_gene_carriers | !Converged)

definitions <- tribble(
  ~Item, ~Definition,
  "Base model", "PGV carrier status modeled using age, sex, and HERPACC recruitment phase.",
  "Expanded individual-FH model", "Base model plus one FH category. FH absence is the reference.",
  "Expanded omnibus model", "Base model plus all eight site-specific FH categories. Any-cancer FH is excluded because it overlaps with the site-specific variables.",
  "Incremental pseudo-R2", "Expanded-model pseudo-R2 minus base-model pseudo-R2, estimated in the same complete-case sample.",
  "Nagelkerke pseudo-R2", "Likelihood-based measure scaled from 0 to 1. It is not a literal proportion of variance explained.",
  "McFadden pseudo-R2", "Alternative likelihood-based pseudo-R2 reported as a robustness measure.",
  "Gene-level reporting rule", paste0("Gene-level estimates with fewer than ", minimum_gene_carriers, " carriers are flagged as sparse and excluded from the reportable sheet."),
  "Likelihood-ratio P value", "Compares the nested base and expanded models."
)

write_xlsx(
  list(
    Primary_overall_pathways = primary_results,
    Gene_level_reportable = gene_reportable,
    Gene_level_sparse = gene_sparse,
    Definitions = definitions
  ),
  output_file
)

cat("Rows in analysis dataset:", nrow(data), "\n")
cat("Genes identified:", length(genes), "\n")
cat("Results written to:", output_file, "\n")
