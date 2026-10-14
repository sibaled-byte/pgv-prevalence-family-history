# Number of qualifying PGVs per participant

library(haven)
library(dplyr)
library(tidyr)
library(writexl)

pv_file <- "familyhistory_pv_hpc.dta"
analysis_file <- "ACC_PV_FH_controls_20260803.dta"

output_file <- "PGV_Carrier_Distribution.xlsx"

pv_data <- read_dta(pv_file)
analysis_data <- read_dta(analysis_file)

# Identify gene-level loss-of-function and ClinVar indicators

lof_vars <- names(pv_data)[
  grepl("_lof$", names(pv_data))
]

clin_vars <- names(pv_data)[
  grepl("_clin$", names(pv_data))
]

genes <- union(
  sub("_lof$", "", lof_vars),
  sub("_clin$", "", clin_vars)
)

cat("\nGenes identified:\n")
print(genes)

cat("\nNumber of genes identified:\n")
print(length(genes))

# Coding:
# 0 = no qualifying PGV
# 1 = qualifying PGV present
# 9 = missing or QC failure

recode_pgv_indicator <- function(x) {
  x_num <- suppressWarnings(
    as.numeric(as.character(x))
  )
  
  case_when(
    x_num == 1 ~ 1L,
    x_num == 0 ~ 0L,
    x_num == 9 ~ NA_integer_,
    TRUE ~ NA_integer_
  )
}

pgv_by_gene <- pv_data %>%
  select(
    acc,
    all_of(lof_vars),
    all_of(clin_vars)
  )

# Create one PGV indicator per gene.
# A gene is positive when either its loss-of-function indicator
# or ClinVar indicator equals 1.

for (gene in genes) {
  
  lof_var <- paste0(gene, "_lof")
  clin_var <- paste0(gene, "_clin")
  
  lof_result <- if (lof_var %in% names(pgv_by_gene)) {
    recode_pgv_indicator(
      pgv_by_gene[[lof_var]]
    )
  } else {
    rep(0L, nrow(pgv_by_gene))
  }
  
  clin_result <- if (clin_var %in% names(pgv_by_gene)) {
    recode_pgv_indicator(
      pgv_by_gene[[clin_var]]
    )
  } else {
    rep(0L, nrow(pgv_by_gene))
  }
  
  pgv_by_gene[[paste0("pgv_", gene)]] <- case_when(
    lof_result == 1L |
      clin_result == 1L ~ 1L,
    
    lof_result == 0L &
      clin_result == 0L ~ 0L,
    
    TRUE ~ NA_integer_
  )
}

pgv_indicator_vars <- names(pgv_by_gene)[
  grepl("^pgv_", names(pgv_by_gene))
]

# Count the number of genes with at least one qualifying PGV.
# Participants with unresolved results for any gene are flagged separately.

pgv_count <- pgv_by_gene %>%
  mutate(
    n_missing_gene_results = rowSums(
      is.na(across(all_of(pgv_indicator_vars)))
    ),
    
    n_genes_with_pgv = rowSums(
      across(
        all_of(pgv_indicator_vars),
        ~ replace_na(.x, 0L)
      )
    )
  ) %>%
  select(
    acc,
    n_genes_with_pgv,
    n_missing_gene_results,
    all_of(pgv_indicator_vars)
  )

# Restrict to the final cancer-free control sample.

controls <- analysis_data %>%
  filter(
    as.character(CaCo_status) %in%
      c("0", "control", "Control")
  )

if (nrow(controls) != 11811) {
  warning(
    "The control count is ",
    nrow(controls),
    "; the expected manuscript count is 11,811."
  )
}

# Merge the PGV counts onto the final cancer-free control sample.

control_pgv_count <- controls %>%
  left_join(
    pgv_count,
    by = "acc"
  )

if (any(is.na(control_pgv_count$n_genes_with_pgv))) {
  warning(
    sum(is.na(control_pgv_count$n_genes_with_pgv)),
    " controls could not be matched to the PGV source dataset."
  )
}

# Distribution of the number of PGV-positive genes among all controls.

distribution_all <- control_pgv_count %>%
  count(
    n_genes_with_pgv,
    name = "N"
  ) %>%
  mutate(
    Percent = round(
      100 * N / sum(N),
      3
    )
  ) %>%
  arrange(n_genes_with_pgv)

# Distribution among participants carrying at least one PGV.

distribution_carriers <- control_pgv_count %>%
  filter(n_genes_with_pgv >= 1) %>%
  count(
    n_genes_with_pgv,
    name = "N"
  ) %>%
  mutate(
    Percent_among_carriers = round(
      100 * N / sum(N),
      2
    )
  ) %>%
  arrange(n_genes_with_pgv)

# Main summary of PGV counts per participant.

summary_results <- control_pgv_count %>%
  summarise(
    N_controls = n(),
    
    N_matched_to_PGV_data = sum(
      !is.na(n_genes_with_pgv)
    ),
    
    N_unmatched_to_PGV_data = sum(
      is.na(n_genes_with_pgv)
    ),
    
    N_with_any_missing_gene_result = sum(
      n_missing_gene_results > 0,
      na.rm = TRUE
    ),
    
    N_with_0_PGVs = sum(
      n_genes_with_pgv == 0,
      na.rm = TRUE
    ),
    
    N_with_1_PGV = sum(
      n_genes_with_pgv == 1,
      na.rm = TRUE
    ),
    
    N_with_2_PGVs = sum(
      n_genes_with_pgv == 2,
      na.rm = TRUE
    ),
    
    N_with_3_or_more_PGVs = sum(
      n_genes_with_pgv >= 3,
      na.rm = TRUE
    ),
    
    N_with_2_or_more_PGVs = sum(
      n_genes_with_pgv >= 2,
      na.rm = TRUE
    ),
    
    Percent_with_2_or_more_among_all = round(
      100 *
        sum(n_genes_with_pgv >= 2, na.rm = TRUE) /
        sum(!is.na(n_genes_with_pgv)),
      3
    ),
    
    Percent_with_2_or_more_among_carriers = round(
      100 *
        sum(n_genes_with_pgv >= 2, na.rm = TRUE) /
        sum(n_genes_with_pgv >= 1, na.rm = TRUE),
      2
    )
  )

# List participants carrying PGVs in two or more genes.

multiple_pgv_participants <- control_pgv_count %>%
  filter(n_genes_with_pgv >= 2) %>%
  pivot_longer(
    cols = all_of(pgv_indicator_vars),
    names_to = "Gene_variable",
    values_to = "PGV_present"
  ) %>%
  filter(PGV_present == 1L) %>%
  mutate(
    Gene = sub("^pgv_", "", Gene_variable)
  ) %>%
  select(
    acc,
    n_genes_with_pgv,
    Gene
  ) %>%
  arrange(
    acc,
    Gene
  )

# Number of carriers for each gene.

gene_carrier_counts <- control_pgv_count %>%
  summarise(
    across(
      all_of(pgv_indicator_vars),
      ~ sum(.x == 1L, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Gene_variable",
    values_to = "Carrier_N"
  ) %>%
  mutate(
    Gene = sub("^pgv_", "", Gene_variable)
  ) %>%
  select(
    Gene,
    Carrier_N
  ) %>%
  arrange(
    desc(Carrier_N),
    Gene
  )

# Check the agreement between the gene-level count and the original
# participant-level carrier indicator.

carrier_agreement <- control_pgv_count %>%
  left_join(
    analysis_data %>%
      transmute(
        acc,
        original_carrier = case_when(
          gene_clinlof == "No_carrier" ~ 0L,
          gene_clinlof == "QC out" ~ NA_integer_,
          is.na(gene_clinlof) ~ NA_integer_,
          TRUE ~ 1L
        )
      ),
    by = "acc"
  ) %>%
  mutate(
    derived_carrier = as.integer(
      n_genes_with_pgv >= 1
    )
  ) %>%
  count(
    original_carrier,
    derived_carrier,
    name = "N"
  )

quality_control <- tibble(
  Check = c(
    "Final cancer-free controls",
    "Controls matched to PGV source data",
    "Controls unmatched to PGV source data",
    "Genes evaluated",
    "Controls with any gene-level code of 9",
    "Participants with at least one PGV-positive gene",
    "Participants with two or more PGV-positive genes"
  ),
  
  N = c(
    nrow(control_pgv_count),
    
    sum(
      !is.na(control_pgv_count$n_genes_with_pgv)
    ),
    
    sum(
      is.na(control_pgv_count$n_genes_with_pgv)
    ),
    
    length(genes),
    
    sum(
      control_pgv_count$n_missing_gene_results > 0,
      na.rm = TRUE
    ),
    
    sum(
      control_pgv_count$n_genes_with_pgv >= 1,
      na.rm = TRUE
    ),
    
    sum(
      control_pgv_count$n_genes_with_pgv >= 2,
      na.rm = TRUE
    )
  )
)

write_xlsx(
  list(
    Summary = summary_results,
    Distribution_all_controls = distribution_all,
    Distribution_carriers = distribution_carriers,
    Multiple_PGV_participants = multiple_pgv_participants,
    Gene_carrier_counts = gene_carrier_counts,
    Carrier_agreement = carrier_agreement,
    Quality_control = quality_control
  ),
  output_file
)

cat("\nSummary of PGV-positive genes per participant:\n")
print(summary_results, width = Inf)

cat("\nDistribution among all controls:\n")
print(distribution_all, n = Inf)

cat("\nDistribution among carriers:\n")
print(distribution_carriers, n = Inf)

cat("\nParticipants carrying PGVs in two or more genes:\n")
print(multiple_pgv_participants, n = Inf)

cat("\nCarrier-status agreement check:\n")
print(carrier_agreement, n = Inf)

cat("\nQuality-control checks:\n")
print(quality_control, n = Inf)

cat(
  "\nResults saved to:\n",
  output_file,
  "\n"
)
