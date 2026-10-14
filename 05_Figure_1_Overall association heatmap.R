# Overall heat map of adjusted associations between family history
# and PGV carrier status
#
# This script reads the overall results from Supplementary Table S3

library(readxl)
library(dplyr)
library(ggplot2)

input_file <- "Supplementary_Table_S3_sex_stratified.xlsx"
output_tiff <- "Heatmap_FH_PGV_overall.tiff"
output_png <- "Heatmap_FH_PGV_overall.png"

overall_heatmap_data <- read_excel(
  input_file,
  sheet = "Overall_full_results"
) %>%
  mutate(
    # Force model estimates to numeric after reading from Excel
    across(
      c(OR, LCL, UCL, P, FDR_q),
      ~ suppressWarnings(as.numeric(.x))
    ),
    
    Family_history = factor(
      Family_history,
      levels = rev(c(
        "Any cancer",
        "Breast cancer",
        "Ovarian cancer",
        "Pancreatic cancer",
        "Prostate cancer",
        "Colorectal cancer",
        "Gastric cancer",
        "Lung cancer",
        "Other cancers"
      ))
    ),
    
    PGV_group = factor(
      PGV_group,
      levels = c("All PGVs", "BRCA1/2", "HRR", "MMR")
    ),
    
    significance = ifelse(
      !is.na(FDR_q) & FDR_q < 0.05,
      "*",
      ""
    ),
    
    # A result is estimable only when the OR and both CI limits are finite
    estimable =
      !is.na(OR) &
      !is.na(LCL) &
      !is.na(UCL) &
      is.finite(OR) &
      is.finite(LCL) &
      is.finite(UCL) &
      OR > 0 &
      LCL > 0 &
      UCL > 0,
    
    log2_OR = ifelse(
      estimable,
      log2(OR),
      NA_real_
    ),
    
    
    label = ifelse(
      estimable,
      paste0(
        sprintf("%.2f", OR),
        significance,
        "\n(",
        sprintf("%.2f", LCL),
        "-",
        sprintf("%.2f", UCL),
        ")"
      ),
      "NE"
    ),
    
    label_face = ifelse(
      significance == "*",
      "bold",
      "plain"
    )
  )

overall_heatmap <- ggplot(
  overall_heatmap_data,
  aes(
    x = PGV_group,
    y = Family_history,
    fill = log2_OR
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label = label,
      fontface = label_face
    ),
    size = 3.2,
    lineheight = 0.9
  ) +
  scale_fill_gradient2(
    low = "#7BAFD4",
    mid = "white",
    high = "#C97A95",
    midpoint = 0,
    na.value = "grey70",
    name = expression(log[2](OR))
  ) +
  coord_fixed(ratio = 0.8) +
  scale_x_discrete(
    labels = c(
      "All PGVs",
      expression(italic("BRCA1/2")),
      "HRR",
      "MMR"
    )
  ) +
  labs(
    title = "Association between FH and PGV carrier status",
    x = "PGV group",
    y = "Family history"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 13
    ),
    axis.text.x = element_text(
      face = "bold",
      size = 9
    ),
    axis.text.y = element_text(
      size = 9,
      color = "grey25"
    ),
    axis.title.x = element_text(
      face = "bold",
      margin = margin(t = 12),
      size = 14
    ),
    axis.title.y = element_text(
      face = "bold",
      margin = margin(r = 18),
      size = 12
    ),
    panel.border = element_rect(
      colour = "grey90",
      fill = NA,
      linewidth = 0.4
    ),
    panel.grid = element_blank(),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )

print(overall_heatmap)

ggsave(
  filename = output_tiff,
  plot = overall_heatmap,
  width = 7.5,
  height = 8.5,
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = output_png,
  plot = overall_heatmap,
  width = 7.5,
  height = 8.5,
  dpi = 300
)
