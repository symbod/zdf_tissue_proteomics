# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Explicitly ensure dplyr functions are prioritized
select <- dplyr::select
filter <- dplyr::filter

# ============================================================================
# CONFIGURATION: Specify comparisons of interest and their labels
# ============================================================================
comparisons_of_interest <- c(
  "diabetic_empty_42-nondiabetic_empty_42" = "Diabetic Untreated vs Non-diabetic Untreated",
  "diabetic_PCL_42-nondiabetic_PCL_42" = "Diabetic PCL Scaffold vs Non-diabetic PCL Scaffold"
)

############################################
# Read data files
# ---
output_dir <- "valid_DEPs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read protein abundance data
abundance_file_AI_unscaled <- "normalization_results_AI/RobNorm_normalized_data.csv"
abundance_file_AS_unscaled <- "normalization_results_AS/RobNorm_normalized_data.csv"
metadata_file_AI <- "normalization_results_AI/meta_data.csv"
metadata_file_AS <- "normalization_results_AS/meta_data.csv"
ai_data_unscaled <- read.csv(abundance_file_AI_unscaled, stringsAsFactors = FALSE, check.names = FALSE)
as_data_unscaled <- read.csv(abundance_file_AS_unscaled, stringsAsFactors = FALSE, check.names = FALSE)
meta_ai <- read.csv(metadata_file_AI, stringsAsFactors = FALSE)
meta_as <- read.csv(metadata_file_AS, stringsAsFactors = FALSE)

# Read meta-analysis proteins
meta_proteins <- read.csv("input/bone_caps_meta_analysis.csv", stringsAsFactors = FALSE)
meta_genes <- unique(toupper(meta_proteins$Gene))

# Read DE results
de_results_AI <- read.csv("de_analysis_results_AI_RobNorm_scaled/de_results_raw_with_orthologs.csv", stringsAsFactors = FALSE)
de_results_AS <- read.csv("de_analysis_results_AS_RobNorm_scaled/de_results_raw_with_orthologs.csv", stringsAsFactors = FALSE)
# Filter DE results to only include specified comparisons
de_results_AI <- de_results_AI[de_results_AI$Comparison %in% names(comparisons_of_interest), ]
de_results_AS <- de_results_AS[de_results_AS$Comparison %in% names(comparisons_of_interest), ]
############################################

############################################
# Create Sample_ID to match AI and AS (same animal, condition, scaffold, intervention, timepoint, tissue) 
# ---
meta_ai$Sample_ID <- paste(meta_ai$Animal, meta_ai$Condition, meta_ai$Scaffold, 
                           meta_ai$Intervention, meta_ai$Timepoint, meta_ai$Tissue, sep="_")
meta_as$Sample_ID <- paste(meta_as$Animal, meta_as$Condition, meta_as$Scaffold, 
                           meta_as$Intervention, meta_as$Timepoint, meta_as$Tissue, sep="_")

# Match AI and AS samples
matched_samples <- inner_join(
  meta_ai %>% select(Sample_ID, Column_AI = Column, 
                     Condition, Scaffold, Intervention, Timepoint, Tissue, Animal),
  meta_as %>% select(Sample_ID, Column_AS = Column),
  by = "Sample_ID"
)

cat("Matched", nrow(matched_samples), "AI/AS sample pairs\n\n")
############################################

############################################
# HELPER FUNCTIONS
############################################

# Per-protein indicator: does ANY of a protein's candidate orthologs (a
# semicolon-separated string, since rat-to-human orthology can be ambiguous)
# appear in a reference gene list? Each protein contributes at most one hit,
# regardless of how many candidate orthologs it lists — the orthologs all
# refer to the same single measured protein, so they must not be counted
# individually (that would let ambiguous, multi-paralog proteins inflate
# overlap counts relative to unambiguous ones).
any_ortholog_in <- function(ortholog_strings, ref_genes) {
  sapply(strsplit(as.character(ortholog_strings), ";"), function(orths) {
    any(toupper(trimws(orths)) %in% ref_genes)
  })
}

# Parse a group name like "diabetic_PCL_42" from a comparison
parse_group <- function(group_str) {
  parts <- strsplit(group_str, "_")[[1]]
  list(
    condition = parts[1],
    scaffold  = parts[length(parts) - 1],
    timepoint = as.integer(parts[length(parts)])
  )
}

# Per-protein summary for an arbitrary set of sample IDs, derived from
# complete_fcs_proteinwise (the single source of truth for per-sample AI/AS data).
#   - n_samples:         how many sample rows contribute
#   - ai_exclusive:      never detected in AS across this sample set
#   - as_exclusive:      never detected in AI across this sample set
#   - dominant_fraction: "AI"/"AS"/"Tie" if all samples agree, else "Inconsistent"
#   - geo_mean_fc:       geometric mean of per-sample AI/AS ratios
#                        (one-fraction-only samples imputed to p99)
summarise_proteins_across_samples <- function(sample_ids, fc_table, p99_value) {
  fc_table %>%
    filter(Sample_ID %in% sample_ids) %>%
    mutate(fc_imputed = ifelse(is.na(Fold_Change), p99_value, Fold_Change)) %>%
    group_by(Protein.IDs) %>%
    summarise(
      n_samples         = n(),
      ai_exclusive      = all(is.na(AS_linear)),
      as_exclusive      = all(is.na(AI_linear)),
      dominant_fraction = if (n_distinct(Dominant_Fraction) == 1)
        first(Dominant_Fraction) else "Inconsistent",
      geo_mean_fc       = exp(mean(log(fc_imputed))),
      .groups           = "drop"
    )
}

############################################
# Build per-sample AI/AS ratios for every protein in every sample pair.
# Retain one-fraction-only detections (NA Fold_Change) so they can be imputed later.
############################################
all_fc <- list()
cat("Processing", nrow(matched_samples), "sample pairs...\n")

for (i in 1:nrow(matched_samples)) {
  ai_col <- matched_samples$Column_AI[i]
  as_col <- matched_samples$Column_AS[i]
  
  ai_sample <- data.frame(
    Protein.IDs = ai_data_unscaled$Protein.IDs,
    AI_abundance = as.numeric(ai_data_unscaled[[ai_col]]),
    stringsAsFactors = FALSE
  )
  as_sample <- data.frame(
    Protein.IDs = as_data_unscaled$Protein.IDs,
    AS_abundance = as.numeric(as_data_unscaled[[as_col]]),
    stringsAsFactors = FALSE
  )
  
  merged <- full_join(ai_sample, as_sample, by = "Protein.IDs") %>%
    filter(!is.na(AI_abundance) | !is.na(AS_abundance))
  
  merged$AI_linear <- 2^merged$AI_abundance
  merged$AS_linear <- 2^merged$AS_abundance
  
  merged$Fold_Change <- ifelse(
    !is.na(merged$AI_linear) & !is.na(merged$AS_linear),
    pmax(merged$AI_linear, merged$AS_linear) / pmin(merged$AI_linear, merged$AS_linear),
    NA_real_
  )
  
  merged$Dominant_Fraction <- case_when(
    is.na(merged$AS_linear) & !is.na(merged$AI_linear) ~ "AI",
    is.na(merged$AI_linear) & !is.na(merged$AS_linear) ~ "AS",
    merged$AI_linear >  merged$AS_linear               ~ "AI",
    merged$AS_linear >  merged$AI_linear               ~ "AS",
    merged$AI_linear == merged$AS_linear               ~ "Tie"
  )
  
  all_fc[[i]] <- data.frame(
    Protein.IDs       = merged$Protein.IDs,
    Sample_ID         = matched_samples$Sample_ID[i],
    Fold_Change       = merged$Fold_Change,
    AI_linear         = merged$AI_linear,
    AS_linear         = merged$AS_linear,
    Dominant_Fraction = merged$Dominant_Fraction,
    stringsAsFactors  = FALSE
  )
}

complete_fcs_proteinwise <- bind_rows(all_fc)

cat("Total sample-protein rows:", nrow(complete_fcs_proteinwise), "\n")
cat("  ...with both fractions measured:",
    sum(!is.na(complete_fcs_proteinwise$Fold_Change)), "\n\n")

############################################
# Percentiles of the background fold-change distribution
############################################
background_fc <- complete_fcs_proteinwise$Fold_Change
background_fc <- background_fc[!is.na(background_fc) & is.finite(background_fc)]

median_fc <- median(background_fc)
p70  <- quantile(background_fc, 0.70)
p80  <- quantile(background_fc, 0.80)
p90  <- quantile(background_fc, 0.90)
p99  <- quantile(background_fc, 0.99)    # imputation value for one-fraction-only samples
p999 <- quantile(background_fc, 0.999)   # plot x-axis display limit

cat("Percentiles (background distribution):\n")
cat("  Median:", round(median_fc, 2), "\n")
cat("  70th:",   round(p70, 2),  "\n")
cat("  80th:",   round(p80, 2),  "\n")
cat("  90th:",   round(p90, 2),  "\n")
cat("  99th:",   round(p99, 2),  "  (imputation value)\n")
cat("  99.9th:", round(p999, 2), "  (plot display limit)\n")
cat("  Max:",    round(max(background_fc), 2), "\n\n")

percentile_data <- data.frame(
  label = c("Median", "70th", "80th", "90th"),
  value = c(median_fc, p70, p80, p90),
  text  = sprintf("%.2f", c(median_fc, p70, p80, p90))
)

plot <- ggplot(data.frame(Fold_Change = background_fc), aes(x = Fold_Change)) +
  geom_histogram(bins = 80, fill = "#2E86AB", alpha = 0.85, color = "white") +
  geom_vline(data = percentile_data, aes(xintercept = value, color = label),
             linewidth = 1.5, linetype = "dashed") +
  scale_color_manual(
    values = c("Median" = "#FFB703", "70th" = "#06A77D", "80th" = "#D62828", "90th" = "#7209B7"),
    breaks = c("Median", "70th", "80th", "90th"),
    labels = setNames(paste0(percentile_data$label, ": ", percentile_data$text),
                      percentile_data$label)
  ) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 500)) +
  coord_cartesian(xlim = c(1, p999)) +
  labs(
    title    = "AI vs AS Fold Change Distribution",
    subtitle = paste0("Every protein in every sample (n = ", length(background_fc),
                      "; top 0.1% extreme ratios not shown)"),
    x        = "Abundance Ratio (Larger / Smaller)",
    y        = "Count",
    color    = "Percentiles"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title    = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 13, color = "grey40"),
    axis.title    = element_text(size = 16),
    axis.text     = element_text(size = 14),
    legend.title  = element_text(size = 15),
    legend.text   = element_text(size = 14),
    legend.position = "right",
    panel.grid    = element_blank()
  )

ggsave(file.path(output_dir, "FC_distribution.png"), plot, width = 12, height = 7, dpi = 300)

extreme_fc <- complete_fcs_proteinwise %>%
  filter(!is.na(Fold_Change) & Fold_Change > 100) %>%
  arrange(desc(Fold_Change))

cat("\nTop 10 extreme fold changes:\n")
for (i in 1:min(10, nrow(extreme_fc))) {
  cat("\nProtein:", extreme_fc$Protein.IDs[i], "\n")
  cat("  Sample:", extreme_fc$Sample_ID[i], "\n")
  cat("  AI intensity:", extreme_fc$AI_linear[i], "\n")
  cat("  AS intensity:", extreme_fc$AS_linear[i], "\n")
  cat("  Fold Change:", extreme_fc$Fold_Change[i], "\n")
  cat("  Dominant:", extreme_fc$Dominant_Fraction[i], "\n")
}

###################################################################
# ============================================================================
# MAIN THRESHOLD TESTING LOOPS
# ============================================================================
###################################################################

# FC thresholds derived from the background distribution computed above — no hardcoding.
fc_thresholds        <- round(c(p70, p80, p90), 2)
# Stability threshold = max allowed ratio between the two groups' geo-mean FCs.
# e.g. 1.2 means "the larger group's AI/AS separation is at most 1.2-fold the smaller's".
stability_thresholds <- c(1.1, 1.2, 1.3)


for (fc_thresh in fc_thresholds) {
  for (stab_thresh in stability_thresholds) {
    
    cat("\n========================================\n")
    cat("Testing FC threshold:", fc_thresh, "| Stability:", stab_thresh, "-fold\n")
    cat("========================================\n")
    
    thresh_dir <- file.path(output_dir, paste0("FC", fc_thresh, "_Stab", stab_thresh))
    dir.create(thresh_dir, showWarnings = FALSE, recursive = TRUE)
    
    all_validated_deps           <- list()
    dominant_removed             <- list()
    tested_second_level_per_comp <- list()
    
    for (comp in names(comparisons_of_interest)) {
      cat("\nProcessing:", comparisons_of_interest[comp], "\n")
      
      comp_parts <- strsplit(comp, "-")[[1]]
      group1 <- parse_group(comp_parts[1])
      group2 <- parse_group(comp_parts[2])
      
      group1_samples <- matched_samples %>%
        filter(Condition == group1$condition, Scaffold == group1$scaffold, Timepoint == group1$timepoint)
      group2_samples <- matched_samples %>%
        filter(Condition == group2$condition, Scaffold == group2$scaffold, Timepoint == group2$timepoint)
      comparison_samples <- rbind(group1_samples, group2_samples)
      
      # Per-protein summaries: pooled (for dominance/exclusivity) + per-group (for stability)
      comp_summary   <- summarise_proteins_across_samples(comparison_samples$Sample_ID, complete_fcs_proteinwise, p99)
      group1_summary <- summarise_proteins_across_samples(group1_samples$Sample_ID,     complete_fcs_proteinwise, p99)
      group2_summary <- summarise_proteins_across_samples(group2_samples$Sample_ID,     complete_fcs_proteinwise, p99)
      
      # Stability: ratio of per-group geo-mean FCs (always >= 1)
      stability_data <- inner_join(
        group1_summary %>% select(Protein.IDs, GeoMeanFC_Group1 = geo_mean_fc),
        group2_summary %>% select(Protein.IDs, GeoMeanFC_Group2 = geo_mean_fc),
        by = "Protein.IDs"
      ) %>%
        mutate(
          GeoMeanFC_Ratio = pmax(GeoMeanFC_Group1, GeoMeanFC_Group2) /
            pmin(GeoMeanFC_Group1, GeoMeanFC_Group2)
        )
      
      # DEPs in each fraction
      deps_ai <- de_results_AI[de_results_AI$Comparison == comp & de_results_AI$Change != "No Change", ]
      deps_as <- de_results_AS[de_results_AS$Comparison == comp & de_results_AS$Change != "No Change", ]
      all_deps <- unique(c(deps_ai$Protein.IDs, deps_as$Protein.IDs))
      
      # ---------------------------------------------------------------------
      # Validation in one pass: joins + derivations + final Validation label.
      # Stability is folded into the Validation case_when (no separate filter).
      # ---------------------------------------------------------------------
      validation_table <- tibble(Protein.IDs = all_deps) %>%
        left_join(deps_ai %>% select(Protein.IDs,
                                     Direction_AI  = Change,
                                     Gene.Names_AI = Gene.Names,
                                     Orthologs_AI  = Orthologs),
                  by = "Protein.IDs") %>%
        left_join(deps_as %>% select(Protein.IDs,
                                     Direction_AS  = Change,
                                     Gene.Names_AS = Gene.Names,
                                     Orthologs_AS  = Orthologs),
                  by = "Protein.IDs") %>%
        left_join(comp_summary %>% select(Protein.IDs,
                                          ai_exclusive, as_exclusive,
                                          dominant_fraction, geo_mean_fc),
                  by = "Protein.IDs") %>%
        left_join(stability_data, by = "Protein.IDs") %>%
        mutate(
          is_dep_in_ai = !is.na(Direction_AI),
          is_dep_in_as = !is.na(Direction_AS),
          single_fraction_dep = case_when(
            is_dep_in_ai & !is_dep_in_as ~ "AI",
            is_dep_in_as & !is_dep_in_ai ~ "AS",
            TRUE                         ~ NA_character_
          ),
          dep_type = case_when(
            is_dep_in_ai & !is_dep_in_as                                ~ "AI_only",
            is_dep_in_as & !is_dep_in_ai                                ~ "AS_only",
            is_dep_in_ai & is_dep_in_as & Direction_AI == Direction_AS  ~ "Both_same",
            is_dep_in_ai & is_dep_in_as & Direction_AI != Direction_AS  ~ "Both_opposite"
          ),
          reached_second_level = (is_dep_in_ai | is_dep_in_as) & !(is_dep_in_ai & is_dep_in_as) &
            !is.na(geo_mean_fc) &
            !(ai_exclusive & is_dep_in_ai) &
            !(as_exclusive & is_dep_in_as),
          qualifies_dominant = reached_second_level & geo_mean_fc >= fc_thresh &
            !is.na(dominant_fraction) & dominant_fraction == single_fraction_dep,
          is_stable_between_groups = !is.na(GeoMeanFC_Ratio) & GeoMeanFC_Ratio <= stab_thresh,
          Validation = case_when(
            ai_exclusive & is_dep_in_ai ~
              "Exclusive to AI fraction",
            as_exclusive & is_dep_in_as ~
              "Exclusive to AS fraction",
            dep_type == "Both_same" ~
              "Significant in both fractions (same direction)",
            qualifies_dominant & is_stable_between_groups ~
              paste0("Dominant and stable in ", single_fraction_dep,
                     " fraction (geo-mean FC=", round(geo_mean_fc, 2), ")"),
            TRUE ~ NA_character_
          )
        )
      
      write.csv(validation_table,
                file.path(thresh_dir, paste0("validation_table_", comp, ".csv")),
                row.names = FALSE)
      
      # Counts
      tested_second_level       <- sum(validation_table$reached_second_level, na.rm = TRUE)
      n_dominant_pre_stability  <- sum(validation_table$qualifies_dominant, na.rm = TRUE)
      n_dominant_post_stability <- sum(grepl("Dominant and stable", validation_table$Validation))
      
      cat("  Second-level candidates: ", tested_second_level,
          "; qualifying before stability: ", n_dominant_pre_stability,
          "; after stability (<= ", stab_thresh, "-fold): ", n_dominant_post_stability,
          "\n", sep = "")
      
      # Filter to validated DEPs, then build output columns
      validated_deps <- validation_table %>%
        filter(!is.na(Validation)) %>%
        mutate(
          Direction  = ifelse(is_dep_in_ai, Direction_AI, Direction_AS),
          Gene.Names = ifelse(is_dep_in_ai, Gene.Names_AI, Gene.Names_AS),
          Orthologs  = ifelse(is_dep_in_ai, Orthologs_AI, Orthologs_AS)
        ) %>%
        select(
          Protein.IDs,
          Gene.Names,
          Orthologs,
          Direction,
          DEP_in_AI            = is_dep_in_ai,
          DEP_in_AS            = is_dep_in_as,
          FC_between_fractions = geo_mean_fc,
          Dominant_Fraction    = dominant_fraction,
          Validation,
          GeoMeanFC_Group1,
          GeoMeanFC_Group2,
          GeoMeanFC_Ratio
        )
      
      if (nrow(validated_deps) > 0) {
        validated_deps$Comparison       <- comp
        validated_deps$Comparison_Label <- comparisons_of_interest[comp]
      }
      
      all_validated_deps[[comp]]           <- validated_deps
      dominant_removed[[comp]]             <- n_dominant_pre_stability - n_dominant_post_stability
      tested_second_level_per_comp[[comp]] <- tested_second_level
      
      cat("  Total unique DEPs (AI+AS):", length(all_deps), "\n")
      cat("  Validated DEPs:", nrow(validated_deps), "\n")
    }
    
    # ========================================================================
    # Save validated DEPs
    # ========================================================================
    validated_deps_df <- bind_rows(all_validated_deps)
    write.csv(validated_deps_df,
              file.path(thresh_dir, "validated_DEPs.csv"),
              row.names = FALSE)
    cat("\nSaved validated DEPs to:", file.path(thresh_dir, "validated_DEPs.csv"), "\n")
    
    # ========================================================================
    # Summary statistics
    # ========================================================================
    summary_stats <- data.frame()
    
    for (comp in names(comparisons_of_interest)) {
      # Bone-healing reference proteins among the proteins tested for DE in either
      # fraction. This is the "successes in the population" term of the downstream
      # hypergeometric test, so it must be counted in the same unit as the
      # population itself (distinct Protein.IDs) and as the draw (validated DEP
      # rows). De-duplicating on the Orthologs *string* instead would collapse
      # distinct proteins that happen to share an ortholog list and undercount it.
      tested_proteins <- rbind(
        de_results_AI[de_results_AI$Comparison == comp, c("Protein.IDs", "Orthologs")],
        de_results_AS[de_results_AS$Comparison == comp, c("Protein.IDs", "Orthologs")]
      )
      tested_proteins <- tested_proteins[!duplicated(tested_proteins$Protein.IDs), ]
      bone_caps_tested <- sum(any_ortholog_in(tested_proteins$Orthologs, meta_genes), na.rm = TRUE)
      
      deps_ai_count <- sum(de_results_AI$Comparison == comp & de_results_AI$Change != "No Change")
      deps_as_count <- sum(de_results_AS$Comparison == comp & de_results_AS$Change != "No Change")
      total_unique_deps <- length(unique(c(
        de_results_AI[de_results_AI$Comparison == comp & de_results_AI$Change != "No Change", "Protein.IDs"],
        de_results_AS[de_results_AS$Comparison == comp & de_results_AS$Change != "No Change", "Protein.IDs"]
      )))
      
      all_dep_orthologs <- unique(c(
        de_results_AI[de_results_AI$Comparison == comp & de_results_AI$Change != "No Change", "Orthologs"],
        de_results_AS[de_results_AS$Comparison == comp & de_results_AS$Change != "No Change", "Orthologs"]
      ))
      all_dep_genes <- unique(toupper(trimws(unlist(strsplit(as.character(all_dep_orthologs), ";")))))
      overlap_total <- sum(any_ortholog_in(all_dep_orthologs, meta_genes), na.rm = TRUE)
      
      validated_count <- nrow(all_validated_deps[[comp]])
      
      if (validated_count > 0) {
        validated <- all_validated_deps[[comp]]
        n_both      <- sum(validated$Validation == "Significant in both fractions (same direction)")
        n_dominant_stable  <- sum(grepl("Dominant and stable in", validated$Validation))
        n_exclusive <- sum(grepl("Exclusive to", validated$Validation))
        
        first_level_proteins <- validated[
          validated$Validation == "Significant in both fractions (same direction)" |
            grepl("Exclusive to", validated$Validation), ]
        overlap_first_level <- sum(any_ortholog_in(first_level_proteins$Orthologs, meta_genes), na.rm = TRUE)

        second_level_proteins <- validated[grepl("Dominant and stable in", validated$Validation), ]
        overlap_second_level <- sum(any_ortholog_in(second_level_proteins$Orthologs, meta_genes), na.rm = TRUE)
      } else {
        n_both <- 0; n_dominant_stable <- 0; n_exclusive <- 0
        overlap_first_level <- 0; overlap_second_level <- 0
      }
      
      summary_stats <- rbind(summary_stats, data.frame(
        Comparison                    = comp,
        Comparison_Label              = comparisons_of_interest[comp],
        DEPs_AI                       = deps_ai_count,
        DEPs_AS                       = deps_as_count,
        Bone_Caps_Tested              = bone_caps_tested,
        Total_Unique_DEPs             = total_unique_deps,
        total_unique_orthologs_count  = length(all_dep_genes),
        Overlap_Total                 = overlap_total,
        Validated_Both_Fractions      = n_both,
        Validated_Exclusive           = n_exclusive,
        Validated_Second_Level            = n_dominant_stable,
        Validated_Total               = validated_count,
        Dominant_Removed_By_Stability = dominant_removed[[comp]],
        DEPs_Tested_Second_Level      = tested_second_level_per_comp[[comp]],
        Overlap_First_Level           = overlap_first_level,
        Overlap_Second_Level          = overlap_second_level,
        FC_Threshold                  = fc_thresh,
        Stability_Threshold           = stab_thresh,
        stringsAsFactors              = FALSE
      ))
    }
    
    write.csv(summary_stats,
              file.path(thresh_dir, "summary_statistics.csv"),
              row.names = FALSE)
    
    # ========================================================================
    # Summary plot
    # ========================================================================
    plot_data <- summary_stats %>%
      mutate(
        First_Level_Total  = Validated_Both_Fractions + Validated_Exclusive,
        Second_Level_Total = Validated_Second_Level
      ) %>%
      select(Comparison_Label, Total_Unique_DEPs, Overlap_Total,
             First_Level_Total, Overlap_First_Level,
             Second_Level_Total, Overlap_Second_Level) %>%
      pivot_longer(cols = -Comparison_Label,
                   names_to = "Category", values_to = "Count")
    
    plot_data$Category <- factor(plot_data$Category,
                                 levels = c("Total_Unique_DEPs", "Overlap_Total",
                                            "First_Level_Total", "Overlap_First_Level",
                                            "Second_Level_Total", "Overlap_Second_Level"))
    
    p <- ggplot(plot_data, aes(x = Comparison_Label, y = Count, fill = Category)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_text(aes(label = Count),
                position = position_dodge(width = 0.9),
                vjust = -0.5, size = 3, fontface = "bold") +
      scale_fill_manual(
        values = c("Total_Unique_DEPs"    = "#E63946",
                   "Overlap_Total"        = "#F5A3AD",
                   "First_Level_Total"    = "#2E86AB",
                   "Overlap_First_Level"  = "#A8DADC",
                   "Second_Level_Total"   = "#F77F00",
                   "Overlap_Second_Level" = "#FCBF49"),
        labels = c("Total Unique DEPs (AI + AS)",
                   "Overlap Total with Meta-analysis",
                   "First-Level Validated DEPs",
                   "Overlap First-Level with Meta-analysis",
                   "Second-Level Validated DEPs",
                   "Overlap Second-Level with Meta-analysis")
      ) +
      theme_minimal() +
      labs(
        title    = "DEPs with Two-Level Validation Strategy",
        subtitle = paste0("FC >= ", fc_thresh, " | Stability ratio <= ", stab_thresh, "-fold"),
        y        = "Number of Proteins",
        x        = "",
        fill     = ""
      ) +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.text.y     = element_text(size = 10),
        legend.position = "top",
        plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle   = element_text(hjust = 0.5, size = 10),
        plot.margin     = margin(20, 20, 20, 60)
      )
    
    ggsave(file.path(thresh_dir, "validated_deps_summary.png"),
           plot = p, width = 14, height = 8, dpi = 300)
    
    cat("\n=== Completed FC", fc_thresh, "Stab", stab_thresh, "===\n")
  }
}

cat("\n\n========================================")
cat("\nALL THRESHOLD COMBINATIONS COMPLETE")
cat("\n========================================\n")