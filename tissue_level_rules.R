# tissue_level_rules.R
#
# The tissue-level validation rules of tissue_level_DEP_collection.R as functions, so that the
# permutation scripts (tissue_level_permutation_check.R, scaffold_label_permutation_test.R)
# apply exactly the same rules. Nothing is read or written here.
#
#   match_sample_pairs(meta_ai, meta_as)
#       the AI/AS sample pairs (same animal, condition, scaffold, intervention, timepoint, tissue)
#   sample_pair_ratios(ai_data, as_data, matched_samples)
#       the AI/AS ratio of every protein in every sample pair (one-fraction-only detections kept
#       with NA ratio, imputed later)
#   background_ratios(fc_table)
#       the background distribution: the ratios of all proteins measured in both fractions
#   summarise_proteins_across_samples(sample_ids, fc_table, p99_value)
#       per-protein summary over a set of samples: exclusivity, dominant fraction, geometric mean ratio
#   share_contrast(sample_groups, fc_table, p99_value)
#       per protein, the sham-corrected contrast applied to the log2 share of the protein in a fraction:
#       the part of the fraction-level effect that comes from a change in fraction distribution
#   validate_tissue_level(deps_ai, deps_as, sample_groups, fc_table, p99_value, fc_thresh, stab_thresh)
#       the validation table of one comparison: every protein that is a DAP in at least one
#       fraction, with its Validation label (NA = not validated); sample_groups = the four groups
#       of the sham-corrected contrast

suppressPackageStartupMessages(library(dplyr))

# Match AI and AS samples
match_sample_pairs <- function(meta_ai, meta_as) {
  meta_ai$Sample_ID <- paste(meta_ai$Animal, meta_ai$Condition, meta_ai$Scaffold,
                             meta_ai$Intervention, meta_ai$Timepoint, meta_ai$Tissue, sep = "_")
  meta_as$Sample_ID <- paste(meta_as$Animal, meta_as$Condition, meta_as$Scaffold,
                             meta_as$Intervention, meta_as$Timepoint, meta_as$Tissue, sep = "_")
  inner_join(
    meta_ai %>% dplyr::select(Sample_ID, Column_AI = Column,
                              Condition, Scaffold, Intervention, Timepoint, Tissue, Animal),
    meta_as %>% dplyr::select(Sample_ID, Column_AS = Column),
    by = "Sample_ID"
  )
}

# Build per-sample AI/AS ratios for every protein in every sample pair.
# Retain one-fraction-only detections (NA Fold_Change) so they can be imputed later.
sample_pair_ratios <- function(ai_data, as_data, matched_samples) {
  all_fc <- list()
  for (i in 1:nrow(matched_samples)) {
    ai_col <- matched_samples$Column_AI[i]
    as_col <- matched_samples$Column_AS[i]

    ai_sample <- data.frame(
      Protein.IDs = ai_data$Protein.IDs,
      AI_abundance = as.numeric(ai_data[[ai_col]]),
      stringsAsFactors = FALSE
    )
    as_sample <- data.frame(
      Protein.IDs = as_data$Protein.IDs,
      AS_abundance = as.numeric(as_data[[as_col]]),
      stringsAsFactors = FALSE
    )

    merged <- full_join(ai_sample, as_sample, by = "Protein.IDs") %>%
      dplyr::filter(!is.na(AI_abundance) | !is.na(AS_abundance))

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
  bind_rows(all_fc)
}

# The background fold-change distribution: every protein in every sample pair with both fractions measured
background_ratios <- function(fc_table) {
  background_fc <- fc_table$Fold_Change
  background_fc[!is.na(background_fc) & is.finite(background_fc)]
}

# Per-protein summary for an arbitrary set of sample IDs, derived from the
# sample-pair ratio table (the single source of truth for per-sample AI/AS data).
#   - n_samples:         how many sample rows contribute
#   - ai_exclusive:      never detected in AS across this sample set
#   - as_exclusive:      never detected in AI across this sample set
#   - dominant_fraction: "AI"/"AS"/"Tie" if all samples agree, else "Inconsistent"
#   - geo_mean_fc:       geometric mean of per-sample AI/AS ratios
#                        (one-fraction-only samples imputed to p99)
summarise_proteins_across_samples <- function(sample_ids, fc_table, p99_value) {
  fc_table %>%
    dplyr::filter(Sample_ID %in% sample_ids) %>%
    mutate(fc_imputed = ifelse(is.na(Fold_Change), p99_value, Fold_Change)) %>%
    group_by(Protein.IDs) %>%
    summarise(
      n_samples         = n(),
      ai_exclusive      = all(is.na(AS_linear)),
      as_exclusive      = all(is.na(AI_linear)),
      dominant_fraction = if (n_distinct(Dominant_Fraction) == 1)
        dplyr::first(Dominant_Fraction) else "Inconsistent",
      geo_mean_fc       = exp(mean(log(fc_imputed))),
      .groups           = "drop"
    )
}

# Geometric-mean share of a protein in each fraction over a set of samples (mean of the log2 shares,
# share = AI / (AI + AS) on the linear scale). Samples in which the protein was quantified in one
# fraction only get the share implied by the imputed ratio p99, as in the dominance rule.
mean_log2_shares <- function(sample_ids, fc_table, p99_value) {
  fc_table %>%
    dplyr::filter(Sample_ID %in% sample_ids) %>%
    mutate(ai_share = case_when(
      !is.na(AI_linear) & !is.na(AS_linear) ~ AI_linear / (AI_linear + AS_linear),
      !is.na(AI_linear)                     ~ p99_value / (p99_value + 1),
      TRUE                                  ~ 1 / (p99_value + 1)
    )) %>%
    group_by(Protein.IDs) %>%
    summarise(log2_ai_share = mean(log2(ai_share)),
              log2_as_share = mean(log2(1 - ai_share)),
              .groups = "drop")
}

# The share contrast. The abundance of a protein measured in a fraction is its tissue abundance times
# its share in that fraction, so on the log2 scale the sham-corrected contrast of the DA,
#     (group1 defect - group1 sham) - (group2 defect - group2 sham),
# applied to the log2 share is the part of the fraction-level effect that comes from a change in
# fraction distribution rather than in tissue abundance. sample_groups: named list of the sample IDs
# of the four groups (group1_defect, group1_sham, group2_defect, group2_sham). Returns per protein the
# geometric-mean AI share of every group and the contrast for the AI share and for the AS share (log2).
share_contrast <- function(sample_groups, fc_table, p99_value) {
  shares <- lapply(sample_groups, mean_log2_shares, fc_table = fc_table, p99_value = p99_value)
  shares$group1_defect %>% dplyr::rename(g1d_ai = log2_ai_share, g1d_as = log2_as_share) %>%
    inner_join(shares$group1_sham   %>% dplyr::rename(g1s_ai = log2_ai_share, g1s_as = log2_as_share), by = "Protein.IDs") %>%
    inner_join(shares$group2_defect %>% dplyr::rename(g2d_ai = log2_ai_share, g2d_as = log2_as_share), by = "Protein.IDs") %>%
    inner_join(shares$group2_sham   %>% dplyr::rename(g2s_ai = log2_ai_share, g2s_as = log2_as_share), by = "Protein.IDs") %>%
    transmute(Protein.IDs,
              AIShare_Group1_Defect = 2^g1d_ai, AIShare_Group1_Sham = 2^g1s_ai,
              AIShare_Group2_Defect = 2^g2d_ai, AIShare_Group2_Sham = 2^g2s_ai,
              share_contrast_AI_log2 = (g1d_ai - g1s_ai) - (g2d_ai - g2s_ai),
              share_contrast_AS_log2 = (g1d_as - g1s_as) - (g2d_as - g2s_as))
}

# The validation table of one comparison. deps_ai / deps_as: the DAPs of the comparison in each
# fraction (columns Protein.IDs, Change, Gene.Names, Orthologs); sample_groups: named list of the
# sample IDs of the four groups of the sham-corrected contrast (group1_defect, group1_sham,
# group2_defect, group2_sham; group1 = diabetic, group2 = non-diabetic); fc_thresh: the AR threshold;
# stab_thresh: the largest fold a change in fraction distribution may contribute to the effect.
validate_tissue_level <- function(deps_ai, deps_as, sample_groups, fc_table, p99_value, fc_thresh, stab_thresh) {
  # Exclusivity and dominance over all samples of the four groups; stability from the share contrast
  comp_summary <- summarise_proteins_across_samples(unlist(sample_groups), fc_table, p99_value)
  shares       <- share_contrast(sample_groups, fc_table, p99_value)

  all_deps <- unique(c(deps_ai$Protein.IDs, deps_as$Protein.IDs))

  # ---------------------------------------------------------------------
  # Validation in one pass: joins + derivations + final Validation label.
  # Stability is folded into the Validation case_when (no separate filter).
  # ---------------------------------------------------------------------
  tibble(Protein.IDs = all_deps) %>%
    left_join(deps_ai %>% dplyr::select(Protein.IDs,
                                        Direction_AI  = Change,
                                        Gene.Names_AI = Gene.Names,
                                        Orthologs_AI  = Orthologs),
              by = "Protein.IDs") %>%
    left_join(deps_as %>% dplyr::select(Protein.IDs,
                                        Direction_AS  = Change,
                                        Gene.Names_AS = Gene.Names,
                                        Orthologs_AS  = Orthologs),
              by = "Protein.IDs") %>%
    left_join(comp_summary %>% dplyr::select(Protein.IDs,
                                             ai_exclusive, as_exclusive,
                                             dominant_fraction, geo_mean_fc),
              by = "Protein.IDs") %>%
    left_join(shares, by = "Protein.IDs") %>%
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
      # Stability: the share contrast of the fraction in which the protein is a DAP must not
      # exceed stab_thresh-fold, i.e. a change in fraction distribution explains at most that much
      ShareContrast_log2 = case_when(
        single_fraction_dep == "AI" ~ share_contrast_AI_log2,
        single_fraction_dep == "AS" ~ share_contrast_AS_log2,
        TRUE                        ~ NA_real_
      ),
      ShareContrast_Fold = 2^abs(ShareContrast_log2),
      is_stable = !is.na(ShareContrast_Fold) & ShareContrast_Fold <= stab_thresh,
      Validation = case_when(
        ai_exclusive & is_dep_in_ai ~
          "Exclusive to AI fraction",
        as_exclusive & is_dep_in_as ~
          "Exclusive to AS fraction",
        dep_type == "Both_same" ~
          "Significant in both fractions (same direction)",
        qualifies_dominant & is_stable ~
          paste0("Dominant and stable in ", single_fraction_dep,
                 " fraction (geo-mean FC=", round(geo_mean_fc, 2), ")"),
        TRUE ~ NA_character_
      )
    )
}
