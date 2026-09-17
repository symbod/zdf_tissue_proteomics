# scaffold_label_permutation_test.R
# ==============================================================================
#
# QUESTION
#   The empty-defect comparison (diabetic vs non-diabetic) yields more tissue-level
#   DEPs than the PCL comparison (148 vs 106), and its tissue-level DEPs overlap
#   more significantly with the bone-healing reference proteins (p = 2.4e-5 vs
#   0.018). If the scaffold made no difference to how diabetes affects the
#   proteome, how often would an analysis like ours produce differences at least
#   this large?
#
# NULL HYPOTHESIS
#   Within each genotype, the labels "empty defect" and "PCL scaffold" are
#   exchangeable: which animals received a scaffold is irrelevant to the diabetic
#   effect.
#
# INPUT   normalization_results_{AI,AS}/scaled/RobNorm_scaled_normalized_data.csv + meta_data.csv
#         normalization_results_{AI,AS}/RobNorm_normalized_data.csv + meta_data.csv
#         input/bone_caps_meta_analysis.csv                                        bone-healing reference proteins
#         valid_DEPs/FC5.59_Stab1.2/summary_statistics.csv                          } self-check only
#         de_analysis_results_{AI,AS}_RobNorm_scaled/de_results_raw_with_orthologs.csv }
# OUTPUT  valid_DEPs/FC5.59_Stab1.2/scaffold_label_permutation_test/
#         permutations.csv   per permutation and comparison: tissue-level DEPs, tested proteins, bone-healing
#                            reference proteins among them, overlap p-value; and the two statistics
#         summary.csv        observed value, null distribution and p-value of both statistics
# ==============================================================================

suppressPackageStartupMessages({ library(PRONE); library(SummarizedExperiment); library(data.table); library(parallel) })

N_PERMUTATIONS            <- 1000
N_CORES                   <- 2        # each worker needs about 1 GB of RAM
SEED                      <- 20260916
ABUNDANCE_RATIO_THRESHOLD <- 5.59     # AR threshold of the canonical tissue-level result (80th percentile of the background)
STABILITY_THRESHOLD       <- 1.2
MIN_VALUES_PER_GROUP      <- 3
COMPARISONS <- list("empty defect" = c(diabetic = "diabetic_empty_42", nondiabetic = "nondiabetic_empty_42"),
                    "PCL scaffold" = c(diabetic = "diabetic_PCL_42",   nondiabetic = "nondiabetic_PCL_42"))
threshold_dir <- sprintf("valid_DEPs/FC%s_Stab%s", ABUNDANCE_RATIO_THRESHOLD, STABILITY_THRESHOLD)
output_dir    <- file.path(threshold_dir, "scaffold_label_permutation_test")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Data --------------------------------------------------------------------
# The DE input as SummarizedExperiment, built like the DE notebook
load_de_input <- function(fraction) {
  directory    <- sprintf("normalization_results_%s/scaled", fraction)
  protein_data <- fread(file.path(directory, "RobNorm_scaled_normalized_data.csv"))
  metadata     <- fread(file.path(directory, "meta_data.csv"))[Column %in% names(protein_data)]
  summarized_experiment <- suppressWarnings(suppressMessages(
    load_data(data = protein_data, md = metadata, protein_column = "Protein.IDs", gene_column = "Gene.Names",
              ref_samples = NULL, batch_column = NULL, condition_column = NULL, label_column = "Label")))
  assays(summarized_experiment)[["log2"]] <- NULL
  names(assays(summarized_experiment))[names(assays(summarized_experiment)) == "raw"] <- "RobNorm_scaled"
  summarized_experiment
}
de_input <- list(AI = load_de_input("AI"), AS = load_de_input("AS"))

# The unscaled intensities for the AI/AS ratios: one matrix per fraction, one column
# per animal, both matrices on the union of all proteins (NA = not measured)
read_intensity_matrix <- function(directory, file_name) {
  protein_data <- fread(file.path(directory, file_name)); metadata <- fread(file.path(directory, "meta_data.csv"))
  intensities <- as.matrix(protein_data[, metadata$Column, with = FALSE])
  rownames(intensities) <- protein_data$Protein.IDs; colnames(intensities) <- metadata$Animal
  intensities
}
unscaled_intensity <- list(AI = read_intensity_matrix("normalization_results_AI", "RobNorm_normalized_data.csv"),
                           AS = read_intensity_matrix("normalization_results_AS", "RobNorm_normalized_data.csv"))
all_proteins          <- union(rownames(unscaled_intensity$AI), rownames(unscaled_intensity$AS))
unscaled_intensity$AI <- unscaled_intensity$AI[match(all_proteins, rownames(unscaled_intensity$AI)), ]
unscaled_intensity$AS <- unscaled_intensity$AS[match(all_proteins, rownames(unscaled_intensity$AS)), colnames(unscaled_intensity$AI)]
rownames(unscaled_intensity$AI) <- rownames(unscaled_intensity$AS) <- all_proteins

sample_metadata    <- fread("normalization_results_AI/scaled/meta_data.csv")           # the 32 defect animals
defect_animals     <- sample_metadata$Animal
genotype_of_animal <- setNames(sample_metadata$Condition, defect_animals)
scaffold_of_animal <- setNames(sample_metadata$Scaffold,  defect_animals)
TIMEPOINT          <- unique(sample_metadata$Timepoint); stopifnot(length(TIMEPOINT) == 1)
group_of_animal    <- function(scaffold_labels) setNames(paste(genotype_of_animal, scaffold_labels, TIMEPOINT, sep = "_"), defect_animals)

background_ratios <- 2^abs(unscaled_intensity$AI - unscaled_intensity$AS)          # all 42 animals, as in tissue_level_DEP_collection.R
stopifnot(round(quantile(background_ratios, 0.80, na.rm = TRUE), 2) == ABUNDANCE_RATIO_THRESHOLD)
RATIO_IF_ONE_FRACTION_ONLY <- quantile(background_ratios, 0.99, na.rm = TRUE)

# Human orthologs of every protein (identical in both fractions) and the bone-healing reference genes
ortholog_table       <- unique(rbind(fread("normalization_results_AI/scaled/RobNorm_scaled_normalized_data.csv")[, .(Protein.IDs, Orthologs)],
                                     fread("normalization_results_AS/scaled/RobNorm_scaled_normalized_data.csv")[, .(Protein.IDs, Orthologs)]))
orthologs_of_protein <- setNames(ortholog_table$Orthologs, ortholog_table$Protein.IDs)
bone_reference_genes <- unique(toupper(fread("input/bone_caps_meta_analysis.csv")$Gene))

# ---- Step 1: differential expression (the DE notebook) ---------------------------
differential_expression <- function(summarized_experiment, groups, diabetic_group, nondiabetic_group) {
  colData(summarized_experiment)$Cond_Scaffold_Time <- unname(groups[colData(summarized_experiment)$Animal])
  intensities     <- as.matrix(assays(summarized_experiment)[["RobNorm_scaled"]])
  group_of_sample <- colData(summarized_experiment)$Cond_Scaffold_Time
  enough_values   <- rowSums(!is.na(intensities[, group_of_sample == diabetic_group]))    >= MIN_VALUES_PER_GROUP &
                     rowSums(!is.na(intensities[, group_of_sample == nondiabetic_group])) >= MIN_VALUES_PER_GROUP
  de_result <- suppressMessages(suppressWarnings(
    run_DE(se = summarized_experiment[enough_values, ], comparisons = paste0(diabetic_group, "-", nondiabetic_group),
           ain = "RobNorm_scaled", condition = "Cond_Scaffold_Time", DE_method = "limma",
           logFC = TRUE, logFC_up = 1, logFC_down = -1, p_adj = TRUE, alpha = 0.05)))
  list(tested_proteins = rowData(summarized_experiment)$Protein.IDs[enough_values],
       deps            = as.data.table(de_result)[Change != "No Change", .(Protein.IDs, Direction = as.character(Change))])
}

# ---- Step 2: tissue-level DEPs (tissue_level_DEP_collection.R) ------------------
# Per-protein summary of the AI/AS relationship over the AI/AS sample pairs of a set of animals
ratio_summary <- function(animals_subset) {
  intensity_AI <- unscaled_intensity$AI[, animals_subset, drop = FALSE]
  intensity_AS <- unscaled_intensity$AS[, animals_subset, drop = FALSE]
  measured  <- !is.na(intensity_AI) | !is.na(intensity_AS)                                  # in at least one fraction
  ai_higher <- is.na(intensity_AS) | (!is.na(intensity_AI) & intensity_AI > intensity_AS)  # only AI measured, or AI higher
  as_higher <- is.na(intensity_AI) | (!is.na(intensity_AS) & intensity_AS > intensity_AI)  # only AS measured, or AS higher
  dominant_per_sample <- ifelse(ai_higher, "AI", ifelse(as_higher, "AS", "Tie")); dominant_per_sample[!measured] <- NA
  ratio_per_sample    <- ifelse(!is.na(intensity_AI) & !is.na(intensity_AS), 2^abs(intensity_AI - intensity_AS), RATIO_IF_ONE_FRACTION_ONLY)
  ratio_per_sample[!measured] <- NA
  data.table(Protein.IDs       = all_proteins,
             ai_exclusive      = rowSums(!is.na(intensity_AS)) == 0,
             as_exclusive      = rowSums(!is.na(intensity_AI)) == 0,
             dominant_fraction = apply(dominant_per_sample, 1, function(higher_fraction_per_animal) {   # one row = one protein
                                 distinct_fractions <- unique(na.omit(higher_fraction_per_animal))
                                 if (length(distinct_fractions) == 1) distinct_fractions else "Inconsistent" }),
             geo_mean_ratio    = exp(rowMeans(log(ratio_per_sample), na.rm = TRUE)))[rowSums(measured) > 0]
}
# The validation rules, applied to every protein that is a DEP in at least one fraction.
# Returns the validated (tissue-level) proteins and the proteins tested in either fraction.
tissue_level_deps <- function(groups, diabetic_group, nondiabetic_group) {
  de_AI <- differential_expression(de_input$AI, groups, diabetic_group, nondiabetic_group)
  de_AS <- differential_expression(de_input$AS, groups, diabetic_group, nondiabetic_group)
  diabetic_animals    <- defect_animals[groups == diabetic_group]
  nondiabetic_animals <- defect_animals[groups == nondiabetic_group]
  candidates <- data.table(Protein.IDs = union(de_AI$deps$Protein.IDs, de_AS$deps$Protein.IDs))
  candidates <- merge(candidates, de_AI$deps[, .(Protein.IDs, Direction_AI = Direction)], all.x = TRUE)
  candidates <- merge(candidates, de_AS$deps[, .(Protein.IDs, Direction_AS = Direction)], all.x = TRUE)
  candidates <- merge(candidates, ratio_summary(c(diabetic_animals, nondiabetic_animals)), all.x = TRUE)   # both groups together
  candidates <- merge(candidates, ratio_summary(diabetic_animals)[, .(Protein.IDs, ratio_diabetic = geo_mean_ratio)], all.x = TRUE)
  candidates <- merge(candidates, ratio_summary(nondiabetic_animals)[, .(Protein.IDs, ratio_nondiabetic = geo_mean_ratio)], all.x = TRUE)
  candidates[, dep_in_AI := !is.na(Direction_AI)]
  candidates[, dep_in_AS := !is.na(Direction_AS)]
  candidates[, single_dep_fraction := fifelse(dep_in_AI & !dep_in_AS, "AI", fifelse(dep_in_AS & !dep_in_AI, "AS", NA_character_))]
  candidates[, stability_ratio := pmax(ratio_diabetic, ratio_nondiabetic) / pmin(ratio_diabetic, ratio_nondiabetic)]
  candidates[, validated := fcase(
    dep_in_AI & ai_exclusive,                             TRUE,   # exclusive to AI
    dep_in_AS & as_exclusive,                             TRUE,   # exclusive to AS
    dep_in_AI & dep_in_AS & Direction_AI == Direction_AS, TRUE,   # significant in both fractions, same direction
    geo_mean_ratio >= ABUNDANCE_RATIO_THRESHOLD & dominant_fraction == single_dep_fraction & stability_ratio <= STABILITY_THRESHOLD, TRUE,   # dominant and stable
    default = FALSE)]
  list(validated_proteins = candidates[validated == TRUE, Protein.IDs],
       tested_proteins    = union(de_AI$tested_proteins, de_AS$tested_proteins))
}

# ---- Step 3: overlap with the bone-healing reference proteins (protein_sets_overlap_bone_CAPs.py) ----
# A protein is a reference protein when any of its human orthologs (semicolon-separated) is in the reference list.
is_bone_reference_protein <- function(protein_ids) {
  vapply(strsplit(orthologs_of_protein[protein_ids], ";"),
         function(orthologs) any(toupper(trimws(orthologs)) %in% bone_reference_genes), logical(1))
}
overlap_with_bone_reference <- function(tested_proteins, validated_proteins) {
  n_tested             <- length(tested_proteins)
  n_reference_tested   <- sum(is_bone_reference_protein(tested_proteins))
  n_validated          <- length(validated_proteins)
  n_reference_overlap  <- sum(is_bone_reference_protein(validated_proteins))
  c(tested_proteins = n_tested, reference_proteins_among_tested = n_reference_tested, reference_protein_overlap = n_reference_overlap,
    overlap_p_value = phyper(n_reference_overlap - 1, n_reference_tested, n_tested - n_reference_tested, n_validated, lower.tail = FALSE))
}

# ---- The two statistics for one labelling of the animals --------------------------
analysis_for_one_labelling <- function(groups) {
  per_comparison <- lapply(COMPARISONS, function(groups_of_comparison) {
    result <- tissue_level_deps(groups, groups_of_comparison["diabetic"], groups_of_comparison["nondiabetic"])
    c(tissue_level_DEPs = length(result$validated_proteins), overlap_with_bone_reference(result$tested_proteins, result$validated_proteins))
  })
  empty <- per_comparison[["empty defect"]]; PCL <- per_comparison[["PCL scaffold"]]
  data.table(t(setNames(empty, paste0(names(empty), "_empty"))), t(setNames(PCL, paste0(names(PCL), "_PCL"))),
             DEP_count_difference            = unname(empty["tissue_level_DEPs"] - PCL["tissue_level_DEPs"]),           # statistic (1)
             overlap_significance_difference = unname(log10(PCL["overlap_p_value"]) - log10(empty["overlap_p_value"])))  # statistic (2)
}

# ---- Observed result and self-check -----------------------------------------------
observed         <- analysis_for_one_labelling(group_of_animal(scaffold_of_animal))
pipeline_summary <- fread(file.path(threshold_dir, "summary_statistics.csv"))
for (comparison_name in names(COMPARISONS)) {
  pipeline_row <- pipeline_summary[Comparison == paste(COMPARISONS[[comparison_name]], collapse = "-")]
  suffix       <- if (comparison_name == "empty defect") "_empty" else "_PCL"
  tested_in_pipeline <- unlist(lapply(c("AI", "AS"), function(fraction)             # the population of the overlap test
    fread(sprintf("de_analysis_results_%s_RobNorm_scaled/de_results_raw_with_orthologs.csv", fraction))[Comparison == pipeline_row$Comparison, Protein.IDs]))
  stopifnot(observed[[paste0("tissue_level_DEPs", suffix)]]               == pipeline_row$Validated_Total,
            observed[[paste0("tested_proteins", suffix)]]                 == uniqueN(tested_in_pipeline),
            observed[[paste0("reference_proteins_among_tested", suffix)]] == pipeline_row$Bone_Caps_Tested,
            observed[[paste0("reference_protein_overlap", suffix)]]       == pipeline_row$Overlap_First_Level + pipeline_row$Overlap_Second_Level)
}
cat(sprintf("Self-check passed. Observed: %d vs %d tissue-level DEPs (difference %d); overlap %d of %d tested (p = %.2g) vs %d of %d (p = %.2g), significance difference %.2f\n",
            observed$tissue_level_DEPs_empty, observed$tissue_level_DEPs_PCL, observed$DEP_count_difference,
            observed$reference_protein_overlap_empty, observed$tested_proteins_empty, observed$overlap_p_value_empty,
            observed$reference_protein_overlap_PCL,   observed$tested_proteins_PCL,   observed$overlap_p_value_PCL, observed$overlap_significance_difference))

# ---- Permutations -----------------------------------------------------------------
set.seed(SEED)
permuted_scaffold_labels <- t(replicate(N_PERMUTATIONS, {      # one relabelling per permutation, one row each
  scaffold_labels <- scaffold_of_animal
  for (genotype in c("diabetic", "nondiabetic")) {
    animals_of_genotype <- genotype_of_animal == genotype
    scaffold_labels[animals_of_genotype] <- sample(scaffold_labels[animals_of_genotype])   # shuffles the 9 "empty" and 7 "PCL" labels within the genotype
  }
  scaffold_labels
}))
cat(sprintf("Running %d permutations on %d cores ...\n", N_PERMUTATIONS, N_CORES))
permutations <- rbindlist(mclapply(seq_len(N_PERMUTATIONS), function(permutation_number)
  cbind(data.table(permutation = permutation_number),
        analysis_for_one_labelling(group_of_animal(permuted_scaffold_labels[permutation_number, ]))), mc.cores = N_CORES))
fwrite(permutations, file.path(output_dir, "permutations.csv"))

# ---- p-values ----------------------------------------------------------------------
summarise_statistic <- function(name) {
  at_least_as_large <- permutations[[name]] >= observed[[name]]
  data.table(statistic = name, observed = observed[[name]],
             null_mean = mean(permutations[[name]]),
             null_2.5_percentile  = quantile(permutations[[name]], 0.025), null_97.5_percentile = quantile(permutations[[name]], 0.975),
             p_value = (sum(at_least_as_large) + 1) / (N_PERMUTATIONS + 1), n_permutations = N_PERMUTATIONS)
}
summary_table <- rbind(summarise_statistic("DEP_count_difference"), summarise_statistic("overlap_significance_difference"))
fwrite(summary_table, file.path(output_dir, "summary.csv"))
print(summary_table)
