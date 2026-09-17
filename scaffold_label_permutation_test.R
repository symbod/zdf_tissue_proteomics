# scaffold_label_permutation_test.R
# ==============================================================================
#
# QUESTION
#   The empty-defect comparison (diabetic vs non-diabetic) yields more tissue-level
#   DEPs than the PCL comparison (148 vs 106). If the scaffold made no difference
#   to how diabetes affects the proteome, how often would an analysis like ours
#   produce a difference at least this large?
#
# NULL HYPOTHESIS
#   Within each genotype, the labels "empty defect" and "PCL scaffold" are
#   exchangeable: which animals received a scaffold is irrelevant to the diabetic
#   effect.
#
# INPUT   normalization_results_{AI,AS}/scaled/RobNorm_scaled_normalized_data.csv + meta_data.csv
#         normalization_results_{AI,AS}/RobNorm_normalized_data.csv + meta_data.csv
#         valid_DEPs/FC5.59_Stab1.2/summary_statistics.csv                          self-check only
# OUTPUT  valid_DEPs/FC5.59_Stab1.2/scaffold_label_permutation_test/
#         permutations.csv   tissue-level DEP counts of both comparisons and their difference, per permutation
#         summary.csv        observed difference, its null distribution and the p-value
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
  as.data.table(de_result)[Change != "No Change", .(Protein.IDs, Direction = as.character(Change))]   # the DEPs
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
# The validation rules, applied to every protein that is a DEP in at least one fraction; returns the number validated
count_tissue_level_deps <- function(groups, diabetic_group, nondiabetic_group) {
  deps_AI <- differential_expression(de_input$AI, groups, diabetic_group, nondiabetic_group)
  deps_AS <- differential_expression(de_input$AS, groups, diabetic_group, nondiabetic_group)
  diabetic_animals    <- defect_animals[groups == diabetic_group]
  nondiabetic_animals <- defect_animals[groups == nondiabetic_group]
  candidates <- data.table(Protein.IDs = union(deps_AI$Protein.IDs, deps_AS$Protein.IDs))
  candidates <- merge(candidates, deps_AI[, .(Protein.IDs, Direction_AI = Direction)], all.x = TRUE)
  candidates <- merge(candidates, deps_AS[, .(Protein.IDs, Direction_AS = Direction)], all.x = TRUE)
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
  sum(candidates$validated)
}

# ---- The statistic for one labelling of the animals ------------------------------
dep_count_difference <- function(groups) {
  n_tissue_level_deps <- sapply(COMPARISONS, function(groups_of_comparison)
    count_tissue_level_deps(groups, groups_of_comparison["diabetic"], groups_of_comparison["nondiabetic"]))
  data.table(tissue_level_DEPs_empty    = n_tissue_level_deps[["empty defect"]],
             tissue_level_DEPs_PCL      = n_tissue_level_deps[["PCL scaffold"]],
             difference_empty_minus_PCL = n_tissue_level_deps[["empty defect"]] - n_tissue_level_deps[["PCL scaffold"]])
}

# ---- Observed result and self-check -----------------------------------------------
observed         <- dep_count_difference(group_of_animal(scaffold_of_animal))
pipeline_summary <- fread(file.path(threshold_dir, "summary_statistics.csv"))
stopifnot(observed$tissue_level_DEPs_empty == pipeline_summary[Comparison == paste(COMPARISONS[["empty defect"]], collapse = "-"), Validated_Total],
          observed$tissue_level_DEPs_PCL   == pipeline_summary[Comparison == paste(COMPARISONS[["PCL scaffold"]], collapse = "-"), Validated_Total])
cat(sprintf("Self-check passed. Observed: %d vs %d tissue-level DEPs, difference %d\n",
            observed$tissue_level_DEPs_empty, observed$tissue_level_DEPs_PCL, observed$difference_empty_minus_PCL))

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
        dep_count_difference(group_of_animal(permuted_scaffold_labels[permutation_number, ]))), mc.cores = N_CORES))
fwrite(permutations, file.path(output_dir, "permutations.csv"))

# ---- p-value -----------------------------------------------------------------------
at_least_as_large <- permutations$difference_empty_minus_PCL >= observed$difference_empty_minus_PCL
summary_table <- data.table(
  statistic = "tissue-level DEPs empty defect minus PCL scaffold",
  observed  = observed$difference_empty_minus_PCL,
  null_mean = mean(permutations$difference_empty_minus_PCL),
  p_value = (sum(at_least_as_large) + 1) / (N_PERMUTATIONS + 1),
  n_permutations = N_PERMUTATIONS)
fwrite(summary_table, file.path(output_dir, "summary.csv"))
print(summary_table)
