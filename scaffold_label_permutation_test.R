# scaffold_label_permutation_test.R
# ==============================================================================
#
# QUESTION
#   The empty-defect comparison (diabetic vs non-diabetic) yields more tissue-level
#   DAPs than the PCL comparison (76 vs 53). If the scaffold made no difference to
#   how diabetes affects the proteome, how often would an analysis like ours produce
#   a difference at least this large?
#
# NULL HYPOTHESIS
#   Within each genotype, the labels "empty defect" and "PCL scaffold" are
#   exchangeable: which animals received a scaffold is irrelevant to the diabetic
#   effect. The sham animals keep their label in every permutation.
#
# INPUT   normalization_results_{AI,AS}/RobNorm_normalized_data.csv + meta_data.csv  all 42 animals per fraction
#         valid_DEPs/FC5.59_Stab1.2/summary_statistics.csv                             self-check only
# OUTPUT  valid_DEPs/FC5.59_Stab1.2/scaffold_label_permutation_test/
#         permutations.csv   per permutation: tissue-level DAPs of both comparisons and their difference
#         summary.csv        observed difference, null distribution and p-value
# ==============================================================================

suppressPackageStartupMessages({ library(limma); library(data.table); library(parallel) })
source("sham_model_differential_abundance.R")   # read_fraction() and fraction_level_test(); its analysis part does not run
source("tissue_level_rules.R")                  # the validation rules of tissue_level_DEP_collection.R

N_PERMUTATIONS            <- 1000
N_CORES                   <- 2
SEED                      <- 20260916
ABUNDANCE_RATIO_THRESHOLD <- 5.59     # AR threshold of the canonical tissue-level result (80th percentile of the background)
STABILITY_THRESHOLD       <- 1.2
threshold_dir <- sprintf("valid_DEPs/FC%s_Stab%s", ABUNDANCE_RATIO_THRESHOLD, STABILITY_THRESHOLD)
output_dir    <- file.path(threshold_dir, "scaffold_label_permutation_test")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Data --------------------------------------------------------------------
# For the fraction-level test: normalized log2 intensities with the animal, genotype and surgery of every sample
fraction_data      <- lapply(setNames(FRACTIONS, FRACTIONS), read_fraction)
defect_samples     <- fraction_data$AI$surgery != "sham"
defect_animals     <- fraction_data$AI$animal[defect_samples]                                  # the 32 bone-resection animals
genotype_of_animal <- setNames(fraction_data$AI$genotype[defect_samples], defect_animals)
scaffold_of_animal <- setNames(fraction_data$AI$surgery[defect_samples],  defect_animals)      # "empty" or "PCL"
sham_animals       <- fraction_data$AI$animal[!defect_samples]                                 # the 10 sham animals; their labels never change
genotype_of_sham   <- setNames(fraction_data$AI$genotype[!defect_samples], sham_animals)

# For the tissue-level rules: the AI/AS ratio of every protein in every sample pair, read and built
# as in tissue_level_DEP_collection.R
matched_samples <- match_sample_pairs(read.csv("normalization_results_AI/meta_data.csv", stringsAsFactors = FALSE),
                                      read.csv("normalization_results_AS/meta_data.csv", stringsAsFactors = FALSE))
ratio_table     <- sample_pair_ratios(read.csv("normalization_results_AI/RobNorm_normalized_data.csv", stringsAsFactors = FALSE, check.names = FALSE),
                                      read.csv("normalization_results_AS/RobNorm_normalized_data.csv", stringsAsFactors = FALSE, check.names = FALSE),
                                      matched_samples)
background      <- background_ratios(ratio_table)
stopifnot(round(quantile(background, 0.80), 2) == ABUNDANCE_RATIO_THRESHOLD)
p99             <- quantile(background, 0.99)     # imputation value for one-fraction-only samples
sample_ids_of   <- function(animal_subset) matched_samples$Sample_ID[matched_samples$Animal %in% animal_subset]

# ---- Tissue-level DAPs of one defect type under a labelling of the defect animals ----
# The defect animals take their surgery label from `scaffold_labels`; the sham animals keep theirs.
tissue_level_daps <- function(scaffold_labels, defect_type) {
  fraction_daps <- lapply(fraction_data, function(data) {
    relabelled   <- unname(scaffold_labels[data$animal])          # NA for the sham animals, because they are not in this vector
    data$surgery <- ifelse(is.na(relabelled), data$surgery, relabelled)
    result <- fraction_level_test(data, defect_type)
    as.data.frame(merge(result[Change != "No Change", .(Protein.IDs, Change)], data$annotation, by = "Protein.IDs"))
  })
  sample_groups <- list(
    group1_defect = sample_ids_of(defect_animals[genotype_of_animal == "diabetic"    & scaffold_labels == defect_type]),
    group1_sham   = sample_ids_of(sham_animals[genotype_of_sham == "diabetic"]),
    group2_defect = sample_ids_of(defect_animals[genotype_of_animal == "nondiabetic" & scaffold_labels == defect_type]),
    group2_sham   = sample_ids_of(sham_animals[genotype_of_sham == "nondiabetic"]))
  validation <- validate_tissue_level(fraction_daps$AI, fraction_daps$AS, sample_groups,
                                      ratio_table, p99, ABUNDANCE_RATIO_THRESHOLD, STABILITY_THRESHOLD)
  sum(!is.na(validation$Validation))
}

# ---- The statistic for one labelling of the defect animals ------------------------
analysis_for_one_labelling <- function(scaffold_labels) {
  counts <- sapply(DEFECT_TYPES, function(defect_type) tissue_level_daps(scaffold_labels, defect_type))
  data.table(tissue_level_DEPs_empty = counts[["empty"]],
             tissue_level_DEPs_PCL   = counts[["PCL"]],
             DEP_count_difference    = counts[["empty"]] - counts[["PCL"]])
}

# ---- Observed result and self-check -----------------------------------------------
observed         <- analysis_for_one_labelling(scaffold_of_animal)
pipeline_summary <- fread(file.path(threshold_dir, "summary_statistics.csv"))
for (defect_type in DEFECT_TYPES) {
  pipeline_row <- pipeline_summary[Comparison == sprintf("diabetic_%s_42-nondiabetic_%s_42", defect_type, defect_type)]
  stopifnot(nrow(pipeline_row) == 1, observed[[paste0("tissue_level_DEPs_", defect_type)]] == pipeline_row$Validated_Total)
}
cat(sprintf("Self-check passed. Observed: %d vs %d tissue-level DAPs (difference %d)\n",
            observed$tissue_level_DEPs_empty, observed$tissue_level_DEPs_PCL, observed$DEP_count_difference))

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
        analysis_for_one_labelling(permuted_scaffold_labels[permutation_number, ])), mc.cores = N_CORES))
fwrite(permutations, file.path(output_dir, "permutations.csv"))

# ---- p-value -----------------------------------------------------------------------
summarise_statistic <- function(name) {
  at_least_as_large <- permutations[[name]] >= observed[[name]]
  data.table(statistic = name, observed = observed[[name]],
             null_mean = mean(permutations[[name]]),
             null_2.5_percentile  = quantile(permutations[[name]], 0.025), null_97.5_percentile = quantile(permutations[[name]], 0.975),
             p_value = mean(at_least_as_large),      # share of permutations at least as large as observed, as in tissue_level_permutation_check.R
             n_permutations = N_PERMUTATIONS)
}
summary_table <- summarise_statistic("DEP_count_difference")
fwrite(summary_table, file.path(output_dir, "summary.csv"))
print(summary_table)
