# tissue_level_permutation_check.R
#
# How many tissue-level DAPs are expected by chance? The fraction-level test uses nominal p-values
# (sham_model_differential_abundance.R), so this check shuffles the genotype labels per animal
# within each surgery group (sham, empty, PCL), with the same shuffle for both fractions because
# they come from the same animals, repeats the fraction-level test and the tissue-level rules of
# tissue_level_DEP_collection.R (canonical thresholds), and counts the tissue-level DAPs per
# shuffle and comparison.
#
# Input : normalization_results_{AI,AS}/RobNorm_normalized_data.csv + meta_data.csv
#         valid_DEPs/FC5.59_Stab1.2/summary_statistics.csv (self-check only)
# Output: valid_DEPs/FC5.59_Stab1.2/tissue_level_permutation_check/
#         permutations.csv  counts per shuffle and comparison: fraction-level DAPs (AI, AS) and
#                           tissue-level DAPs by rule (both fractions, exclusive, second level)
#         summary.csv       observed count of tissue-level DAPs, null distribution and share of
#                           shuffles with at least the observed count
# Runtime: about 15 min on 4 cores (1000 shuffles).

suppressPackageStartupMessages({ library(limma); library(data.table); library(parallel) })
source("sham_model_differential_abundance.R")   # read_fraction() and fraction_level_test(); its analysis part does not run
source("tissue_level_rules.R")                  # the validation rules of tissue_level_DEP_collection.R

N_SHUFFLES                <- 1000
N_CORES                   <- 4
SEED                      <- 1
ABUNDANCE_RATIO_THRESHOLD <- 5.59     # canonical thresholds of tissue_level_DEP_collection.R
STABILITY_THRESHOLD       <- 1.2
threshold_dir <- sprintf("valid_DEPs/FC%s_Stab%s", ABUNDANCE_RATIO_THRESHOLD, STABILITY_THRESHOLD)
output_dir    <- file.path(threshold_dir, "tissue_level_permutation_check")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Data --------------------------------------------------------------------
# For the fraction-level test: log2 intensities with the animal, genotype and surgery of every sample
fraction_data      <- lapply(setNames(FRACTIONS, FRACTIONS), read_fraction)
animals            <- unique(fraction_data$AI$animal)
genotype_of_animal <- fraction_data$AI$genotype[match(animals, fraction_data$AI$animal)]
surgery_of_animal  <- fraction_data$AI$surgery[match(animals, fraction_data$AI$animal)]

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

# ---- The DAPs of one defect type under a genotype labelling of the animals ----
# Fraction-level DAPs of both fractions and tissue-level DAPs by validation rule
dap_counts <- function(genotype, defect_type) {
  fraction_daps <- lapply(fraction_data, function(data) {
    result <- fraction_level_test(data, defect_type, genotype[match(data$animal, animals)])
    as.data.frame(merge(result[Change != "No Change", .(Protein.IDs, Change)], data$annotation, by = "Protein.IDs"))
  })
  samples_in <- function(genotype_label, surgery) sample_ids_of(animals[genotype == genotype_label & surgery_of_animal == surgery])
  sample_groups <- list(group1_defect = samples_in("diabetic",    defect_type), group1_sham = samples_in("diabetic",    "sham"),
                        group2_defect = samples_in("nondiabetic", defect_type), group2_sham = samples_in("nondiabetic", "sham"))
  validation <- validate_tissue_level(fraction_daps$AI, fraction_daps$AS, sample_groups,
                                      ratio_table, p99, ABUNDANCE_RATIO_THRESHOLD, STABILITY_THRESHOLD)
  counts <- c(DAPs_AI        = nrow(fraction_daps$AI),
              DAPs_AS        = nrow(fraction_daps$AS),
              both_fractions = sum(validation$Validation == "Significant in both fractions (same direction)", na.rm = TRUE),
              exclusive      = sum(grepl("Exclusive to", validation$Validation)),
              second_level   = sum(grepl("Dominant and stable", validation$Validation)))
  setNames(counts, paste(defect_type, names(counts), sep = "_"))
}
dap_counts_both_comparisons <- function(genotype) unlist(lapply(DEFECT_TYPES, function(defect_type) dap_counts(genotype, defect_type)))

tissue_level_count <- function(counts, defect_type) {
  counts[[paste0(defect_type, "_both_fractions")]] + counts[[paste0(defect_type, "_exclusive")]] + counts[[paste0(defect_type, "_second_level")]]
}

# ---- Self-check: the real labels reproduce summary_statistics.csv ------------
observed           <- as.list(dap_counts_both_comparisons(genotype_of_animal))
summary_statistics <- fread(file.path(threshold_dir, "summary_statistics.csv"))
for (defect_type in DEFECT_TYPES) {
  expected <- summary_statistics[Comparison == sprintf("diabetic_%s_42-nondiabetic_%s_42", defect_type, defect_type)]
  stopifnot(nrow(expected) == 1,
            observed[[paste0(defect_type, "_DAPs_AI")]]        == expected$DEPs_AI,
            observed[[paste0(defect_type, "_DAPs_AS")]]        == expected$DEPs_AS,
            observed[[paste0(defect_type, "_both_fractions")]] == expected$Validated_Both_Fractions,
            observed[[paste0(defect_type, "_exclusive")]]      == expected$Validated_Exclusive,
            observed[[paste0(defect_type, "_second_level")]]   == expected$Validated_Second_Level)
}
cat("Self-check passed: the real labels reproduce summary_statistics.csv\n")

# ---- Shuffles ----------------------------------------------------------------
set.seed(SEED)
shuffled_genotypes <- replicate(N_SHUFFLES, {
  genotype <- genotype_of_animal
  for (surgery in unique(surgery_of_animal)) {
    animals_of_surgery <- which(surgery_of_animal == surgery)
    genotype[animals_of_surgery] <- sample(genotype[animals_of_surgery])
  }
  genotype
}, simplify = FALSE)
null_counts <- rbindlist(mclapply(shuffled_genotypes, function(genotype) as.list(dap_counts_both_comparisons(genotype)), mc.cores = N_CORES))
null_counts[, shuffle := .I]
setcolorder(null_counts, "shuffle")
fwrite(null_counts, file.path(output_dir, "permutations.csv"))

# ---- Summary -----------------------------------------------------------------
summary <- rbindlist(lapply(DEFECT_TYPES, function(defect_type) {
  observed_count <- tissue_level_count(observed, defect_type)
  null_count     <- tissue_level_count(null_counts, defect_type)
  data.table(defect_type = defect_type, observed = observed_count,
             null_mean = mean(null_count), null_median = median(null_count),
             null_95th_percentile = unname(quantile(null_count, 0.95)), null_max = max(null_count),
             share_of_shuffles_at_least_observed = mean(null_count >= observed_count), shuffles = N_SHUFFLES)
}))
fwrite(summary, file.path(output_dir, "summary.csv"))
for (i in seq_len(nrow(summary))) with(summary[i], cat(sprintf(
  "%s defect: observed %3d tissue-level DAPs | by chance: mean %5.1f, 95th percentile %3g, maximum %3d (%d shuffles) | share of shuffles with at least the observed number: %.3f\n",
  defect_type, observed, null_mean, null_95th_percentile, null_max, shuffles, share_of_shuffles_at_least_observed)))
