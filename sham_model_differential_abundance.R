# sham_model_differential_abundance.R
#
# Differential abundance per fraction with the sham animals inside the limma model.
# Here all six groups (genotype x empty / PCL / sham) are in one model per fraction, and the
# diabetic effect is the contrast
#     (diabetic_defect - diabetic_sham) - (nondiabetic_defect - nondiabetic_sham)
#
# Significance at the fraction level is nominal (p < 0.05 and |log2FC| >= 1, no FDR adjustment):
# with 5 sham animals per genotype the FDR-adjusted analysis yields no protein. False positives
# are controlled at the tissue level (tissue_level_DEP_collection.R); tissue_level_permutation_check.R
# quantifies how many tissue-level DAPs are expected by chance.
#
# Input : normalization_results_{AI,AS}/RobNorm_normalized_data.csv + meta_data.csv;
#         LoessF_normalized_data.csv is read only for the robustness check below.
# Output: de_analysis_results_{AI,AS}_sham_model/{de_results_raw, de_results_significant,
#         de_results_raw_with_orthologs}.csv, in the format of the previous step. The comparison
#         labels are kept (e.g. "diabetic_empty_42-nondiabetic_empty_42") so that the downstream
#         scripts run unchanged; they now denote the sham-baseline-corrected comparison.
#         normalization_robustness_jaccard.csv: agreement of the DAP sets with those of the LoessF
#         normalization (the test is repeated on the LoessF data; only this table is kept of it).

suppressMessages({library(limma); library(data.table)})

MIN_SAMPLES_PER_GROUP <- 3      # as in the previous pipeline
LOG2FC_THRESHOLD      <- 1      # as in the previous pipeline
P_VALUE_THRESHOLD     <- 0.05   # nominal p-value, see header
DEFECT_TYPES          <- c("empty", "PCL")
FRACTIONS             <- c("AI", "AS")

# Read one fraction: log2 intensities (proteins x samples), the animal, genotype and surgery
# (sham / empty / PCL) of every sample, and the protein annotation. RobNorm is the primary
# normalization; the LoessF data are read only for the robustness check.
read_fraction <- function(fraction, normalization = "RobNorm") {
  data <- fread(sprintf("normalization_results_%s/%s_normalized_data.csv", fraction, normalization))
  meta <- fread(sprintf("normalization_results_%s/meta_data.csv", fraction))
  intensities <- as.matrix(data[, meta$Column, with = FALSE])
  rownames(intensities) <- data$Protein.IDs
  list(intensities = intensities,
       animal      = meta$Animal,
       genotype    = meta$Condition,
       surgery     = ifelse(meta$Tissue == "ctrl_42", "sham", meta$Scaffold),
       annotation  = data[, .(Protein.IDs, Gene.Names, Orthologs)])
}

# limma with all six groups of one fraction. Returns the baseline-corrected diabetic effect of
# every protein with enough values, for one defect type. The genotype vector can be replaced
# (used by tissue_level_permutation_check.R to shuffle the labels).
fraction_level_test <- function(fraction_data, defect_type, genotype = fraction_data$genotype) {
  intensities <- fraction_data$intensities
  group       <- paste(genotype, fraction_data$surgery, sep = "_")
  defect_groups <- paste(c("diabetic", "nondiabetic"), defect_type, sep = "_")
  sham_groups   <- c("diabetic_sham", "nondiabetic_sham")

  enough_values <- function(g) rowSums(!is.na(intensities[, group == g, drop = FALSE])) >= MIN_SAMPLES_PER_GROUP
  keep <- enough_values(defect_groups[1]) & enough_values(defect_groups[2]) &
          enough_values(sham_groups[1])   & enough_values(sham_groups[2])

  design <- model.matrix(~ 0 + factor(group))
  colnames(design) <- levels(factor(group))
  contrast <- makeContrasts(contrasts = sprintf("(%s-%s)-(%s-%s)", defect_groups[1], sham_groups[1],
                                                defect_groups[2], sham_groups[2]), levels = design)
  fit    <- eBayes(contrasts.fit(lmFit(intensities[keep, ], design), contrast), trend = TRUE, robust = TRUE)
  result <- topTable(fit, number = Inf, sort.by = "none")

  significant <- result$P.Value < P_VALUE_THRESHOLD & abs(result$logFC) >= LOG2FC_THRESHOLD
  data.table(Protein.IDs = rownames(result),
             logFC = result$logFC, P.Value = result$P.Value, adj.P.Val = result$adj.P.Val,
             Change = ifelse(!significant, "No Change", ifelse(result$logFC > 0, "Up Regulated", "Down Regulated")),
             Comparison = sprintf("diabetic_%s_42-nondiabetic_%s_42", defect_type, defect_type))
}

# Robustness check: the same test on the LoessF-normalized data of the fraction, compared with the
# RobNorm results per comparison (Jaccard index of the DAP sets, proteins called in opposite directions).
normalization_agreement <- function(robnorm_data, loessf_data, robnorm_results) {
  loessf_data$intensities <- loessf_data$intensities[rownames(robnorm_data$intensities), ]
  stopifnot(identical(is.na(loessf_data$intensities), is.na(robnorm_data$intensities)))   # same proteins pass the filter
  rbindlist(lapply(DEFECT_TYPES, function(defect_type) {
    loessf_results <- fraction_level_test(loessf_data, defect_type)
    comparison   <- unique(loessf_results$Comparison)
    daps_robnorm <- robnorm_results[Comparison == comparison & Change != "No Change"]
    daps_loessf  <- loessf_results[Change != "No Change"]
    shared <- merge(daps_robnorm, daps_loessf, by = "Protein.IDs", suffixes = c("_robnorm", "_loessf"))
    data.table(Comparison = comparison,
               DAPs_RobNorm = nrow(daps_robnorm), DAPs_LoessF = nrow(daps_loessf), DAPs_shared = nrow(shared),
               Jaccard = round(nrow(shared) / length(union(daps_robnorm$Protein.IDs, daps_loessf$Protein.IDs)), 3),
               opposite_direction = sum(shared$Change_robnorm != shared$Change_loessf))
  }))
}

# Both comparisons of one fraction, written in the format of the previous step, plus the robustness check.
run_differential_abundance <- function(fraction) {
  fraction_data <- read_fraction(fraction)
  results <- rbindlist(lapply(DEFECT_TYPES, function(defect_type) fraction_level_test(fraction_data, defect_type)))
  agreement <- normalization_agreement(fraction_data, read_fraction(fraction, "LoessF"), results)
  results <- merge(results, fraction_data$annotation, by = "Protein.IDs")
  results[, Assay := "RobNorm"]
  results[, IDs := as.character(.I)]

  out_dir <- sprintf("de_analysis_results_%s_sham_model", fraction)
  dir.create(out_dir, showWarnings = FALSE)
  raw <- results[, .(Protein.IDs, logFC, P.Value, adj.P.Val, Change, Comparison, Assay, Gene.Names, IDs)]
  # write.csv as in the previous step: every text cell quoted, so that cells such as "SMAD2;SMAD3" survive a
  # spreadsheet import that also splits on semicolons
  write.csv(raw, file.path(out_dir, "de_results_raw.csv"), row.names = FALSE)
  write.csv(raw[Change != "No Change"], file.path(out_dir, "de_results_significant.csv"), row.names = FALSE)
  write.csv(results[, .(Protein.IDs, Gene.Names, logFC, P.Value, adj.P.Val, Change, Comparison, Assay, IDs, Orthologs)],
            file.path(out_dir, "de_results_raw_with_orthologs.csv"), row.names = FALSE)
  write.csv(agreement, file.path(out_dir, "normalization_robustness_jaccard.csv"), row.names = FALSE)

  for (comparison in unique(results$Comparison)) {
    r <- results[Comparison == comparison]
    a <- agreement[Comparison == comparison]
    cat(sprintf("%s %-40s tested %4d | changed (p < %.2f, |log2FC| >= %g): %3d (down %3d, up %3d) | adj.P < 0.05: %d | LoessF agreement: Jaccard %.2f, opposite direction %d\n",
                fraction, comparison, nrow(r), P_VALUE_THRESHOLD, LOG2FC_THRESHOLD, sum(r$Change != "No Change"),
                sum(r$Change == "Down Regulated"), sum(r$Change == "Up Regulated"), sum(r$adj.P.Val < 0.05),
                a$Jaccard, a$opposite_direction))
  }
}

# Run the analysis only when this file is executed directly, not when its functions are sourced
# by tissue_level_permutation_check.R and scaffold_label_permutation_test.R.
if (sys.nframe() == 0) {
  for (fraction in FRACTIONS) run_differential_abundance(fraction)
}
