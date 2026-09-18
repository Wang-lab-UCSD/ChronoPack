
# cell 1

NOTEBOOK_DIR <- normalizePath(getwd())
REPO_ROOT <- normalizePath(file.path(NOTEBOOK_DIR, ".."))
source(file.path(REPO_ROOT, "chronopack", "R", "grn_inference_helpers.R"))


# cell 2

VALID_TARGETS_PATH <- "/new-stg/home/alex/tmp/grn_inference_v2/reproduction/inputs/external/ref/genome/hg38/annotations/protein_coding_genes.csv"
GENE_TABLES <- list(
  K562 = "/new-stg/home/alex/tmp/grn_inference_v2/reproduction/inputs/external/projects/chrono_seq/analyses/k562-tnfa.hg38.rep/intmd/sc-K562-1e4_to_tradeSeq/asso_df.1trajectory.lfdr.csv",
  THP1 = "/new-stg/home/alex/tmp/grn_inference_v2/reproduction/inputs/external/projects/chrono_seq/analyses/thp1-k562-lps.hg38.rep/intmd/sc-THP1_to_tradeSeq/asso_df.1trajectory.lfdr.csv"
)
OUTPUT_DIR <- file.path(NOTEBOOK_DIR, "outputs", "tables")
REFERENCE_DIR <- "/new-stg/home/alex/tmp/grn_inference_v2/reproduction/benchmark_outputs/tables"


# cell 3

result <- calculate_lfdr_recalibration(VALID_TARGETS_PATH, GENE_TABLES)
paths <- write_lfdr_recalibration_outputs(result, OUTPUT_DIR)
print(paths)


# cell 4

reference_paths <- list(
  scores = file.path(REFERENCE_DIR, "thp1_k562_lfdr_recalibration_scores.csv"),
  summary = file.path(REFERENCE_DIR, "thp1_k562_lfdr_recalibration_summary.csv"),
  params = file.path(REFERENCE_DIR, "thp1_k562_lfdr_recalibration_params.csv")
)
compare_csv <- function(label) {
  generated <- read.csv(paths[[label]], check.names = FALSE, stringsAsFactors = FALSE)
  reference <- read.csv(reference_paths[[label]], check.names = FALSE, stringsAsFactors = FALSE)
  same_shape <- identical(dim(generated), dim(reference))
  same_names <- identical(names(generated), names(reference))
  same_values <- isTRUE(all.equal(generated, reference, tolerance = 0, check.attributes = FALSE))
  data.frame(output = label, same_shape = same_shape, same_names = same_names, same_values = same_values)
}
do.call(rbind, lapply(names(paths), compare_csv))
