#!/usr/bin/env Rscript

# logsum.R
#
# For each per-fragment CpG methylation file, this script computes
# fragment-level log-likelihood sums against a reference methylation matrix.
#
# The script is designed for SLURM array execution:
#   - each array task processes one input sample
#   - each fragment is matched to reference CpG loci
#   - methylated / unmethylated observations are converted to log-probabilities
#   - log-probabilities are summed across CpGs within each fragment
#
# Input
# -----
# 1. Per-fragment CpG methylation calls:
#      *.perRead_cpg.tsv.gz
# 2. Reference methylation matrix:
#      rows   = CpG loci / probe IDs
#      cols   = reference profiles / tissues / components
#      values = methylation probabilities in [0, 1]
#
# Output
# ------
# One gzipped TSV per sample containing fragment-level metadata and
# summed log-probabilities for each reference column.
#
# Notes
# -----
# - CpG loci are mapped to reference IDs using genomic coordinates.
# - Fragments with fewer than 3 matched CpGs are discarded.
# - Reference values are clipped away from 0 and 1 to avoid log(0).
# - The script assumes execution in an HPC environment with
#   SLURM_ARRAY_TASK_ID set.

suppressMessages({
  library(data.table)
})

# ------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------

# Directory containing per-fragment CpG methylation files.
input_dir <- "/proj/nobackup/sens2024549/perRead_methylation/per_read_per_cpg/"

# Output directory for fragment-level log-probability summaries.
out_dir <- "/proj/nobackup/sens2024549/perRead_methylation/20260319_logsum_cfD"

# Reference methylation matrix.
ref_file <- "/proj/nobackup/sens2024549/perRead_methylation/20260319_logsum_cfD/dat_methy_merged_df.tsv.gz"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# Discover input files
# ------------------------------------------------------------------

perread_files <- list.files(
  input_dir,
  pattern = "\\.perRead_cpg\\.tsv\\.gz$",
  # pattern = ".*Curve.*\\.perRead_cpg\\.tsv\\.gz$",
  full.names = TRUE
)

# Optional exclusion of specific samples.
# perread_files <- perread_files[!grepl("Curve", perread_files)]

sample_names <- sub("\\.perRead_cpg\\.tsv\\.gz$", "", basename(perread_files))

# ------------------------------------------------------------------
# Resolve SLURM array task
# ------------------------------------------------------------------

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA))

if (is.na(task_id)) {
  stop("SLURM_ARRAY_TASK_ID is not set.")
}
if (task_id < 1L || task_id > length(perread_files)) {
  stop("SLURM_ARRAY_TASK_ID out of range: ", task_id)
}

infile <- perread_files[task_id]
sample_name <- sample_names[task_id]
outfile <- file.path(out_dir, paste0(sample_name, "_prob_sum.tsv.gz"))

cat("Task:", task_id, "/", length(perread_files), "\n")
cat("Sample:", sample_name, "\n")
cat("Input :", infile, "\n")
cat("Output:", outfile, "\n")

# Skip completed outputs to support restartable array execution.
if (file.exists(outfile) && file.info(outfile)$size > 0) {
  cat("Output exists, skipping.\n")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------
# Load reference methylation matrix
# ------------------------------------------------------------------

df_ref <- fread(ref_file)
setDT(df_ref)

# First column is assumed to contain CpG / probe IDs.
# Remaining columns are reference methylation probabilities.
ref_mat <- as.matrix(df_ref[, -1, with = FALSE])
rownames(ref_mat) <- df_ref$id

# Prevent numerical issues in downstream log-probability calculations:
#   log(0)     is undefined
#   log1p(-1)  is undefined
ref_mat[ref_mat <= 0] <- 1e-12
ref_mat[ref_mat >= 1] <- 1 - 1e-12

# Optional global methylation decay adjustment.
# This was kept commented to preserve original behavior.
# decay_rate <- 0.05
# ref_mat[, -1] <- ref_mat[, -1] * (1 - decay_rate)

# ------------------------------------------------------------------
# Load per-fragment CpG methylation calls
# ------------------------------------------------------------------

df_cpg <- fread(
  infile,
  select = c("fragment_id", "chrom", "ref_pos", "strand", "call", "fragment_span")
)

# Prefix fragment IDs with chromosome to reduce the chance of collisions
# if fragment IDs were assigned independently within chromosomes.
df_cpg[, fragment_id := paste(chrom, fragment_id)]

# Convert methylation state to binary encoding:
#   m -> 1
#   u -> 0
df_cpg[, call := as.integer(call == "m")]

# ------------------------------------------------------------------
# Map CpG observations to reference IDs
# ------------------------------------------------------------------

# The reference ID convention assumes:
#   - forward-strand CpG uses position = ref_pos
#   - reverse-strand CpG uses position = ref_pos - 1
#
# This effectively maps both strands of the same CpG dyad to a shared ID.
df_cpg[, id := paste(chrom, ref_pos + fifelse(strand == "-", -1L, 0L))]

# ------------------------------------------------------------------
# Filter to CpGs present in the reference matrix
# ------------------------------------------------------------------

ref_idx <- match(df_cpg$id, df_ref$id)
keep <- !is.na(ref_idx)

if (!any(keep)) {
  fwrite(data.table(), outfile, sep = "\t")
  quit(save = "no", status = 0)
}

df_cpg <- df_cpg[keep]
ref_idx <- ref_idx[keep]

# Compute fragment-level summary statistics before applying the minimum
# CpG-count filter.
df_cpg[, `:=`(
  n_probe = .N,
  n_probe_meth = sum(call == 1L)
), by = fragment_id]

# Retain only fragments supported by at least 3 matched CpGs.
df_cpg <- df_cpg[n_probe >= 3]

if (nrow(df_cpg) == 0L) {
  fwrite(data.table(), outfile, sep = "\t")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------
# Extract reference rows corresponding to retained CpGs
# ------------------------------------------------------------------

ref_idx <- match(df_cpg$id, rownames(ref_mat))
keep <- !is.na(ref_idx)

if (!all(keep)) {
  df_cpg <- df_cpg[keep]
  ref_idx <- ref_idx[keep]
}

if (nrow(df_cpg) == 0L) {
  fwrite(data.table(), outfile, sep = "\t")
  quit(save = "no", status = 0)
}

ref_sub <- ref_mat[ref_idx, , drop = FALSE]

# ------------------------------------------------------------------
# Compute per-CpG log-probabilities
# ------------------------------------------------------------------

# For methylated observations:
#   log P(call = methylated) = log(p)
#
# For unmethylated observations:
#   log P(call = unmethylated) = log(1 - p)
#
# where p is the reference methylation probability at that CpG.
log_prob <- log(ref_sub)

idx_unmeth <- df_cpg$call == 0L
log_prob[idx_unmeth, ] <- log1p(-ref_sub[idx_unmeth, , drop = FALSE])

# ------------------------------------------------------------------
# Sum log-probabilities across CpGs within each fragment
# ------------------------------------------------------------------

group <- df_cpg$fragment_id

# rowsum performs efficient grouped summation across all reference columns.
# reorder = FALSE preserves first-seen fragment order.
log_prob_sum <- rowsum(log_prob, group = group, reorder = FALSE)

# Identify the first row for each fragment so fragment metadata can be
# carried forward into the output table.
first_idx <- !duplicated(group)

out <- data.table(
  fragment_id   = df_cpg$fragment_id[first_idx],
  pos_id        = df_cpg$id[first_idx],
  fragment_span = df_cpg$fragment_span[first_idx],
  n_probe       = df_cpg$n_probe[first_idx],
  n_probe_meth  = df_cpg$n_probe_meth[first_idx]
)

# Append fragment-level summed log-probabilities.
out <- cbind(out, as.data.table(round(log_prob_sum, 5)))

# Preserve original output behavior by dropping fragment_id.
out[, fragment_id := NULL]

# ------------------------------------------------------------------
# Write output
# ------------------------------------------------------------------

fwrite(out, outfile, sep = "\t")

cat("Done.\n")
