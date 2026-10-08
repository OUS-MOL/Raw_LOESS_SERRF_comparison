# ============================================================
# Raw vs LOESS vs SERRF normalization comparison
# (untargeted LC-MS, Compound Discoverer 3.4 exports)
#
# Conditions compared:  Uncorrected | LOESS | SERRF
#   - LOESS: per-batch LOESS drift correction fitted on the PQC injections
#   - SERRF: per-batch random-forest drift correction fitted on the PQC
#            injections (Systematic Error Removal using Random Forest)
#   - Both are followed by PQC-anchored between-batch median scaling.
# QC reference materials (QC/PQC, Oslo-RM/LTQC, NIST) are used only for
# assessment (CV per feature) — only the PQC drives the correction.
#
# Usage (from the repository root):
#   Rscript normalization_comparison.R <compounds.xlsx> <sequence.xlsx> <output_dir>
# or edit the defaults below and source() the script.
#
# Outputs (<output_dir>/figures, <output_dir>/tables):
#   *_for_publication  — A: PCA (coloured by batch, 95% ellipses per batch)
#                        B: median log10 intensity vs injection order
#                        C: CV violins per QC type          (PDF + TIFF)
#   *_v1_violins, *_v2_cv_scatter, *_v3_batch_stability, *_v4_summary,
#   *_pca_labeled      — additional views of the same comparison
#   CV_summary.csv, CV_per_feature.csv, and corrected data tables per condition
#
# Internal sample-type codes are PQC / LTQC / NIST / Sample; figures show them
# as QC / Oslo-RM / NIST / Sample via `display_names` below.
# ============================================================

# ============================================================
# CONFIGURATION — modify these settings before running
# ============================================================

# Input file paths (overridable from the command line, see Usage)
cli_args       <- commandArgs(trailingOnly = TRUE)
cd_export_path <- if (length(cli_args) >= 1) cli_args[1] else "data/compounds.xlsx"
sequence_path  <- if (length(cli_args) >= 2) cli_args[2] else "data/sequence.xlsx"

# Output directory (created automatically if it doesn't exist)
output_dir     <- if (length(cli_args) >= 3) cli_args[3] else "output"

# Dataset label (used in file names and plot titles)
dataset_label    <- "Lipidomics_POS"

# Correction methods to run
run_loess  <- TRUE   # condition 2: per-batch LOESS + PQC between-batch scaling
run_serrf  <- TRUE   # condition 3: per-batch SERRF + PQC between-batch scaling

# LOESS parameters
loess_span_fixed       <- NULL  # NULL = adaptive per-batch CV; numeric = fixed for all
loess_span_sample_size <- 200   # features sampled per batch for span estimation

# SERRF parameters
serrf_num_corr_features <- 10
serrf_n_trees           <- 200
serrf_n_cores_ratio     <- 0.5

# CV thresholds for reporting (also shown as y-axis tick marks on violin plots)
cv_threshold_1 <- 0.15   # lower threshold (default 15%)
cv_threshold_2 <- 0.30   # upper threshold (default 30%)

# Colour palette (sample types) — used across injection order, CV, and PCA fill
# Keys stay as the internal sample_type values (PQC/LTQC/NIST); display text
# is remapped separately via `display_names` below (PQC -> "QC", LTQC -> "Oslo-RM").
sample_type_colors <- c(
  "Sample" = "#999999",
  "PQC"    = "#00A087",  # displayed as "QC" — green
  "LTQC"   = "#E64B35",  # displayed as "Oslo-RM" — red
  "NIST"   = "#4DBBD5"   # blue
)
# Batch colours are generated dynamically from the data (ggsci NPG palette)

# Display names for figure text — internal sample_type values are unchanged
# (so CSV exports / CV computations stay identical to prior runs)
display_names <- c("Sample" = "Sample", "PQC" = "QC", "LTQC" = "Oslo-RM", "NIST" = "NIST")

# ============================================================
# SECTION 1: PACKAGE INSTALLATION AND LOADING
# ============================================================

cat("Loading packages...\n")

cran_packages <- c(
  "readxl", "tidyverse", "dplyr", "tidyr", "purrr", "stringr", "lubridate",
  "ggplot2", "ggsci", "patchwork", "cowplot", "scales",
  "ranger", "parallel", "data.table", "ggrepel"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cran.r-project.org")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("pcaMethods", quietly = TRUE)) BiocManager::install("pcaMethods")
suppressPackageStartupMessages(library(pcaMethods))

cat("All packages loaded.\n\n")

# ============================================================
# SECTION 2: DATA IMPORT AND METADATA CONSTRUCTION
# ============================================================

cat("======== SECTION 2: DATA IMPORT ========\n")

dir.create(file.path(output_dir, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

# ---- 2.1 Load sequence file ----

cat("Loading sequence file...\n")

seq_raw <- readxl::read_xlsx(sequence_path, sheet = "InputFiles")

batch_col_name <- grep("^Batch(\\.\\.\\.\\d+)?$", colnames(seq_raw), value = TRUE)[1]
if (is.na(batch_col_name)) stop("Could not find a 'Batch' column in the sequence file.")
cat(sprintf("  Using batch column: '%s'\n", batch_col_name))

seq_raw$creation_date <- lubridate::dmy_hms(seq_raw[["Creation Date"]])

seq_filtered <- seq_raw %>%
  dplyr::filter(as.numeric(.data[[batch_col_name]]) >= 0, `Sample Type` != "Blank")

seq_filtered$file_id <- as.integer(gsub("^F", "", seq_filtered[["Study File ID"]]))

seq_filtered$sample_type <- dplyr::case_when(
  seq_filtered[["Sample Type"]] == "Quality Control" ~ "PQC",
  seq_filtered[["Group"]] == "LTQC"                  ~ "LTQC",
  seq_filtered[["Group"]] == "NIST"                  ~ "NIST",
  TRUE                                                ~ "Sample"
)
seq_filtered$sample_type <- factor(seq_filtered$sample_type,
                                   levels = c("Sample", "PQC", "LTQC", "NIST"))

seq_filtered <- seq_filtered %>%
  dplyr::arrange(creation_date) %>%
  dplyr::mutate(injection_order = seq_len(dplyr::n()))

metadata <- seq_filtered %>%
  dplyr::transmute(
    file_id         = file_id,
    file_name_base  = basename(`File Name`),
    creation_date   = creation_date,
    batch           = as.integer(.data[[batch_col_name]]),
    sample_type     = sample_type,
    injection_order = injection_order
  )

cat(sprintf("Metadata loaded: %d samples after filtering.\n", nrow(metadata)))

# ---- 2.2 Load CD export ----

cat("Loading CD export...\n")

cd_raw <- readxl::read_xlsx(cd_export_path, sheet = "Compounds")

feature_meta_cols_present <- intersect(c("Name", "Formula", "RT [min]", "m/z", "Calc. MW"),
                                        colnames(cd_raw))
feature_meta <- cd_raw[, feature_meta_cols_present, drop = FALSE]

area_col_mask <- grepl("^Area:.*\\.raw", colnames(cd_raw)) &
                 !grepl("^Group Area:", colnames(cd_raw))
area_cols     <- colnames(cd_raw)[area_col_mask]

if (length(area_cols) == 0)
  stop("No area columns found matching 'Area: *.raw'. Check CD export column names.")

extract_file_id <- function(col_name) {
  m <- regmatches(col_name, regexpr("\\(F(\\d+)\\)", col_name))
  if (length(m) == 0) return(NA_integer_)
  as.integer(sub("\\(F(\\d+)\\)", "\\1", m))
}

area_file_ids   <- sapply(area_cols, extract_file_id)
area_matrix_raw <- as.matrix(cd_raw[, area_cols, drop = FALSE])
colnames(area_matrix_raw) <- as.character(area_file_ids)

cat(sprintf("CD export loaded: %d features, %d area columns.\n",
            nrow(area_matrix_raw), ncol(area_matrix_raw)))

# ---- 2.3 Join and order ----

cd_ids_in_meta  <- which(colnames(area_matrix_raw) %in% as.character(metadata$file_id))
cd_ids_not_meta <- which(!(colnames(area_matrix_raw) %in% as.character(metadata$file_id)))

if (length(cd_ids_not_meta) > 0)
  cat(sprintf("WARNING: %d CD area columns unmatched (dropped): IDs = %s\n",
              length(cd_ids_not_meta),
              paste(colnames(area_matrix_raw)[cd_ids_not_meta], collapse = ", ")))

area_matrix_joined <- area_matrix_raw[, cd_ids_in_meta, drop = FALSE]

meta_ordered <- metadata[order(metadata$injection_order), ]
common_ids   <- intersect(as.character(meta_ordered$file_id), colnames(area_matrix_joined))
meta_ordered <- meta_ordered[meta_ordered$file_id %in% as.integer(common_ids), ]
area_matrix  <- area_matrix_joined[, as.character(meta_ordered$file_id), drop = FALSE]
metadata     <- meta_ordered

cat(sprintf("After join: %d features x %d samples.\n", nrow(area_matrix), ncol(area_matrix)))

# ---- 2.4 Duplicate feature name handling ----

if (!"Name" %in% colnames(feature_meta)) feature_meta$Name <- NA_character_
feature_meta$Name_original <- feature_meta$Name

unknown_mask <- is.na(feature_meta$Name) | trimws(feature_meta$Name) == ""
if (any(unknown_mask)) {
  rt_col <- intersect(c("RT [min]", "RT"), colnames(feature_meta))[1]
  mz_col <- intersect(c("m/z", "mz"),      colnames(feature_meta))[1]

  feature_meta$Name[unknown_mask] <- vapply(which(unknown_mask), function(i) {
    rt_str <- if (!is.na(rt_col) && !is.na(feature_meta[[rt_col]][i]))
      sprintf("RT%.2f", as.numeric(feature_meta[[rt_col]][i])) else sprintf("row%d", i)
    mz_str <- if (!is.na(mz_col) && !is.na(feature_meta[[mz_col]][i]))
      sprintf("_mz%.4f", as.numeric(feature_meta[[mz_col]][i])) else ""
    paste0("Unk_", rt_str, mz_str)
  }, character(1))
  cat(sprintf("  %d unnamed features labelled as Unk_RT<rt>_mz<mz>.\n", sum(unknown_mask)))
}

name_counts <- table(feature_meta$Name)
for (dn in names(name_counts[name_counts > 1])) {
  idx <- which(feature_meta$Name == dn)
  feature_meta$Name[idx] <- paste0(dn, "(", letters[seq_along(idx)], ")")
}
rownames(area_matrix) <- feature_meta$Name

# ---- 2.5 Data validation printout ----

ltqc_counts <- table(metadata$batch[metadata$sample_type == "LTQC"])
cat(sprintf(
  "\nDataset: %s\nFeatures: %d\nTotal samples (incl. QCs): %d\n  - Biological samples: %d\n  - PQC: %d (across %d batches)\n  - LTQC: %d (across %d batches, %s per batch)\n  - NIST: %d (across %d batches)\nBatches: %s\nInjection order range: %s to %s\n\n",
  dataset_label, nrow(feature_meta), nrow(metadata),
  sum(metadata$sample_type == "Sample"),
  sum(metadata$sample_type == "PQC"),
  length(unique(metadata$batch[metadata$sample_type == "PQC"])),
  sum(metadata$sample_type == "LTQC"), length(ltqc_counts),
  paste(ltqc_counts, collapse = "/"),
  sum(metadata$sample_type == "NIST"),
  length(unique(metadata$batch[metadata$sample_type == "NIST"])),
  paste(sort(unique(metadata$batch)), collapse = ", "),
  format(min(metadata$creation_date), "%Y-%m-%d %H:%M"),
  format(max(metadata$creation_date), "%Y-%m-%d %H:%M")
))

# ============================================================
# SECTION 3: LOG TRANSFORMATION
# ============================================================

cat("======== SECTION 3: LOG TRANSFORMATION ========\n")

area_matrix <- apply(area_matrix, 2, as.numeric)
rownames(area_matrix) <- feature_meta$Name

n_replaced <- sum(area_matrix <= 0, na.rm = TRUE)
area_matrix[area_matrix <= 0] <- NA
cat(sprintf("Replaced %d zero/negative values with NA.\n", n_replaced))

log2_matrix <- log2(area_matrix)
cat("Log2 transformation complete.\n\n")

# ============================================================
# SECTION 4: LOESS CORRECTION FUNCTION
# ============================================================

cat("======== SECTION 4: LOESS CORRECTION FUNCTION ========\n")

select_loess_span <- function(x, y, span_vals = seq(0.25, 1, by = 0.05), folds = 5) {
  valid <- !is.na(y)
  x_v <- x[valid]; y_v <- y[valid]
  if (length(x_v) < 4) return(0.75)
  set.seed(42)
  fold_ids <- sample(rep(seq_len(folds), length.out = length(x_v)))
  mae_per_span <- vapply(span_vals, function(sp) {
    mean(vapply(seq_len(folds), function(f) {
      tryCatch({
        fit  <- loess(y_v[fold_ids != f] ~ x_v[fold_ids != f], span = sp,
                      control = loess.control(surface = "direct"))
        pred <- predict(fit, newdata = data.frame(
          `x_v[fold_ids != f]` = x_v[fold_ids == f]))
        mean(abs(y_v[fold_ids == f] - pred), na.rm = TRUE)
      }, error = function(e) NA_real_)
    }, numeric(1)), na.rm = TRUE)
  }, numeric(1))
  if (all(is.na(mae_per_span))) return(0.75)
  span_vals[which.min(mae_per_span)]
}

apply_loess_correction <- function(log2_mat, metadata, per_batch = TRUE,
                                   fixed_span = NULL,
                                   span_sample_size = loess_span_sample_size) {
  corrected  <- log2_mat
  n_features <- nrow(log2_mat)
  if (!per_batch) metadata$batch <- 0L
  batches <- sort(unique(metadata$batch))

  for (b in batches) {
    t_batch   <- proc.time()["elapsed"]
    batch_idx <- which(metadata$batch == b)
    pqc_idx   <- which(metadata$batch == b & metadata$sample_type == "PQC")

    if (length(pqc_idx) < 3) {
      cat(sprintf("  Batch %s: fewer than 3 PQCs — skipping LOESS.\n", b)); next
    }

    if (per_batch) {
      rank_map    <- rank(metadata$injection_order[batch_idx])
      names(rank_map) <- as.character(metadata$injection_order[batch_idx])
      all_x_local <- rank_map
      pqc_x_local <- rank_map[as.character(metadata$injection_order[pqc_idx])]
    } else {
      all_x_local <- metadata$injection_order[batch_idx]
      pqc_x_local <- metadata$injection_order[pqc_idx]
    }

    if (!is.null(fixed_span)) {
      batch_span <- fixed_span
    } else {
      cands <- which(rowSums(!is.na(log2_mat[, pqc_idx, drop = FALSE])) >= 3)
      s_idx <- sample(cands, min(span_sample_size, length(cands)))
      cat(sprintf("  Batch %s: estimating span from %d features...\n", b, length(s_idx)))
      spans      <- vapply(s_idx, function(f)
        select_loess_span(pqc_x_local, log2_mat[f, pqc_idx]), numeric(1))
      batch_span <- median(spans, na.rm = TRUE)
      if (is.na(batch_span)) batch_span <- 0.75
      cat(sprintf("  Batch %s: span = %.2f\n", b, batch_span))
    }

    for (f in seq_len(n_features)) {
      tryCatch({
        y_pqc <- log2_mat[f, pqc_idx]
        if (sum(!is.na(y_pqc)) < 3) next
        fit  <- loess(y_pqc ~ pqc_x_local, span = batch_span,
                      control = loess.control(surface = "direct"))
        pred <- predict(fit, newdata = data.frame(pqc_x_local = all_x_local))
        if (any(is.na(pred))) {
          nn <- which(!is.na(pred))
          if (!length(nn)) next
          for (i in which(is.na(pred)))
            pred[i] <- pred[nn[which.min(abs(nn - i))]]
        }
        ref_med <- median(log2_mat[f, pqc_idx], na.rm = TRUE)
        corrected[f, batch_idx] <- log2_mat[f, batch_idx] - (pred - ref_med)
      }, error = function(e) NULL)
    }
    cat(sprintf("  LOESS batch %s done in %.1f s.\n", b,
                round(proc.time()["elapsed"] - t_batch, 1)))
  }
  corrected
}

# ============================================================
# SECTION 5: SERRF CORRECTION FUNCTION
# ============================================================

cat("======== SECTION 5: SERRF CORRECTION FUNCTION ========\n")

#' Core SERRF computation — explicit arguments, no closure/ambient state
serrf_compute_one <- function(f, pred_idx, pqc_b, batch_b,
                               pqc_inj, batch_inj, n_tr) {
  tryCatch({
    y_pqc <- pqc_b[f, ]
    if (sum(!is.na(y_pqc)) < 5) return(batch_b[f, ])
    train_df <- as.data.frame(t(pqc_b[pred_idx, , drop = FALSE]))
    colnames(train_df) <- paste0("V", seq_len(ncol(train_df)))
    train_df$inj_order <- pqc_inj
    train_df$y         <- y_pqc
    ok <- complete.cases(train_df)
    if (sum(ok) < 5) return(batch_b[f, ])
    train_df <- train_df[ok, ]
    rf <- ranger::ranger(y ~ ., data = train_df,
                         num.trees = n_tr, seed = 42, verbose = FALSE)
    all_df <- as.data.frame(t(batch_b[pred_idx, , drop = FALSE]))
    colnames(all_df) <- paste0("V", seq_len(ncol(all_df)))
    all_df$inj_order <- batch_inj
    pred    <- predict(rf, data = all_df)$predictions
    ref_med <- median(pqc_b[f, ], na.rm = TRUE)
    batch_b[f, ] - (pred - ref_med)
  }, error = function(e) batch_b[f, ])
}

#' Parallel entry point: reads .sw_* worker globals
serrf_worker_parallel <- function(item) {
  serrf_compute_one(item$f, item$pred_idx,
                    .sw_pqc, .sw_batch, .sw_pqc_inj, .sw_batch_inj, .sw_n_trees)
}

apply_serrf_correction <- function(log2_mat, metadata, per_batch = TRUE,
                                   num_features = 10, n_trees = 200,
                                   n_cores_ratio = 0.5) {
  set.seed(42)
  if (!per_batch) metadata$batch <- 0L
  batches    <- sort(unique(metadata$batch))
  n_features <- nrow(log2_mat)
  corrected  <- log2_mat

  n_cores <- max(1L, floor(parallel::detectCores(logical = FALSE) * n_cores_ratio))
  cat(sprintf("  SERRF using %d core(s), %d trees.\n", n_cores, n_trees))

  cl <- if (n_cores > 1) {
    cl_tmp <- parallel::makeCluster(n_cores)
    parallel::clusterEvalQ(cl_tmp, library(ranger))
    parallel::clusterExport(cl_tmp,
      varlist = c("serrf_compute_one", "serrf_worker_parallel"), envir = .GlobalEnv)
    cl_tmp
  } else NULL
  on.exit(try(if (!is.null(cl)) parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  for (b in batches) {
    t_batch   <- proc.time()["elapsed"]
    batch_idx <- which(metadata$batch == b)
    pqc_idx   <- which(metadata$batch == b & metadata$sample_type == "PQC")

    if (length(pqc_idx) < 5) {
      cat(sprintf("  Batch %s: <5 PQCs — median scaling fallback.\n", b))
      gm <- apply(log2_mat, 1, median, na.rm = TRUE)
      bm <- apply(log2_mat[, pqc_idx, drop = FALSE], 1, median, na.rm = TRUE)
      corrected[, batch_idx] <- log2_mat[, batch_idx] + (gm - bm); next
    }

    inj_batch <- metadata$injection_order[batch_idx]
    inj_pqc   <- metadata$injection_order[pqc_idx]
    pqc_mat_b <- log2_mat[, pqc_idx,   drop = FALSE]
    log2_b    <- log2_mat[, batch_idx, drop = FALSE]

    cat(sprintf("  Batch %s: computing correlations (%d × %d)...\n",
                b, n_features, n_features))
    t_cor   <- proc.time()["elapsed"]
    pqc_imp <- pqc_mat_b
    rmeans  <- rowMeans(pqc_imp, na.rm = TRUE)
    for (j in seq_len(ncol(pqc_imp))) {
      na_j <- is.na(pqc_imp[, j])
      pqc_imp[na_j, j] <- rmeans[na_j]
    }
    cor_mat <- cor(t(pqc_imp))
    diag(cor_mat) <- 0  # zero diagonal before abs() ordering

    # Build top-correlated-feature index matrix.
    # For each feature fi, cor_f[fi] is set to -Inf so that feature fi
    # cannot rank as its own predictor (self-prediction/overfitting prevention).
    # Using abs() on the remaining values selects strongest absolute correlation.
    k_sel   <- min(num_features, n_features - 1L)
    top_idx <- matrix(0L, nrow = n_features, ncol = k_sel)
    for (fi in seq_len(n_features)) {
      cor_f       <- abs(cor_mat[fi, ])
      cor_f[fi]   <- -Inf  # explicit self-exclusion: -Inf < 0 guarantees last rank
      top_idx[fi, ] <- order(cor_f, decreasing = TRUE)[seq_len(k_sel)]
    }
    rm(cor_mat, pqc_imp); gc(verbose = FALSE)
    cat(sprintf("  Batch %s: correlations done in %.1f s.\n", b,
                round(proc.time()["elapsed"] - t_cor, 1)))

    items <- lapply(seq_len(n_features), function(f)
      list(f = f, pred_idx = top_idx[f, ]))

    if (!is.null(cl)) {
      .sw_pqc <- pqc_mat_b; .sw_batch <- log2_b
      .sw_pqc_inj <- inj_pqc; .sw_batch_inj <- inj_batch; .sw_n_trees <- n_trees
      parallel::clusterExport(cl,
        varlist = c(".sw_pqc", ".sw_batch", ".sw_pqc_inj", ".sw_batch_inj", ".sw_n_trees"),
        envir = environment())
      rm(.sw_pqc, .sw_batch, .sw_pqc_inj, .sw_batch_inj, .sw_n_trees)
    }

    results <- if (!is.null(cl))
      parallel::parLapply(cl, items, serrf_worker_parallel)
    else
      lapply(items, function(item)
        serrf_compute_one(item$f, item$pred_idx,
                          pqc_mat_b, log2_b, inj_pqc, inj_batch, n_trees))

    for (f in seq_len(n_features)) corrected[f, batch_idx] <- results[[f]]
    cat(sprintf("  SERRF batch %s done in %.1f s.\n", b,
                round(proc.time()["elapsed"] - t_batch, 1)))
  }

  # Part B — SERRF diagnostics
  cat("\n--- SERRF diagnostics ---\n")
  cat(sprintf("Features processed:    %d\n", n_features))
  cat(sprintf("Batches processed:     %d\n", length(batches)))
  cat(sprintf("RF trees used:         %d\n", n_trees))
  cat(sprintf("Corr. features used:   %d\n", num_features))
  cat(sprintf("Self-exclusion active: YES (cor_f[fi] <- -Inf before abs() ordering)\n"))
  cat(sprintf("Non-NA corrected:      %d / %d (%.1f%%)\n",
              sum(!is.na(corrected)), length(corrected),
              100 * mean(!is.na(corrected))))

  corrected
}

# ============================================================
# SECTION 6: PQC BETWEEN-BATCH SCALING FUNCTION
# ============================================================

cat("======== SECTION 6: PQC BATCH SCALING FUNCTION ========\n")

#' Apply PQC-based between-batch median scaling
#'
#' Shifts each batch so its PQC median equals the grand PQC median across all batches.
#' Equivalent to the between-batch step performed by CD 3.4 SERRF.
#'
#' @param log2_mat Features x samples log2 matrix (within-batch corrected)
#' @param metadata Data frame with batch, sample_type
#' @param min_pqc Minimum PQC injections required per batch (default 3)
#' @return PQC-scaled log2 matrix
apply_pqc_batch_scaling <- function(log2_mat, metadata, min_pqc = 3) {
  corrected   <- log2_mat
  batches     <- sort(unique(metadata$batch))
  n_features  <- nrow(log2_mat)
  all_pqc_idx <- which(metadata$sample_type == "PQC")

  if (length(all_pqc_idx) == 0) {
    warning("No PQC samples — skipping PQC batch scaling."); return(corrected)
  }

  for (f in seq_len(n_features)) {
    grand_median <- median(log2_mat[f, all_pqc_idx], na.rm = TRUE)
    for (b in batches) {
      batch_idx <- which(metadata$batch == b)
      pqc_b_idx <- which(metadata$batch == b & metadata$sample_type == "PQC")
      if (length(pqc_b_idx) < min_pqc) {
        cat(sprintf("  WARNING: Batch %s has %d PQC samples (min = %d) — skipping.\n",
                    b, length(pqc_b_idx), min_pqc)); next
      }
      batch_median <- median(log2_mat[f, pqc_b_idx], na.rm = TRUE)
      corrected[f, batch_idx] <- log2_mat[f, batch_idx] + (grand_median - batch_median)
    }
  }
  corrected
}

# ============================================================
# SECTION 7: APPLY THREE CONDITIONS
# ============================================================

cat("======== SECTION 7: APPLYING CORRECTIONS ========\n")

corrected_list <- list()

# 1. Uncorrected
corrected_list[["Uncorrected"]] <- log2_matrix
cat("Condition 1 (Uncorrected): stored.\n")

# 2. LOESS: per-batch drift correction + PQC between-batch scaling
if (run_loess) {
  cat("Running LOESS correction (per-batch + PQC scaling)...\n")
  loess_pb <- apply_loess_correction(log2_matrix, metadata,
                                     per_batch  = TRUE,
                                     fixed_span = loess_span_fixed)
  corrected_list[["LOESS"]] <- apply_pqc_batch_scaling(loess_pb, metadata)
  cat("Condition 2 (LOESS): done.\n")
}

# 3. SERRF: per-batch drift correction + PQC between-batch scaling (CD 3.4 equivalent)
if (run_serrf) {
  cat("Running SERRF correction (per-batch + PQC scaling) — may take several minutes...\n")
  serrf_pb <- apply_serrf_correction(log2_matrix, metadata,
                                     per_batch     = TRUE,
                                     num_features  = serrf_num_corr_features,
                                     n_trees       = serrf_n_trees,
                                     n_cores_ratio = serrf_n_cores_ratio)
  corrected_list[["SERRF"]] <- apply_pqc_batch_scaling(serrf_pb, metadata)
  cat("Condition 3 (SERRF): done.\n")
}

cat("\nAll conditions applied.\n\n")

# ============================================================
# SECTION 8: EXPORT CORRECTED DATA TABLES
# ============================================================

cat("======== SECTION 8: EXPORTING TABLES ========\n")

for (cond in names(corrected_list)) {
  raw_mat   <- 2^corrected_list[[cond]]
  meta_cols <- feature_meta[, intersect(c("Name", "Formula", "RT [min]", "m/z", "Calc. MW"),
                                         colnames(feature_meta)), drop = FALSE]
  export_wide <- cbind(meta_cols, as.data.frame(raw_mat))
  colnames(export_wide)[seq(ncol(meta_cols) + 1, ncol(export_wide))] <- metadata$file_name_base

  cond_safe <- gsub("\\+", "plus", gsub("[^A-Za-z0-9_]", "_", cond))
  data.table::fwrite(export_wide,
    file.path(output_dir, "tables", paste0(dataset_label, "_", cond_safe, "_corrected.csv")))

  long_df <- as.data.frame(raw_mat) %>%
    tibble::rownames_to_column("feature_name") %>%
    tidyr::pivot_longer(-feature_name, names_to = "file_id_char",
                        values_to = "corrected_area") %>%
    dplyr::left_join(
      metadata %>% dplyr::mutate(file_id_char = as.character(file_id)),
      by = "file_id_char") %>%
    dplyr::select(feature_name, file_id, file_name_base, batch, sample_type,
                  injection_order, corrected_area)
  data.table::fwrite(long_df,
    file.path(output_dir, "tables", paste0(dataset_label, "_", cond_safe, "_long.csv")))
  cat(sprintf("  Exported: %s\n", cond_safe))
}
cat("\n")

# ============================================================
# SECTION 9: CV CALCULATION
# ============================================================

cat("======== SECTION 9: CV CALCULATION ========\n")

compute_cv <- function(log2_mat, metadata, qc_type) {
  idx <- which(metadata$sample_type == qc_type)
  if (length(idx) < 2) return(rep(NA_real_, nrow(log2_mat)))
  apply(2^log2_mat[, idx, drop = FALSE], 1, function(x) {
    x <- x[!is.na(x)]; if (length(x) < 2) NA_real_ else sd(x) / mean(x)
  })
}

qc_types   <- c("PQC", "LTQC", "NIST")
cv_results <- list()
for (cond in names(corrected_list)) {
  cv_results[[cond]] <- setNames(
    lapply(qc_types, function(qc) compute_cv(corrected_list[[cond]], metadata, qc)),
    qc_types)
}

cv_summary <- do.call(rbind, lapply(names(corrected_list), function(cond)
  do.call(rbind, lapply(qc_types, function(qc) {
    cvs <- cv_results[[cond]][[qc]]; cvs <- cvs[!is.na(cvs)]
    data.frame(condition = cond, qc_type = qc,
               median_cv    = if (length(cvs)) median(cvs) else NA_real_,
               pct_below_20 = if (length(cvs)) mean(cvs < cv_threshold_1) * 100 else NA_real_,
               pct_below_30 = if (length(cvs)) mean(cvs < cv_threshold_2) * 100 else NA_real_,
               n_features   = length(cvs), stringsAsFactors = FALSE)
  }))))
rownames(cv_summary) <- NULL

data.table::fwrite(cv_summary,
  file.path(output_dir, "tables", paste0(dataset_label, "_CV_summary.csv")))
cv_per_feature <- do.call(rbind, lapply(names(corrected_list), function(cond)
  do.call(rbind, lapply(qc_types, function(qc)
    data.frame(feature = feature_meta$Name, condition = cond, qc_type = qc,
               cv = cv_results[[cond]][[qc]], stringsAsFactors = FALSE)))))
data.table::fwrite(cv_per_feature,
  file.path(output_dir, "tables", paste0(dataset_label, "_CV_per_feature.csv")))
cat("CV tables exported.\n\n")

# ============================================================
# SECTION 10: PCA CALCULATION
# ============================================================

cat("======== SECTION 10: PCA CALCULATION ========\n")

compute_pca <- function(log2_mat, metadata) {
  mat_t     <- t(log2_mat)
  col_means <- colMeans(mat_t, na.rm = TRUE)
  for (j in seq_len(ncol(mat_t)))
    mat_t[is.na(mat_t[, j]), j] <- col_means[j]
  sds   <- apply(mat_t, 2, sd, na.rm = TRUE)
  mat_t <- mat_t[, sds > 0, drop = FALSE]
  if (ncol(mat_t) < 2) { warning("PCA skipped."); return(NULL) }
  mat_scaled <- apply(mat_t, 2, function(x) {
    m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
    if (s == 0) x - m else (x - m) / sqrt(s)
  })
  pca_res <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)
  var_exp <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
  list(
    scores  = data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
                          batch = as.character(metadata$batch),
                          sample_type = as.character(metadata$sample_type),
                          injection_order = metadata$injection_order,
                          sample_name = metadata$file_name_base,
                          stringsAsFactors = FALSE),
    var_exp = var_exp,
    pc1_var = round(var_exp[1], 1),
    pc2_var = round(var_exp[2], 1)
  )
}

pca_results <- lapply(names(corrected_list), function(cond) {
  r <- compute_pca(corrected_list[[cond]], metadata)
  cat(sprintf("  PCA: %s\n", cond)); r
})
names(pca_results) <- names(corrected_list)
cat("\n")

# ============================================================
# SECTION 11: VISUALISATION
# ============================================================

cat("======== SECTION 11: VISUALISATION ========\n")

batch_ids    <- sort(unique(metadata$batch))
batch_colors <- setNames(ggsci::pal_npg("nrc")(min(length(batch_ids), 10)),
                          as.character(batch_ids))

# Oslo-RM (LTQC) = solid square (15); NIST = filled inverted triangle (25,
# requires colour/border + fill/interior); QC (PQC) unchanged triangle (17).
sample_type_shapes <- c("Sample" = 16, "PQC" = 17, "LTQC" = 15, "NIST" = 25)
sample_type_sizes  <- c("Sample" = 2,  "PQC" = 3,  "LTQC" = 3,  "NIST" = 3)

batch_ends <- if (length(batch_ids) > 1) {
  vapply(batch_ids[-length(batch_ids)], function(b)
    max(metadata$injection_order[metadata$batch == b]), numeric(1))
} else numeric(0)

conditions <- names(corrected_list)

# ---- A4: Shared y-axis for injection order plots (log10 units, matches panel B) ----
all_median_intensities <- lapply(conditions, function(cond)
  apply(corrected_list[[cond]], 2, median, na.rm = TRUE) * log10(2))
inj_y_min <- min(unlist(all_median_intensities), na.rm = TRUE)
inj_y_max <- max(unlist(all_median_intensities), na.rm = TRUE)
inj_y_pad <- (inj_y_max - inj_y_min) * 0.05

# ---- A3: Shared y-cap for CV panels ----
all_cv_vals <- unlist(lapply(conditions, function(cond)
  lapply(c("PQC", "LTQC", "NIST"), function(qc) cv_results[[cond]][[qc]])))
all_cv_vals <- all_cv_vals[!is.na(all_cv_vals)]
shared_cv_cap <- min(quantile(all_cv_vals, 0.99, na.rm = TRUE), 1.5)
shared_cv_cap <- max(shared_cv_cap, cv_threshold_2 * 2)

# ---- 11.2 Injection order plots (A1, A4, A5) ----
# LTQC is a passive observer — no connecting line in any condition.
make_injection_plot <- function(log2_mat, metadata, cond_name,
                                 y_lim = NULL) {
  plot_df <- metadata %>%
    dplyr::mutate(median_log10 = colMeans(log2_mat, na.rm = TRUE) * log10(2),
                  sample_type = as.character(sample_type))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = injection_order, y = median_log10)) +
    ggplot2::geom_point(ggplot2::aes(colour = sample_type, fill = sample_type,
                                      shape = sample_type, size = sample_type)) +
    ggplot2::scale_colour_manual(values = sample_type_colors, name = "Sample Type",
                                  labels = display_names) +
    ggplot2::scale_fill_manual(values   = sample_type_colors, name = "Sample Type",
                                labels = display_names) +
    ggplot2::scale_shape_manual(values  = sample_type_shapes, name = "Sample Type",
                                 labels = display_names) +
    ggplot2::scale_size_manual(values   = sample_type_sizes,  name = "Sample Type",
                                labels = display_names) +
    ggplot2::labs(title = cond_name, x = "Injection Order",
                  y = expression("Median log"[10]*" Intensity")) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(size = 9, face = "bold"))

  for (be in batch_ends)
    p <- p + ggplot2::geom_vline(xintercept = be + 0.5, linetype = "dashed",
                                  colour = "grey40", linewidth = 0.4)

  if (!is.null(y_lim))
    p <- p + ggplot2::coord_cartesian(ylim = y_lim)
  p
}

# ---- 11.3 PCA plots (A1, A5) ----
# show_legend = TRUE used for extracting the shared legend only.
make_pca_plot <- function(pca_res, cond_name, show_legend = FALSE) {
  if (is.null(pca_res))
    return(ggplot2::ggplot() +
             ggplot2::labs(title = paste(cond_name, "(PCA unavailable)")) +
             ggplot2::theme_void())

  scores_df <- pca_res$scores
  scores_df$sample_type <- factor(scores_df$sample_type,
                                   levels = c("Sample", "PQC", "LTQC", "NIST"))
  sample_only <- scores_df[scores_df$sample_type == "Sample", ]

  pca_size_map <- c("Sample" = 2, "PQC" = 3, "LTQC" = 3.5, "NIST" = 3.5)

  # PCA encoding: colour = batch, fill = batch (makes LTQC ▽ fully batch-coloured,
  # no separate green border). Shape encodes sample type. Fill legend suppressed
  # to avoid duplicating the batch legend.
  p <- ggplot2::ggplot(scores_df, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(ggplot2::aes(colour = batch, fill = batch,
                                      shape = sample_type, size = sample_type)) +
    ggplot2::scale_colour_manual(values = batch_colors, name = "Batch") +
    ggplot2::scale_fill_manual(values   = batch_colors, guide = "none") +
    ggplot2::scale_shape_manual(values  = sample_type_shapes, name = "Sample Type",
                                 labels = display_names,
                                 guide = ggplot2::guide_legend(
                                   override.aes = list(colour = "#333333",
                                                       fill   = "#333333",
                                                       size   = 3))) +
    ggplot2::scale_size_manual(values   = pca_size_map, guide = "none") +
    ggplot2::labs(title = cond_name,
                  x = sprintf("PC1 (%.1f%%)", pca_res$pc1_var),
                  y = sprintf("PC2 (%.1f%%)", pca_res$pc2_var)) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none",
                   plot.title = ggplot2::element_text(size = 9, face = "bold"))

  if (nrow(sample_only) > 0)
    p <- p + ggplot2::stat_ellipse(data = sample_only,
                                    ggplot2::aes(colour = batch),
                                    type = "norm", level = 0.95,
                                    linewidth = 0.4, linetype = "dashed",
                                    show.legend = FALSE)
  p
}

# ---- 11.3b Labeled PCA: sample names via ggrepel ----
make_pca_labeled_plot <- function(pca_res, cond_name) {
  if (is.null(pca_res))
    return(ggplot2::ggplot() +
             ggplot2::labs(title = paste(cond_name, "(PCA unavailable)")) +
             ggplot2::theme_void())

  scores_df <- pca_res$scores
  scores_df$sample_type <- factor(scores_df$sample_type,
                                   levels = c("Sample", "PQC", "LTQC", "NIST"))
  sample_only <- scores_df[scores_df$sample_type == "Sample", ]

  pca_size_map <- c("Sample" = 2, "PQC" = 3, "LTQC" = 3.5, "NIST" = 3.5)

  p <- ggplot2::ggplot(scores_df, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(ggplot2::aes(colour = batch, fill = batch,
                                      shape = sample_type, size = sample_type)) +
    ggrepel::geom_text_repel(ggplot2::aes(label = sample_name, colour = batch),
                              size = 1.8, max.overlaps = 50,
                              show.legend = FALSE,
                              segment.size = 0.2, segment.colour = "grey60") +
    ggplot2::scale_colour_manual(values = batch_colors, name = "Batch") +
    ggplot2::scale_fill_manual(values   = batch_colors, guide = "none") +
    ggplot2::scale_shape_manual(values  = sample_type_shapes, name = "Sample Type",
                                 labels = display_names,
                                 guide = ggplot2::guide_legend(
                                   override.aes = list(colour = "#333333",
                                                       fill   = "#333333",
                                                       size   = 3))) +
    ggplot2::scale_size_manual(values   = pca_size_map, guide = "none") +
    ggplot2::labs(title = cond_name,
                  x = sprintf("PC1 (%.1f%%)", pca_res$pc1_var),
                  y = sprintf("PC2 (%.1f%%)", pca_res$pc2_var)) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(size = 9, face = "bold"))

  if (nrow(sample_only) > 0)
    p <- p + ggplot2::stat_ellipse(data = sample_only,
                                    ggplot2::aes(colour = batch),
                                    type = "norm", level = 0.95,
                                    linewidth = 0.4, linetype = "dashed",
                                    show.legend = FALSE)
  p
}

# ---- 11.4 CV plots: one panel per correction condition, x = QC type ----
# Each panel shows PQC / LTQC / NIST violins within a single condition.
# Fill colour = sample_type_colors, matching injection order plots.
# Shared y-cap across all panels (computed in A3 block above).
make_cv_condition_plot <- function(cv_results, condition,
                                    qc_types = c("PQC", "LTQC", "NIST"),
                                    y_cap    = shared_cv_cap) {
  cv_df <- do.call(rbind, lapply(qc_types, function(qc)
    data.frame(qc_type = qc, cv = cv_results[[condition]][[qc]],
               stringsAsFactors = FALSE)))
  cv_df <- cv_df[!is.na(cv_df$cv), ]
  cv_df$qc_type <- factor(cv_df$qc_type, levels = qc_types)

  med_df <- cv_df %>%
    dplyr::group_by(qc_type) %>%
    dplyr::summarise(med_cv = median(cv, na.rm = TRUE), .groups = "drop")

  fill_vals <- sample_type_colors[qc_types]

  ggplot2::ggplot(cv_df, ggplot2::aes(x = qc_type, y = cv, fill = qc_type)) +
    ggplot2::geom_violin(scale = "width", bw = 0.03, trim = TRUE,
                          alpha = 0.7, colour = NA) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA,
                           fill = "white", colour = "grey30", alpha = 0.8) +
    ggplot2::scale_fill_manual(values = fill_vals, guide = "none") +
    ggplot2::scale_x_discrete(labels = display_names[qc_types]) +
    ggplot2::geom_hline(yintercept = cv_threshold_1, linetype = "dashed",
                         colour = "firebrick", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = cv_threshold_2, linetype = "dashed",
                         colour = "darkorange", linewidth = 0.5) +
    ggplot2::annotate("text", x = 0.6, y = cv_threshold_1,
                       label = paste0(round(cv_threshold_1 * 100), "%"),
                       colour = "firebrick", size = 2.5, hjust = 0, vjust = -0.4) +
    ggplot2::annotate("text", x = 0.6, y = cv_threshold_2,
                       label = paste0(round(cv_threshold_2 * 100), "%"),
                       colour = "darkorange", size = 2.5, hjust = 0, vjust = -0.4) +
    ggplot2::geom_text(data = med_df,
                        ggplot2::aes(x = qc_type, y = med_cv,
                                      label = sprintf("%.1f%%", med_cv * 100)),
                        vjust = -0.8, size = 2.5, colour = "black",
                        inherit.aes = FALSE) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      breaks = function(limits) {
        # Include standard pretty breaks plus the two threshold lines
        sort(unique(c(pretty(c(0, limits[2]), n = 5),
                      cv_threshold_1, cv_threshold_2)))
      }
    ) +
    ggplot2::labs(title = condition, x = NULL, y = "CV") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 30, hjust = 1, size = 9),
      plot.title      = ggplot2::element_text(size = 10, face = "bold"),
      legend.position = "none"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, y_cap))
}

# ---- 11.5 CV scatter: LTQC vs NIST per feature (PDF 2) ----
make_cv_scatter <- function(cv_ltqc, cv_nist, condition_name,
                             cv_threshold_1, cv_threshold_2) {
  df <- data.frame(ltqc_cv = cv_ltqc, nist_cv = cv_nist)
  df <- df[!is.na(df$ltqc_cv) & !is.na(df$nist_cv), ]
  pct_better <- mean(df$ltqc_cv < df$nist_cv, na.rm = TRUE) * 100
  df$better  <- ifelse(df$ltqc_cv < df$nist_cv, "Oslo-RM more stable", "NIST more stable")
  axis_cap   <- min(quantile(c(df$ltqc_cv, df$nist_cv), 0.99, na.rm = TRUE), 1.5)

  ggplot2::ggplot(df, ggplot2::aes(x = ltqc_cv, y = nist_cv, colour = better)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.2) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "black",
                          linewidth = 0.5, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = cv_threshold_1, colour = "#E69F00",
                         linewidth = 0.4, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = cv_threshold_1, colour = "#E69F00",
                         linewidth = 0.4, linetype = "dotted") +
    ggplot2::scale_colour_manual(
      values = c("Oslo-RM more stable" = "#E64B35", "NIST more stable" = "#4DBBD5")) +
    ggplot2::coord_cartesian(xlim = c(0, axis_cap), ylim = c(0, axis_cap)) +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = "Oslo-RM CV (per feature)", y = "NIST CV (per feature)",
                  title = condition_name, colour = NULL) +
    ggplot2::annotate("text", x = axis_cap * 0.05, y = axis_cap * 0.95,
                       label = sprintf("%.0f%% of features:\nOslo-RM more stable", pct_better),
                       hjust = 0, vjust = 1, size = 3, colour = "#E64B35") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "none", aspect.ratio = 1)
}

# ---- 11.6 Per-batch median trajectory (PDF 3) ----
make_batch_stability_plot <- function(log2_mat, metadata, qc_type,
                                       condition_name) {
  idx     <- which(metadata$sample_type == qc_type)
  batches <- sort(unique(metadata$batch[idx]))
  if (length(batches) == 0)
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = sprintf("%s — %s (no data)", qc_type, condition_name)))

  df <- do.call(rbind, lapply(batches, function(b) {
    b_idx   <- which(metadata$batch == b & metadata$sample_type == qc_type)
    raw_mat <- 2^log2_mat[, b_idx, drop = FALSE]
    data.frame(
      batch         = factor(b),
      sample_median = apply(raw_mat, 2, median, na.rm = TRUE),
      inj_order     = metadata$injection_order[b_idx]
    )
  }))

  df_sum <- df %>%
    dplyr::group_by(batch) %>%
    dplyr::summarise(
      med = median(sample_median),
      q25 = quantile(sample_median, 0.25),
      q75 = quantile(sample_median, 0.75),
      .groups = "drop")

  ggplot2::ggplot(df_sum, ggplot2::aes(x = batch, y = med, colour = batch, group = 1)) +
    ggplot2::geom_line(colour = sample_type_colors[[qc_type]],
                        linewidth = 0.7, linetype = "solid") +
    ggplot2::geom_point(ggplot2::aes(colour = batch), size = 3.5, shape = 16) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = q25, ymax = q75, colour = batch),
                            width = 0.25, linewidth = 0.6) +
    ggplot2::scale_colour_manual(values = batch_colors) +
    ggplot2::scale_y_continuous(labels = scales::scientific) +
    ggplot2::labs(x = "Batch", y = "Median raw intensity",
                  title = sprintf("%s — %s", display_names[[qc_type]], condition_name),
                  colour = "Batch") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "none")
}

# ---- 11.7 CV summary bar chart (PDF 4) ----
make_cv_summary_bars <- function(cv_results,
                                  qc_types   = c("LTQC", "NIST"),
                                  conditions = names(cv_results)) {
  df <- do.call(rbind, lapply(conditions, function(cond) {
    do.call(rbind, lapply(qc_types, function(qt) {
      data.frame(
        condition = factor(cond, levels = conditions),
        qc_type   = qt,
        median_cv = median(cv_results[[cond]][[qt]], na.rm = TRUE)
      )
    }))
  }))

  ggplot2::ggplot(df, ggplot2::aes(x = condition, y = median_cv, fill = qc_type)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                       width = 0.6, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = 0.20, colour = "#E69F00",
                         linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0.30, colour = "#E64B35",
                         linetype = "dashed", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = c("LTQC" = sample_type_colors[["LTQC"]],
                 "NIST" = sample_type_colors[["NIST"]]),
      labels = display_names[c("LTQC", "NIST")]) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 0.5)) +
    ggplot2::labs(x = NULL, y = "Median CV",
                  fill  = "Reference\nMaterial",
                  title = "Reference Material CV Summary") +
    ggplot2::theme_bw(base_size = 9)
}

# ============================================================
# SECTION 11.8: LEGENDS AND PRE-BUILT PLOT OBJECTS
# ============================================================

cat("Building row-specific legends...\n")

legend_theme <- ggplot2::theme(
  legend.position  = "right",
  legend.box       = "vertical",
  legend.title     = ggplot2::element_text(size = 9, face = "bold"),
  legend.text      = ggplot2::element_text(size = 8),
  legend.key.size  = ggplot2::unit(0.5, "cm"),
  legend.spacing.y = ggplot2::unit(0.2, "cm")
)

# --- Legend 1 (PCA row): batch colours + all sample type shapes/fills ---
legend_pca <- cowplot::get_legend(
  make_pca_plot(pca_results[["Uncorrected"]], "Uncorrected", show_legend = TRUE) +
    legend_theme)

# --- Legend 2 (injection row): QC type colour + shape only (no "Sample") ---
inj_legend_df <- data.frame(
  qc_type = factor(c("PQC", "LTQC", "NIST"), levels = c("PQC", "LTQC", "NIST")),
  x = 1L, y = 1:3,
  stringsAsFactors = FALSE)

legend_inj <- cowplot::get_legend(
  ggplot2::ggplot(inj_legend_df,
                   ggplot2::aes(x = x, y = y,
                                 colour = qc_type, fill = qc_type, shape = qc_type)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(
      values = sample_type_colors[c("PQC", "LTQC", "NIST")], name = "Sample Type",
      labels = display_names[c("PQC", "LTQC", "NIST")]) +
    ggplot2::scale_fill_manual(
      values = sample_type_colors[c("PQC", "LTQC", "NIST")], name = "Sample Type",
      labels = display_names[c("PQC", "LTQC", "NIST")]) +
    ggplot2::scale_shape_manual(
      values = sample_type_shapes[c("PQC", "LTQC", "NIST")], name = "Sample Type",
      labels = display_names[c("PQC", "LTQC", "NIST")]) +
    ggplot2::theme_void() + legend_theme)

# --- Legend 3 (CV row): fill colour swatches only, no shapes ---
cv_legend_df <- data.frame(
  qc_type = factor(c("PQC", "LTQC", "NIST"), levels = c("PQC", "LTQC", "NIST")),
  x = 1L, y = 1:3,
  stringsAsFactors = FALSE)

legend_cv <- cowplot::get_legend(
  ggplot2::ggplot(cv_legend_df, ggplot2::aes(x = x, y = y, fill = qc_type)) +
    ggplot2::geom_point(shape = 22, size = 5) +
    ggplot2::scale_fill_manual(
      values = sample_type_colors[c("PQC", "LTQC", "NIST")], name = "CV",
      labels = display_names[c("PQC", "LTQC", "NIST")]) +
    ggplot2::theme_void() + legend_theme)

# --- Combined 3-row legend: each section occupies equal height to match its row ---
# rel_heights = c(1, 1, 1) mirrors plot_layout(heights = c(2, 2, 2)) in the main panel
combined_legend_3row <- cowplot::plot_grid(
  legend_pca,
  legend_inj,
  legend_cv,
  ncol = 1, rel_heights = c(1, 1, 1))

# --- 2-row legend for PDFs with PCA + injection rows but no CV row ---
# Third cell is empty (ggdraw()) to keep alignment with the scatter row
combined_legend_2row <- cowplot::plot_grid(
  legend_pca,
  legend_inj,
  cowplot::ggdraw(),    # empty placeholder for scatter row
  ncol = 1, rel_heights = c(1, 1, 1))

# Pre-build plot objects reused across multiple PDFs
pca_plots <- lapply(conditions, function(cond)
  make_pca_plot(pca_results[[cond]], cond))

inj_plots <- lapply(conditions, function(cond)
  make_injection_plot(corrected_list[[cond]], metadata, cond,
                       y_lim = c(inj_y_min - inj_y_pad, inj_y_max + inj_y_pad)))

# CV plots: one panel per condition, x = QC type (PQC / LTQC / NIST)
cv_cond_plots <- lapply(conditions, function(cond)
  make_cv_condition_plot(cv_results, cond))

p_scatter <- lapply(conditions, function(cond)
  make_cv_scatter(cv_results[[cond]][["LTQC"]], cv_results[[cond]][["NIST"]],
                   cond, cv_threshold_1, cv_threshold_2))

p_ltqc_stab <- lapply(conditions, function(cond)
  make_batch_stability_plot(corrected_list[[cond]], metadata, "LTQC", cond))
p_nist_stab <- lapply(conditions, function(cond)
  make_batch_stability_plot(corrected_list[[cond]], metadata, "NIST", cond))

# Determine best condition for PDF 4 (lowest combined LTQC + NIST median CV)
best_cond <- conditions[which.min(vapply(conditions, function(cond) {
  median(cv_results[[cond]][["LTQC"]], na.rm = TRUE) +
  median(cv_results[[cond]][["NIST"]], na.rm = TRUE)
}, numeric(1)))]
cat(sprintf("Best condition (lowest combined LTQC+NIST CV): %s\n", best_cond))

p_cv_bars <- make_cv_summary_bars(cv_results, conditions = conditions)

# Helper: attach a legend grob to the right of a patchwork panel
attach_legend <- function(panel, legend_grob = combined_legend_3row,
                           rel_widths = c(12, 1.2)) {
  cowplot::plot_grid(panel, legend_grob, ncol = 2, rel_widths = rel_widths)
}

# Helper: save a plot as both PDF and TIFF (LZW compressed)
save_figure <- function(plot, path_stem, width, height, units = "mm", dpi = 300) {
  ggplot2::ggsave(paste0(path_stem, ".pdf"),  plot = plot,
                  width = width, height = height, units = units,
                  device = "pdf", dpi = dpi)
  ggplot2::ggsave(paste0(path_stem, ".tiff"), plot = plot,
                  width = width, height = height, units = units,
                  device = "tiff", dpi = dpi, compression = "lzw")
  cat(sprintf("  Saved: %s.pdf + .tiff\n", basename(path_stem)))
}

# Pre-build labeled PCA plots (for separate export in Section 11.13)
labeled_pca_plots <- lapply(conditions, function(cond)
  make_pca_labeled_plot(pca_results[[cond]], cond))

# ============================================================
# SECTION 11.9: PDF 1 — ORIGINAL CONCEPT CLEANED UP (v1_violins)
# ============================================================

cat("Assembling PDF 1 (v1_violins)...\n")
# Row 3: one CV panel per condition, comparing PQC / LTQC / NIST within each condition

row1_v1 <- patchwork::wrap_plots(pca_plots,     nrow = 1)
row2_v1 <- patchwork::wrap_plots(inj_plots,     nrow = 1)
row3_v1 <- patchwork::wrap_plots(cv_cond_plots, nrow = 1)

panel_v1 <- (row1_v1 / row2_v1 / row3_v1) +
  patchwork::plot_layout(heights = c(2, 2, 2))

save_figure(
  plot      = attach_legend(panel_v1, combined_legend_3row),
  path_stem = file.path(output_dir, "figures",
                         paste0(dataset_label, "_ref_material_v1_violins")),
  width = 420, height = 297, units = "mm", dpi = 300)
cat("PDF 1 saved.\n\n")

# ============================================================
# SECTION 11.10: PDF 2 — CV SCATTER COMPARISON (v2_cv_scatter)
# ============================================================

cat("Assembling PDF 2 (v2_cv_scatter)...\n")

row3_v2 <- patchwork::wrap_plots(p_scatter, nrow = 1)

panel_v2 <- (row1_v1 / row2_v1 / row3_v2) +
  patchwork::plot_layout(heights = c(2, 2, 2))

save_figure(
  plot      = attach_legend(panel_v2, combined_legend_2row),
  path_stem = file.path(output_dir, "figures",
                         paste0(dataset_label, "_ref_material_v2_cv_scatter")),
  width = 420, height = 297, units = "mm", dpi = 300)
cat("PDF 2 saved.\n\n")

# ============================================================
# SECTION 11.11: PDF 3 — BATCH STABILITY (v3_batch_stability)
# ============================================================

cat("Assembling PDF 3 (v3_batch_stability)...\n")

row_ltqc_v3 <- patchwork::wrap_plots(p_ltqc_stab, nrow = 1)
row_nist_v3 <- patchwork::wrap_plots(p_nist_stab, nrow = 1)

panel_v3 <- (row_ltqc_v3 / row_nist_v3) +
  patchwork::plot_layout(heights = c(1, 1))

save_figure(
  plot      = attach_legend(panel_v3, legend_pca),
  path_stem = file.path(output_dir, "figures",
                         paste0(dataset_label, "_ref_material_v3_batch_stability")),
  width = 420, height = 297, units = "mm", dpi = 300)
cat("PDF 3 saved.\n\n")

# ============================================================
# SECTION 11.12: PDF 4 — COMPACT SUMMARY (v4_summary)
# ============================================================

cat(sprintf("Assembling PDF 4 (v4_summary) — best condition: %s...\n", best_cond))

p_best_ltqc <- make_batch_stability_plot(corrected_list[[best_cond]], metadata,
                                          "LTQC", best_cond)
p_best_nist <- make_batch_stability_plot(corrected_list[[best_cond]], metadata,
                                          "NIST", best_cond)

row1_v4 <- patchwork::wrap_plots(p_scatter, nrow = 1)
row2_v4 <- (p_best_ltqc | p_best_nist | p_cv_bars)

panel_v4 <- (row1_v4 / row2_v4) +
  patchwork::plot_layout(heights = c(1, 1))

save_figure(
  plot      = attach_legend(panel_v4, legend_pca),
  path_stem = file.path(output_dir, "figures",
                         paste0(dataset_label, "_ref_material_v4_summary")),
  width = 297, height = 210, units = "mm", dpi = 300)
cat("PDF 4 saved.\n\n")

# ============================================================
# SECTION 11.13: LABELED PCA FIGURE (separate export)
# ============================================================

cat("Assembling labeled PCA figure...\n")

panel_pca_labeled <- patchwork::wrap_plots(labeled_pca_plots, nrow = 1)

save_figure(
  plot      = attach_legend(panel_pca_labeled, legend_pca),
  path_stem = file.path(output_dir, "figures",
                         paste0(dataset_label, "_pca_labeled")),
  width = 420, height = 140, units = "mm", dpi = 300)
cat("Labeled PCA figure saved.\n\n")

# ============================================================
# SECTION 11.14: PUBLICATION FIGURE — A/B/C row labels (for_publication)
# ============================================================

cat("Assembling publication figure (A/B/C labels)...\n")

# Rebuild the same 3-row panel as v1 (reuses pre-built plot lists)
panel_pub <- (
  patchwork::wrap_plots(pca_plots,     nrow = 1) /
  patchwork::wrap_plots(inj_plots,     nrow = 1) /
  patchwork::wrap_plots(cv_cond_plots, nrow = 1)
) + patchwork::plot_layout(heights = c(1, 1, 1))

pub_base <- attach_legend(panel_pub, combined_legend_3row)

# Overlay A, B, C at the top-left of each row using cowplot::draw_label.
# With 3 equal-height rows and the full figure height = 1.0:
#   Row 1 top ≈ 0.99, Row 2 top ≈ 0.66, Row 3 top ≈ 0.33
# x = 0.005 sits inside the left-most PCA panel, just clear of the y-axis text.
pub_figure <- cowplot::ggdraw(pub_base) +
  cowplot::draw_label("A", x = 0.005, y = 0.99, hjust = 0, vjust = 1,
                       fontface = "bold", size = 16) +
  cowplot::draw_label("B", x = 0.005, y = 0.66, hjust = 0, vjust = 1,
                       fontface = "bold", size = 16) +
  cowplot::draw_label("C", x = 0.005, y = 0.33, hjust = 0, vjust = 1,
                       fontface = "bold", size = 16)

save_figure(
  plot      = pub_figure,
  path_stem = file.path(output_dir, "figures",
                         paste0(dataset_label, "_for_publication")),
  width = 420, height = 297, units = "mm", dpi = 300)
cat("Publication figure saved.\n\n")

# ============================================================
# SECTION 13: CONSOLE SUMMARY REPORT
# ============================================================

cat("======== SECTION 13: SUMMARY REPORT ========\n\n")

get_cv_stat <- function(cond, qc, stat = "median") {
  cvs <- cv_results[[cond]][[qc]]; cvs <- cvs[!is.na(cvs)]
  if (!length(cvs)) return(NA_real_)
  if (stat == "median") median(cvs) else mean(cvs < cv_threshold_1) * 100
}
fmt_med <- function(x) if (is.na(x)) "  N/A  " else sprintf("%5.1f%%", x * 100)
fmt_pct <- function(x) if (is.na(x)) "  N/A  " else sprintf("%5.1f%%", x)

header <- sprintf("%-16s | %14s | %6s | %15s | %6s | %15s | %6s",
                   "Condition", "PQC CV (med)", "%<20%",
                   "LTQC CV (med)", "%<20%", "NIST CV (med)", "%<20%")
sep    <- paste(rep("-", nchar(header)), collapse = "")

cat(paste(rep("=", 54), collapse = ""), "\n")
cat(sprintf(" LTQC REFERENCE MATERIAL SUMMARY — %s\n", dataset_label))
cat(paste(rep("=", 54), collapse = ""), "\n\n")
cat(header, "\n", sep, "\n")

for (cond in names(corrected_list)) {
  cat(sprintf("%-16s | %14s | %6s | %15s | %6s | %15s | %6s\n",
              cond,
              fmt_med(get_cv_stat(cond, "PQC",  "median")),
              fmt_pct(get_cv_stat(cond, "PQC",  "pct")),
              fmt_med(get_cv_stat(cond, "LTQC", "median")),
              fmt_pct(get_cv_stat(cond, "LTQC", "pct")),
              fmt_med(get_cv_stat(cond, "NIST", "median")),
              fmt_pct(get_cv_stat(cond, "NIST", "pct"))))
}
cat(paste(rep("=", 54), collapse = ""), "\n")
cat(sprintf("\nBest condition (lowest LTQC+NIST CV): %s\n", best_cond))
cat("\nOutputs saved to:", output_dir, "\n")
cat("PDFs: v1_violins | v2_cv_scatter | v3_batch_stability | v4_summary\n")
