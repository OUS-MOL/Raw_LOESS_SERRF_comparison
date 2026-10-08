# ============================================================
# Generate a fully synthetic example dataset in the input format expected by
# normalization_comparison.R (Compound Discoverer "Compounds" export + sequence).
# No real data are used: feature names, intensities, drift and batch effects are
# all simulated with a fixed seed.
#
# Usage:  Rscript make_example_data.R          # writes example_data/*.xlsx
# ============================================================

if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")

set.seed(2026)
out_dir <- "example_data"
dir.create(out_dir, showWarnings = FALSE)

n_features  <- 150
n_batches   <- 6
n_study     <- 10     # study samples per batch
n_pqc       <- 12     # pooled QC injections per batch
n_ltqc      <- 3      # long-term QC (Oslo-RM) injections per batch
n_nist      <- 3      # NIST injections per batch

# ---- feature table -------------------------------------------------------
feature_mean <- rnorm(n_features, mean = 20, sd = 1.5)       # log2 area
feat <- data.frame(
  Name        = sprintf("Example lipid %03d", seq_len(n_features)),
  Formula     = sprintf("C%dH%dO2", sample(14:40, n_features, TRUE), sample(20:80, n_features, TRUE)),
  `RT [min]`  = round(runif(n_features, 1, 15), 2),
  `m/z`       = round(runif(n_features, 200, 1000), 4),
  `Calc. MW`  = round(runif(n_features, 200, 1000), 4),
  check.names = FALSE, stringsAsFactors = FALSE)

# ---- sequence: one block of injections per batch ---------------------------
rows <- list(); fid <- 0L; t0 <- as.POSIXct("2025-01-06 08:00:00", tz = "UTC")
for (b in 0:(n_batches - 1)) {
  types <- c(rep("PQC", n_pqc), rep("LTQC", n_ltqc), rep("NIST", n_nist),
             rep("Study", n_study), "Blank")
  # spread QCs through the batch: shuffle study/LTQC/NIST, keep PQC evenly spaced
  body  <- sample(c(rep("Study", n_study), rep("LTQC", n_ltqc), rep("NIST", n_nist)))
  pos   <- round(seq(1, length(body) + n_pqc, length.out = n_pqc))
  seqb  <- rep(NA_character_, length(body) + n_pqc)
  seqb[pos] <- "PQC"; seqb[is.na(seqb)] <- body
  seqb  <- c(seqb, "Blank")
  for (k in seq_along(seqb)) {
    fid <- fid + 1L
    ty  <- seqb[k]
    nm  <- switch(ty,
      PQC   = sprintf("PQC_%d_B%d_fullMS_POS.raw", k, b),
      LTQC  = sprintf("LTQC_%d_B%d_fullMS_POS.raw", k, b),
      NIST  = sprintf("NIST_%d_B%d_fullMS_POS.raw", k, b),
      Study = sprintf("S%03d_B%d_fullMS_POS.raw", fid, b),
      Blank = sprintf("Blank_%d_B%d_fullMS_POS.raw", k, b))
    rows[[length(rows) + 1]] <- data.frame(
      `Study File ID` = paste0("F", fid),
      `File Name`     = paste0("D:\\example\\", nm),
      `Creation Date` = format(t0 + (fid - 1) * 1800, "%d.%m.%Y %H:%M:%S"),
      Batch           = b,
      `Sample Type`   = switch(ty, PQC = "Quality Control", Blank = "Blank", "Sample"),
      Group           = switch(ty, LTQC = "LTQC", NIST = "NIST", Study = "Ctrl", "n/a"),
      Type = ty, FileBase = nm, fid = fid, k = k, nk = length(seqb),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  t0 <- t0 + 6 * 3600   # gap between batches
}
seqdf <- do.call(rbind, rows)

# ---- two identification-only (MS/MS) runs with batch -1: should be excluded ----
ident <- data.frame(
  `Study File ID` = c("F900", "F901"),
  `File Name`     = paste0("D:\\example\\", c("PQC_ID_1_MSMS_POS.raw", "PQC_ID_2_MSMS_POS.raw")),
  `Creation Date` = c("05.01.2025 08:00:00", "05.01.2025 09:00:00"),
  Batch = -1, `Sample Type` = "Identification Only", Group = "n/a",
  check.names = FALSE, stringsAsFactors = FALSE)

# ---- simulate log2 areas ----------------------------------------------------
n <- nrow(seqdf)
loading   <- rnorm(n_features, 1, 0.5)                         # feature sensitivity to drift
batch_eff <- rnorm(n_batches, 0, 0.8)[seqdf$Batch + 1]         # batch offset
pos_in_b  <- (seqdf$k - 1) / (seqdf$nk - 1)                    # 0..1 within batch
drift     <- 0.9 * (pos_in_b - 0.5) + 0.5 * sin(2 * pi * pos_in_b)  # within-batch trend
type_shift <- c(PQC = 0, Study = 0, LTQC = 0.15, NIST = -0.2, Blank = -6)[seqdf$Type]
bio_sd     <- c(PQC = 0.05, Study = 0.6, LTQC = 0.05, NIST = 0.05, Blank = 0.1)[seqdf$Type]

# study samples also carry a feature-specific biological deviation
bio <- matrix(rnorm(n * n_features, 0, 1), n, n_features) * bio_sd
log2_area <- t(sapply(seq_len(n), function(i)
  feature_mean + loading * (batch_eff[i] + drift[i]) + type_shift[i] + bio[i, ] +
  rnorm(n_features, 0, 0.18)))                                  # measurement noise
area <- 2^log2_area

# ---- assemble Compound Discoverer-style table ---------------------------------
area_df <- as.data.frame(t(area))
colnames(area_df) <- sprintf("Area: %s (%s)", seqdf$FileBase, seqdf$`Study File ID`)
compounds <- data.frame(feat, area_df, check.names = FALSE)

seq_out <- seqdf[, c("Study File ID", "File Name", "Creation Date", "Batch",
                     "Sample Type", "Group")]
seq_out <- rbind(seq_out, ident)

writexl::write_xlsx(list(Compounds  = compounds), file.path(out_dir, "compounds_example.xlsx"))
writexl::write_xlsx(list(InputFiles = seq_out),   file.path(out_dir, "sequence_example.xlsx"))
cat(sprintf("Wrote %d features x %d injections to %s/\n", n_features, n, out_dir))
