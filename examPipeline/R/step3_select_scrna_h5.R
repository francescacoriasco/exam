#' @title Step 3 — Select scRNA-seq 10x HDF5 files by content (Gene Expression/RNA)
#' @description
#' Scans a directory for 10x Genomics \code{.h5} files and keeps only those that contain
#' \dQuote{Gene Expression}/\dQuote{RNA} features (via \code{rowData()} columns \code{Type}
#' or \code{feature_type}). The helper used to inspect file content is defined *inside*
#' this function (not exported, not separate).
#' @param data_dir Directory to scan for \code{.h5} files.
#' @param recursive If TRUE, scan recursively.
#' @param pattern Regex for \code{list.files()} (default \code{"\\\\.h5$"}).
#' @param bpparam \code{BiocParallelParam} for \code{DropletUtils::read10xCounts()}.
#' @param rna_types Feature types treated as scRNA (default \code{c("Gene Expression","RNA")}).
#' @param assume_rna_if_missing If TRUE and type metadata missing, assume RNA-only.
#' @param sort_files If TRUE, return deterministically sorted paths.
#' @param verbose If TRUE, prints progress.
#' @return Character vector of kept \code{.h5} paths.
#' @export
step3_select_scrna_h5 <- function(data_dir,
                                  recursive = TRUE,
                                  pattern = "\\.h5$",
                                  bpparam = BiocParallel::SerialParam(),
                                  rna_types = c("Gene Expression", "RNA"),
                                  assume_rna_if_missing = TRUE,
                                  sort_files = TRUE,
                                  verbose = TRUE) {
  stopifnot(is.character(data_dir), length(data_dir) == 1L, dir.exists(data_dir))
  stopifnot(is.logical(recursive), length(recursive) == 1L)
  stopifnot(is.character(pattern), length(pattern) == 1L)
  stopifnot(is.logical(assume_rna_if_missing), length(assume_rna_if_missing) == 1L)
  stopifnot(is.logical(sort_files), length(sort_files) == 1L)
  stopifnot(is.logical(verbose), length(verbose) == 1L)

  # --- helper INSIDE step3 (as requested) ---
  has_gene_expression <- function(f) {
    sce <- DropletUtils::read10xCounts(
      samples   = f,
      type      = "HDF5",
      col.names = TRUE,
      BPPARAM   = bpparam
    )
    rd <- SummarizedExperiment::rowData(sce)

    ft <- NULL
    if ("Type" %in% colnames(rd)) {
      ft <- rd$Type
    } else if ("feature_type" %in% colnames(rd)) {
      ft <- rd$feature_type
    }

    ok <- if (is.null(ft)) {
      isTRUE(assume_rna_if_missing)
    } else {
      any(ft %in% rna_types, na.rm = TRUE)
    }

    rm(sce, rd); gc()
    ok
  }
  # --- end helper ---

  h5_files <- list.files(data_dir, pattern = pattern, full.names = TRUE, recursive = recursive)
  if (length(h5_files) == 0L) stop("No .h5 files found in: ", data_dir)
  if (isTRUE(sort_files)) h5_files <- sort(h5_files)

  keep <- logical(length(h5_files))
  for (i in seq_along(h5_files)) {
    f <- h5_files[i]
    if (verbose) message("Checking content: ", basename(f), " ...")
    keep[i] <- tryCatch(has_gene_expression(f), error = function(e) {
      if (verbose) message("  -> ERROR reading file; skipping: ", conditionMessage(e))
      FALSE
    })
    if (verbose) message("  -> Gene Expression/RNA present? ", keep[i])
  }

  scrna_files <- h5_files[keep]
  if (length(scrna_files) == 0L) stop("No scRNA .h5 files detected in: ", data_dir)
  if (isTRUE(sort_files)) scrna_files <- sort(scrna_files)
  scrna_files
}
