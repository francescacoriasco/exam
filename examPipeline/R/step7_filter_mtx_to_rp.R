#' @title Step 7 — Filter Step 5 combined scRNA MTX to RP-only features
#' @description
#' Loads a 10x-like MTX bundle (\code{matrix.mtx}, \code{features.tsv}, \code{barcodes.tsv})
#' and keeps only RP features using the RP Ensembl \code{gene_id} list from Step 6
#' (or falls back to \code{startsWith("RP")} if features are symbols).
#' Supports gz inputs (\code{.gz}) for convenience.
#' @param in_dir Input directory containing the MTX bundle (from Step 5).
#' @param rp_ids_rds Path to RP gene_id RDS (from Step 6).
#' @param out_dir Output directory for RP-only MTX bundle.
#' @param ensembl_fraction_threshold If fraction of IDs matching \code{^ENSG} is > threshold, treat as Ensembl IDs.
#' @param verbose If TRUE, prints summary.
#' @return Invisible list with \code{out_dir}, \code{dims}, \code{kept_features}.
#' @export
step7_filter_mtx_to_rp <- function(in_dir,
                                   rp_ids_rds,
                                   out_dir,
                                   ensembl_fraction_threshold = 0.5,
                                   verbose = TRUE) {
  stopifnot(is.character(in_dir), length(in_dir) == 1L, dir.exists(in_dir))
  stopifnot(is.character(rp_ids_rds), length(rp_ids_rds) == 1L, file.exists(rp_ids_rds))
  stopifnot(is.character(out_dir), length(out_dir) == 1L)
  stopifnot(is.numeric(ensembl_fraction_threshold), length(ensembl_fraction_threshold) == 1L)
  stopifnot(is.logical(verbose), length(verbose) == 1L)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  pick_existing <- function(paths) {
    p <- paths[file.exists(paths)][1]
    if (is.na(p)) stop("Missing file. Tried:\n", paste(paths, collapse = "\n"))
    p
  }
  read_lines_maybe_gz <- function(path) {
    if (grepl("\\.gz$", path)) {
      con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
      readLines(con, warn = FALSE)
    } else {
      readLines(path, warn = FALSE)
    }
  }
  read_tsv_maybe_gz <- function(path) {
    if (grepl("\\.gz$", path)) {
      con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
      utils::read.delim(con, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    } else {
      utils::read.delim(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    }
  }
  read_mtx_maybe_gz <- function(path) {
    if (grepl("\\.gz$", path)) {
      con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
      Matrix::readMM(con)
    } else {
      Matrix::readMM(path)
    }
  }

  mtx_in <- pick_existing(file.path(in_dir, c("matrix.mtx", "matrix.mtx.gz")))
  fea_in <- pick_existing(file.path(in_dir, c("features.tsv", "features.tsv.gz")))
  bar_in <- pick_existing(file.path(in_dir, c("barcodes.tsv", "barcodes.tsv.gz")))

  mat      <- read_mtx_maybe_gz(mtx_in)
  features <- read_tsv_maybe_gz(fea_in)
  barcodes <- read_lines_maybe_gz(bar_in)

  stopifnot(nrow(features) == nrow(mat))
  stopifnot(length(barcodes) == ncol(mat))

  rp_ids <- readRDS(rp_ids_rds)
  if (!is.character(rp_ids) || length(rp_ids) == 0) stop("rp_ids_rds must contain a non-empty character vector.")
  rp_ids <- unique(rp_ids)

  feat1 <- features[[1]]
  is_ensembl <- mean(grepl("^ENSG", feat1)) > ensembl_fraction_threshold

  if (is_ensembl) {
    keep <- feat1 %in% rp_ids
    features_out <- cbind(feat1[keep], feat1[keep], rep("Gene Expression", sum(keep)))
  } else {
    keep <- startsWith(feat1, "RP")
    features_out <- cbind(feat1[keep], feat1[keep], rep("Gene Expression", sum(keep)))
  }

  if (sum(keep) == 0) stop("RP filter kept 0 features. Check features.tsv and RP list.")
  mat_rp <- mat[keep, , drop = FALSE]

  Matrix::writeMM(mat_rp, file = file.path(out_dir, "matrix.mtx"))
  writeLines(barcodes, con = file.path(out_dir, "barcodes.tsv"))
  utils::write.table(features_out,
                     file = file.path(out_dir, "features.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  if (verbose) {
    message("Step 7 DONE: ", out_dir)
    message("Dims RP-only (features x cells): ", nrow(mat_rp), " x ", ncol(mat_rp))
    message("Kept RP features: ", sum(keep), " / ", length(keep))
  }

  invisible(list(
    out_dir = out_dir,
    dims = c(features = nrow(mat_rp), cells = ncol(mat_rp)),
    kept_features = sum(keep)
  ))
}