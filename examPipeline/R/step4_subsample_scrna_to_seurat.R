#' @title Step 4 — Subsample scRNA-seq cells reproducibly and build Seurat objects
#' @description
#' For each scRNA \code{.h5} file, reads counts, keeps only \dQuote{Gene Expression}/\dQuote{RNA}
#' features when annotated, subsamples \code{sample_n} cells reproducibly, and returns one
#' Seurat object per dataset.
#' @param scrna_files Character vector of \code{.h5} paths (typically output of Step 3).
#' @param sample_n Number of cells per dataset (keeps all if fewer).
#' @param seed Random seed for reproducible sampling.
#' @param bpparam \code{BiocParallelParam} for \code{DropletUtils::read10xCounts()}.
#' @param rna_types Feature types treated as scRNA.
#' @param prefix_cells If TRUE, prefixes cell IDs with sample_id to avoid collisions.
#' @param sample_id_fun Function mapping file path -> sample ID.
#' @param verbose If TRUE, prints progress.
#' @return Named list of Seurat objects.
#' @export
step4_subsample_scrna_to_seurat <- function(scrna_files,
                                            sample_n = 500L,
                                            seed = 1L,
                                            bpparam = BiocParallel::SerialParam(),
                                            rna_types = c("Gene Expression", "RNA"),
                                            prefix_cells = TRUE,
                                            sample_id_fun = function(path) tools::file_path_sans_ext(basename(path)),
                                            verbose = TRUE) {
  stopifnot(is.character(scrna_files), length(scrna_files) >= 1L, all(file.exists(scrna_files)))
  stopifnot(is.numeric(sample_n), length(sample_n) == 1L, sample_n >= 1)
  stopifnot(is.numeric(seed), length(seed) == 1L)
  stopifnot(is.logical(prefix_cells), length(prefix_cells) == 1L)
  stopifnot(is.function(sample_id_fun))
  stopifnot(is.logical(verbose), length(verbose) == 1L)

  scrna_files <- sort(scrna_files)
  set.seed(as.integer(seed))

  objs <- list()

  for (f in scrna_files) {
    sample_id <- sample_id_fun(f)
    if (verbose) message("Reading: ", sample_id)

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

    if (!is.null(ft) && any(ft %in% rna_types, na.rm = TRUE)) {
      sce <- sce[ft %in% rna_types, , drop = FALSE]
    }

    counts <- SummarizedExperiment::assay(sce, "counts")
    if (ncol(counts) == 0L) {
      if (verbose) message("  -> 0 cells; skipping.")
      rm(sce, rd, counts); gc()
      next
    }

    k <- min(as.integer(sample_n), ncol(counts))
    idx <- if (k == ncol(counts)) seq_len(ncol(counts)) else sample.int(ncol(counts), k)

    mat_k <- as(counts[, idx, drop = FALSE], "dgCMatrix")
    rownames(mat_k) <- rownames(sce)

    if (isTRUE(prefix_cells)) {
      colnames(mat_k) <- paste0(sample_id, "_", colnames(mat_k))
    }

    objs[[sample_id]] <- Seurat::CreateSeuratObject(counts = mat_k, project = sample_id)

    rm(sce, rd, counts, mat_k); gc()
  }

  if (length(objs) == 0L) stop("No Seurat objects created (all datasets empty or skipped).")
  objs
}