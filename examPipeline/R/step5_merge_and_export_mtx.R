#' @title Step 5 — Merge subsampled scRNA Seurat objects and export 10x-like MTX bundle
#' @description
#' Merges Seurat objects, consolidates counts layers when available (Seurat v5),
#' and writes \code{matrix.mtx}, \code{barcodes.tsv}, \code{features.tsv}.
#' @param objs Named list of Seurat objects (output of Step 4).
#' @param out_dir Output directory.
#' @param prefix Optional filename prefix (NULL writes standard 10x names).
#' @param overwrite If FALSE, stops when output exists.
#' @param assay Assay to export (default \code{"RNA"}).
#' @param layer Layer/slot name to export (default \code{"counts"}).
#' @param feature_type Third column in \code{features.tsv}.
#' @param join_layers If TRUE and available, calls \code{SeuratObject::JoinLayers()}.
#' @param verbose If TRUE, prints summary.
#' @return Invisible list with \code{combined}, \code{dims}, and \code{out_files}.
#' @export
step5_merge_and_export_mtx <- function(objs,
                                       out_dir,
                                       prefix = NULL,
                                       overwrite = TRUE,
                                       assay = "RNA",
                                       layer = "counts",
                                       feature_type = "Gene Expression",
                                       join_layers = TRUE,
                                       verbose = TRUE) {
  stopifnot(is.list(objs), length(objs) >= 1L)
  stopifnot(is.character(out_dir), length(out_dir) == 1L)
  stopifnot(is.null(prefix) || (is.character(prefix) && length(prefix) == 1L && nzchar(prefix)))
  stopifnot(is.logical(overwrite), length(overwrite) == 1L)
  stopifnot(is.logical(join_layers), length(join_layers) == 1L)
  stopifnot(is.logical(verbose), length(verbose) == 1L)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  fname <- function(stem) if (is.null(prefix)) stem else paste0(prefix, "_", stem)
  f_matrix   <- file.path(out_dir, fname("matrix.mtx"))
  f_barcodes <- file.path(out_dir, fname("barcodes.tsv"))
  f_features <- file.path(out_dir, fname("features.tsv"))

  if (!overwrite) {
    existing <- c(f_matrix, f_barcodes, f_features)
    existing <- existing[file.exists(existing)]
    if (length(existing) > 0L) stop("Output file(s) exist and overwrite=FALSE:\n", paste(existing, collapse = "\n"))
  }

  if (!is.null(names(objs)) && all(nzchar(names(objs)))) {
    objs <- objs[sort(names(objs))]
  }

  combined <- objs[[1]]
  if (length(objs) > 1L) {
    for (i in 2:length(objs)) combined <- merge(combined, y = objs[[i]])
  }
  Seurat::DefaultAssay(combined) <- assay

  has_joinlayers <- requireNamespace("SeuratObject", quietly = TRUE) &&
    "JoinLayers" %in% getNamespaceExports("SeuratObject")

  if (isTRUE(join_layers) && has_joinlayers) {
    combined <- SeuratObject::JoinLayers(object = combined, assay = assay, layers = layer)
    mat <- Seurat::GetAssayData(combined, assay = assay, layer = layer)
  } else {
    mat <- Seurat::GetAssayData(combined, assay = assay, slot = layer)
  }

  Matrix::writeMM(mat, file = f_matrix)
  writeLines(colnames(mat), con = f_barcodes)
  utils::write.table(
    cbind(rownames(mat), rownames(mat), rep(feature_type, nrow(mat))),
    file = f_features,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  if (verbose) {
    message("DONE: ", out_dir)
    message("Combined dims (genes x cells): ", nrow(mat), " x ", ncol(mat))
    message("Datasets merged: ", length(objs))
  }

  invisible(list(
    combined  = combined,
    dims      = c(genes = nrow(mat), cells = ncol(mat)),
    out_files = list(matrix = f_matrix, barcodes = f_barcodes, features = f_features)
  ))
}
