#' @title Step 8 — Seurat clustering + UMAP, plot colored by PDX (RP scRNA)
#' @description
#' Loads the RP-only MTX bundle from Step 7, builds a Seurat object, derives PDX labels
#' from cell names, runs a standard scRNA pipeline (PCA/Neighbors/Clusters/UMAP),
#' and saves a UMAP plot colored by PDX.
#' @param in_dir Input directory containing Step 7 MTX bundle.
#' @param out_dir Output directory for plots.
#' @param project Seurat project name.
#' @param pdx_fun Function mapping cell name -> PDX label (default: prefix before first underscore).
#' @param npcs Number of PCs (NULL = auto).
#' @param dims_use PCA dims used (NULL = auto).
#' @param k_param Neighbors k (NULL = auto).
#' @param resolution Clustering resolution.
#' @param nfeatures_var Variable features count (capped by number of genes).
#' @param save_plots If TRUE, saves png+pdf.
#' @param verbose If TRUE, prints summary.
#' @return Invisible list with \code{obj}, \code{plot}, and \code{plot_files}.
#' @export
step8_seurat_umap_by_pdx <- function(in_dir,
                                     out_dir,
                                     project = "RP_only",
                                     pdx_fun = function(cell) sub("_.*$", "", cell),
                                     npcs = NULL,
                                     dims_use = NULL,
                                     k_param = NULL,
                                     resolution = 0.5,
                                     nfeatures_var = 30L,
                                     save_plots = TRUE,
                                     verbose = TRUE) {
  stopifnot(is.character(in_dir), length(in_dir) == 1L, dir.exists(in_dir))
  stopifnot(is.character(out_dir), length(out_dir) == 1L)
  stopifnot(is.function(pdx_fun))
  stopifnot(is.logical(save_plots), length(save_plots) == 1L)
  stopifnot(is.logical(verbose), length(verbose) == 1L)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  mat <- Matrix::readMM(file.path(in_dir, "matrix.mtx"))
  features <- utils::read.delim(file.path(in_dir, "features.tsv"),
                                header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  barcodes <- readLines(file.path(in_dir, "barcodes.tsv"))
  stopifnot(nrow(features) == nrow(mat), length(barcodes) == ncol(mat))

  rownames(mat) <- features[[1]]
  colnames(mat) <- barcodes

  pdx <- vapply(colnames(mat), pdx_fun, character(1))

  obj <- Seurat::CreateSeuratObject(counts = mat, project = project)
  obj$pdx <- pdx

  nfeat <- nrow(obj)
  if (is.null(npcs)) npcs <- max(2, min(20, nfeat - 1))
  if (is.null(dims_use)) dims_use <- 1:min(10, npcs)
  if (is.null(k_param))  k_param <- min(20, max(10, floor(sqrt(ncol(obj)))))

  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst",
                                      nfeatures = min(as.integer(nfeatures_var), nfeat),
                                      verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, npcs = as.integer(npcs), verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, dims = dims_use, k.param = as.integer(k_param), verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = resolution, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, dims = dims_use, verbose = FALSE)

  p_pdx <- Seurat::DimPlot(obj, reduction = "umap", group.by = "pdx", label = TRUE)

  plot_files <- character()
  if (isTRUE(save_plots)) {
    png_path <- file.path(out_dir, "umap_by_pdx.png")
    pdf_path <- file.path(out_dir, "umap_by_pdx.pdf")
    grDevices::png(png_path, width = 1200, height = 900, res = 150)
    print(p_pdx); grDevices::dev.off()
    grDevices::pdf(pdf_path, width = 8, height = 6)
    print(p_pdx); grDevices::dev.off()
    plot_files <- c(png = png_path, pdf = pdf_path)
  }

  if (verbose) {
    message("Step 8 DONE")
    message("Features: ", nrow(obj), " Cells: ", ncol(obj))
    message("PCs: ", npcs, " dims: ", paste(dims_use, collapse = ","), " k: ", k_param)
    message("Saved plots to: ", out_dir)
  }

  invisible(list(obj = obj, plot = p_pdx, plot_files = plot_files))
}
