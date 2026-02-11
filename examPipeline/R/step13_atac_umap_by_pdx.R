#' @title Step 13 — ATAC clustering + UMAP, plot colored by PDX
#' @description
#' Loads the combined ATAC object (Step 12) or accepts it directly, runs TF-IDF/LSI,
#' clustering and UMAP, then saves a UMAP plot colored by PDX.
#' @param atac Either a Seurat object or a path to \code{.rds}.
#' @param out_dir Output directory for plots.
#' @param pdx_col Metadata column name to use for PDX (if NULL, tries common candidates).
#' @param pdx_candidates Candidate metadata columns if \code{pdx_col} is NULL.
#' @param resolution Clustering resolution.
#' @param k_param Neighbor k (NULL = auto).
#' @param max_lsi Max LSI components to consider.
#' @param dims_use LSI dims to use (NULL = auto, default skips dim 1).
#' @param save_plots If TRUE, saves png+pdf.
#' @param verbose If TRUE, prints summary.
#' @return Invisible list with \code{atac}, \code{plot}, \code{plot_files}.
#' @export
step13_atac_umap_by_pdx <- function(atac,
                                    out_dir,
                                    pdx_col = NULL,
                                    pdx_candidates = c("pdx", "PDX", "pdx_scRNA", "pdx_of_origin"),
                                    resolution = 0.5,
                                    k_param = NULL,
                                    max_lsi = 30L,
                                    dims_use = NULL,
                                    save_plots = TRUE,
                                    verbose = TRUE) {
  stopifnot(is.character(out_dir), length(out_dir) == 1L)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.character(atac)) {
    stopifnot(length(atac) == 1L, file.exists(atac))
    atac <- readRDS(atac)
  }
  stopifnot(inherits(atac, "Seurat"))

  md <- atac[[]]
  if (is.null(pdx_col)) {
    hit <- pdx_candidates[pdx_candidates %in% colnames(md)][1]
    if (is.na(hit)) stop("No PDX column found. Provide pdx_col or add metadata.")
    pdx_col <- hit
  }
  atac$pdx <- md[[pdx_col]]

  Seurat::DefaultAssay(atac) <- "ATAC"

  atac <- Signac::RunTFIDF(atac)
  atac <- Signac::FindTopFeatures(atac, min.cutoff = "q0")
  atac <- Signac::RunSVD(atac)

  ndim <- min(as.integer(max_lsi), ncol(Seurat::Embeddings(atac, "lsi")))
  if (ndim < 3) stop("Too few LSI dims: ", ndim)

  if (is.null(dims_use)) {
    dims_use <- 2:min(15, ndim)
  }
  if (is.null(k_param)) {
    k_param <- min(20, ncol(atac) - 1)
  }

  atac <- Seurat::FindNeighbors(atac, reduction = "lsi", dims = dims_use, k.param = as.integer(k_param), verbose = FALSE)
  atac <- Seurat::FindClusters(atac, resolution = resolution, verbose = FALSE)
  atac <- Seurat::RunUMAP(atac, reduction = "lsi", dims = dims_use, verbose = FALSE)

  p_pdx <- Seurat::DimPlot(atac, reduction = "umap", group.by = "pdx", label = TRUE)

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
    message("Step 13 DONE")
    message("Saved plots to: ", out_dir)
    message("Cells: ", ncol(atac), " Features: ", nrow(atac))
  }

  invisible(list(atac = atac, plot = p_pdx, plot_files = plot_files))
}
