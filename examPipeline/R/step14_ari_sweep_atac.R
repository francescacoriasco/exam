#' @title Step 14 — Sweep ATAC clustering parameters and score PDX separation with ARI
#' @description
#' Runs TF-IDF/LSI (optionally) and sweeps LSI dims, neighbors k, and clustering resolution.
#' Computes ARI between clusters and PDX labels.
#' @param atac Either a Seurat object or a path to \code{.rds} (Step 12 output).
#' @param out_dir Output directory to save CSV/RDS.
#' @param pdx_col Metadata column containing PDX labels (default \code{"pdx"}).
#' @param do_preprocess If TRUE, runs TF-IDF/TopFeatures/SVD before sweep.
#' @param grid_ndims Number of LSI dims to use (from a pool that skips dim 1).
#' @param grid_k Neighbors k grid.
#' @param grid_res Resolution grid.
#' @param seed Seed for reproducibility (set each iteration).
#' @param verbose If TRUE, prints summary.
#' @return Invisible list with \code{results}, \code{best}, \code{out_files}.
#' @export
step14_ari_sweep_atac <- function(atac,
                                  out_dir,
                                  pdx_col = "pdx",
                                  do_preprocess = TRUE,
                                  grid_ndims = c(10, 20, 30),
                                  grid_k     = c(10, 20, 30),
                                  grid_res   = c(0.2, 0.5, 0.8, 1.2),
                                  seed = 1L,
                                  verbose = TRUE) {
  stopifnot(is.character(out_dir), length(out_dir) == 1L)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.character(atac)) {
    stopifnot(length(atac) == 1L, file.exists(atac))
    atac <- readRDS(atac)
  }
  stopifnot(inherits(atac, "Seurat"))

  Seurat::DefaultAssay(atac) <- "ATAC"
  if (!(pdx_col %in% colnames(atac[[]]))) stop("No PDX labels found in atac[[pdx_col]].")
  pdx_labels <- as.integer(factor(atac[[pdx_col]][, 1]))

  if (isTRUE(do_preprocess)) {
    atac <- Signac::RunTFIDF(atac)
    atac <- Signac::FindTopFeatures(atac, min.cutoff = "q0")
    atac <- Signac::RunSVD(atac)
  }
  if (!"lsi" %in% Seurat::Reductions(atac)) stop("Missing LSI reduction. Run Step 13 or set do_preprocess=TRUE.")

  avail_lsi <- ncol(Seurat::Embeddings(atac, "lsi"))
  if (avail_lsi < 5) stop("Too few LSI components: ", avail_lsi)

  max_dim <- min(30L, avail_lsi)
  dims_pool <- 2:max_dim

  grid_ndims <- sort(unique(pmin(as.integer(grid_ndims), length(dims_pool))))
  grid_k     <- sort(unique(pmin(as.integer(grid_k), ncol(atac) - 1)))

  results <- vector("list", length(grid_ndims) * length(grid_k) * length(grid_res))
  ii <- 1L

  for (nd in grid_ndims) {
    dims_use <- dims_pool[1:nd]
    for (k in grid_k) {
      for (res in grid_res) {
        set.seed(as.integer(seed))
        tmp <- atac
        tmp <- Seurat::FindNeighbors(tmp, reduction = "lsi", dims = dims_use, k.param = k, verbose = FALSE)
        tmp <- Seurat::FindClusters(tmp, resolution = res, verbose = FALSE)

        clu_labels <- as.integer(factor(Seurat::Idents(tmp)))
        ari <- mclust::adjustedRandIndex(clu_labels, pdx_labels)

        results[[ii]] <- data.frame(
          ndims = nd,
          lsi_dims = paste0(min(dims_use), ":", max(dims_use)),
          k = k,
          resolution = res,
          ARI = ari
        )
        ii <- ii + 1L

        rm(tmp); gc()
      }
    }
  }

  res_df <- do.call(rbind, results)
  res_df <- res_df[order(-res_df$ARI), , drop = FALSE]
  best <- res_df[1, , drop = FALSE]

  f_rds <- file.path(out_dir, "step14_ARI_grid_results.rds")
  f_csv <- file.path(out_dir, "step14_ARI_grid_results.csv")
  f_best_rds <- file.path(out_dir, "best_params_step14.rds")
  f_best_csv <- file.path(out_dir, "best_params_step14.csv")

  saveRDS(res_df, f_rds)
  utils::write.csv(res_df, f_csv, row.names = FALSE)
  saveRDS(best, f_best_rds)
  utils::write.csv(best, f_best_csv, row.names = FALSE)

  if (verbose) {
    message("Step 14 DONE: ", out_dir)
    message("Best: ndims=", best$ndims, " k=", best$k, " res=", best$resolution, " ARI=", signif(best$ARI, 4))
  }

  invisible(list(
    results = res_df,
    best = best,
    out_files = c(results_rds = f_rds, results_csv = f_csv, best_rds = f_best_rds, best_csv = f_best_csv)
  ))
}
