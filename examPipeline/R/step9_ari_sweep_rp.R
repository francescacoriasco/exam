#' @title Step 9 — Sweep clustering parameters and score PDX separation with ARI (RP scRNA)
#' @description
#' Sweeps \code{npc}, \code{k.param}, and \code{resolution}; computes ARI between clusters
#' and PDX labels.
#' @param obj Seurat object (output of Step 8) containing PCA and \code{obj[[pdx_col]]}.
#' @param out_dir Output directory to save CSV/RDS.
#' @param pdx_col Metadata column containing PDX labels (default \code{"pdx"}).
#' @param grid_npc Vector of PC counts to test.
#' @param grid_k Vector of neighbor k to test.
#' @param grid_res Vector of resolutions to test.
#' @param seed Seed for reproducibility of graph construction (set each iteration).
#' @param verbose If TRUE, prints summary.
#' @return Invisible list with \code{results}, \code{best}, \code{out_files}.
#' @export
step9_ari_sweep_rp <- function(obj,
                               out_dir,
                               pdx_col = "pdx",
                               grid_npc = c(5, 10, 15, 20),
                               grid_k   = c(10, 20, 30),
                               grid_res = c(0.2, 0.5, 0.8, 1.2),
                               seed = 1L,
                               verbose = TRUE) {
  stopifnot(inherits(obj, "Seurat"))
  stopifnot(is.character(out_dir), length(out_dir) == 1L)
  stopifnot(is.character(pdx_col), length(pdx_col) == 1L)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (!"pca" %in% Seurat::Reductions(obj)) stop("obj has no PCA reduction. Run Step 8 first.")
  if (!(pdx_col %in% colnames(obj[[]]))) stop("Missing PDX labels in obj[[pdx_col]].")

  avail_pcs <- ncol(Seurat::Embeddings(obj, "pca"))
  if (avail_pcs < 2) stop("Too few PCs available: ", avail_pcs)

  grid_npc <- sort(unique(pmin(as.integer(grid_npc), avail_pcs)))
  grid_k   <- sort(unique(pmin(as.integer(grid_k), ncol(obj) - 1)))

  pdx_labels <- as.integer(factor(obj[[pdx_col]][, 1]))

  results <- vector("list", length(grid_npc) * length(grid_k) * length(grid_res))
  ii <- 1L

  for (npc in grid_npc) {
    dims_use <- 1:npc
    for (k in grid_k) {
      for (res in grid_res) {
        set.seed(as.integer(seed))
        tmp <- obj
        tmp <- Seurat::FindNeighbors(tmp, dims = dims_use, k.param = k, verbose = FALSE)
        tmp <- Seurat::FindClusters(tmp, resolution = res, verbose = FALSE)

        clu_labels <- as.integer(factor(Seurat::Idents(tmp)))
        ari <- mclust::adjustedRandIndex(clu_labels, pdx_labels)

        results[[ii]] <- data.frame(npc = npc, k = k, resolution = res, ARI = ari)
        ii <- ii + 1L

        rm(tmp); gc()
      }
    }
  }

  res_df <- do.call(rbind, results)
  res_df <- res_df[order(-res_df$ARI), , drop = FALSE]
  best <- res_df[1, , drop = FALSE]

  f_rds <- file.path(out_dir, "step9_ARI_grid_results.rds")
  f_csv <- file.path(out_dir, "step9_ARI_grid_results.csv")
  f_best_rds <- file.path(out_dir, "best_params_step9.rds")
  f_best_csv <- file.path(out_dir, "best_params_step9.csv")

  saveRDS(res_df, f_rds)
  utils::write.csv(res_df, f_csv, row.names = FALSE)
  saveRDS(best, f_best_rds)
  utils::write.csv(best, f_best_csv, row.names = FALSE)

  if (verbose) {
    message("Step 9 DONE: ", out_dir)
    message("Best: npc=", best$npc, " k=", best$k, " res=", best$resolution, " ARI=", signif(best$ARI, 4))
  }

  invisible(list(
    results = res_df,
    best = best,
    out_files = c(results_rds = f_rds, results_csv = f_csv, best_rds = f_best_rds, best_csv = f_best_csv)
  ))
}