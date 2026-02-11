#' @title Step 12 — Merge per-PDX ATAC objects into a single combined dataset
#' @description
#' Loads per-PDX ATAC objects from Step 11, prefixes cell names for uniqueness, intersects
#' shared features, column-binds counts, and constructs a single combined ATAC Seurat object.
#' @param in_dir Directory containing Step 11 per-PDX \code{.rds}.
#' @param out_rds Output \code{.rds} path for the combined object.
#' @param genome Genome string for \code{CreateChromatinAssay} (can be NULL).
#' @param verbose If TRUE, prints summary.
#' @return Invisible list with \code{atac_combined} and \code{out_rds}.
#' @export
step12_merge_atac_objects <- function(in_dir,
                                      out_rds,
                                      genome = "GRCh38",
                                      verbose = TRUE) {
  stopifnot(is.character(in_dir), length(in_dir) == 1L, dir.exists(in_dir))
  stopifnot(is.character(out_rds), length(out_rds) == 1L)
  stopifnot(is.logical(verbose), length(verbose) == 1L)

  rds_files <- list.files(in_dir, pattern = "^ATAC_RP_only_.*\\.rds$", full.names = TRUE)
  if (length(rds_files) < 2) stop("Need >=2 ATAC objects in: ", in_dir)

  objs <- lapply(rds_files, readRDS)

  get_counts <- function(obj) {
  Seurat::DefaultAssay(obj) <- "ATAC"

  # SeuratObject v5: layer=
  m <- tryCatch(
    Seurat::GetAssayData(obj, assay = "ATAC", layer = "counts"),
    error = function(e) NULL
  )
  if (!is.null(m)) return(m)

  # Fallback SeuratObject < 5: slot=
  Seurat::GetAssayData(obj, assay = "ATAC", slot = "counts")
}


  counts_list <- vector("list", length(objs))
  pdx_list    <- vector("list", length(objs))

  for (i in seq_along(objs)) {
    obj <- objs[[i]]

    pdx <- NA_character_
    md <- obj@meta.data
    if ("pdx" %in% colnames(md)) pdx <- unique(md$pdx)[1]
    if ("pdx_scRNA" %in% colnames(md)) pdx <- unique(md$pdx_scRNA)[1]
    if (is.na(pdx) || pdx == "") {
      pdx <- sub("^ATAC_RP_only_PDX_", "", sub("\\.rds$", "", basename(rds_files[i])))
    }

    m <- get_counts(obj)

    new_cells <- paste0(pdx, "_", colnames(m))
    colnames(m) <- new_cells

    counts_list[[i]] <- m
    pdx_list[[i]] <- data.frame(cell = new_cells, pdx = pdx, stringsAsFactors = FALSE)

    rm(obj, m); gc()
  }

  common_features <- Reduce(intersect, lapply(counts_list, rownames))
  if (!length(common_features)) stop("No shared features across ATAC objects.")
  counts_list <- lapply(counts_list, function(m) m[common_features, , drop = FALSE])

  counts_merged <- do.call(cbind, counts_list)

  assay_merged <- Signac::CreateChromatinAssay(counts = counts_merged, genome = genome)
  atac_combined <- Seurat::CreateSeuratObject(counts = assay_merged, assay = "ATAC", project = "ATAC_RP_combined")

  pdx_df <- do.call(rbind, pdx_list)
  pdx_df <- pdx_df[match(colnames(atac_combined), pdx_df$cell), ]
  atac_combined$pdx <- pdx_df$pdx

  dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
  saveRDS(atac_combined, out_rds)

  if (verbose) {
    message("Step 12 DONE")
    message("Features x Cells: ", nrow(atac_combined), " x ", ncol(atac_combined))
    message("Saved: ", out_rds)
  }

  invisible(list(atac_combined = atac_combined, out_rds = out_rds))
}