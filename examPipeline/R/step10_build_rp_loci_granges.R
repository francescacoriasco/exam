#' @title Step 10 — Build GRanges of RP gene loci from RP-only GTF
#' @description
#' Imports an RP-only GTF (\code{.gtf.gz}) and returns loci for \code{type == "gene"} as \code{GRanges}.
#' Optionally deduplicates by \code{gene_id} and saves to RDS.
#' @param gtf_rp_gz Path to RP-only GTF (\code{.gtf.gz}).
#' @param out_rds Optional output RDS path (if NULL, saved alongside the GTF).
#' @param unique_by_gene_id If TRUE, keeps one locus per \code{gene_id}.
#' @param verbose If TRUE, prints summary.
#' @return Invisible list with \code{gr_gene} and \code{out_rds}.
#' @export
step10_build_rp_loci_granges <- function(gtf_rp_gz,
                                         out_rds = NULL,
                                         unique_by_gene_id = TRUE,
                                         verbose = TRUE) {
  stopifnot(is.character(gtf_rp_gz), length(gtf_rp_gz) == 1L, file.exists(gtf_rp_gz))
  stopifnot(is.logical(unique_by_gene_id), length(unique_by_gene_id) == 1L)
  stopifnot(is.logical(verbose), length(verbose) == 1L)

  if (is.null(out_rds)) {
    out_rds <- file.path(dirname(gtf_rp_gz), "rp_gene_loci_GRanges.rds")
  }

  gr_all  <- rtracklayer::import(gtf_rp_gz)
  gr_gene <- gr_all[gr_all$type == "gene"]

  if (isTRUE(unique_by_gene_id) && "gene_id" %in% names(S4Vectors::mcols(gr_gene))) {
    gr_gene <- gr_gene[!duplicated(S4Vectors::mcols(gr_gene)$gene_id)]
  }

  saveRDS(gr_gene, out_rds)

  if (verbose) {
    message("Step 10 DONE")
    message("RP gene loci: ", length(gr_gene))
    message("Saved: ", out_rds)
  }

  invisible(list(gr_gene = gr_gene, out_rds = out_rds))
}