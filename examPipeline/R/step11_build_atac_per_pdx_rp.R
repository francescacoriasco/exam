#' @title Step 11 — Create one ATAC object per PDX from fragments restricted to RP loci
#' @description
#' For each PDX, selects the same \code{required_cells} cells used in Step 5 (from \code{barcodes.tsv}),
#' finds the corresponding ATAC fragments file, ensures BGZF + tabix index, harmonizes contig naming,
#' computes a FeatureMatrix on RP loci, and saves one Seurat ATAC object per PDX as RDS.
#' @param data_dir Directory containing fragments files (searched recursively by \code{fragments_finder}).
#' @param barcodes_tsv Path to Step 5 \code{barcodes.tsv} (merged scRNA barcodes).
#' @param rp_loci_rds Path to Step 10 RP loci RDS (GRanges).
#' @param out_dir Output directory for per-PDX ATAC objects.
#' @param required_cells Cells required per PDX (default 500).
#' @param strict If TRUE, stops if not exactly \code{required_cells} are present in fragments.
#' @param barcode_parser Function mapping a merged barcode -> list(pdx=..., bc=...).
#' @param fragments_finder Function(pdx_id, data_dir) -> fragments path (.bgz or .gz).
#' @param ensure_bgzf If TRUE, converts \code{.gz} -> \code{.bgz} using bgzip.
#' @param chunk Chunk size when scanning fragments for existing barcodes.
#' @param verbose If TRUE, prints progress.
#' @return Invisible list with \code{created}, \code{paths}, \code{out_dir}.
#' @export
step11_build_atac_per_pdx_rp <- function(data_dir,
                                         barcodes_tsv,
                                         rp_loci_rds,
                                         out_dir,
                                         required_cells = 500L,
                                         strict = TRUE,
                                         barcode_parser = function(x) {
                                           raw_bc <- sub(".*_", "", x)
                                           m <- regexpr("PDX_([^_]+)", x)
                                           pdx_id <- if (m[1] == -1) NA_character_ else sub("^PDX_", "", regmatches(x, m))
                                           list(pdx = pdx_id, bc = raw_bc)
                                         },
                                         fragments_finder = function(pdx_id, data_dir) {
                                           patt_bgz <- paste0("_PDX_", pdx_id, "_atac_fragments\\.tsv\\.bgz$")
                                           patt_gz  <- paste0("_PDX_", pdx_id, "_atac_fragments\\.tsv\\.gz$")
                                           hit_bgz <- list.files(data_dir, pattern = patt_bgz, full.names = TRUE, recursive = TRUE)
                                           if (length(hit_bgz)) return(hit_bgz[1])
                                           hit_gz <- list.files(data_dir, pattern = patt_gz, full.names = TRUE, recursive = TRUE)
                                           if (length(hit_gz)) return(hit_gz[1])
                                           NA_character_
                                         },
                                         ensure_bgzf = TRUE,
                                         chunk = 200000L,
                                         verbose = TRUE) {
  stopifnot(is.character(data_dir), length(data_dir) == 1L, dir.exists(data_dir))
  stopifnot(is.character(barcodes_tsv), length(barcodes_tsv) == 1L, file.exists(barcodes_tsv))
  stopifnot(is.character(rp_loci_rds), length(rp_loci_rds) == 1L, file.exists(rp_loci_rds))
  stopifnot(is.character(out_dir), length(out_dir) == 1L)
  stopifnot(is.numeric(required_cells), length(required_cells) == 1L, required_cells >= 1)
  stopifnot(is.logical(strict), length(strict) == 1L)
  stopifnot(is.function(barcode_parser))
  stopifnot(is.function(fragments_finder))
  stopifnot(is.logical(ensure_bgzf), length(ensure_bgzf) == 1L)
  stopifnot(is.logical(verbose), length(verbose) == 1L)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  rp_loci <- readRDS(rp_loci_rds)
  stopifnot(inherits(rp_loci, "GRanges"))

  ensure_bgzip_file <- function(frag_path) {
    if (grepl("\\.bgz$", frag_path)) return(frag_path)
    if (!isTRUE(ensure_bgzf)) stop("Fragments are not .bgz and ensure_bgzf=FALSE: ", frag_path)

    bgzip <- Sys.which("bgzip")
    zcat  <- Sys.which("zcat")
    if (zcat == "") zcat <- Sys.which("gzcat")
    if (bgzip == "") stop("bgzip not found (install tabix package in OS).")
    if (zcat  == "") stop("zcat/gzcat not found (install gzip in OS).")

    out_path <- sub("\\.gz$", ".bgz", frag_path)
    if (file.exists(out_path) && file.size(out_path) > 0) return(out_path)

    if (verbose) message("Converting to BGZF (.bgz): ", basename(frag_path), " -> ", basename(out_path))
    cmd <- sprintf("%s %s | %s -c > %s",
                   shQuote(zcat), shQuote(frag_path),
                   shQuote(bgzip), shQuote(out_path))
    status <- system(cmd)
    if (!identical(status, 0L) || !file.exists(out_path) || file.size(out_path) == 0) {
      stop("Failed to create BGZF file: ", out_path)
    }
    out_path
  }

  ensure_tabix <- function(frag_path) {
    tabix <- Sys.which("tabix")
    if (tabix == "") stop("tabix not found (install tabix in OS).")
    tbi <- paste0(frag_path, ".tbi")

    if (file.exists(tbi) && file.mtime(tbi) < file.mtime(frag_path)) file.remove(tbi)
    if (file.exists(tbi)) return(invisible(TRUE))

    out <- system2(tabix, c("-p", "bed", frag_path), stdout = TRUE, stderr = TRUE)
    if (!file.exists(tbi)) stop("Failed to create .tbi for: ", frag_path, "\n", paste(out, collapse = "\n"))
    invisible(TRUE)
  }

  get_contigs <- function(frag_path) {
    tabix <- Sys.which("tabix")
    if (tabix != "") {
      ctg <- system2(tabix, c("-l", frag_path), stdout = TRUE, stderr = TRUE)
      ctg <- ctg[nzchar(ctg)]
      if (length(ctg)) return(ctg)
    }
    con <- gzfile(frag_path, "rt"); on.exit(close(con), add = TRUE)
    x <- readLines(con, n = 50000, warn = FALSE)
    if (!length(x)) return(character())
    chr <- vapply(strsplit(x, "\t", fixed = TRUE), `[[`, character(1), 1)
    unique(chr[nzchar(chr)])
  }

  harmonize_loci <- function(loci, contigs) {
    loci2 <- loci
    cont_has_chr <- any(startsWith(contigs, "chr"))
    loci_has_chr <- any(startsWith(GenomeInfoDb::seqlevels(loci2), "chr"))

    if (cont_has_chr && !loci_has_chr) {
      old <- GenomeInfoDb::seqlevels(loci2)
      loci2 <- GenomeInfoDb::renameSeqlevels(loci2, setNames(paste0("chr", old), old))
    } else if (!cont_has_chr && loci_has_chr) {
      old <- GenomeInfoDb::seqlevels(loci2)
      loci2 <- GenomeInfoDb::renameSeqlevels(loci2, setNames(sub("^chr", "", old), old))
    }

    common <- intersect(GenomeInfoDb::seqlevels(loci2), contigs)
    if (!length(common)) stop("No matching chromosomes between loci and fragments after harmonization.")
    GenomeInfoDb::keepSeqlevels(loci2, common, pruning.mode = "coarse")
  }

  filter_to_cells_present <- function(target_cells, frag_path, chunk_n = 200000L) {
    target_cells <- unique(target_cells)
    found <- character(0)

    con <- gzfile(frag_path, "rt"); on.exit(close(con), add = TRUE)
    repeat {
      x <- readLines(con, n = as.integer(chunk_n), warn = FALSE)
      if (!length(x)) break
      spl <- strsplit(x, "\t", fixed = TRUE)
      bc  <- vapply(spl, function(z) if (length(z) >= 4) z[[4]] else NA_character_, character(1))
      bc  <- unique(bc[!is.na(bc)])
      hits <- intersect(bc, target_cells)
      if (length(hits)) {
        found <- unique(c(found, hits))
        if (length(found) >= length(target_cells)) break
      }
    }
    target_cells[target_cells %in% found]
  }

  merged_bc <- readLines(barcodes_tsv)
  if (!length(merged_bc)) stop("barcodes_tsv is empty: ", barcodes_tsv)

  parsed <- lapply(merged_bc, barcode_parser)
  pdx_vec <- vapply(parsed, `[[`, character(1), "pdx")
  bc_vec  <- vapply(parsed, `[[`, character(1), "bc")

  pdx_levels <- sort(unique(pdx_vec[!is.na(pdx_vec)]))
  if (!length(pdx_levels)) stop("No PDX labels parsed from barcodes.tsv. Check barcode_parser().")
  if (verbose) message("PDX found in Step5 barcodes: ", paste(pdx_levels, collapse = ", "))

  created <- 0L
  paths <- character()

  for (pdx in pdx_levels) {
    if (verbose) message("\n--- PDX: ", pdx, " ---")

    frag_path0 <- fragments_finder(pdx, data_dir)
    if (is.na(frag_path0) || !file.exists(frag_path0)) {
      if (verbose) message("Fragments not found for PDX ", pdx, " -> SKIP")
      next
    }

    frag_path <- ensure_bgzip_file(frag_path0)
    ensure_tabix(frag_path)
    if (verbose) message("Fragments (bgzf): ", frag_path)

    target_cells <- unique(bc_vec[pdx_vec == pdx])
    target_cells <- target_cells[1:min(as.integer(required_cells), length(target_cells))]
    if (verbose) message("Target cells from Step5 = ", length(target_cells))

    contigs <- get_contigs(frag_path)
    rp_use  <- harmonize_loci(rp_loci, contigs)

    cells_use <- filter_to_cells_present(target_cells, frag_path, chunk_n = chunk)
    if (verbose) message("Cells present in ATAC fragments = ", length(cells_use))

    if (isTRUE(strict) && length(cells_use) != as.integer(required_cells)) {
      stop("PDX ", pdx, ": expected ", required_cells, " cells present, got ", length(cells_use), ".")
    }
    if (length(cells_use) == 0L) next

    frags <- Signac::CreateFragmentObject(path = frag_path, cells = cells_use, validate.fragments = FALSE)
    if (verbose) message("Computing FeatureMatrix on RP loci...")
    counts_rp <- Signac::FeatureMatrix(fragments = frags, features = rp_use, cells = cells_use)

    atac_assay <- Signac::CreateChromatinAssay(counts = counts_rp, fragments = frags, genome = NULL)
    atac_obj   <- Seurat::CreateSeuratObject(counts = atac_assay, assay = "ATAC", project = paste0("ATAC_PDX_", pdx))
    atac_obj$pdx <- pdx

    save_path <- file.path(out_dir, paste0("ATAC_RP_only_PDX_", pdx, ".rds"))
    saveRDS(atac_obj, save_path)
    paths <- c(paths, save_path)
    created <- created + 1L

    rm(atac_obj, atac_assay, counts_rp, frags); gc()
  }

  if (verbose) {
    message("\nStep 11 DONE")
    message("ATAC objects created: ", created)
    message("Output folder: ", out_dir)
  }

  invisible(list(created = created, paths = paths, out_dir = out_dir))
}