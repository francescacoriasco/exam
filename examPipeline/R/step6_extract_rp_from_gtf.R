#' @title Step 6 — Extract ribosomal protein (RP) genes from an Ensembl GTF and write RP-only GTF
#' @description
#' Downloads (optionally) an Ensembl GTF, extracts gene_ids for genes whose gene_name starts
#' with \code{rp_prefix} (default \code{"RP"}), saves the gene_id list to RDS, and writes a
#' filtered RP-only GTF (\code{.gtf.gz}).
#' @param refs_dir Output/reference directory.
#' @param gtf_url Optional explicit URL. If NULL, it is constructed from \code{species}, \code{assembly}, \code{release}.
#' @param species Ensembl species folder (default \code{"homo_sapiens"}).
#' @param assembly Assembly string used in filename (default \code{"GRCh38"}).
#' @param release Ensembl release integer (default 115).
#' @param gtf_gz Optional local path to full GTF (.gtf.gz). If NULL, uses \code{refs_dir}.
#' @param rp_prefix Prefix for \code{gene_name} to define RP genes (default \code{"RP"}).
#' @param rp_ids_rds Path to save RP gene_ids RDS (default in \code{refs_dir}).
#' @param gtf_rp_gz Path to write RP-only GTF (default in \code{refs_dir}).
#' @param download_if_missing If TRUE, downloads full GTF if not found.
#' @param chunk Number of lines to read per chunk (streaming; low RAM).
#' @param verbose If TRUE, prints progress.
#' @return Invisible list with \code{rp_gene_ids}, \code{gtf_gz}, \code{rp_ids_rds}, \code{gtf_rp_gz}.
#' @export
step6_extract_rp_from_gtf <- function(refs_dir,
                                      gtf_url = NULL,
                                      species = "homo_sapiens",
                                      assembly = "GRCh38",
                                      release = 115L,
                                      gtf_gz = NULL,
                                      rp_prefix = "RP",
                                      rp_ids_rds = NULL,
                                      gtf_rp_gz = NULL,
                                      download_if_missing = TRUE,
                                      chunk = 20000L,
                                      verbose = TRUE) {
  stopifnot(is.character(refs_dir), length(refs_dir) == 1L)
  dir.create(refs_dir, showWarnings = FALSE, recursive = TRUE)

  rel <- as.integer(release)

  if (is.null(gtf_url)) {
    gtf_url <- sprintf(
      "https://ftp.ensembl.org/pub/release-%d/gtf/%s/%s.%s.%d.gtf.gz",
      rel, species,
      tools::toTitleCase(gsub("_", " ", species)) |> gsub(" ", "_", x = _), # safe-ish
      assembly, rel
    )
    # For homo_sapiens this yields Homo_sapiens.GRCh38.115.gtf.gz (matches your notebook).
    # If you pass a custom gtf_url, this construction is ignored.
  }

  default_gtf_name <- sprintf("Homo_sapiens.%s.%d.gtf.gz", assembly, rel)
  if (is.null(gtf_gz)) gtf_gz <- file.path(refs_dir, default_gtf_name)

  if (is.null(rp_ids_rds)) rp_ids_rds <- file.path(refs_dir, sprintf("rp_gene_ids_release%d.rds", rel))
  if (is.null(gtf_rp_gz))  gtf_rp_gz  <- file.path(refs_dir, sprintf("Homo_sapiens.%s.%d.RP_only.gtf.gz", assembly, rel))

  extract_gene_id   <- function(attr) sub('.*gene_id "([^"]+)".*', "\\1", attr)
  extract_gene_name <- function(attr) sub('.*gene_name "([^"]+)".*', "\\1", attr)

  extract_rp_gene_ids <- function(gtf_gz_in, prefix = "RP", chunk_n = 20000L) {
    con <- gzfile(gtf_gz_in, open = "rt")
    on.exit(close(con), add = TRUE)
    rp_ids <- character()

    repeat {
      x <- readLines(con, n = as.integer(chunk_n), warn = FALSE)
      if (!length(x)) break
      x <- x[!startsWith(x, "#")]
      if (!length(x)) next

      spl <- strsplit(x, "\t", fixed = TRUE)
      is_gene <- vapply(spl, function(z) length(z) >= 9 && z[[3]] == "gene", logical(1))
      if (!any(is_gene)) next

      attr <- vapply(spl[is_gene], `[[`, character(1), 9)
      gnm  <- extract_gene_name(attr)
      ok   <- startsWith(gnm, prefix)

      if (any(ok)) rp_ids <- c(rp_ids, extract_gene_id(attr[ok]))
    }
    unique(rp_ids)
  }

  write_rp_only_gtf <- function(in_gz, out_gz, gene_ids, chunk_n = 20000L) {
    gene_set <- unique(gene_ids)

    con_in  <- gzfile(in_gz,  open = "rt")
    on.exit(close(con_in), add = TRUE)
    con_out <- gzfile(out_gz, open = "wt")
    on.exit(close(con_out), add = TRUE)

    wrote_header <- FALSE

    repeat {
      x <- readLines(con_in, n = as.integer(chunk_n), warn = FALSE)
      if (!length(x)) break

      if (!wrote_header) {
        hdr <- x[startsWith(x, "#")]
        if (length(hdr)) writeLines(hdr, con_out)
        wrote_header <- TRUE
      }

      y <- x[!startsWith(x, "#")]
      if (!length(y)) next

      spl <- strsplit(y, "\t", fixed = TRUE)
      ok9 <- vapply(spl, function(z) length(z) >= 9, logical(1))
      if (!any(ok9)) next

      y2   <- y[ok9]
      spl2 <- spl[ok9]
      attr <- vapply(spl2, `[[`, character(1), 9)

      gid  <- extract_gene_id(attr)
      keep <- gid %in% gene_set
      if (any(keep)) writeLines(y2[keep], con_out)
    }

    invisible(out_gz)
  }

  # 0) Download full GTF if needed
  if (!file.exists(gtf_gz)) {
    if (!isTRUE(download_if_missing)) stop("Missing full GTF: ", gtf_gz, " and download_if_missing=FALSE")
    if (verbose) message("Downloading GTF...")
    utils::download.file(gtf_url, destfile = gtf_gz, mode = "wb", quiet = !verbose)
  }
  stopifnot(file.exists(gtf_gz))
  if (verbose) message("GTF present: ", gtf_gz)

  # 1) RP gene_ids
  if (!file.exists(rp_ids_rds)) {
    rp_gene_ids <- extract_rp_gene_ids(gtf_gz, prefix = rp_prefix, chunk_n = chunk)
    saveRDS(rp_gene_ids, rp_ids_rds)
  } else {
    rp_gene_ids <- readRDS(rp_ids_rds)
  }
  stopifnot(is.character(rp_gene_ids), length(rp_gene_ids) > 0)
  if (verbose) message("RP gene_ids found: ", length(rp_gene_ids))

  # 2) RP-only GTF
  if (!file.exists(gtf_rp_gz)) {
    write_rp_only_gtf(gtf_gz, gtf_rp_gz, rp_gene_ids, chunk_n = chunk)
  }
  stopifnot(file.exists(gtf_rp_gz))

  # quick sample check
  con <- gzfile(gtf_rp_gz, "rt")
  x <- readLines(con, n = 5000, warn = FALSE); close(con)
  x <- x[!startsWith(x, "#")]
  if (length(x)) {
    attr <- vapply(strsplit(x, "\t", fixed = TRUE), `[[`, character(1), 9)
    gnm  <- extract_gene_name(attr)
    stopifnot(all(startsWith(gnm, rp_prefix)))
  }

  if (verbose) {
    message("RP-only GTF written to: ", gtf_rp_gz)
    message("Check OK (sampled lines).")
  }

  invisible(list(
    rp_gene_ids = rp_gene_ids,
    gtf_gz = gtf_gz,
    rp_ids_rds = rp_ids_rds,
    gtf_rp_gz = gtf_rp_gz
  ))
}