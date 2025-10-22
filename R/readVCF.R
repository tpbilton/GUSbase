#' Read a VCF file and convert to RA object
#'
#' This function reads a VCF (Variant Call Format) file, extracts reference and alternate allele counts for each SNP and individual
#' using a fast Rcpp implementation, and returns an RA class object suitable for downstream analysis in GUSbase.
#'
#' @param vcf_file Path to the VCF file to read.
#' @param infilename Optional filename to use for reporting (defaults to \code{vcf_file}).
#' @param min_ind_with_depth Minimum number of individuals per SNP with non-zero depth required for retention (default = 2).
#'
#' @return An object of class \code{RA} containing genotype calls, allele counts, SNP and sample metadata, and summary information.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' ra_obj <- readVCF("mydata.vcf")
#' }
#'
#' @export

# -- Use the RA class definition from GUSbase/R/RA.R --

readVCF <- function(vcf_file, infilename = NULL, min_ind_depth = 2) {
  # Argument checks
  if (!is.character(vcf_file) || length(vcf_file) != 1) {
    stop("vcf_file must be a single character string specifying the path to the VCF file.")
  }
  if (!file.exists(vcf_file)) {
    stop("The specified vcf_file does not exist: ", vcf_file)
  }
  if (!is.null(infilename) && (!is.character(infilename) || length(infilename) != 1)) {
    stop("infilename must be a single character string if provided (or NULL).")
  }
  if (!is.numeric(min_ind_depth) || length(min_ind_depth) != 1 || min_ind_depth < 1 || min_ind_depth != as.integer(min_ind_depth)) {
    stop("min_ind_depth must be a single positive integer.")
  }

  vcf <- data.table::fread(vcf_file, skip = "#CHROM", header = TRUE, sep = "\t", showProgress = TRUE, data.table = FALSE)
  gt_cols <- colnames(vcf)[10:ncol(vcf)]
  indID <- gt_cols
  nInd <- length(indID)
  chrom <- as.character(vcf$`#CHROM`)
  pos <- as.integer(vcf$POS)
  SNP_Names <- if ("ID" %in% colnames(vcf)) as.character(vcf$ID) else paste0(chrom, "_", pos)
  nSnps <- length(SNP_Names)
  infilename <- if(is.null(infilename)) vcf_file else infilename

  # Get FORMAT and sample fields per SNP (row)
  FORMAT <- matrix(vcf$FORMAT, ncol=1)
  DAT <- as.matrix(vcf[, gt_cols])

  # Allocate allele count matrices
  ref <- matrix(NA_real_, ncol = nSnps, nrow = nInd)
  alt <- matrix(NA_real_, ncol = nSnps, nrow = nInd)

  # Use wrapper now
  extractAlleleCounts_wrapper(FORMAT, DAT, ref, alt)

  # Genon matrix: simple genotype call from counts
  genon = NULL
  #genon <- matrix(NA_integer_, ncol=nSnps, nrow=nInd)
  #genon[ref > 0 & alt == 0] <- 0
  #genon[ref > 0 & alt > 0] <- 1
  #genon[alt > 0 & ref == 0] <- 2
  #genon[ref == 0 & alt == 0] <- NA_integer_

  # --- REMOVE SNPs with only reference or only alternate counts across all individuals ---
  only_ref <- colSums((is.na(alt) | alt == 0)) == nrow(alt)   # all alt missing/0
  only_alt <- colSums((is.na(ref) | ref == 0)) == nrow(ref)   # all ref missing/0
  keep <- !(only_ref | only_alt)

  # --- ADDITIONAL SNP FILTER: retain only SNPs with >= min_ind_with_depth individuals with non-zero depth ---
  depth = ref + alt
  nonzero_depth <- depth > 0 
  ind_with_depth_per_snp <- colSums(nonzero_depth)
  keep <- keep & (ind_with_depth_per_snp >= min_ind_depth)

  ref <- ref[, keep, drop = FALSE]
  alt <- alt[, keep,  drop = FALSE]
  #genon <- genon[, keep, drop = FALSE]
  chrom <- chrom[keep]
  pos <- pos[keep]
  SNP_Names <- SNP_Names[keep]
  nSnps <- sum(keep)
  gform <- "reference"

  # AFrq: mean alt/(ref+alt) per SNP (if ref+alt>0)
  AFrq <- NULL

  # summaryInfo
  summaryInfo = list(
    header = "Data Summary:\n",
    file = paste0("Data file:\t\t", infilename, "\n"),
    meandepth = paste0("Mean Depth:\t\t", round(mean(depth, na.rm=TRUE), 2), "\n"),
    callrate = paste0("Mean Call Rate:\t", round(sum(depth != 0, na.rm=TRUE) / length(depth), 2), "\n"),
    num = "Number of...\n",
    samples = paste0("  Samples:\t\t", nInd, "\n"),
    snps = paste0("  SNPs:\t\t", nSnps, "\n"),
    reads = paste0("  Reads:\t\t", sum(depth, na.rm=TRUE), "\n")
  )

  List <- list(
    genon = genon,
    ref = ref,
    alt = alt,
    chrom = chrom,
    pos = pos,
    SNP_Names = SNP_Names,
    indID = indID,
    nSnps = nSnps,
    nInd = nInd,
    gform = gform,
    AFrq = AFrq,
    infilename = infilename,
    summaryInfo = summaryInfo
  )
  RAobj <- RA$new(List)
  return(RAobj)
}
