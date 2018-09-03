##########################################################################
# Genotyping Uncertainty with Sequencing data (GUSbase)
# Copyright 2017-2018 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################

#' Read an Reference/Alternate (RA) file.
#'
#' Function which processes an RA file into an RA object.
#'
#' RA format is a tab-delimited with columns, CHROM, POS, SAMPLES
#' where SAMPLES consists of sampleIDs, which typically consist of a colon-delimited sampleID, flowcellID, lane, seqlibID.
#' e.g.,
#' \tabular{llll}{
#' CHROM \tab  POS  \tab   999220:C4TWKACXX:7:56 \tab  999204:C4TWKACXX:7:56 \cr
#' 1     \tab  415  \tab   5,0                   \tab  0,3                   \cr
#' 1     \tab  443  \tab   1,0                   \tab  4,4                   \cr
#' 1     \tab  448  \tab   0,0                   \tab  0,2
#' }
#'
#' @param RAfile Character string giving the path to the RA file to be read into R. Typically the required string is
#' returned from the VCFtoRA function when the VCF file is converted to RA format.
# #' @param gform Character string specifying whether the SNPs in the RA data have been called using Uneak (\code{gform="uneak"})
# #' or using an reference based assembly (\code{gform="reference"}).
#' @param sampthres A numeric value giving the filtering threshold for which individual samples are removed. Default is 0.01
#' which means that samples with an average number of reads per SNP that is less than 0.01 are removed.
#' @param excsamp A character vector of the sample IDs that are to be excluded (or discarded). Note that the sample IDs must correspond
#' to those given in the RA file that is to be processed.
#' @return An R6 object of class RA.
#' @author Timothy P. Bilton
#' @examples
#' file <- simDS()
#' RAfile <- VCFtoRA(file$vcf)
#' simdata <- readRA(RAfile)
#' @export
#### Function for reading in RA data and converting to genon and depth matrices.
readRA <- function(genofile, sampthres = 0.01, excsamp = NULL, ...){

  if(!is.character(genofile) || !is.vector(genofile) || length(genofile) != 1)
    stop("File name of RA data set is not a string of length one")
  #if(!is.character(gform) || length(gform) != 1 || !(gform %in% c("reference","uneak")))
  #  stop("gform argument must be either 'reference' or 'uneak'")
  if(missing(gform))
    gform = "reference"

  ## separate character between reference and alternate allele count
  gsep <- switch(gform, uneak = "|", reference = ",")
  ## Process the individuals info
  ghead <- scan(genofile, what = "", nlines = 1, sep = "\t")

  ## Read in the data
  # If reference based
  if (gform == "reference"){
    genosin <- scan(genofile, skip = 1, sep = "\t", what = c(list(chrom = "", coord = 0), rep(list(""), length(ghead) - 2)))
    chrom <- genosin[[1]]
    pos <- genosin[[2]]
    SNP_Names <- paste(genosin[[1]],genosin[[2]],sep="_")
    indID <- ghead[3:length(ghead)]
    AFrq <- NULL
  }
  else if (gform == "uneak"){
    genosin <- scan(genofile, skip = 1, sep = "\t", what = c(list(chrom = ""), rep(list(""), length(ghead) - 6), list(hetc1 = 0, hetc2 = 0, acount1 = 0, acount2 = 0, p = 0)))
    SNP_Names <- genosin[[1]]
    indID <- ghead[2:(length(ghead)-5)]
    AFrq <- genosin[[length(genosin)]]
    chrom <- pos <- NULL
  }

  ## compute dimensions
  nSnps <- length(SNP_Names)
  nInd <- length(ghead) - switch(gform, reference=2, uneak=6)

  ## generate the genon and depth matrices
  ref <- alt <- matrix(as.integer(0), nrow = nInd, ncol = nSnps)
  start.ind <- switch(gform, uneak=1, reference=2)
  for (i in 1:nInd){
    depths <- strsplit(genosin[[start.ind+i]], split = gsep, fixed = TRUE)
    ref[i, ] <- as.integer(unlist(lapply(depths,function(z) z[1])))
    alt[i, ] <- as.integer(unlist(lapply(depths,function(z) z[2])))
  }
  genon <- (ref > 0) + (alt == 0)
  genon[which(ref == 0 & alt == 0)] <- NA

  ## Check that the samples meet the minimum sample treshold
  sampDepth <- rowMeans(ref + alt)
  badSamp <- which(sampDepth < sampthres)
  if(length(badSamp) > 0){
    cat("Samples removed due to having a minimum sample threshold below ",sampthres,":\n",sep="")
    cat(paste0(indID[badSamp],collapse = "\n"),"\n\n")
    excsamp <- unique(c(excsamp,indID[badSamp]))
  }
  ## Remove any sample which we don't want
  if(!is.null(excsamp)){
    toRemove <- which(indID %in% excsamp)
    if(length(excsamp) > 0){
      ref <- ref[-toRemove,]
      alt <- alt[-toRemove,]
      genon <- genon[-toRemove,]
      indID <- indID[-toRemove]
      nInd <- length(indID)
    }
  }

  ## Create the R6 object
  obj <- RA$new(
    list(genon = genon, ref = ref, alt = alt, chrom = chrom, pos = pos,
         SNP_Names = SNP_Names, indID = indID, nSnps = as.integer(nSnps), nInd = as.integer(nInd),
         gform = gform, AFrq = AFrq, infilename = genofile)
  )

  return(obj)
}

