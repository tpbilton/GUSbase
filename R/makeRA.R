##########################################################################
# Genotyping Uncertainty with Sequencing data (GUSbase)
# Copyright 2017-2024 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
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

#' Make an Reference/Alternate (RA) object.
#'
#' Function which creates an RA object from inputs of count matrices of reference and alternate alleles.
#'
#' @param ref Matrix of non-negative counts of reference alleles. Number of rows correspond to the number of samples. 
#' Number of columns correspond to the number of SNPs.
#' @param alt Matrix of non-negative counts of alternate alleles. Number of rows correspond to the number of samples. 
#' Number of columns correspond to the number of SNPs.
#' @param indID Vector of IDs for samples. Needs to have the same length as the number of rows in `ref` and `alt` arguments. 
#' Note: Duplicate IDs will be merged into one sample.
#' @param chrom Vector of chromosome numbers. Needs to have the same length as the number of columns in `ref` and `alt` arguments. 
#' If `NULL`, then each SNP is given a random unique chromsome number. 
#' @param pos Vector of SNP positions. Needs to have the same length as the number of columns in `ref` and `alt` arguments. 
#' If `NULL`, then each SNP is given a random unique positi number. 
#' @param sampthres A numeric value giving the filtering threshold for which individual samples are removed. Default is 0.01
#' which means that samples with an average number of reads per SNP that is less than 0.01 are removed.
#' @return An R6 object of class RA.
#' @author Timothy P. Bilton
#' @examples
#' ## Load data from package
#' datloc = system.file("extdata", "simdata.Rdata", package = "GUSbase")
#' load(datloc)
#' 
#' ## Create RA object from the data
#' makeRA(ref, alt, indID, chrom, pos)
#'
#' @export
#### Function for creating RA object from matrix of reference and alternative allele counts.
makeRA = function(ref, alt, indID, chrom = NULL, pos = NULL, sampthres = 0.01){
  
  ## Check the inputs:
  if(!is.matrix(ref) || !is.numeric(ref) || any(ref < 0) || !all(ref == floor(ref)))
    stop("Incorrect input for argument `ref`. An integer matrix of non-negative values expected.")
  if(!is.matrix(alt) || !is.numeric(alt) || any(alt < 0) || !all(alt == floor(alt)))
    stop("Incorrect input for argument `alt`. An integer matrix of non-negative values expected.")
  if(ncol(ref) != ncol(alt))
    stop("Number of columns in the reference and alternate count allele matrices do not match")
  if(nrow(ref) != nrow(alt))
    stop("Number of rows in the reference and alternate allele count matrices do not match")
  nSnps = ncol(ref)
  nInd = nrow(ref)
  ## Check for missing values:
  if(any(is.na(ref))){
    message("Missing values found in reference allele count matrix. Setting to zero")
    ref[which(is.na(ref))] = 0
  }
  if(any(is.na(alt))){
    message("Missing values found in alternate allele count matrix. Setting to zero")
    ref[which(is.na(alt))] = 0
  }
  
  ## Check individual IDs input
  if(!is.vector(indID) || length(indID) != nInd)
    stop(paste0("Incorrect input for argument `indID`. A vector of length ",nInd," is expected"))
  ## check for duplicate individuals 
  if(any(duplicated(indID))){
    ref = rowsum(ref, indID)
    alt = rowsum(alt, indID)
    message(paste0("Samples with duplicate IDs found. ", sum(duplicated(indID))," samples have been merged to other samples."))
    indID = rownames(ref)
  }
  
  ## Check SNP information
  if(is.null(chrom) & is.null(pos)){
    chrom = paste0("UNK", 1:nSnps)
    pos = 1:nSnps
  } else{
    if(!is.vector(chrom) || length(chrom) != nSnps)
      stop(paste0("Incorrect input for argument `chrom`. A vector of length ",nSnps," is expected"))
    if(any(is.na(chrom)))
      stop(paste0("Input for argument `chrom` contains missing values")) 
    if(!is.vector(pos) || length(pos) != nSnps)
      stop(paste0("Incorrect input for argument `pos`. A vector of length ",nSnps," is expected"))
    if(any(is.na(pos)))
      stop(paste0("Input for argument `pos` contains missing values")) 
  }
  
  ## compute some variables:
  genon <- (ref > 0) + (alt == 0)
  genon[ref == 0 & alt == 0] <- NA
  SNP_Names <- make.unique(paste(chrom, pos, sep="_"))
  
  ## Check that the samples meet the minimum sample threshold
  sampDepth <- rowMeans(ref + alt)
  badSamp <- which(sampDepth < sampthres)
  if(length(badSamp) > 0){
    cat("Samples removed due to having a minimum sample threshold below ",sampthres,":\n",sep="")
    cat(paste0(indID[badSamp],collapse = "\n"),"\n\n")
    excsamp <- unique(c(excsamp,indID[badSamp]))
  }
  
  ## Compute summary information
  temp <- ref + alt
  summaryInfo = list(
    header="Data Summary:\n",
    #file=paste0("Data file:\t\t",rafile,"\n"),
    meandepth=paste0("Mean Depth:\t\t", round(mean(temp),2),"\n"),
    callrate=paste0("Mean Call Rate:\t",round(sum(temp!=0)/length(temp),2),"\n"),
    num="Number of...\n",
    samples=paste0("  Samples:\t\t",nInd,"\n"),
    snps=paste0("  SNPs:\t\t",nSnps,"\n"),
    reads=paste0("  Reads:\t\t",sum(temp),"\n"))
  
  ## Create the R6 object
  obj <- RA$new(
    list(genon = genon, ref = ref, alt = alt, chrom = chrom, pos = pos,
         SNP_Names = SNP_Names, indID = indID, nSnps = as.integer(nSnps), nInd = as.integer(nInd),
         gform = NULL, AFrq = NULL, infilename = NULL, summaryInfo = summaryInfo)
  )
}
