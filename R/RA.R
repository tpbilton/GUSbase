##########################################################################
# Genotyping Uncertainty with Sequencing data - Base package (GUSbase)
# Copyright 2017-2020 Timothy P. Bilton <timothy.bilton@agresearch.co.nz>
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
#' RA Object
#'
#' Class for storing reference/alternative (RA) data and methods for performing exporartory data analysis
#' on sequencing data.
#'
#' @usage NULL
#' @section Usage:
#' \preformatted{
#' ## Create RA object
#' RAobj <- readRA(rafile, snpsubset=NULL, sampthres = 0.01, excsamp = NULL)
#'
#' ## RA Functions (Methods)
#' RAobj$extractVar(nameList)
#' RAobj$mergeSamples(samID, useID=TRUE)
#' RAobj$writRA(snpsubset=NULL, indsubset=NULL, file="GUSbase")
#' RAobj$writeVCF(snpsubset=NULL, indsubset=NULL, file="GUSbase", IDuse=NULL)
#' }
#'
#' @section Details:
#' An RA object is returned from the \code{\link{readRA}} function and contains the RA data, various
#' statistics of the dataset that have been computed, and functions (or methods) for analyzing the data.
#'
#' @section Methods (Functions):
#'
#' A list of the methods that are available to an RA object:
#' \describe{
# #'   \item{\code{\link{$cometPlot}}}{Function for create a comet plot}
#'     \item{\code{\link{$extractVar}}}{Extract private variables stored in an RA object}
#'     \item{\code{\link{$mergeSamples}}}{Merge read counts for different samples}
#'     \item{\code{\link{$writeRA}}}{Convert the data in the RA object back to an RA file}
#'     \item{\code{\link{$writeVCF}}}{Convert the data in the RA object back to VCF format}
#' }
#' @format NULL
#' @author Timothy P. Bilton
#' @seealso \code{\link{readRA}}
#' @name RA
#' @export


### R6 class for data aligned to reference assembly
RA <- R6::R6Class("RA",
              public = list(
                initialize = function(List){
                  private$genon      <- List$genon
                  private$ref        <- List$ref
                  private$alt        <- List$alt
                  private$chrom      <- List$chrom
                  private$pos        <- List$pos
                  private$SNP_Names  <- List$SNP_Names
                  private$indID      <- List$indID
                  private$nSnps      <- List$nSnps
                  private$nInd       <- List$nInd
                  private$gform      <- List$gform
                  private$AFrq       <- List$AFrq
                  private$infilename <- List$infilename
                  private$summaryInfo <- List$summaryInfo
                },
                print = function(...){
                  cat(unlist(private$summaryInfo))
                },
                ## Function for merging samples
                mergeSamples = function(samID, useID=TRUE){

                  ## check input
                  if(!is.vector(samID) || !is.character(samID) || length(samID) != private$nInd)
                    stop("Argument `samID` is invalid")
                  if(!is.logical(useID) || !is.vector(useID) || length(useID) != 1 || is.na(useID))
                    stop("Argument `useID` is invalid")

                  ## determine groupings
                  group = as.numeric(factor(samID, levels=unique(samID), labels=unique(samID)))

                  ## check which IDs to use after the merge
                  if(!useID) newID = unique(samID)
                  else newID = private$indID[which(!duplicated(samID))]

                  ## Now merge information
                  private$ref = rowsum(private$ref, group)
                  private$alt = rowsum(private$alt, group)
                  private$genon <- (private$ref > 0) + (private$alt == 0)
                  private$genon[private$ref == 0 & private$alt == 0] <- NA
                  private$indID = newID
                  private$nInd = length(private$indID)

                  ## update summary information
                  private$computeSummary()
                  return(invisible(NULL))
                },
                extractVar = function(nameList){
                  res <- NULL
                  if(!is.character(nameList) || !is.vector(nameList) || any(is.na(nameList)))
                    stop("Argument `nameList` is not a character vector.")
                  for(name in nameList){
                    res <- c(res,list(private[[name]]))
                  }
                  names(res) <- nameList
                  return(res)
                },
                #### Diagonostic plots ####
                # Ratio of alleles for heterozygous genotype calls (observed vs expected)
                cometPlot = function(ploid=2, filename=NULL, cex=1, maxdepth=500, maxSNPs=1e6, res=300, ind=FALSE, ncores=1, ...){
                  cometPlot(private$ref, private$alt, ploid=ploid, file=filename, cex=cex, maxdepth=maxdepth, maxSNPs=maxSNPs, res=res,
                            ind=ind, ncores=ncores, indID=private$indID, ...)
                },
                # Counts of reference and alternate reads scaled by square root of the expected counts
                rocketPlot = function(ploid=2, filename=NULL, cex=1, maxdepth=500, maxSNPs=1e5, res=300, scaled=TRUE, ...){
                  rocketPlot(private$ref, private$alt, ploid=ploid, file=filename, cex=cex, maxdepth=maxdepth, maxSNPs=maxSNPs, res=res, scaled=scaled, ...)
                },
                # Ratio of alleles for heterozygous genotype calls (observed vs expected)
                RDDPlot = function(ploid=2, filename=NULL, maxdepth=500, maxSNPs=1e5, ...){
                  RDDPlot(private$ref, private$alt, ploid=ploid, file=filename, maxdepth=maxdepth, maxSNPs=maxSNPs, ...)
                },
                ###############################
                writeVCF = function(snpsubset=NULL, indsubset=NULL, file="GUSbase", IDuse=NULL){

                  ## Do some checks
                  if(is.null(snpsubset)) snpsubset <- 1:private$nSnps
                  else if(checkVector(snpsubset, type="pos_integer", minv=1, maxv=private$nSnps))
                    stop("Invalid for snpsubset argument.")
                  if(is.null(indsubset)) indsubset <- 1:private$nInd
                  else if(checkVector(indsubset, type="pos_integer", minv=1, maxv=private$nInd))
                    stop("Invalid for indsubset argument.")
                  if(is.null(IDuse)) IDuse <- private$indID[indsubset]
                  else if(!is.vector(IDuse) || length(IDuse) != length(indsubset))
                    stop("Invalid input for alternative individual ID argument.")
                  if(!is.vector(file) || length(file) != 1 || !is.character(file))
                    stop("Invalid input for file name.")
                  outfilename <- paste0(tail(strsplit(file,split=.Platform$file.sep)[[1]],1),".vcf")
                  outpath <- dts(normalizePath("./", winslash=.Platform$file.sep, mustWork=T))

                  ## Subset the data
                  ref <- private$ref[indsubset,snpsubset]
                  alt <- private$alt[indsubset,snpsubset]
                  ## create genotype matrix
                  genon <- 2*(ref > 0 & alt == 0) + (ref > 0 & alt > 0) - (ref == 0 & alt == 0)
                  ## create header information and initize VCF file
                  metaInfo <- paste('##fileformat=VCFv4.3',paste0("##fileDate=",Sys.Date()),"##source=GUSbase",
                                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                                    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele Read Counts">\n',sep="\n")
                  cat(metaInfo, file=outfilename)
                  ## colnames:
                  cat(c("#CHROM", "POS", "ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT", IDuse, "\n"), file=outfilename, append=TRUE, sep="\t")
                  ## set up the matrix to be written
                  out <- matrix(nrow=length(snpsubset),ncol=9+length(indsubset))
                  out[,1] <- private$chrom[snpsubset]
                  out[,2] <- private$pos[snpsubset]
                  out[,3] <- rep(".", length(snpsubset))
                  out[,4] <- rep("C", length(snpsubset))
                  out[,5] <- rep("G", length(snpsubset))
                  out[,6] <- rep(".", length(snpsubset))
                  out[,7] <- rep(".", length(snpsubset))
                  out[,8] <- rep(".", length(snpsubset))
                  out[,9] <- rep("GT:AD", length(snpsubset))
                  ## Create the data part
                  #temp <- options()$scipen
                  #options(scipen=10)  #needed for formating
                  out[,-c(1:9)] <- matrix(paste(sapply(as.vector(genon), function(x) switch(x+2,"./.","1/1","0/1","0/0")),
                                                paste(ref,alt,sep=","), sep=":"), nrow=length(snpsubset), ncol=length(indsubset), byrow=TRUE)
                  #options(scipen=temp)
                  ## write data to vcf file
                  data.table::fwrite(split(t(out), 1:(length(indsubset)+9)), file=outfilename, sep="\t", append=TRUE, nThread = 1)
                  cat("Name of VCF file:    \t",outfilename,"\n")
                  cat("Location of VCF file:\t",outpath,"/\n\n",sep="")
                  return(invisible(NULL))
                },
                writeRA = function(snpsubset=NULL, indsubset=NULL, file="GUSbase"){
                  ## Do some checks
                  if(is.null(snpsubset)) snpsubset <- 1:private$nSnps
                  else if(checkVector(snpsubset, type="pos_integer", minv=1, maxv=private$nSnps))
                    stop("Invalid for snpsubset argument.")
                  if(is.null(indsubset)) indsubset <- 1:private$nInd
                  else if(checkVector(indsubset, type="pos_integer", minv=1, maxv=private$nInd))
                    stop("Invalid for indsubset argument.")
                  if(!is.vector(file) || length(file) != 1 || !is.character(file))
                    stop("Invalid input for file name.")
                  outfilename <- paste0(tail(strsplit(file,split=.Platform$file.sep)[[1]],1),".ra.tab")
                  outpath <- dts(normalizePath("./", winslash=.Platform$file.sep, mustWork=T))

                  ## Subset the data
                  ref <- t(private$ref[indsubset,snpsubset])
                  alt <- t(private$alt[indsubset,snpsubset])

                  out <- matrix(nrow=length(snpsubset),ncol=2+length(indsubset))
                  out[,1] <- private$chrom[snpsubset]
                  out[,2] <- private$pos[snpsubset]
                  out[,-c(1:2)] <- paste(ref,alt, sep = ",")
                  colnames(out) <- c("#CHROM", "POS", private$indID[indsubset])
                  ## write out the data
                  write.table(x = out, file = outfilename, sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
                  cat("Name of RA file:    \t",outfilename,"\n")
                  cat("Location of RA file:\t",outpath,"/\n\n", sep = "")
                  return(invisible(NULL))
                }
              ),
              private = list(
                genon = NULL,
                ref = NULL,
                alt = NULL,
                chrom = NULL,
                pos = NULL,
                SNP_Names = NULL,
                indID = NULL,
                nSnps = NULL,
                nInd = NULL,
                gform = NULL,
                AFrq = NULL,
                infilename = NULL,
                summaryInfo = NULL,
                updatePrivate = function(List){
                  for(elem in names(List)){
                    private[[elem]] <- List[[elem]]
                  }
                },
                computeSummary = function(){
                  temp <- private$ref + private$alt
                  private$summaryInfo = list(
                    header="Data Summary:\n",
                    file="Data file:\t\t",private$infilename,"\n",
                    meandepth="Mean Depth:\t\t", round(mean(temp),2),"\n",
                    callrate="Mean Call Rate:\t",round(sum(temp!=0)/length(temp),2),"\n",
                    num="Number of...\n",
                    samples="  Samples:\t\t",private$nInd,"\n",
                    snps="  SNPs:\t\t",private$nSnps,"\n",
                    reads="  Reads:\t\t",sum(temp),"\n")
                }
              )
)
