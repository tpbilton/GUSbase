##########################################################################
# Genotyping Uncertainty with Sequencing data - Base package (GUSbase)
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

#' Make an unrelated (UR) population
#'
#' Create an UR object from an RA object and perform standard filtering and compute statistics
#' specific to unrelated populations.
#'
#' If \code{mafEst=TRUE}, then the major allele frequency and sequencing error rate for each SNP is estimated based on optimizing the likelihood
#' \deqn{P(Y=a) = \sum_{G} P(Y=a|G)P(G)}
#' where \eqn{P(G)} are genotype probabilities under Hardy Weinberg Equilibrium (HWE) and \eqn{P(Y=a|G)} are the probilities given in Equation (5) of
#' \insertCite{bilton2018genetics2;textual}{GUSbase}. Otherwise, the allele
#' frequencies are taken as the mean of the allele ratio (defined as the number of reference reads divided by the total number of reads) and
#' the sequencing error rate is assumed to be zero.
#'
#' The filtering criteria currently implemented are
#' \itemize{
#' \item{Minor allele frequency (MAF): }{SNPs are discarded if their MAF is less than the threshold (default is 0.01)}
#' \item{Proportion of missing data (MISS): }{SNPs are discarded if the proportion of individuals with no reads
#' (e.g. missing genotype) is greater than the threshold value (default is 0.5)}
#' \item{Hardy Weinberg Distance (HW): }{SNPs are discarded if their Hardy Weinberg distance is less than the first threshold
#' value (default=\code{-0.05}) or if their Hardy Weinberg distance is greater than the second threshold value (default=\code{Inf}).
#' This filtering criteria has been taken from the KGD software (\url{https://github.com/AgResearch/KGD}).}
#' \item{Maximum average SNP read depth (MAXDEPTH): }{SNPs are discarded if the average read depth for the SNP
#' is larger than the threshold (default is 500)}
#' }
#' If \code{filter = NULL}, then no filtering is performed.
#'
#' Estimation of the allele frequencies when \code{mafEst=TRUE} is parallelized using openMP in compiled C code, where the
#' number of threads used in the parallelization is specified by the argument \code{nThreads}.
#'
#' @param RAobj Object of class RA created via the \code{\link{readRA}} function.
#' @param indsubset Integer vector specifying which samples of the RA dataset to retain in the UR
#' population.
#' @param ploid An integer number specifying the ploidy level of the population. Currently, only
#' a ploidy level of two (diploid) is implemented.
#' @param filter Named list of thresholds for various criteria used to fiter SNPs.
#' See below for details.
#' @param mafEst Logical value indicating whether the allele frequences and sequencing
#' error parameters are to estimated for each SNP (see details).
#' @param nThreads Integer vector specifying the number of clusters to use in the foreach loop. Only used in the estimation of
#' allele frequencies when \code{mafEst=TRUE}.
#' @return An R6 object of class UR.
#' @author Timothy P. Bilton and Ken G. Dodds
#' @references
#' \insertRef{bilton2018genetics2}{GUSbase}
#' @export makeUR
#' @examples
#' file <- simDS()
#' RAfile <- VCFtoRA(file$vcf)
#' simdata <- readRA(RAfile)
#'
#' ## make unrelated population
#' urpop <- makeUR(simdata)

#### Make an unrelated population
makeUR <- function(RAobj, ploid = 2, indsubset=NULL, filter=list(MAF=0.01, MISS=0.5, HW=c(-0.05,Inf), MAXDEPTH=500),
                   mafEst=TRUE, nThreads=2){

  ## Do some checks
  if(!all(class(RAobj) %in% c("RA","R6")))
    stop("The `RAobj` argument supplied is not of class 'R6' and 'RA'")
  if(!is.null(filter)){
    if(is.null(filter$MAF)) filter$MAF <- 0.01
    else if( length(filter$MAF) != 1 || !is.numeric(filter$MAF) || filter$MAF<0 || filter$MAF>1)
      stop("Minor allele frequency filter is invalid")
    if(is.null(filter$MISS)) filter$MISS <- 0.5
    else if( length(filter$MISS) != 1 || !is.numeric(filter$MISS) || filter$MISS<0 || filter$MISS>1 )
      stop("Proportion of missing data filter is invalid")
    if(is.null(filter$MAXDEPTH)) filter$MAXDEPTH <- 500
    else if(checkVector(filter$MAXDEPTH, type="pos_numeric", minv=0, equal=FALSE) || length(filter$MAXDEPTH) != 1)
      stop("Maximum mean SNP depth filter is invalid.")
    if(is.null(filter$HW)) filter$HW=c(-0.05,Inf)
    else if(!is.vector(filter$HW) || any(!is.numeric(filter$HW)) || length(filter$HW) != 2 || any(is.na(filter$HW)) || filter$HW[1] >= filter$HW[2])
      stop("Hardy Weinberg (HW) filter is invalid.")
  }
  #if(is.null(filter$PVALUE)) filter$PVALUE <- 1e-6
  #else if( length(filter$PVALUE) != 1 || !is.numeric(filter$PVALUE) || filter$PVALUE<0 || filter$PVALUE>1 )
   #stop("P-value for Hardy-Weinberg equilibrium filter is invalid.")
  if(!is.vector(ploid) || !is.numeric(ploid) || length(ploid) != 1 || ploid != 2) #  round(ploid/2) != ploid/2)
    stop("Argument for ploid level is invalid.")
  if(!is.numeric(nThreads) || length(nThreads) != 1 || nThreads < 0 || round(nThreads) != nThreads)
    stop("Argument for the number of cores for the parallelization is invalid")
  if(is.null(indsubset)) indsubset <- 1:RAobj$.__enclos_env__$private$nInd
  else if(!is.vector(indsubset) || !is.character(indsubset) || any(is.na(indsubset)) ||
          any(!(indsubset %in% RAobj$.__enclos_env__$private$indID)))
    stop("Error in `indsubset` argument. At least one sample name was not found.")
  else indsubset <- match(indsubset, RAobj$.__enclos_env__$private$indID)
  if(checkVector(mafEst, type="one_logical"))
    stop("The `mafEst` argument needs to be a single logical value.")
  #if(checkVector(err, type="one_logical"))
  #  stop("The `err` argument needs to be a single logical value.")

  ## initalize the UR object
  URobj <- UR$new(RAobj, ploid)

  cat("-------------\n")
  cat("Processing Data.\n\n")

  if(!is.null(filter)){
    cat("Filtering criteria for removing SNPs:\n")
    cat("Minor allele frequency (MAF) < ", filter$MAF,"\n",sep="")
    cat("Percentage of missing genotypes > ", filter$MISS*100,"%\n",sep="")
    cat("Maximum average SNP read depth > ", filter$MAXDEPTH, "\n", sep="")
  } else cat("No Filtering performed\n")
  #cat("Hardy-Weinberg equilibrium: < ", filter$HWdis[1]," and > ",filter$HWdis[2],"\n\n",sep="")

  ## Extract the private variables we want
  indsubset <- sort(unique(indsubset))
  indID <- URobj$.__enclos_env__$private$indID[indsubset]
  nSnps <- URobj$.__enclos_env__$private$nSnps
  genon <- URobj$.__enclos_env__$private$genon[indsubset,]

  ## Calculate the proportion of missing data and mean depths
  if(!is.null(filter)){
    miss <- apply(genon,2, function(x) sum(is.na(x))/length(x))
    mdepth <- colMeans(URobj$.__enclos_env__$private$ref[indsubset,] + URobj$.__enclos_env__$private$alt[indsubset,])
    ## Calculate HWE
    naa <- colSums(genon == 2, na.rm = TRUE)
    nab <- colSums(genon == 1, na.rm = TRUE)
    nbb <- colSums(genon == 0, na.rm = TRUE)
    n1 <- 2 * naa + nab
    n2 <- nab + 2 * nbb
    n <- n1 + n2  #n alleles
    p1 <- n1/n
    p2 <- 1 - p1
    HWdis <- naa/(naa + nab + nbb) - p1 * p1
    ## do prelimiary filtering
    snpsubset <- which(miss < filter$MISS & mdepth < filter$MAXDEPTH & HWdis > filter$HW[1] & HWdis < filter$HW[2])
  } else snpsubset <- 1:nSnps

  ## estimate allele frequencies and sequencing error parameters
  if(mafEst){
    temp <- URobj$.__enclos_env__$private$p_est(snpsubset=snpsubset, indsubset=indsubset, nThreads=nThreads)
    pfreq <- temp$p
    ep <- temp$ep
  }
  else{
    ratio <- URobj$.__enclos_env__$private$ref[indsubset,snpsubset]/(URobj$.__enclos_env__$private$ref[indsubset,snpsubset]+URobj$.__enclos_env__$private$alt[indsubset,snpsubset])
    pfreq <- colMeans(ratio, na.rm=T)
    ep <- rep(0, length(snpsubset))
  }
  if(!is.null(filter)){
    ## Compute MAF
    maf <- pmin(pfreq,1-pfreq)
    ## subset SNPs
    maf_indx <- which(maf > filter$MAF)
    indx <- snpsubset[maf_indx]
  } else indx <- maf_indx <- 1:nSnps

  ## check that there are still some SNPs left
  if(length(indx) == 0)
    stop("No SNPs remaining after filtering. Consider using a different set of filtering criteria.")

  ## Update the data in the R6 object
  genon <- genon[,indx]
  ref <- URobj$.__enclos_env__$private$ref[indsubset,indx]
  alt <- URobj$.__enclos_env__$private$alt[indsubset,indx]
  SNP_Names <- URobj$.__enclos_env__$private$SNP_Names[indx]
  nSnps <- length(indx)
  pfreq <- pfreq[maf_indx]
  ep <- ep[maf_indx]
  if(URobj$.__enclos_env__$private$gform == "reference"){
    chrom = URobj$.__enclos_env__$private$chrom[indx]
    pos = URobj$.__enclos_env__$private$pos[indx]
    AFrq = NULL
  }
  else if(URobj$.__enclos_env__$private$gform == "uneak"){
    chrom = pos = NULL
    AFrq = URobj$.__enclos_env__$private$AFrq[indx]
  }

  ## Update the R6 objective
  URobj$.__enclos_env__$private$updatePrivate(list(
    genon = genon, ref = ref, alt = alt, chrom = chrom, pos = pos, nInd = length(indID),
    SNP_Names = SNP_Names, nSnps = nSnps, AFrq = AFrq, pfreq = pfreq, ep = ep)
  )

  return(URobj)
}
