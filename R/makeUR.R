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
#' Create an UR object from an RA object and performs standard filtering and computes statistics specific to unrelated populations.
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
#' \item{Minor allele frequency (MAF): }{SNPs are discarded if their MAF is less than the threshold (default is 0.05)}
#' \item{Proportion of missing data (MISS): }{SNPs are discarded if the proportion of individuals with no reads (e.g. missing genotype) is greater than the threshold value (default is 0.5)}
#' }
#'
#' Estimation of the allele frequencies when \code{mafEst=TRUE} is parallelized using the \code{\link[foreach]{foreach}} function, where the
#' number of cores to use in the parallelization is specified by the argument \code{nClust}. Note: Do not
#' set the number of cores to be more than what is available on your computer (or bad things will happen!).
#'
#' @param RAobj Object of class RA created via the \code{\link{readRA}} function.
#' @param filter Named list of of thresholds for various criteria of fitering SNPs.
#' See below for details.
#' @param ploid An integer number specifying the ploidy level of the population.
#' @param mafEst Logical value indicating whether the allele frequences and sequencing
#' error parameters are to estimated for each SNP (see details).
#' @param nClust Integer vector specifying the number of clusters to use in the foreach loop. Only used in the estimation of
#' allele frequencies when \code{mafEst=TRUE}.
#' @return An R6 object of class UR.
#' @author Timothy P. Bilton
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
makeUR <- function(RAobj, filter=list(MAF=0.05, MISS=0.5), ploid=2, mafEst=TRUE, nClust=3){

  ## Do some checks
  if(!all(class(RAobj) %in% c("RA","R6")))
    stop("First argument supplied is not of class 'R6' and 'RA'")
  if(is.null(filter$MAF)) filter$MAF <- 0.05
  else if( length(filter$MAF) != 1 || !is.numeric(filter$MAF) || filter$MAF<0 || filter$MAF>1)
    stop("Minor allele frequency filter is invalid")
  if(is.null(filter$MISS)) filter$MISS <- 0.2
  else if( length(filter$MISS) != 1 || !is.numeric(filter$MISS) || filter$MISS<0 || filter$MISS>1 )
    stop("Proportion of missing data filter is invalid")
  if(!is.vector(ploid) || !is.numeric(ploid) || length(ploid) != 1 || round(ploid/2) != ploid/2)
    stop("Argument for ploid level is invalid.")
  if(!is.numeric(nClust) || length(nClust) != 1 || nClust < 0 || round(nClust) != nClust)
    stop("Argument for the number of cores for the parallelization is invalid")
  #if(is.null(filter$HWdis)) filter$HWdis <- c(-0.05, 1)
  #else if(!is.vector(filter$HWdis) || length(filter$HWdis) != 2 || !is.numeric(filter$HWdis) || filter$HWdis<0 || filter$HWdis >1)
  #  stop("Hardy Weinberg Equilibrium (HWE) filter is invalid")

  ## initalize the UR object
  URobj <- UR$new(RAobj, ploid)

  cat("-------------\n")
  cat("Processing Data.\n\n")

  cat("Filtering criteria for removing SNPs :\n")
  cat("Minor allele frequency (MAF) < ", filter$MAF,"\n")
  cat("Percentage of missing genotypes > ", filter$MISS*100,"%\n",sep="")
  #cat("Hardy-Weinberg equilibrium: < ", filter$HWdis[1]," and > ",filter$HWdis[2],"\n\n",sep="")

  ## Extract the private variables we want
  indID <- URobj$.__enclos_env__$private$indID
  nSnps <- URobj$.__enclos_env__$private$nSnps
  genon <- URobj$.__enclos_env__$private$genon
  ## Calculate the MAF
  if(mafEst){
    temp <- URobj$.__enclos_env__$private$p_est(nClust=nClust)
    pfreq <- unname(temp[1,])
    ep <- unname(temp[2,])
    ll_HWE <- unname(temp[3,])
  }
  else{
    ratio <- URobj$.__enclos_env__$private$ref/(URobj$.__enclos_env__$private$ref+URobj$.__enclos_env__$private$alt)
    pfreq <- colMeans(ratio, na.rm=T)
    ep <- rep(0, nSnps)
  }

  ## Calculate the proportion of missing data
  miss <- apply(genon,2, function(x) sum(is.na(x))/length(x))

  # ## Compute the HWE distance
  # naa <- colSums(genon == 2, na.rm = TRUE)
  # nab <- colSums(genon == 1, na.rm = TRUE)
  # nbb <- colSums(genon == 0, na.rm = TRUE)
  # n1 <- 2 * naa + nab
  # n2 <- nab + 2 * nbb
  # n <- n1 + n2  #n alleles
  # p1 <- n1/n
  # p2 <- 1 - p1
  # HWdis <- naa/(naa + nab + nbb) - p1 * p1

  ## Indx the filtered SNPs
  maf <- pmin(pfreq,1-pfreq)
  indx <- (maf > filter$MAF) & (miss < filter$MISS) #& (HWdis > filter$HWdis[1]) & (HWdis < filter$HWdis[2])

  ## Update the data in the R6 object
  genon <- genon[,indx]
  ref <- URobj$.__enclos_env__$private$ref[,indx]
  alt <- URobj$.__enclos_env__$private$alt[,indx]
  SNP_Names <- URobj$.__enclos_env__$private$SNP_Names[indx]
  nSnps <- sum(indx)
  pfreq <- pfreq[indx]
  ep <- ep[indx]
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
    genon = genon, ref = ref, alt = alt, chrom = chrom, pos = pos,
    SNP_Names = SNP_Names, nSnps = nSnps, AFrq = AFrq, pfreq = pfreq, ep = ep)
  )

  return(URobj)
}
