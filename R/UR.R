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
#' @export UR



### R6 class for creating a data format for unrelated individuals
UR <- R6Class("UR",
              inherit = RA,
              public = list(
                ## variables
                ## initialize function
                initialize = function(R6obj,ploid){
                  private$genon     <- R6obj$.__enclos_env__$private$genon
                  private$ref       <- R6obj$.__enclos_env__$private$ref
                  private$alt       <- R6obj$.__enclos_env__$private$alt
                  private$chrom     <- R6obj$.__enclos_env__$private$chrom
                  private$pos       <- R6obj$.__enclos_env__$private$pos
                  private$SNP_Names <- R6obj$.__enclos_env__$private$SNP_Names
                  private$indID     <- R6obj$.__enclos_env__$private$indID
                  private$nSnps     <- R6obj$.__enclos_env__$private$nSnps
                  private$nInd      <- R6obj$.__enclos_env__$private$nInd
                  private$gform     <- R6obj$.__enclos_env__$private$gform
                  private$AFrq      <- R6obj$.__enclos_env__$private$AFrq
                  private$infilename<- R6obj$.__enclos_env__$private$infilename
                  private$ploid     <- as.integer(ploid)
                  private$pfreq     <- NULL
                  private$ep        <- NULL
                }
              ),
              private = list(
                ploid = NULL,
                pfreq = NULL,
                ep    = NULL,
                ## Function for updating the private variables
                updatePrivate = function(List){
                  if(!is.null(List$genon))
                    private$genon        = List$genon
                  if(!is.null(List$ref))
                    private$ref          = List$ref
                  if(!is.null(List$alt))
                    private$alt          = List$alt
                  if(!is.null(List$chrom))
                    private$chrom        = List$chrom
                  if(!is.null(List$pos))
                    private$pos          = List$pos
                  if(!is.null(List$SNP_Names))
                    private$SNP_Names    = List$SNP_Names
                  if(!is.null(List$indID))
                    private$indID        = List$indID
                  if(!is.null(List$nSnps))
                    private$nSnps        = List$nSnps
                  if(!is.null(List$nInd))
                    private$nInd         = List$nInd
                  if(!is.null(List$masked))
                    private$masked       = List$masked
                  if(!is.null(List$AFrq))
                    private$AFrq         = List$AFrq
                  if(!is.null(List$pfreq))
                    private$pfreq        = List$pfreq
                  if(!is.null(List$ep))
                    private$ep        = List$ep
                },
                p_est = function(snpsubset=NULL,indsubset=NULL, nClust=3, para=NULL, multerr=TRUE){
                  ## Do some checks
                  if(is.null(snpsubset)) snpsubset <- 1:private$nSnps
                  else if(!is.vector(snpsubset) || !is.numeric(snpsubset) || min(snpsubset) < 0 || max(snpsubset) > private$nSnps)
                    stop("Index for SNPs is invalid")
                  if(is.null(indsubset)) indsubset <- 1:private$nInd
                  else if(!is.vector(indsubset) || !is.numeric(indsubset) || min(indsubset) < 0 || max(indsubset) > private$nInd)
                    stop("Index for individuals is invalid")
                  ref <- private$ref[indsubset,snpsubset]
                  alt <- private$alt[indsubset,snpsubset]
                  ratio <- ref/(ref+alt)
                  nSnps <- length(snpsubset)
                  nInd <- length(indsubset)
                  if(!is.numeric(nClust) || nClust < 0 || round(nClust) != nClust)
                    stop("Argument for the number of cores for the parallelization is invalid")
                  ## inital value
                  if(!is.null(para)){
                    if(!is.list(para) || length(para) != 2)
                      stop("Starting values for the parameters are invalid")
                    if(is.null(para$p)) {
                      pinit <- colMeans(ratio, na.rm=T)
                      pinit[which(pinit > 0.99)] <- 0.99
                      pinit[which(pinit < 0.01)] <- 0.01
                    }
                    else if(!vector(para$p) || !is.numeric(para$p) || length(para$p) != length(snpsubset) || any(para$p < 0.01) || any(para$p > 0.99) )
                      stop("Starting values for the allele frequency parameters are invalid")
                    else pinit <- para$p
                    if(is.null(para$ep)) epinit <- rep(0.01, nSnps)
                    else if(!vector(para$ep) || !is.numeric(para$ep) || length(para$ep) != length(snpsubset) || any(para$ep <= 0) || any(para$ep >= 0.5) )
                      stop("Starting value for the error parameter parameter is invalid")
                    else epinit <- para$ep
                  }
                  else{
                    pinit <- colMeans(ratio, na.rm=T)
                    pinit[which(pinit > 0.99)] <- 0.99
                    pinit[which(pinit < 0.01)] <- 0.01
                    epinit <- rep(0.01, nSnps)
                  }
                  ## perform the estimation

                  # Set up the Clusters
                  if(multerr){
                    ploid = private$ploid
                    cl <- makeCluster(nClust)
                    registerDoSNOW(cl)
                    res <- foreach(snp = iter(snpsubset), .combine="cbind") %dopar% {
                      #parscale <- c(logit(p),logit2(ep))/10
                      #parscale[which(abs(inv.logit(para[-(nSnps+1)]) - 0.5) < 0.0001)] <- 0.0001
                      MLE <- optim(par = c(logit(pinit[snp]), logit2(epinit[snp])), fn=ll_pest, gr=score_pest, method="BFGS",
                                   v=ploid, ref=ref[,snp], alt=alt[,snp], nInd=nInd, nSnps=as.integer(1))
                      return(c(inv.logit(MLE$par[1]), inv.logit2(MLE$par[2])))
                    }
                    stopCluster(cl)
                  }
                  else{
                    stop("Yet to be implemented")
                  }
                  return(res)
                }
              )
)

#### Make an unrelated population
makePop.UR <- function(R6obj, filter=list(MAF=0.05, MISS=0.2), mafEst=TRUE){

  ## Do some checks
  if(is.null(filter$MAF)) filter$MAF <- 0.05
  else if( length(filter$MAF) != 1 || !is.numeric(filter$MAF) || filter$MAF<0 || filter$MAF>1)
    stop("Minor allele frequency filter is invalid")
  if(is.null(filter$MISS)) filter$MISS <- 0.2
  else if( length(filter$MISS) != 1 || !is.numeric(filter$MISS) || filter$MISS<0 || filter$MISS>1 )
    stop("Proportion of missing data filter is invalid")
  #if(is.null(filter$HWdis)) filter$HWdis <- c(-0.05, 1)
  #else if(!is.vector(filter$HWdis) || length(filter$HWdis) != 2 || !is.numeric(filter$HWdis) || filter$HWdis<0 || filter$HWdis >1)
  #  stop("Hardy Weinberg Equilibrium (HWE) filter is invalid")

  cat("-------------\n")
  cat("Processing Data.\n\n")

  cat("Filtering criteria for removing SNPs :\n")
  cat("Minor allele frequency (MAF) < ", filter$MAF,"\n")
  cat("Percentage of missing genotypes > ", filter$MISS*100,"%\n",sep="")
  cat("Hardy-Weinberg equilibrium: < ", filter$HWdis[1]," and > ",filter$HWdis[2],"\n\n",sep="")

  ## Extract the private variables we want
  indID <- R6obj$.__enclos_env__$private$indID
  nSnps <- R6obj$.__enclos_env__$private$nSnps
  genon <- R6obj$.__enclos_env__$private$genon
  ## Calculate the MAF
  if(mafEst){
    temp <- R6obj$.__enclos_env__$private$p_est()
    pfreq <- unname(temp[1,])
    ep <- unname(temp[2,])
  }
  else{
    ratio <- R6obj$.__enclos_env__$private$ref/(R6obj$.__enclos_env__$private$ref+R6obj$.__enclos_env__$private$alt)
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
  ref <- R6obj$.__enclos_env__$private$ref[,indx]
  alt <- R6obj$.__enclos_env__$private$alt[,indx]
  SNP_Names <- R6obj$.__enclos_env__$private$SNP_Names[indx]
  nSnps <- sum(indx)
  pfreq <- pfreq[indx]
  ep <- ep[indx]
  if(R6obj$.__enclos_env__$private$gform == "reference"){
    chrom = R6obj$.__enclos_env__$private$chrom[indx]
    pos = R6obj$.__enclos_env__$private$pos[indx]
    AFrq = NULL
  }
  else if(R6obj$.__enclos_env__$private$gform == "uneak"){
    chrom = pos = NULL
    AFrq = R6obj$.__enclos_env__$private$AFrq[indx]
  }

  ## Update the R6 objective
  R6obj$.__enclos_env__$private$updatePrivate(list(
    genon = genon, ref = ref, alt = alt, chrom = chrom, pos = pos,
    SNP_Names = SNP_Names, nSnps = nSnps, AFrq = AFrq, pfreq = pfreq, ep = ep)
  )

  return(R6obj)
}
