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

#' UR object
#'
#' Class for storing RA data and associated functions for analyzing of unrelated populations.
#'
#' An UR object is created from the \code{\link{makeUR}} function and contains RA data,
#' various statistics of the dataset that have been computed, and functions (or methods)
#' for analyzing the data. Information in an UR are specific to unrelated populations (or
#' populations with no known relationships).
#' @usage
#' URobj <- makeUR()
#' @format NULL
#' @author Timothy P. Bilton
#' @seealso \code{\link{makeUR}}
#' @name UR
#' @export
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
                p_est = function(snpsubset=NULL,indsubset=NULL, nClust=3, para=NULL, multerr=TRUE, err=TRUE){
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
                  if(!is.numeric(nClust) || length(nClust) != 1 || nClust < 0 || round(nClust) != nClust)
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

                  ploid = private$ploid
                  # Set up the Clusters
                  if(multerr){
                    cl <- makeCluster(nClust)
                    registerDoSNOW(cl)
                    if(err){
                      res <- foreach(snp = iter(snpsubset), .combine="cbind") %dopar% {
                        #parscale <- c(logit(p),logit2(ep))/10
                        #parscale[which(abs(inv.logit(para[-(nSnps+1)]) - 0.5) < 0.0001)] <- 0.0001
                        MLE <- optim(par = c(logit(pinit[snp]), logit2(epinit[snp])), fn=ll_pest, gr=score_pest, method="BFGS",
                                     v=ploid, ref=ref[,snp], alt=alt[,snp], nInd=nInd, nSnps=as.integer(1))
                        return(c(inv.logit(MLE$par[1]), inv.logit2(MLE$par[2]), -MLE$value))
                      }
                    }
                    else{
                      res <- foreach(snp = iter(snpsubset), .combine="cbind") %dopar% {
                        MLE <- optim(par = c(logit(pinit[snp])), fn=ll_pest, gr=score_pest, method="BFGS",
                                     v=ploid, ref=ref[,snp], alt=alt[,snp], nInd=nInd, nSnps=as.integer(1),
                                     seqErr=F, extra=epinit[snp])
                        return(c(inv.logit(MLE$par[1]), epinit[snp], -MLE$value))
                      }
                    }
                    stopCluster(cl)
                  }
                  else{
                    stop("Yet to be implemented")
                  }
                  return(res)
                },
                g_est = function(snpsubset=NULL,indsubset=NULL, nClust=3, para=NULL, multerr=TRUE, err=T){
                  ploid = private$ploid
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
                  if(!is.numeric(nClust) || length(nClust) != 1 || nClust < 0 || round(nClust) != nClust)
                    stop("Argument for the number of cores for the parallelization is invalid")
                  ## inital value
                  if(!is.null(para)){
                    if(!is.list(para))
                      stop("Starting values for the parameters are invalid")
                    if(is.null(para$g)) {
                      ginit <- rep(1/(ploid+1),ploid)
                    }
                    else if(!vector(para$g) || !is.numeric(para$g) || length(para$g) != length(snpsubset) ||
                            any(para$g < 0.01) || any(para$g > 0.99) || sum(para$g) >= 1)
                      stop("Starting values for the genotype frequency parameters are invalid")
                    else ginit <- para$g
                    if(is.null(para$ep)) epinit <- rep(0.01, nSnps)
                    else if(!vector(para$ep) || !is.numeric(para$ep) || length(para$ep) != length(snpsubset) || any(para$ep <= 0) || any(para$ep >= 0.5) )
                      stop("Starting value for the error parameter parameter is invalid")
                    else epinit <- para$ep
                  }
                  else{
                    ginit <- rep(1/(ploid+1),ploid)
                    epinit <- rep(0.01, nSnps)
                  }
                  ## perform the estimation

                  cl <- makeCluster(nClust)
                  registerDoSNOW(cl)
                  # Set up the Clusters
                  if(multerr){
                    if(err){
                      res <- foreach(snp = iter(snpsubset), .combine="cbind") %dopar% {
                        MLE <- optim(par = c(mlogit(ginit), logit2(epinit[snp])), fn=ll_gest, gr=score_pest, method="BFGS",
                                     v=ploid, ref=ref[,snp], alt=alt[,snp], nInd=nInd, nSnps=as.integer(1))
                        return(c(inv.mlogit(MLE$par[1:ploid], 2), inv.logit2(MLE$par[ploid]), -MLE$value))
                      }
                    }
                    else{
                      res <- foreach(snp = iter(snpsubset), .combine="cbind") %dopar% {
                        MLE <- optim(par = c(mlogit(ginit)), fn=ll_gest, gr=score_pest, method="BFGS",
                                     v=ploid, ref=ref[,snp], alt=alt[,snp], nInd=nInd, nSnps=as.integer(1),
                                     seqErr=F, extra=epinit[snp])
                        return(c(inv.mlogit(MLE$par[1:ploid], 2), epinit[snp], -MLE$value))
                      }
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

