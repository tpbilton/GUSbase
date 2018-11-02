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
#' Class for storing RA data and associated functions for analysis of unrelated populations.
#'
#' @usage NULL
#' @section Usage:
#' \preformatted{
#' URobj <- makeUR(RAobj, ploid = 2,indsubset = NULL, filter = list(MAF = 0.01,
#'   MISS = 0.5, HW = c(-0.05, Inf), MAXDEPTH = 500), mafEst = TRUE,
#'   nThreads = 2)
#' }
#'
#' @section Details:
#' An UR object is created from the \code{\link{makeUR}} function and contains RA data,
#' various statistics of the dataset that have been computed, and functions (or methods)
#' for analyzing the data. Information in an UR object are specific to unrelated populations (or
#' populations with no known relationships).
#'
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
                p_est = function(snpsubset=NULL,indsubset=NULL, nThreads=2, para=NULL, control=NULL
                                 #, method="EM", newStarts=c(0.001,0.01,0.1,0.2)
                                 ){
                  ## Do some checks
                  if(is.null(snpsubset)) snpsubset <- 1:private$nSnps
                  else if(checkVector(snpsubset, type="pos_integer", minv=1, maxv = private$nSnps))
                    stop("Index for SNPs is invalid")
                  if(is.null(indsubset)) indsubset <- 1:private$nInd
                  else if(checkVector(indsubset, type="pos_integer", minv=1, maxv = private$nInd))
                    stop("Index for individuals is invalid")
                  ref <- private$ref[indsubset,snpsubset]
                  alt <- private$alt[indsubset,snpsubset]
                  ratio <- ref/(ref+alt)
                  nSnps <- length(snpsubset)
                  nInd <- length(indsubset)
                  if(!is.numeric(nThreads) || length(nThreads) != 1 || nThreads < 0 || round(nThreads) != nThreads)
                    stop("Argument for the number of cores for the parallelization is invalid")
                  ## inital value
                  if(!is.null(para)){
                    if(!is.list(para))
                      stop("Starting values for the parameters are invalid")
                    if(is.null(para$p)) {
                        pinit <- colMeans(ratio, na.rm=T)
                        pinit[which(pinit > 0.99)] <- 0.99
                        pinit[which(pinit < 0.01)] <- 0.01
                    }
                    else if(checkVector(para$p, type="pos_numeric", minv=0.01, maxv=0.99, equal=FALSE) || length(para$p) != nSnps)
                      stop("Starting values for the allele frequency parameters are invalid")
                    else pinit <- para$p
                    ## check error parameters
                    if(is.null(para$ep)) epinit <- rep(0.01, nSnps)
                    else if(checkVector(para$ep, type="pos_numeric", minv=0, maxv=0.5, equal=FALSE))
                      stop("Starting value for the error parameter parameter is invalid")
                    else{
                      if(length(para$ep) == 1)
                        epinit <- rep(para$ep, nSnps)
                      else if(length(para$ep) != nSnps)
                        stop("The number of error parameters does not equal the number of SNPs")
                      else epinit <- para$ep
                    }
                  }
                  else{
                    pinit <- colMeans(ratio, na.rm=T)
                    pinit[which(pinit > 0.99)] <- 0.99
                    pinit[which(pinit < 0.01)] <- 0.01
                    epinit <- rep(0.01, nSnps)
                  }
                  ## perform the estimation

                  # if(method == "optim"){
                  #   if(is.null(control))
                  #     control <- list(maxit = 200, reltol=1e-10)
                  #   else if(!is.list(control))
                  #     stop("Argument `control` must be a list object")
                  #   else {
                  #     if(is.null(control$maxit)) control$maxit = 200
                  #     if(is.null(control$reltol)) control$reltol = 1e-10
                  #   }
                  #   ploid <- private$ploid
                  #   cl <- parallel::makeCluster(nClust)
                  #   doParallel::registerDoParallel(cl)
                  #   res_temp <- foreach::foreach(snp = 1:nSnps, .combine="cbind") %dopar% {
                  #     MLE <- stats::optim(par = c(logit(pinit[snp]), logit2(epinit[snp])), fn=ll_pest, gr=score_pest, method="BFGS",
                  #                         v=ploid, ref=ref[,snp], alt=alt[,snp], nInd=nInd, nSnps=as.integer(1), control=control)
                  #     ## Check for badly behaved estimates
                  #     if(MLE$par[2] > logit2(0.48)){
                  #       MLE.list <- vector(mode="list", length=length(newStarts))
                  #       for(i in 1:length(newStarts)){
                  #         MLE.list[[i]] <- stats::optim(par = c(logit(pinit[snp]), logit2(newStarts[i])), fn=ll_pest, gr=score_pest, method="BFGS",
                  #                                       v=ploid, ref=ref[,snp], alt=alt[,snp], nInd=nInd, nSnps=as.integer(1), control=control)
                  #       }
                  #       MLE <- MLE.list[[which.min(lapply(MLE.list, function(x) x$value))]]
                  #     }
                  #     ## Return the MLEs
                  #     return(c(inv.logit(MLE$par[1]), inv.logit2(MLE$par[2]), -MLE$value))
                  #   }
                  #   parallel::stopCluster(cl)
                  #
                  #   res <- list(unname(res_temp[3,]), unname(res_temp[1,]), unname(res_temp[2,]))
                  # } else{
                    if(is.null(control))
                      control <- list(maxit = 200, reltol=1e-10)
                    else if(!is.list(control))
                      stop("Argument `control` must be a list object")
                    else {
                      if(is.null(control$maxit)) control$maxit = 200
                      if(is.null(control$reltol)) control$reltol = 1e-10
                    }
                    res <- .Call("pest_em_c", pinit=pinit, epinit=epinit, ref=ref, alt=alt, nInd=nInd,
                                nSnps=nSnps, nThreads=nThreads, EMpara=c(as.numeric(control$maxit), as.numeric(control$reltol)))
                  # }
                  names(res) <- c("loglik", "p", "ep")
                  return(res)
                }
               )
)

