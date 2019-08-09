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

#' Produce a reference ratio density (RRD) plot
#' #'
#' #' Some discription
#' #'
#' @usage
#' RAobj$RRDPlot(model="random", alpha=NULL, filename=NULL, cex=1, maxdepth=500, ...)
#'
#' @name $RRDPlot
#'
#'
#'
#' @export
RRDPlot <- function(ref, alt, ploid=2, gfreq=NULL, file=NULL, maxdepth=500, maxSNPs=1e5, ...){

  ## Do some checks
  if(!is.matrix(ref) || GUSbase::checkVector(as.vector(ref)))
    stop("Matrix of reference allele counts is invalid")
  if(!is.matrix(alt) || GUSbase::checkVector(as.vector(alt)))
    stop("Matrix of alternate allele counts is invalid")
  if(GUSbase::checkVector(ploid, minv=2))
    stop("Matrix of alternate allele counts is invalid")
  if (!is.null(gfreq)) {
    if (!is.matrix(gfreq) || !is.numeric(gfreq) || any(is.na(gfreq)) ||
        any(gfreq < 0) || any(gfreq > 1) || nrow(gfreq) !=
        (ploid + 1) || ncol(gfreq) != ncol(gfreq) ||
        any(sapply(colSums(gfreq), function(x) isTRUE(!all.equal(x, 1, tolerance = 1e-04)))))
      stop("Matrix of genotype frequencies is invalid")
  }
  if(!is.null(file)){
    if(!is.vector(file) || !is.character(file) || length(file) != 1)
      stop("Filename input is invalid")
    filename <- paste0(tail(strsplit(file,split=.Platform$file.sep)[[1]],1),"_RRD.png")
    outfile = file.path(outpath,filename)
    if(!file.create(outfile,showWarnings = F))
      stop("Unable to create output file.")
  }
  if(GUSbase::checkVector(maxdepth, minv = 2))
    stop("Argument `maxdepth` is invalid")
  if(GUSbase::checkVector(maxSNPs, minv = 2))
    stop("Argument `maxSNPs` is invalid")
  if(GUSbase::checkVector(res, type="pos_numeric", minv = 2))
    stop("Argument `res` is invalid")
  if(GUSbase::checkVector(width, type="pos_numeric"))
    stop("Argument `width` is invalid")
  if(GUSbase::checkVector(height, type="pos_numeric"))
    stop("Argument `height` is invalid")

  ### check if there are too many SNPs
  if(ncol(ref) > maxSNPs){
    SNPsam <- sample(1:ncol(ref), size=maxSNPs)
    ref <- ref[,SNPsam]
    alt <- alt[,SNPsam]
  }

  ## Remove SNPs with really large read depths
  depth <- ref+alt
  toRemove = which(depth>maxdepth)
  ref[toRemove] <- 0
  alt[toRemove] <- 0
  depth[toRemove] <- 0

  ## compute the observed depths
  heteCall <- which(ref > 0 | alt > 0)
  maxCount <- min(max(50,max(depth)), max(ref,alt)+2)
  ref_sub <- factor(ref[heteCall],levels=0:maxCount,labels=0:maxCount)
  alt_sub <- factor(alt[heteCall],levels=0:maxCount,labels=0:maxCount)
  obsTab <- table(ref_sub,alt_sub)
  ## varibles for histrogram
  counts_obs <- as.vector(obsTab)[-1]
  counts_obs <- counts_obs/sum(counts_obs)
  temp <- seq(maxCount+1)-1
  tempref <- matrix(temp, nrow=maxCount+1,ncol=maxCount+1, byrow=F)
  tempd <- matrix(temp, nrow=maxCount+1,ncol=maxCount+1, byrow=T) + tempref
  ratio_obs <- tempref/tempd
  ratio_obs <- as.vector(ratio_obs)[-1]

  rm(list=c("ref_sub","alt_sub"))
  gc()

  ## Compute the expected counts
  if(is.null(gfreq)){ # if genotype frequencies are not specified
    p <- matrix(colMeans(ref/(ref+alt),na.rm=T), nrow=nrow(ref), ncol=ncol(ref), byrow=T)
    ## data for producing density plots
    expCounts <- sapply(1:maxCount, function(x){
      indx <- which(depth == x)
      if(length(indx) == 0) return(numeric(x+1))
      else{
        pvec <- sort(unique(round(p[indx],2)))
        countp <- tabulate(round(p[indx]*100))
        countp <- countp[which(countp > 0)]
        out <- sapply(1:length(pvec), function(z){
          return(countp[z]*sapply(0:x, function(y) sum(dbinom(y, x, prob=0:ploid/ploid)*dbinom(0:ploid, size=ploid, prob=pvec[z]))))
        })
        return(rev(rowSums(out)))
      }
    })
  } else{
    pdist <- unique(gfreq, MARGIN=2)
    indx_p <- sapply(1:ncol(pdist), function(x) which(apply(gfreq,2, function(y) identical(y,pdist[,x]))), simplify = F)
    pvec <- numeric(ncol(ref))
    for(i in 1:length(indx_p))
      pvec[indx_p[[i]]] <- i
    expCounts <- sapply(1:maxCount, function(x){
      indx <- which(depth == x, arr.ind = T)[,2]
      if(length(indx) == 0) return(numeric(x+1))
      else{
        countp <- tabulate(pvec[indx])
        nzp <- which(countp > 0)
        countp <- countp[nzp]
        out <- sapply(nzp, function(z){
          return(countp[z]*sapply(0:x, function(y) sum(dbinom(y, x, prob=0:ploid/ploid)*pdist[,z])))
        })
        return(rev(rowSums(out)))
      }
    })
  }

  expCounts <- c(0,expCounts)
  tempMat <- matrix(nrow=maxCount+1, ncol=maxCount+1)
  expCounts_order <-  sapply(1:(maxCount+1), function(x)
    cbind(unlist(lapply(expCounts[x:(maxCount+1)], function(y) y[x])), 0:(maxCount+1-x)+1,x))
  temp <- do.call("rbind",expCounts_order)
  indx_entry <- temp[,c(2:3)]
  tempMat[indx_entry] <- temp[,1]
  counts_exp <- as.vector(tempMat)[-1]
  counts_exp[which(is.na(counts_exp))] <- 0
  counts_exp <- counts_exp/sum(counts_exp)
  temp_df <- data.frame(counts=c(counts_exp,counts_obs),ratio=rep(ratio_obs,2)*ploid,
                        method=rep(c("Expected","Observed"),rep(length(ratio_obs),2)))

  ## Produce RRD plot
  pp <- ggplot2::ggplot(temp_df, ggplot2::aes(x=ratio, fill=method)) +
    ggplot2::geom_density(ggplot2::aes(weight=counts), alpha=0.2) + ggplot2::theme_bw() +
    ggplot2::xlab("Dosage of reference allele") + ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position="top")
  if(is.null(file)) print(pp)
  else  ggplot2::ggsave(filename=paste0(file,"_RRD.png"), plot=pp,...)
  return(invisible(NULL))
}
