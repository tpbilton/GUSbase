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

#' @name $RRDPlot
#' @export
RRDPlot <- function(ref, alt, ploid=2, freq=NULL, filename=NULL, maxdepth=500, maxSNPs=1e5, ...){

  ## Do some checks
  if(!is.null(freq)){
    if(!is.matrix(freq) || !is.numeric(freq) || any(is.na(freq)) || any(freq < 0) || any(freq > 1) ||
       nrow(freq) != (ploid+1) || ncol(freq) != ncol(freq) ||
       any(sapply(colSums(freq), function(x) isTRUE(!all.equal(x,1,tolerance=1e-4)))))
      stop("Vector of allele frequencies (argument freq) is invalid.")
    else
      p <- freq
  }

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
  maxCount <- min(max(50,max(depth)),max(ref,alt))
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

  ## set up the count matrix
  countMat <- matrix(nrow=maxCount+2,ncol=maxCount+2)
  countMat[upper.tri(countMat)] <- log(obsTab[-which(lower.tri(obsTab))])

  rm(list=c("ref_sub","alt_sub"))
  gc()

  ## Compute the expected counts
  if(is.null(freq)){ # allele frequencies assuming hardy-weinberg population
    p <- matrix(colMeans(ref/(ref+alt),na.rm=T), nrow=nrow(ref), ncol=ncol(ref), byrow=T)
    expCounts <- sapply(1:maxCount, function(x){
      indx <- which(depth == x)
      if(length(indx) == 0) return(numeric(x+1))
      else{
        pvec <- sort(unique(round(p[indx],2)))
        countp <- tabulate(round(p[indx]*100))
        countp <- countp[which(countp > 0)]
        out <- sapply(1:length(pvec), function(z){
          temp <-  countp[z]*sapply(0:x, function(y) sum(stats::dbinom(y, x, prob=0:ploid/ploid)*stats::dbinom(0:ploid, size=ploid, prob=(1-pvec[z]))))
        })
        return(rowSums(out))
      }
    })
  } else{ # if genotype frequencies are given
    expCounts <- sapply(1:maxCount, function(x){
      indx <- which(depth == x, arr.ind = T)[,2]
      if(length(indx) == 0) return(numeric(x+1))
      else{
        pvec <- matrix(unique(round(p[,unique(indx)]*100), MARGIN=2),nrow=ploid+1)
        pvec <- matrix(pvec[,sort(apply(pvec,2,paste0,collapse="_"), index.return=T)$ix],nrow=ploid+1)
        countp <- table(apply(matrix(round(p[,indx]*100),nrow=ploid+1),2,paste0,collapse="_"))
        out <- sapply(1:ncol(pvec), function(z){
          temp <-  countp[z]*sapply(0:x, function(y) sum(stats::dbinom(y, x, prob=0:ploid/ploid)*pvec[,z]/100))
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

  ## Produce histrograms
  p <- ggplot2::ggplot(temp_df, ggplot2::aes(x=ratio, fill=method)) +
    ggplot2::geom_density(ggplot2::aes(weight=counts), alpha=0.2) + ggplot2::theme_bw() +
    ggplot2::xlab("Dosage of reference allele") + ggplot2::theme(legend.title = ggplot2::element_blank(),legend.position="top")
  if(!is.null(filename)) ggsave(paste0(filename,".png"), ...)
  else print(p)

  return(invisible(p))
}

