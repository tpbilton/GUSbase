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

#' @name $cometPlot
#' @export

cometPlot <- function(ref, alt, ploid=2, freq=NULL, filename=NULL, cex=1, maxdepth=500, maxSNPs=1e5, ...){

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
  ## table for cometplot
  obsTab[upper.tri(obsTab)] <- obsTab[upper.tri(obsTab)] + t(obsTab)[upper.tri(t(obsTab))]

  ## set up the count matrix
  countMat <- matrix(nrow=maxCount+2,ncol=maxCount+2)
  countMat[upper.tri(countMat)] <- log(obsTab[-which(lower.tri(obsTab))])

  rm(list=c("ref_sub","alt_sub"))
  gc()

  ## Compute the expected counts
  if(is.null(freq)){
    p <- matrix(colMeans(ref/(ref+alt),na.rm=T), nrow=nrow(ref), ncol=ncol(ref), byrow=T)
    uni_depth <- 2:maxCount
    expCounts <- sapply(uni_depth, function(x){
      indx <- which(depth == x)
      if(length(indx) == 0) return(numeric(floor(x/2)+1))
      else{
        pvec <- sort(unique(round(p[indx],2)))
        countp <- tabulate(round(p[indx]*100))
        countp <- countp[which(countp > 0)]
        out <- sapply(1:length(pvec), function(z){
          temp <-  countp[z]*sapply(0:x, function(y) sum(dbinom(y, x, prob=0:ploid/ploid)*dbinom(0:ploid, size=ploid, prob=(1-pvec[z]))))
          if((x %% 2) == 0)
            return(temp[0:(x/2)+1] + c(temp[x:(x/2+1)+1],0))
          else
            return(temp[0:floor(x/2)+1] + temp[x:ceiling(x/2)+1])
        })
        return(rowSums(out))
      }
    })
  } else{ # if genotype frequencies are given
    expCounts <- sapply(2:maxCount, function(x){
      indx <- which(depth == x, arr.ind = T)[,2]
      if(length(indx) == 0) return(numeric(floor(x/2)+1))
      else{
        pvec <- matrix(unique(round(p[,unique(indx)]*100), MARGIN=2),nrow=ploid+1)
        pvec <- matrix(pvec[,sort(apply(pvec,2,paste0,collapse="_"), index.return=T)$ix],nrow=ploid+1)
        countp <- table(apply(matrix(round(p[,indx]*100),nrow=ploid+1),2,paste0,collapse="_"))
        out <- sapply(1:ncol(pvec), function(z){
          temp <-  countp[z]*sapply(0:x, function(y) sum(stats::dbinom(y, x, prob=0:ploid/ploid)*pvec[,z]/100))
          if((x %% 2) == 0)
            return(temp[0:(x/2)+1] + c(temp[x:(x/2+1)+1],0))
          else
            return(temp[0:floor(x/2)+1] + temp[x:ceiling(x/2)+1])
        })
        return(rowSums(out))
      }
    })
  }

  expCounts <- c(0,obsTab[1,2],expCounts)
  maxCol = length(utils::tail(expCounts,1)[[1]])
  expCounts_order <- c(list(cbind(unlist(lapply(expCounts, function(x) x[1])), 0:maxCount+1, 1)),
                       sapply(2:maxCol, function(x)
                         cbind(unlist(lapply(expCounts[(x*2-1):(maxCount+1)], function(y) y[x])), (x-1):(maxCount+1-x)+1,x)))

  temp <- do.call("rbind",expCounts_order)
  indx_entry <- temp[,c(2:3)]
  indx_entry[,1] <- indx_entry[,1]+1
  countMat[indx_entry] <- log(round(temp[,1]))

  ### Produce the Comet Plot
  if(!is.null(filename))
    grDevices::png(paste0(filename,".png"), width=max(maxCount*3,640)+sqrt(cex)*maxCount,height=max(maxCount*3,640)+sqrt(cex)*maxCount)
  graphics::par(mar = c(5.1,5.1,4.1,3.1)*sqrt(cex), ...)
  newCol <- grDevices::colorRampPalette(c("red","orange","yellow","green","cyan","blue"))(200)
  graphics::image(0:nrow(countMat),0:ncol(countMat),countMat, xlab="No. Reads (Allele 1)",ylab="No. Reads (Allele 2)",col=newCol,zlim=c(0,max(countMat,na.rm=T)),
                  cex.lab=cex,cex.axis=cex, xlim=c(0,maxCount), ylim=c(0,maxCount), mgp=c(3*sqrt(cex),sqrt(cex),0))
  graphics::abline(0,1)
  if(ploid > 2){
    rat1 = 1:(ceiling(ploid/2)-1)
    rat2 = (ploid-1):(floor(ploid/2)+1)
    junk <- sapply(1:length(rat1), function(x) graphics::abline(0,rat1[x]/rat2[x], lty=3))
    junk <- sapply(1:length(rat1), function(x) graphics::abline(0,rat2[x]/rat1[x], lty=3))
  }
  graphics::mtext("Observed",side=3,cex=cex, font=2, line=1*sqrt(cex))
  graphics::mtext("Expected (Binomial Model)", side=4,cex=cex, font=2,line=1*sqrt(cex))
  legend_image <- grDevices::as.raster(matrix(rev(newCol), ncol = 1))
  adj = maxCount*0.05
  graphics::rasterImage(legend_image,xleft = maxCount-adj*2, ybottom = adj, ytop = 5*adj, xright = maxCount-adj )
  graphics::text(x = maxCount-adj*2, y = adj + seq(0,4*adj,length.out=6), labels = format(exp(seq(0,max(countMat,na.rm=T),length.out=6)),digits=1),cex=cex, pos=2)
  if(!is.null(filename))
    grDevices::dev.off()

  return(invisible())
}
