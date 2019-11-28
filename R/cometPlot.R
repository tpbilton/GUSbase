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

#' #' Produce a cometPlot
#' #'
#' #' Some discription
#' #'
#' #' @usage
#' #' RAobj$comotPlot(model="random", alpha=NULL, filename=NULL, cex=1, maxdepth=500, ...)
#' #'
#' #' @name $cometPlot
#' #'
#' #'
#'
#' @export
cometPlot <- function(ref, alt, ploid=2, gfreq=NULL, file=NULL, cex=1, maxdepth=500, maxSNPs=1e5, res=300, ...){

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
    filename <- paste0(tail(strsplit(file,split=.Platform$file.sep)[[1]],1),"_comet.png")
    if(!file.create(filename,showWarnings = F))
      stop("Unable to create output file.")
  }
  if(GUSbase::checkVector(maxdepth, minv = 2))
    stop("Argument `maxdepth` is invalid")
  if(GUSbase::checkVector(maxSNPs, minv = 2))
    stop("Argument `maxSNPs` is invalid")
  if(GUSbase::checkVector(res, type="pos_numeric", minv = 2))
    stop("Argument `res` is invalid")
  if(GUSbase::checkVector(cex, type="pos_numeric"))
    stop("Argument `cex` is invalid")

  ### check if there are too many SNPs
  if(ncol(ref) > maxSNPs){
    SNPsam <- sample(1:ncol(ref), size=maxSNPs)
    ref <- ref[,SNPsam]
    alt <- alt[,SNPsam]
  }

  ## Remove SNPs with really large read depths
  depth <- ref+alt
  toRemove = which(depth>2*maxdepth)
  ref[toRemove] <- 0
  alt[toRemove] <- 0
  depth[toRemove] <- 0

  ## compute the observed depths
  heteCall <- which(ref > 0 | alt > 0)
  maxCount <- max(100,2*maxdepth)
  ref_sub <- factor(ref[heteCall],levels=0:maxCount,labels=0:maxCount)
  alt_sub <- factor(alt[heteCall],levels=0:maxCount,labels=0:maxCount)
  obsTab <- table(ref_sub,alt_sub)
  ## table for cometplot
  obsTab[upper.tri(obsTab)] <- obsTab[upper.tri(obsTab)] + t(obsTab)[upper.tri(t(obsTab))]

  ## set up the count matrix
  countMat <- matrix(nrow=maxCount+2,ncol=maxCount+2)
  countMat[upper.tri(countMat)] <- log(obsTab[-which(lower.tri(obsTab))])

  rm(list=c("ref_sub","alt_sub"))
  gc()

  ## Compute the expected counts
  if(is.null(gfreq)){ # if genotype frequencies are not specified
    p <- matrix(colMeans(ref/(ref+alt),na.rm=T), nrow=nrow(ref), ncol=ncol(ref), byrow=T)
    uni_depth <- 2:maxCount
    expCounts <- sapply(uni_depth, function(x){
      indx <- which(depth == x)
      if(length(indx) == 0) return(numeric(floor(x/2)+1))
      else{
        pvec <- sort(unique(round(p[indx],2)))
        countp <- tabulate(round(p[indx]*100))
        if(any(pvec == 0)) countp <- c(sum(round(p[indx]*100) == 0), countp)
        countp <- countp[which(countp > 0)]
        out <- as.matrix(sapply(1:length(pvec), function(z){
          temp <-  countp[z]*sapply(0:x, function(y) sum(dbinom(y, x, prob=0:ploid/ploid)*dbinom(0:ploid, size=ploid, prob=(1-pvec[z]))))
          if((x %% 2) == 0)
            return(temp[0:(x/2)+1] + c(temp[x:(x/2+1)+1],0))
          else
            return(temp[0:floor(x/2)+1] + temp[x:ceiling(x/2)+1])
        }))
        return(rowSums(out))
      }
    })
  } else{
    pdist <- unique(gfreq, MARGIN=2)
    indx_p <- sapply(1:ncol(pdist), function(x) which(apply(gfreq,2, function(y) identical(y,pdist[,x]))), simplify = F)
    pvec <- numeric(ncol(ref))
    for(i in 1:length(indx_p))
      pvec[indx_p[[i]]] <- i
    expCounts <- sapply(2:maxCount, function(x){
      indx <- which(depth == x, arr.ind = T)[,2]
      if(length(indx) == 0) return(numeric(x+1))
      else{
        countp <- tabulate(pvec[indx])
        nzp <- which(countp > 0)
        countp <- countp[nzp]
        out <- as.matrix(sapply(1:length(countp), function(z){
          temp <- countp[z]*sapply(0:x, function(y) sum(dbinom(y, x, prob=0:ploid/ploid)*pdist[,z]))
          if((x %% 2) == 0)
            return(temp[0:(x/2)+1] + c(temp[x:(x/2+1)+1],0))
          else
            return(temp[0:floor(x/2)+1] + temp[x:ceiling(x/2)+1])
        }))
        return(rowSums(out))
      }
    })
  }
  expCounts <- c(0,obsTab[1,2],expCounts)
  maxCol = length(tail(expCounts,1)[[1]])
  expCounts_order <- c(list(cbind(unlist(lapply(expCounts, function(x) x[1])), 0:(maxCount)+1, 1)),
                       sapply(2:maxCol, function(x)
                         cbind(unlist(lapply(expCounts[(x*2-1):(maxCount+1)], function(y) y[x])), (x-1):(maxCount+1-x)+1,x)))

  temp <- do.call("rbind",expCounts_order)
  indx_entry <- temp[,c(2:3)]
  indx_entry[,1] <- indx_entry[,1]+1
  countMat[indx_entry] <- log(temp[,1]+1)

  ### Produce the plot
  maxplot <- min(max(ref,alt),maxdepth)
  if(!is.null(file))
    png(filename, width=max(maxplot*3,640)*4+sqrt(cex)*maxplot,height=max(maxplot*3,640)*4+sqrt(cex)*maxplot, res=res)
  par(mar = c(5.1,5.1,5.1,5.1)*sqrt(cex), ...)
  newCol <- colorRampPalette(c("white","red","orange","yellow","green","cyan","blue"))(200)
  grid_add <- function(){
    xaxp <- par("xaxp")
    ticks <- seq(par()$xaxp[1],par()$xaxp[2],par()$xaxp[2]/par()$xaxp[3])
    for(i in ticks){
      lines(x=rep(i-0.5,2),y=c(-2,i-0.5), lty=3, col="grey", cex=cex)
      lines(x=c(i-1.5,maxplot),y=rep(i,2)-1.5, lty=3, col="grey", cex=cex)
      lines(y=rep(i-0.5,2),x=c(-2,i-0.5), lty=3, col="grey", cex=cex)
      lines(y=c(i-1.5,maxplot),x=rep(i,2)-1.5, lty=3, col="grey", cex=cex)
    }
    return(invisible(NULL))
  }
  image(0:nrow(countMat)-2,0:ncol(countMat)-2,countMat, xlab="Major Read Depth",ylab="Major Read Depth",col=newCol,zlim=c(0,max(countMat,na.rm=T)),
        cex.lab=cex,cex.axis=cex, xlim=c(-2,maxplot), ylim=c(-2,maxplot), mgp=c(3*sqrt(cex),sqrt(cex),0),
        xaxt = "n", yaxt="n")
  box()
  grid_add()
  ticks <- seq(par()$xaxp[1],par()$xaxp[2],par()$xaxp[2]/par()$xaxp[3])
  axis(side = 1, at = ticks-0.5,labels = ticks, cex.axis=cex)
  axis(side = 2, at = ticks-0.5,labels = ticks, cex.axis=cex)
  axis(side = 3, at = ticks-1.5,labels = ticks, cex.axis=cex)
  axis(side = 4, at = ticks-1.5,labels = ticks, cex.axis=cex)
  abline(0,1, cex=cex)
  if(ploid > 2){
    rat1 = 1:(ceiling(ploid/2)-1)
    rat2 = (ploid-1):(floor(ploid/2)+1)
    junk <- sapply(1:length(rat1), function(x) abline(-1.5,rat1[x]/rat2[x], lty=3, cex=cex))
    junk <- sapply(1:length(rat1), function(x) abline(rat2[x]/rat1[x]*1.5,rat2[x]/rat1[x], lty=3, cex=cex))
  }
  #mtext("Observed",side=3,cex=cex, font=2, line=1*sqrt(cex))
  #mtext("Expected (Binomial Model)", side=4,cex=cex, font=2,line=1*sqrt(cex))
  mtext("(Observed)\nMinor Read Depth",side=3,cex=cex, font=1, line=2.5*sqrt(cex))
  mtext("Minor Read Depth\n(Expected)", side=4,cex=cex, font=1,line=3.5*sqrt(cex))
  legend_image <- as.raster(matrix(rev(newCol), ncol = 1))
  adj = min(max(ref,alt),maxplot)*0.05
  rasterImage(legend_image,xleft = min(max(ref,alt),maxplot)-adj*2, ybottom = adj, ytop = 5*adj, xright = min(max(ref,alt),maxplot)-adj)
  lines(x = c(rep(min(max(ref,alt),maxplot)-adj*2,2),rep(min(max(ref,alt),maxplot)-adj,2),min(max(ref,alt),maxplot)-adj*2), y=c(adj,rep(5*adj,2),adj,adj),cex=cex)
  text(x = min(max(ref,alt),maxplot)-adj*2-1, y = adj + seq(0,4*adj,length.out=6), labels = format(exp(seq(0,max(countMat,na.rm=T),length.out=6))-1,digits=1, scientific=F),cex=cex, pos=2)
  junk <- sapply(adj + seq(0,4*adj,length.out=6), function(y) lines(min(max(ref,alt),maxplot) - adj*2 + c(0,-1), y=rep(y,2), cex=cex))
  if(!is.null(file))
    dev.off()
  return(invisible())
}
