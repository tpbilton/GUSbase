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
cometPlot <- function(ref, alt, model="random", alpha=NULL, filename=NULL, cex=1, maxdepth=500, ...){

  ## Run some tests on the input variables
  if(model == "beta-binomal" && (alpha < 0 | !is.numeric(alpha)) )
    stop("Parameter for the beta-binomal model is not a positive numeric value.")

  ## Remove SNPs with really large read depths
  toRemove = which(ref>maxdepth | alt>maxdepth)
  ref[toRemove] <- 0
  alt[toRemove] <- 0

  ## Compute the observed depth count
  heteCall <- which(ref != 0 & alt != 0)
  hete <- cbind(ref[heteCall], alt[heteCall])
  hete <- cbind(apply(hete,1,min),apply(hete,1,max))
  ## compute the expected depth count
  heteDepth <- rowSums(hete)
  expHete <- matrix(nrow=length(heteDepth), ncol=2)
  ### simulated the expected ratios
  for(i in 1:length(heteDepth)){
    if(model == "beta-binom"){
      newCall <- 0
      while((newCall %% heteDepth[i]) == 0){
        p <- rbeta(1, alpha, alpha)
        newCall <- rbinom(1, heteDepth[i], p)
      }
      expHete[i,] <- sort(c(newCall,heteDepth[i]-newCall))
    }
    else {
      newCall <- c(0,0)
      while(any(newCall == 0))
        newCall <- table(  factor(rbinom(heteDepth[i], size=1, prob=0.5),levels=0:1,labels=c("0","1")))
      expHete[i,] <- sort(newCall)
    }
  }
  maxCount <- max(hete,expHete,50)
  obsHeteTab <- table(factor(hete[,1],levels=1:maxCount,labels=1:maxCount),
                      factor(hete[,2],levels=1:maxCount,labels=1:maxCount))
  obsHeteTab <- log(obsHeteTab)
  expHeteTab <- table(factor(expHete[,1],levels=1:maxCount,labels=1:maxCount),
                      factor(expHete[,2],levels=1:maxCount,labels=1:maxCount))
  expHeteTab <- log(expHeteTab)

  HeteMat <- matrix(nrow=maxCount+1,ncol=maxCount+1)
  HeteMat[lower.tri(HeteMat)] <- t(expHeteTab)[-which(upper.tri(t(expHeteTab)))]
  HeteMat[upper.tri(HeteMat)] <- obsHeteTab[-which(lower.tri(t(expHeteTab)))]
  HeteMat[which(is.infinite(HeteMat))] <- NA

  ### Produce the plot
  if(!is.null(filename))
    png(paste0(filename,".png"), width=max(maxCount*3,640)+sqrt(cex)*maxCount,height=max(maxCount*3,640)+sqrt(cex)*maxCount)
  par(mar = c(5.1,5.1,4.1,3.1)*sqrt(cex), ...)
  newCol <- colorRampPalette(c("red","orange","yellow","green","cyan","blue"))(200)
  image(1:nrow(HeteMat),1:ncol(HeteMat),HeteMat, xlab="No. Reads (Allele 1)",ylab="No. Reads (Allele 2)",col=newCol,zlim=c(0,max(HeteMat,na.rm=T)),
        cex.lab=cex,cex.axis=cex, xlim=c(1,maxCount), ylim=c(1,maxCount), mgp=c(3*sqrt(cex),sqrt(cex),0))
  abline(0,1)
  mtext("Observed",side=3,cex=cex, font=2, line=1*sqrt(cex))
  mtext(paste0("Expected (",ifelse(model=="beta-binom","Beta-Binomial","Random")," model)"), side=4,cex=cex, font=2,line=1*sqrt(cex))
  legend_image <- as.raster(matrix(rev(newCol), ncol = 1))
  adj = maxCount*0.05
  rasterImage(legend_image,xleft = maxCount-adj*2, ybottom = adj, ytop = 5*adj, xright = maxCount-adj )
  text(x = maxCount-adj*2, y = adj + seq(0,4*adj,length.out=6), labels = format(exp(seq(0,max(HeteMat,na.rm=T),length.out=6)),digits=1),cex=cex, pos=2)
  if(!is.null(filename))
    dev.off()
  return(invisible())
}
