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

#' @export




#### Make an unrelated population
makePop_UR <- function(RAobj, filter=list(MAF=0.05, MISS=0.2), mafEst=TRUE){

  ## Do some checks
  if(!all(class(RAobj) %in% c("RA","R6")))
    stop("First argument supplied is not of class 'R6' and 'RA'")
  if(is.null(filter$MAF)) filter$MAF <- 0.05
  else if( length(filter$MAF) != 1 || !is.numeric(filter$MAF) || filter$MAF<0 || filter$MAF>1)
    stop("Minor allele frequency filter is invalid")
  if(is.null(filter$MISS)) filter$MISS <- 0.2
  else if( length(filter$MISS) != 1 || !is.numeric(filter$MISS) || filter$MISS<0 || filter$MISS>1 )
    stop("Proportion of missing data filter is invalid")
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
  cat("Hardy-Weinberg equilibrium: < ", filter$HWdis[1]," and > ",filter$HWdis[2],"\n\n",sep="")

  ## Extract the private variables we want
  indID <- URobj$.__enclos_env__$private$indID
  nSnps <- URobj$.__enclos_env__$private$nSnps
  genon <- URobj$.__enclos_env__$private$genon
  ## Calculate the MAF
  if(mafEst){
    temp <- URobj$.__enclos_env__$private$p_est()
    pfreq <- unname(temp[1,])
    ep <- unname(temp[2,])
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
