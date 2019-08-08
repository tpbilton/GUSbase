
#' @keywords internal
#' @export p_est_em

p_est_em = function(ref, alt, ploid=2, snpsubset=NULL, indsubset=NULL, nThreads=1, para=NULL, EMpara=NULL){
  ## Do some checks
  if(all(dim(ref) != dim(alt)))
    stop("The matrix of reference alleles has different dimensions to matrix of alternate alleles")
  else{
    nSnps = ncol(ref)
    nInd  = nrow(ref)
  }
  if(is.null(snpsubset)) snpsubset <- 1:nSnps
  else if(checkVector(snpsubset, type="pos_integer", minv=0, maxv=nSnps))
    stop("Index for SNPs is invalid")
  if(is.null(indsubset)) indsubset <- 1:nInd
  else if(checkVector(indsubset, type="pos_integer", minv=0, maxv=nInd))
    stop("Index for individuals is invalid")
  if(is.null(EMpara))
    EMpara <- list()
  else if(!is.list(EMpara))
    stop("The argument `EMpara` for specifying the optimization parameters of the EM algorithm is not a list.")
  if(is.null(EMpara$nIter))
    EMpara$nIter = 1000
  if(is.null(EMpara$delta))
    EMpara$delta = 1e-10
  ref <- ref[indsubset,snpsubset]
  alt <- alt[indsubset,snpsubset]
  ratio <- ref/(ref+alt)
  ## Re comupte the number of snps and individuals
  nSnps <- length(snpsubset)
  nInd <- length(indsubset)
  ## check for large read depths
  badReads <- which(ref > 500 | alt > 500)
  ref[badReads] <- round(ratio[badReads]*500)
  ref <- matrix(as.integer(ref), nrow=nInd, ncol=nSnps)
  alt[badReads] <- round((1-ratio[badReads])*500)
  alt <- matrix(as.integer(alt), nrow=nInd, ncol=nSnps)
  if(checkVector(nThreads, type="pos_integer"))
    stop("Argument `nThreads` which controls the number of threads in the parallelization is invalid.")
  ## inital value
  if(!is.null(para)){
    if(!is.list(para))
      stop("Starting values for the parameters are invalid")
    if(is.null(para$p)) {
      pinit <- colMeans(ratio, na.rm=T)
      pinit[which(pinit > 0.99)] <- 0.99
      pinit[which(pinit < 0.01)] <- 0.01
    }
    else if(!is.vector(para$p) || !is.numeric(para$p) || length(para$p) != length(snpsubset) || any(para$p < 0.01) || any(para$p > 0.99) )
      stop("Starting values for the allele frequency parameters are invalid")
    else pinit <- para$p
    ## check error parameters
    if(is.null(para$ep)) epinit <- rep(0, nSnps)
    else if(!is.vector(para$ep) || !is.numeric(para$ep) || any(para$ep < 0) || any(para$ep >= 0.5) )
      stop("Starting value for the error parameter parameter is invalid")
    else{
      if(length(para$ep) == 1)
        epinit <- rep(para$ep, length(snpsubset))
      else if(length(para$ep) != length(snpsubset))
        stop("The number of error parameters does not equal the number of SNPs")
      else epinit <- para$ep
    }
  }
  else{
    pinit <- colMeans(ratio, na.rm=T)
    pinit[which(pinit > 0.99)] <- 0.99
    pinit[which(pinit < 0.01)] <- 0.01
    epinit <- rep(0, nSnps)
  }

  ## compute the estimates
  res <- .Call("pest_em_c", pinit=pinit, ep=epinit, ploid=ploid,
               ref=ref, alt=alt, nInd=nInd,
               nSnps=nSnps, nThreads=nThreads,
               EMpara=c(as.numeric(EMpara$nIter), as.numeric(EMpara$delta)))
  names(res) <- c("loglik", "p")
  return(res)
}

#' @keywords internal
#' @export g_est_em

g_est_em = function(ref, alt, ploid, snpsubset=NULL,indsubset=NULL, nThreads=1, para=NULL, EMpara=NULL){
  ## Do some checks
  if(all(dim(ref) != dim(alt)))
    stop("The matrix of reference alleles has different dimensions to matrix of alternate alleles")
  else{
    nSnps = ncol(ref)
    nInd  = nrow(ref)
  }
  if(is.null(snpsubset)) snpsubset <- 1:nSnps
  else if(checkVector(snpsubset, type="pos_integer", minv=0, maxv=nSnps))
    stop("Index for SNPs is invalid")
  if(is.null(indsubset)) indsubset <- 1:nInd
  else if(checkVector(indsubset, type="pos_integer", minv=0, maxv=nInd))
    stop("Index for individuals is invalid")
  if(is.null(EMpara))
    EMpara <- list()
  else if(!is.list(EMpara))
    stop("The argument `EMpara` for specifying the optimization parameters of the EM algorithm is not a list.")
  if(is.null(EMpara$nIter))
    EMpara$nIter = 1000
  if(is.null(EMpara$delta))
    EMpara$delta = 1e-10
  ref <- ref[indsubset,snpsubset]
  alt <- alt[indsubset,snpsubset]
  ratio <- ref/(ref+alt)
  ## Re comupte the number of snps and individuals
  nSnps <- length(snpsubset)
  nInd <- length(indsubset)
  ## check for large read depths
  badReads <- which(ref > 500 | alt > 500)
  ref[badReads] <- round(ratio[badReads]*500)
  ref <- matrix(as.integer(ref), nrow=nInd, ncol=nSnps)
  alt[badReads] <- round((1-ratio[badReads])*500)
  alt <- matrix(as.integer(alt), nrow=nInd, ncol=nSnps)
  if(checkVector(nThreads, type="pos_integer"))
    stop("Argument `nThreads` which controls the number of threads in the parallelization is invalid.")
  ## inital value
  if(!is.null(para)){
    if(!is.list(para))
      stop("Starting values for the parameters are invalid")
    if(is.null(para$g)) {
      interval <- seq(1e-5,1-1e-5, length.out=ploid)
      ginit <- as.vector(apply(ratio,2, function(x) {
        temp <- factor(findInterval(x, interval), levels=0:ploid)
        res <- table(temp)/sum(!is.na(temp))
        return(res[-length(res)])
      }))
    }
    else if(!is.vector(para$g) || !is.numeric(para$g) || length(para$g) != length(snpsubset) ||
            any(para$g < 0.01) || any(para$g > 0.99) || sum(para$g) >= 1)
      stop("Starting values for the genotype frequency parameters are invalid")
    else ginit <- para$g
    if(is.null(para$ep)) epinit <- rep(0, nSnps)
    else if(!is.vector(para$ep) || !is.numeric(para$ep) ||any(para$ep < 0) || any(para$ep >= 0.5) )
      stop("Starting value for the error parameter parameter is invalid")
    else{
      if(length(para$ep) == 1)
        epinit <- rep(para$ep, length(snpsubset))
      else if(length(para$ep) != length(snpsubset))
        stop("The number of error parameters does not equal the number of SNPs")
      else epinit <- para$ep
    }
  }
  else{
    interval <- seq(1e-5,1-1e-5, length.out=ploid)
    ginit <- as.vector(apply(ratio,2, function(x) {
      temp <- factor(findInterval(x, interval), levels=0:ploid)
      res <- table(temp)/sum(!is.na(temp))
      return(res[-length(res)])
    }))
    epinit <- rep(0, nSnps)
  }

  ## compute the estimates
  res <- .Call("gest_em_c", ginit=ginit, ep=epinit, ploid=ploid,
               ref=ref, alt=alt, nInd=nInd,
               nSnps=nSnps, nThreads=nThreads,
               EMpara=c(as.numeric(EMpara$nIter), as.numeric(EMpara$delta)))
  res[[2]] <- matrix(res[[2]], nrow=ploid)

  names(res) <- c("loglik", "g")
  return(res)
}
