


## likelihood and score function for the allele and error estimation assuming HWE
ll_pest <- function(para, v=v, ref=ref, alt=alt, nInd=nInd, nSnps=nSnps){
  p = inv.logit(para[1:nSnps])
  ep = inv.logit2(para[nSnps+1])
  out <- .Call("pest_c", p=p, ep=ep, v=v, ref=ref, alt=alt, nInd=nInd, nSnps=nSnps)
  assign(".score", -out[[2]], envir = parent.frame(3))
  return(-out[[1]])
}

score_pest <- function(para, ...){
  return(get(".score", envir = parent.frame(3)))
}
