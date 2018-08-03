

## function needed for foreach loop
comb_mat <- function(...){
  mapply('rbind',...,SIMPLIFY=FALSE)
}

comb_vec <- function(...){
  mapply('c',...,SIMPLIFY=FALSE)
}
