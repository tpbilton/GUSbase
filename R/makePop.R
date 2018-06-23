#' @export createPop

#### Function for creating a particular population structure
createPop <- function(R6obj, pop = "unrelated", ...){

  ## Make new R6 object depending on the family type.
  if(pop == "full-sib")
    newObj <- FS$new(R6obj)
  else if(pop == "unrelated")
    newObj <- UR$new(R6obj)
  else
    stop(paste("Population structure",pop,"has not yet be implemented\n"))

  ## Create the population
  return(makePop(newObj,...))
}



### Generic method for creating a population
makePop <- function(obj, ...){
  UseMethod("makePop")
}

