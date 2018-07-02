#' @export createPop

#### Function for creating a particular population structure
createPop <- function(R6obj, pop = "unrelated", ploid=2, ...){

  # Do some checks
  if(!all(class(R6obj) %in% c("RA","R6")))
    stop("First argumented supplied is not of class 'R6' and 'RA'")
  if(!is.vector(pop) || !is.character(pop) || length(pop) != 1 ||
     all(pop !=  c("full-sib","unrelated")))
     stop("Population specified is invalid. Must be one of 'unrelated' or 'full-sib'")
  if(!is.vector(ploid) || !is.numeric(ploid) || length(ploid) != 1 ||
     ploid < 1 || ploid > 100 || (ploid %% 2) != 0 )
    stop("Ploidy level is invalid.")


  ## Make new R6 object depending on the family type.
  if(pop == "full-sib")
    newObj <- FS$new(R6obj)
  else if(pop == "unrelated")
    newObj <- UR$new(R6obj, ploid)
  else
    stop(paste("Population structure",pop,"has not yet be implemented\n"))

  ## Create the population
  return(makePop(newObj,...))
}



### Generic method for creating a population
makePop <- function(obj, ...){
  UseMethod("makePop")
}

