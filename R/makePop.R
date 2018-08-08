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

#' Create a population
#'
#' Function which creates a object containing information for specific population
#'
#' @param R6obj An object of classes 'RA' and 'R6' that was created from the \link{readRA} function
#' @param pop A character string specifing the type of population to be created. Currently, only \code{unrelated} and \code{full-sib}
#' are implemented.
#' @param ploid Integer specifying the ploid level in polyploids populations. Only used in \code{unrelated} populations.
#' @param ... Arguments passed to ?
#'
#' @return An \code{R6} object of classes \code{RA} and \code{UR} (for unrelated populations) or \code{FS} (for full-sib populations)
#' @author Timothy P. Bilton.




#' @export createPop

#### Function for creating a particular population structure
createPop <- function(R6obj, pop = "unrelated", ploid=2, ...){

  # Do some checks
  if(!all(class(R6obj) %in% c("RA","R6")))
    stop("First argument supplied is not of class 'R6' and 'RA'")
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

