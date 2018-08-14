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

#' Transformation functions.
#'
#' Functions used in the GUS universe to transform parameters estimation onto an unbounded
#' interval. These functions are used in likelihood computations optimised using direct
#' maximumization routines.
#'
#' Parameter transformations used in the GUS universe:
#'
#' 1. The logit transformation:
#' \deqn{logit(p) = log(\frac{p}{1-p})}
#' 2. A modified logit transformation (called logit2):
#' \deqn{logit2(p) = log(\frac{2p}{1-2p}).}
#'
#' @param p Probability value
#' @param logit_p Logit probability value
#' @name transformations
#' @export logit
#' @export logit2
#' @export inv.logit
#' @export inv.logit2

## Functions required for transforming the recombination fraction parameters on the interval [0,1]
## to the interval [-inf,inf]
#' @rdname transformations
logit <- function(p) qlogis(p)
#' @rdname transformations
inv.logit <- function(logit_p) plogis(logit_p)

## Functions requred for transforming the recombination fraction parameters on the interval [0,0.5]
## to the interval [-inf,inf]
#' @rdname transformations
logit2 <- function(p) log(2*p/(1-2*p))
#' @rdname transformations
inv.logit2 <- function(logit_p) {
  p <- 1/(2*(1+1/exp(logit_p)))
  p.na <- is.na(p)
  if(sum(p.na)!=0)
    p[which(p.na)] <- 0.5
  return(p)
}



# Alternatives
#logit3 <- function(p) qlogis(2*p)
#inv.logit3 <- function(logit_p) plogis(logit_p)/2
