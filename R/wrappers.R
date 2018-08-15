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
#' @useDynLib GUSbase


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
