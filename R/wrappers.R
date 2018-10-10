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
#' @useDynLib GUSbase, .registration = TRUE


## likelihood and score function for the allele and error estimation assuming HWE
ll_pest <- function(para, v=v, ref=ref, alt=alt, nInd=nInd, nSnps=nSnps, seqErr=T, extra=0){
  p = inv.logit(para[1:nSnps])
  if(seqErr)
    ep = inv.logit2(para[nSnps+1])
  else
    ep = extra
  out <- .Call("pest_c", p=p, ep=ep, v=v, ref=ref, alt=alt, nInd=nInd, nSnps=nSnps, PACKAGE="GUSbase")
  if(seqErr)
    assign(".score", -out[[2]], envir = parent.frame(3))
  else
    assign(".score", -out[[2]][-(nSnps+1)], envir = parent.frame(3))
  return(-out[[1]])
}

score_pest <- function(para, ...){
  return(get(".score", envir = parent.frame(3)))
}

## likelihood and score function for the allele and error estimation assuming HWE
# ll_gest <- function(para, v=v, ref=ref, alt=alt, nInd=nInd, nSnps=nSnps, seqErr=T, extra=0){
#   g = matrix(inv.mlogit(para[1:(nSnps*v)], n=v), ncol=nSnps, nrow=v)
#   g = rbind(g,1-colSums(g))
#   if(seqErr)
#     ep = inv.logit2(para[nSnps*v+1])
#   else
#     ep = extra
#   out <- .Call("gest_c", geno=g, ep=ep, v=v, ref=ref, alt=alt, nInd=nInd, nSnps=nSnps)
#   if(seqErr)
#     assign(".score", -out[[2]], envir = parent.frame(3))
#   else
#     assign(".score", -out[[2]][-(nSnps*v+1)], envir = parent.frame(3))
#   return(-out[[1]])
# }

