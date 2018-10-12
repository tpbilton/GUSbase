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
#' RA Method: Extract variables from an RA Object
#'
#' Method for extracting the private variables of an RA object.
#'
#' @section Usage:
#' \preformatted{
#' RAobj$extractVar(nameList)
#' }
#'
#' @section Arguments:
#' \describe{
#' \item{nameList}{A list of the variable names to be extracted from the RA object.}
#' }
#' @name $extractVar
#' @author Timothy P. Bilton
#' @seealso \code{\link{RA}}
#' @examples
#' vcffile <- simDS()
#' rafile <- VCFtoRA(vcffile$vcf)
#' RAdata <- readRA(rafile)
#'
#' ## extract the depth matrix
#' depthMat <- RAdata$extractVar(list("ref","alt"))
#' str(depthMat)
#'
#' ## extract the chromosome and positions
#' assemblyInfo <- RAdata$extractVar(list("chrom","pos"))
#' str(assemblyInfo)
NULL
