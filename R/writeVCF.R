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
#' RA Method: Write genotype data to VCF format
#'
#' Method for converting allele count data in an RA object back to VCF format.
#'
#' @section Usage:
#' \preformatted{
#' RAobj$writeVCF(snpsubset=NULL, indsubset=NULL, file="GUSbase", IDuse=NULL)
#' }
#'
#' @section Arguments:
#' \describe{
#' \item{snpsubset}{Integer vector giving the indices of the SNPs to retain in the VCF file}
#' \item{indsubset}{Integer vector giving the indices of the samples to retain in the VCF file}
#' \item{file}{Character giving the name of the VCF file to be written}
#' \item{IDuse}{Character vector specifying alternative samples names. Useful for anonymizing sample IDs.}
#' }
#'
#' @name $writeVCF
#' @author Timothy P. Bilton
#' @seealso \code{\link{RA}}
#' @examples
#' file <- simDS()
#' RAfile <- VCFtoRA(file$vcf)
#' subset <- readRA(RAfile, snpsubset = 10:30)
#'
#' ## write the subset of the data back to VCF format
#' subset$writeVCF(file = "subset")
#'
NULL
