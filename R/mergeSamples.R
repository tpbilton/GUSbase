##########################################################################
# Genotyping Uncertainty with Sequencing data - Base package (GUSbase)
# Copyright 2017-2020 Timothy P. Bilton <timothy.bilton@agresearch.co.nz>
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
#' RA Method: Merge samples in an RA object
#'
#' Method for merging read count information from multiple samples in an RA object
#'
#' The main argument is \code{samID} which expects a character vector (length equal to the number of individuals in
#' the dataset) where the enteries correspond to the sample IDs. Enteries in \code{samID} that are identical
#' result in the corresponding samples being merged. See the examples for an illustration.
#'
#' @section Usage:
#' \preformatted{
#' RAobj$mergeSamples(samID, useID=TRUE)
#' }
#'
#' @section Arguments:
#' \describe{
#' \item{samID}{Character vector:  }
#' \item{useID}{Logical: If \code{TRUE}, the values in \code{samID} are used. Otherwise, the original sample
#' IDs in the RA object are used.}
#' }
#'
#' @name $mergeSamples
#' @author Timothy P. Bilton
#' @seealso \code{\link{RA}}
#' @examples
#' ## Load simulated dataset with GUSbase into R
#' file <- simDS()
#' RAfile <- VCFtoRA(file$vcf)
#' subset <- readRA(RAfile, snpsubset = 10:30)
#'
#' ## create new sample ID vector
#' samID = substr(ra$.__enclos_env__$private$indID, 1, 2) ## sample IDs from RA object
#' samID[5:104] = paste0(samID[5:104],"_", 1:100)
#' head(samID) # will merge first row with third and second row with fourth row
#'
#' ## merge samples
#' subset$mergeSamples(samID)
#'
NULL
