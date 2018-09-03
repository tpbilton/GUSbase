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
#' Simulated dataset and pedigree file
#'
#' Function for extracting the path to the (VCF) file of a simulated dataset and pedigree file.
#'
#' The data consists of a simulated F1 family of 100 progeny, 1000 SNPs and three chromosomes. The data is in VCF format (see this \href{https://samtools.github.io/hts-specs/VCFv4.2.pdf}{page}
#' for specification of VCF format). The pedigree file contains five columns, namely
#' \itemize{
#' \item SampleID: A unique character string of the sample ID. These correspond to those found in the VCF file
#' \item IndividualID: A character giving the ID number of the individual for which the sample corresponds to.
#' Note that some samples can be from the same individual.
#' \item Mother: The ID of the mother as given in the IndividualID. Note, if the mother is unknown then this should be left blank.
#' \item Father: The ID of the father as given in the IndividualID. Note, if the father is unknown then this should be left blank.
#' \item Family: The name of the Family for a group of progeny with the same parents. Note that this is not necessary but if
#' given must be the same for all the progeny.
#' }
#'
#' @return A character string of the complete path to the simulated dataset contained in the package
#' and a pedigree file that goes with the data.
#' @author Timothy P. Bilton
#' @examples
#' ## extract the name of the vcf file and the pedigree file
#' simDS()
#'
#' @export
## Wrapper function for reading in the F1 simulate data set
simDS <- function(){
  return(list(vcf=system.file("extdata", "simDS.vcf", package="GUSbase"),
              ped=system.file("extdata", "simDS_ped.csv", package="GUSbase")))
}
