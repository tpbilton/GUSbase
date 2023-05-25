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
#' Function for extracting the path to the Variant Call Format (VCF) file of a simulated dataset and associated pedigree file.
#'
#' The data set consists of a simulated F1 family of 100 progeny, 1000 SNPs and three chromosomes. The data is in VCF format (see this \href{https://samtools.github.io/hts-specs/VCFv4.2.pdf}{page}
#' for specification of VCF format). The pedigree file contains five columns, namely
#' \itemize{
#' \item SampleID: A unique character string of the sample ID. These correspond to those found in the VCF file
#' \item IndividualID: A character giving the ID number of the individual for which the sample corresponds to.
#' \item Mother: The ID of the mother as given in the IndividualID. Blank fields means the mother is unknown.
#' \item Father: The ID of the father as given in the IndividualID. Blank fields means the father is unknown.
#' \item Family: The name of the Family for a group of progeny with the same parents.
#' }
#'
#' @return A character string of the complete path to the simulated dataset contained in the package
#' and a pedigree file that goes with the data.
#' @author Timothy P. Bilton
#' @examples
#' ## extract the name of the vcf file and the pedigree file
#' simDS()
#' @export
## Wrapper function for reading in the F1 simulate data set
simDS <- function(){
  return(list(vcf=system.file("extdata", "simDS.vcf", package="GUSbase"),
              ped=system.file("extdata", "simDS_ped.csv", package="GUSbase"),
              meta=system.file("extdata", "simDS_metaInfo.csv", package="GUSbase")))
}
