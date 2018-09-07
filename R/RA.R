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
#' RA Object
#'
#' Class for storing reference/alternative (RA) data and methods for performing exporartory data analysis.
#'
#' @usage
#' ## Create RA object
#' RAobj <- readRA(genofile, gform, sampthres = 0.01, excsamp = NULL)
#'
#' ## RA Functions (Methods)
# #' RAobj$cometPlot()
#' RAobj$extractVar(nameList)
#'
#' @section Details:
#' An RA object is returned from the \code{\link{readRA}} function and contains the RA data, various
#' statistics of the dataset that have been computed, and functions (or methods) for analyzing the data.
#'
#' @section Functions(Methods):
#' \describe{
# #'   \item{\code{\link{$cometPlot}}}{Function for create a comet plot}
#'     \item{\code{\link{$extractVar}}}{Extract variables stored in an RA object}
#' }
#' @format NULL
#' @author Timothy P. Bilton
#' @seealso \code{\link{readRA}}
#' @name RA
#' @export


### R6 class for data aligned to reference assembly
RA <- R6Class("RA",
              public = list(
                initialize = function(List){
                  private$genon      <- List$genon
                  private$ref        <- List$ref
                  private$alt        <- List$alt
                  private$chrom      <- List$chrom
                  private$pos        <- List$pos
                  private$SNP_Names  <- List$SNP_Names
                  private$indID      <- List$indID
                  private$nSnps      <- List$nSnps
                  private$nInd       <- List$nInd
                  private$gform      <- List$gform
                  private$AFrq       <- List$AFrq
                  private$infilename <- List$infilename
                },
                print = function(...){
                  cat("Data Summary:\n")
                  cat("Data file:\t",private$infilename,"\n")
                  temp <- private$ref + private$alt
                  cat("Mean Depth:\t", mean(temp),"\n")
                  cat("Mean Call Rate:\t",sum(temp!=0)/length(temp),"\n")
                  cat("Number of...\n")
                  cat("  Individuals:\t",private$nInd,"\n")
                  cat("  SNPs:\t\t",private$nSnps,"\n")
                  cat("  Reads:\t",sum(temp),"\n")
                },
                extractVar = function(nameList){
                  if(!is.character(nameList) || !is.vector(nameList) || any(is.na(nameList)))
                    res <- NULL
                  for(name in nameList){
                    res <- c(res,list(private[[name]]))
                  }
                  names(res) <- nameList
                  return(res)
                }
                #### Diagonostic functions ####
                ## Ratio of alleles for heterozygous genotype calls (observed vs expected)
                #cometPlot = function(model="random", alpha=NULL, filename=NULL, cex=1, maxdepth=500, ...){
                #  cometPlot(private$ref, private$alt, model=model, alpha=alpha, filename=filename, cex=cex, maxdepth=maxdepth, ...)
                #}
                ###############################
              ),
              private = list(
                genon = NULL,
                ref = NULL,
                alt = NULL,
                chrom = NULL,
                pos = NULL,
                SNP_Names = NULL,
                indID = NULL,
                nSnps = NULL,
                nInd = NULL,
                gform = NULL,
                AFrq = NULL,
                infilename = NULL,
                updatePrivate = function(List){
                  for(elem in names(List)){
                    private[[elem]] <- List[[elem]]
                  }
                }
              )
)
