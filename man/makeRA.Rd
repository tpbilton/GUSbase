% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeRA.R
\name{makeRA}
\alias{makeRA}
\title{Make an Reference/Alternate (RA) object.}
\usage{
makeRA(ref, alt, indID, chrom = NULL, pos = NULL, sampthres = 0.01)
}
\arguments{
\item{ref}{Matrix of non-negative counts of reference alleles. Number of rows correspond to the number of samples. 
Number of columns correspond to the number of SNPs.}

\item{alt}{Matrix of non-negative counts of alternate alleles. Number of rows correspond to the number of samples. 
Number of columns correspond to the number of SNPs.}

\item{indID}{Vector of IDs for samples. Needs to have the same length as the number of rows in `ref` and `alt` arguments. 
Note: Duplicate IDs will be merged into one sample.}

\item{chrom}{Vector of chromosome numbers. Needs to have the same length as the number of columns in `ref` and `alt` arguments. 
If `NULL`, then each SNP is given a random unique chromsome number.}

\item{pos}{Vector of SNP positions. Needs to have the same length as the number of columns in `ref` and `alt` arguments. 
If `NULL`, then each SNP is given a random unique positi number.}

\item{sampthres}{A numeric value giving the filtering threshold for which individual samples are removed. Default is 0.01
which means that samples with an average number of reads per SNP that is less than 0.01 are removed.}
}
\value{
An R6 object of class RA.
}
\description{
Function which creates an RA object from inputs of count matrices of reference and alternate alleles.
}
\examples{
## Load data from package
datloc = system.file("extdata", "simdata.Rdata", package = "GUSbase")
load(datloc)

## Create RA object from the data
makeRA(ref, alt, indID, chrom, pos)

}
\author{
Timothy P. Bilton
}
