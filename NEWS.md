
# GUSbase 0.1.1

* Added `$writeVCF` function to RA object. Converts the allele count data in the RA object back to a VCF file.
* Added `$writeRA` function to RA object. Converts the allele count data in the RA object back to an RA file.
* `readRA` function now has argument `snpsubset` which allows only a subset of the SNPs in the RA file to be read in.
* `makeUR` function now has argument `indsubset` which allows only a subset of the individuals to be retained in the creation of the UR population.
* The `$p_est` function in the makeUR was fixed up and is now working properly. In addition, the estimation of allele frequencies and sequencing errors have been implemented using the EM algorithm which is much faster and considerably more reliable.

* A few bugs were also fixed.

# GUSbase 0.1.0

This is the first release of GUSbase. This package forms a base for reading in data into the GUSverse.

## Main Functions:

* Converting VCF files into reference/alternate (RA) count format. 
* Reading RA format into R.
* Creating Unrelated populations
* Transformation functions used in the GUSverse

## Comments:

* Additionally features are due to be added on to this package as the other packages of the GUSverse are developed and released to CRAN.


