
# GUSbase 0.2.3

* Add function makeRA to add user to create a RA object from count matrices of reference and alternate alleles.
* Added another data file for an example of the makeRA object.

# GUSbase 0.2.1

* Comet plot can now be produced for each individual. A few bugs in the comet plot function have also been fixed
* Add a function to the RA object `$mergeSamples` which can merge samples in the RA object.
* Summary information in a RA object is computed once and stored to same time in recomputing

Various other bugs have been fixed


# GUSbase 0.2.0

Version released for Thesis by Bilton (2020). 

New plot functions added:
* Comet plot which is produced using `$cometPlot` in an RA object
* rocket plot which is produced using `$rocketPlot` in an RA object
* RDD (read ratio density) plot which is produced using `$RDDPlot` in an RA object

Various other bugs have also been fixed.

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


# References

* Bilton, T.P. (2020). "Developing statistical methods for genetic analysis of genotypes from genotyping-by-sequencing data" (Doctoral Thesis, University of Otago, Dunedin, New Zealand).
