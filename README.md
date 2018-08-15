# GUSbase v0.1.0
Genotyping Uncertainty with Sequencing data - Base package

An R package for the basis of performing analysis on high-throughput sequencing data for a number of R packages, namely:
- [GUSMap](https://github.com/tpbilton/GUSMap)
- [GUSLD](https://github.com/AgResearch/GUS-LD)

The novelty of these packages is the methods account for errors associated with low sequencing depth and miscalled bases, so that filtering with respect to read depth is not required. 

[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)

### Installation:

The easiest way to install GUSbase in R is using the devtools package.

```
install.packages("devtools")
library(devtools)
install_github("tpbilton/GUSbase")
```
Note: Some of the functions are coded in C and therefore an appropriate C compiler is needed for the package to work. For windows OS, Rtools (https://cran.r-project.org/bin/windows/Rtools/) provides a compiler. 


### License:

GUSbase is Copyright 2017-2018 Timothy P. Bilton, released under the GNU General Public License version 3.

### Funding:
The initial development of this package was partially funed by the Ministry of Business, Innovation and Employment via its funding of the “Genomics for Production & Security in a Biological Economy” programme (Contract ID C10X1306).
