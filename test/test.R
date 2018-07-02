
library(GUSbase)
file <- GUSMap:::Manuka11()$vcf
out <- VCFtoRA(file, direct = "test")
RAdata <- readRA(out, gform="reference")
URdata <- createPop(RAdata)

debug(GUSbase:::makePop.UR)
