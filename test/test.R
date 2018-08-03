
library(GUSbase)
file <- GUSMap:::Manuka11()$vcf
out <- VCFtoRA(file, direct = "test")
RAdata <- readRA(out, gform="reference")
URdata <- createPop(RAdata)

debug(GUSbase:::makePop.UR)

tt <- GUSLD(URdata, SNPsets=list(12:22))
tt <- GUSLD(URdata, SNPsets=list(12:22,1:10))
undebug(GUSLD)
