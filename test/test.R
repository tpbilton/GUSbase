
library(GUSbase)
file <- GUSMap:::Manuka11()$vcf
out <- VCFtoRA(file, direct = "test")
RAdata <- readRA(out, gform="reference", excsamp = c("P1_S1","P1_S2","P2_S1","P2_S2"))
URdata <- makeUR(RAdata, ploid=2)

tt <- GUSLD(URdata)
tt <- GUSLD(URdata, SNPsets=list(12:22,1:10))
undebug(GUSLD)


file <- simDS()
out <- VCFtoRA(file$vcf, direct = "test")
RAdata <- readRA(out, excsamp = c("P1_1","P2_1","P1_2","P2_2"))
URdata <- makeUR(RAdata)
