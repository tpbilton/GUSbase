








data <- readRA("simDS.vcf.ra.tab")
ur <- makeUR(data, mafEst = T)
system.time({tt <- ur$.__enclos_env__$private$p_est(method="EM")})
system.time({tt2 <- ur$.__enclos_env__$private$p_est()})
