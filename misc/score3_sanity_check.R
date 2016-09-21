
#############################
## SANITY CHECK -- SCORE 3 ##
#############################

phen.rec <- res$dat$phen.rec


plot.phen(tree, phen.nodes=phen.rec, snp.nodes = snp.rec)
title("PERFECT SNP  (score3 = 7.13)", line=0)

snp <- phen.rec[1:100]

snp.rec <- asr(var = snp, tree = tree, type = "parsimony")
snp.rec <- snp.rec$var.rec

foo <- subsequent.test(snps.rec, phen.rec, tree)
foo # 7.13

snps.rec <- cbind(snp.rec, snp.rec)
head(snps.rec)


snps.sim <- res$dat$snps.sim

names(res$res$subsequent)

score <- res$vals$subsequent
str(score)

corr.sim <- score$corr.sim
max(corr.sim) # 6.54
which.max(corr.sim)

snp.sim <- snps.sim[,which.max(corr.sim)] # 56864


snp.sim.rec <- asr(var = snp.sim, tree = tree, type = "parsimony")
snp.sim.rec <- snp.sim.rec$var.rec

snps.sim.rec.ori <- snps.sim.rec

temp <- snps.sim.rec
temp <- temp+1
temp <- replace(temp, which(temp == 2), 0)
# temp <- replace(temp, which(temp == 1), 1)

# snps.sim.rec <- cbind(snp.sim.rec, snp.sim.rec)
snps.sim.rec <- temp

plot.phen(tree, phen.nodes=phen.rec, snp.nodes = snps.sim.rec[,1])
title("NULL SNP  (score3 = 6.54)", line=0)


