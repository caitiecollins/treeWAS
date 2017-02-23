
######################
## Troubleshooting! ##
######################

## sample sim data for troubleshooting...

snps <- get(load("E:/treeWAS_Sims/set1/n.subs_1/set1_1_snps.Rdata"))
phen <- get(load("E:/treeWAS_Sims/set1/n.subs_1/set1_1_phen.Rdata"))
snps.assoc <- get(load("E:/treeWAS_Sims/set1/n.subs_1/set1_1_performance.Rdata"))
snps.assoc <- snps.assoc$snps.assoc
tree <- get(load("E:/treeWAS_Sims/set1/n.subs_1/set1_1_tree.Rdata"))


working.dir <- "E:/treeWAS_Sims/set"
set.number <- 1


#################################################
## Correcting for Pop. Strat. w/ PCA & DAPC ?? ##
#################################################

#########
## PCA ##
#########

## STEPS: ##
## (1) Run PCA
## (2) Select sig. number of PCs (--> HOW??)
## (3) Regress snps along sig residuals
## (4) Run assoc test on adjusted snps dataset.

## (Below: taken from Glasgow/practical/practical-GWAS_before_cuts.Rnw ~ practical-GWAS_day4.pdf)

## get snps:
# snps.ori <- snps
snps <- snps.ori
if(is.null(colnames(snps))) colnames(snps) <- c(1:ncol(snps))

## run PCA:
# pca1 <- dudi.pca(snps, scale=FALSE)

## Keep only n.PCs significant axes:
n.PCs <- 5
pca1 <- dudi.pca(snps, scale=FALSE, scannf=FALSE, nf=n.PCs)
# pca1.ori <- pca1

## OPTIONAL PCA STUFF: ####################################################################

## Plot eigenvalues (--> identify (if selected) n.signifiant PCs (seems reasonable)!!)
barplot(pca1$eig, main="PCA eigenvalues", col=c(rep("black", n.PCs), rep("grey", 1000)))

## get distance matrix:
D <- dist(pca1$li[,1:n.PCs])^2

## get hclust tree:
clust <- hclust(D, method="complete")

## plot tree:
plot(clust, main=paste("Clustering (complete linkage) based on the first ", n.PCs, " PCs", sep=""), cex=.4)
## NOTE-- hclust actually does a pretty good job of recreating our coalescent tree!
## (though it shortens the bottom third branches, extending the middle third down instead...)

## identify main pop clusters (for plotting, mainly):
# pop <- factor(cutree(clust, k=5))
# head(pop, 20)
# table(pop)

## OR -- determine n.clusters w BIC?
# grp <- find.clusters(snps, max.n.clust=50)
grp <- find.clusters(snps, n.pca=n.PCs, choose.n.clust=FALSE) # pca.select="percVar", perc.pca=60,
pop <- grp$grp # gives same result as cutree(clust, k=6)
table(pop)
# barplot(grp$Kstat)

# pop2 <- factor(cutree(clust, h=max(clust$height)/3))
## TO CHECK -- IS max.height/3 ALWAYS/USUALLY A GOOD CHOICE TO GET THE MAIN/GOLDILOCKS N.CLUSTERS?
# identical(pop, pop2) # TRUE
## TO CHECK-- WILL IT BE A BETTER IDEA TO AUTOMATE BY HEIGHT (H) OR BY N.CLUSTERS (K)? (Do we even need to do this automatically?)


## Plot PCA showing groups from hclust tree:
## Plot PCs 1 and 2:
## NOTE -- Ellipses may not show if clusters are small/tight and narrow/elongated (can increase cellipse from 1.5 to <= 10 if desired... )
s.class(pca1$li, fac=pop, col=transp(funky(5)), cpoint=2, cellipse=10,
        sub="PCA - axes 1 and 2")
add.scatter.eig(pca1$eig, n.PCs, 1, 2, ratio=.24, posi="bottomright")

## We can do the same for PCs 3 and 4:
s.class(pca1$li, xax=3, yax=4, fac=pop, col=transp(funky(5)), cellipse=8,
        cpoint=2, sub="PCA - axes 3 and 4")
add.scatter.eig(pca1$eig, n.PCs, 3, 4, ratio=.24, posi="bottomright")


## Cross-tabulate phen and pop:
table(phen, pop)

## Check for correlation btw phen and pop clusters --> Would/do we even want to correct for pop structure here?
chisq.test(table(phen, pop), simulate=TRUE)

###################################################################################################### ## end (optional) PCA stuff

############################################
## Correct for pop strat w PCA residuals: ##
############################################

## NOTE-- if any MISSING DATA contained in SNPs dataset, must REPLACE it here (eg. with the mean for that snps column).

## get formula (get string up to n.PCs): lm(e ~ pca1$li[,1] + pca1$li[,2] + pca1$li[,3] + pca1$li[,4] + pca1$li[,5])
PC.string <- sapply(c(1:n.PCs), function(e) paste("pca1$li[, ", e, "]", sep=""))
PC.string <- paste0(PC.string, collapse=" + ")
# print(PC.string)
var.string <- paste("e ~ ", PC.string)

## Correct SNPs w PCA!
snps.corrected <- apply(snps, 2, function(e) residuals(do.call(lm, list(as.formula(var.string))))) # may take a minute
# snps.corrected <- apply(snps, 2, function(e) residuals(lm(e ~ pca1$li[,1] + pca1$li[,2] + pca1$li[,3] + pca1$li[,4] + pca1$li[,5]))) # may take a minute
# snps.corrected.ori <- snps.corrected

# snps.corrected.pca.n.PCs.3 <- snps.corrected

## CHECK snps.corrected??
# snps.corrected[1:10,1:10]
#
# pca2 <- dudi.pca(snps.corrected, scale=FALSE, scannf=FALSE, nf=5)
#
# s.class(pca2$li, fac=pop, col=transp(funky(5)), cpoint=2,
#         sub="PCA - axes 1 and 2")
# add.scatter.eig(pca2$eig,5,1,2, ratio=.24, posi="topleft")
#
# ## PCA axes 3 and 4
# s.class(pca2$li, xax=3, yax=4, fac=pop, col=transp(funky(5)),
#         cpoint=2, sub="PCA - axes 3 and 4")
# add.scatter.eig(pca2$eig,5,3,4, ratio=.24, posi="topleft")

## First make sure PHEN is in BINARY form (0, 1) only!
levs <- levels(as.factor(phen))
phen.ini <- phen
if(any(!levs %in% c(0,1))){
  phen <- as.character(phen)
  phen <- replace(phen, which(phen == levs[1]), 0)
  phen <- replace(phen, which(phen == levs[2]), 1)
  phen <- as.numeric(phen)
  names(phen) <- names(phen.ini)
} # end make phen binary..


## Run association test with ANOVA (~ Fisher, but cannot run Fisher on non-binary data)

############
## TO DO !!! ----- If running this for real in simTest, MAKE SURE TO ADD STEP REDUCING SNPS TO ONLY ITS UNIQUE COLUMNS!!!!!!!!!    ####    ####    ####
############

## SLOW STEP..!
pval2 <- numeric(0)
system.time( # 120.78
for(i in 1:ncol(snps.corrected)){
  foo <- suppressWarnings(glm(phen ~ snps.corrected[,i], family="binomial"))
  ANOVA <- anova(foo, test="Chisq")
  pval2[i] <- ANOVA$"Pr(>Chi)"[2]
} # end for loop
)

## Store pvals and snps.corrected from pca-corrected ANOVA association test:
pval.pca <- pval2
snps.corrected.pca <- snps.corrected

# snps.corrected.pca.n.PCs.5 <- snps.corrected.pca
identical(snps.corrected.pca.n.PCs.3, snps.corrected.pca.n.PCs.5)

# pval2.ori <- pval2
# results:
length(which(pval2 < 0.05))
length(which(pval2 < 0.01))
# w bonf
length(which(pval2 < 0.05/ncol(snps)))
length(which(pval2 < 0.01/ncol(snps)))

# ##  CHECK/compare to ANOVA values for SNPS without correction?
# pval3 <- numeric(0)
# system.time(
#   for(i in 1:ncol(snps)){
#     foo <- suppressWarnings(glm(phen ~ snps[,i], family="binomial"))
#     ANOVA <- anova(foo, test="Chisq")
#     pval3[i] <- ANOVA$"Pr(>Chi)"[2]
#   } # end for loop
# )



##########################################################################################################################################################

##########################################################################################################################################################


##########
## DAPC ##
##########

## USE SAME CLUSTERS AS ABOVE (w BIC:)
grp <- find.clusters(snps, n.pca=n.PCs, choose.n.clust=FALSE, max.n.clust=6) # pca.select="percVar", perc.pca=60,
pop <- grp$grp # gives same result as cutree(clust, k=6)
n.grp <- length(levels(pop))
# table(pop)


## Using our pop clusters as the group factor in DAPC,
## we can generate a new DAPC object ...
## after performing cross-validation to optimise the discrimination between these subpopulations:

# # profiling of both time and memory
# Rprof("E:/treeWAS_Sims/Rprof_simTest", memory.profiling=T)
# xval.pop <- xvalDapc(snps, pop)
# Rprof(NULL)
# summaryRprof("E:/treeWAS_Sims/Rprof_simTest") # , memory="both"

xval.pop <- xvalDapc(snps, pop)
# str(xval.pop)
# xval.pop[2:6]

## store DAPC object:
dapc.pop <- xval.pop$DAPC

###################################
## Correct for pop strat w DAPC: ##
###################################

## As we did when correcting with PCA, we regress along the axes of DAPC.
## When correcting with the DAPC approach, we do not need to determine
## how many axes are ``significant'': we will always correct with (k - 1) axes.
## NOTE that (k - 1) is not necessarily = n.PCs
## ... though we may want to limit the K we work with depending on how system.time scales with K... (?)

## get formula (get string up to n.PCs): lm(e ~ dapc.pop$ind.coord[,1] + dapc.pop$ind.coord[,2] + dapc.pop$ind.coord[,3] + dapc.pop$ind.coord[,4])
# DAPC.string <- sapply(c(1:n.PCs), function(e) paste("dapc.pop$ind.coord[, ", e, "]", sep=""))
if(n.PCs == 5){
  DAPC.string <- sapply(c(1:4), function(e) paste("dapc.pop$ind.coord[, ", e, "]", sep=""))
}else{
  DAPC.string <- sapply(c(1:5), function(e) paste("dapc.pop$ind.coord[, ", e, "]", sep=""))
}
print(DAPC.string)
# DAPC.string <- sapply(c(1:(n.grp - 1)), function(e) paste("dapc.pop$ind.coord[, ", e, "]", sep=""))

DAPC.string <- paste0(DAPC.string, collapse=" + ")
var.string <- paste("e ~ ", DAPC.string)

## Correct SNPs w DAPC!
snps.corrected <- apply(snps, 2, function(e) residuals(do.call(lm, list(as.formula(var.string))))) # may take a minute

## Run association test on DAPC-corrected snps:
pval3 <- numeric(0)
system.time(
for(i in 1:ncol(snps.corrected)){
  foo <- suppressWarnings(glm(phen ~ snps.corrected[,i], family="binomial"))
  ANOVA <- anova(foo, test="Chisq")
  pval3[i] <- ANOVA$"Pr(>Chi)"[2]
} # end for loop
)

## store pvals and snps.corrected for DAPC-corrected assoc test:
pval.dapc <- pval3
snps.corrected.dapc <- snps.corrected

# snps.corrected.dapc.e4 <- snps.corrected
# pval3.new.e4 <- pval3

# snps.corrected.dapc.e5 <- snps.corrected
# pval3.new.e5 <- pval3
snps.corrected.dapc.e5.new <- snps.corrected
pval3.new.e5.new <- pval3

# identical(snps.corrected.pca.n.PCs.5, snps.corrected.dapc.e5)

# pval3.ori <- pval3
# pval3.new <- pval3

# results:
length(which(pval3 < 0.05))
length(which(pval3 < 0.01))
# w bonf
length(which(pval3 < 0.05/ncol(snps)))
length(which(pval3 < 0.01/ncol(snps)))

#




#





#






#








#
