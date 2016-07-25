

#############
## EXAMPLE ##
#############

#######################
## Clear environment ##
#######################
## NOTE TO USER: his step will delete all variables from your environment.
## You may want to save unsaved variables or skip this step.
rm(list=ls())

##############################
## Load sample distribution ##
##############################
## (currently using the ClonalFrame Saureus output
## just so we can see what happens when we
## simulate data based on this distribution
## AND then use it to inform treeWAS
## (as compared to treeWAS's performance with dist=NULL))
data(dist)

################################
## Simulate a coalescent tree ##
################################
tree <- coalescent.tree.sim(n.ind = 100, seed = 1)

#######################################################
## Simulate a phenotype for individuals in this tree ##
#######################################################
## get list of phenotype simulation output
phen.output <- phen.sim(tree, n.subs = 15)

## get phenotype for terminal nodes only
phen <- phen.output$phen

## get phenotype for all nodes,
## terminal and internal
phen.nodes <- phen.output$phen.nodes

## get the indices of phen.subs (ie. branches)
phen.loci <- phen.output$phen.loci

#################################
## Plot Tree showing Phenotype ##
#################################
phen.plot.colours <- plot.phen(tree = tree,
                               phen.nodes = phen.nodes,
                               plot = TRUE)

###################################################################
## Simulate genetic data (SNPs) that fit this tree and phenotype ##
###################################################################
snps.output <- snp.sim(n.snps = 10000, n.subs=dist,
                       n.snps.assoc = 10, assoc.prob = 90,
                       tree = tree,
                       phen.loci = phen.loci,
                       heatmap = FALSE, reconstruct = FALSE,
                       dist.dna.model="JC69",
                       seed = 1)
snps <- snps.output$snps
snps.assoc <- snps.output$snps.assoc
snps.names <- colnames(snps)
snps.indices <- c(1:ncol(snps))

################################################################################
## Note that all previous steps can be performed with this combined function: ##
################################################################################
# sim.output <- coalescent.sim(n.ind=100,
#                       n.snps=10000, n.subs=1,
#                       n.snps.assoc=10, assoc.prob=90,
#                       n.phen.subs=15, phen=NULL,
#                       plot=TRUE,
#                       heatmap=FALSE, reconstruct=FALSE,
#                       seed=1)
# snps <- sim.output$snps
# tree <- sim.output$tree
# phen <- sim.output$phen
# snps.assoc <- sim.output$snps.assoc


#################
## Run treeWAS ##
#################

## First, we'll try treeWAS with dist=NULL
## (so it will use the default Poisson with parameter 1 to
## get the number of substitutions per site to simulate)

treeWAS.output <- treeWAS(snps, phen, n.subs = 1,
                          tree = tree,
                          dist.dna.model = NULL, plot.tree = FALSE,
                          test = "score",
                          p.value = 0.001, p.value.correct = "bonf", p.value.by = "count",
                          sim.n.snps = 10000, n.reps = 1,
                          plot.null.dist = TRUE, plot.dist = FALSE)

str(treeWAS.output)

# out <- treeWAS.output
# corr.dat <- out$corr.dat
# corr.sim <- out$corr.sim


##############
## EVALUATE ##
##############
test.positive <- treeWAS.output$sig.snps$SNP.locus
test.negative <- snps.indices[-which(snps.indices %in% test.positive)]
## get true positives
snps.not <- snps.names[-which(snps.indices %in% snps.assoc)]
true.positive <- test.positive[which(test.positive %in% snps.assoc)]
TP <- length(true.positive)
## get true negatives
true.negative <- test.negative[which(test.negative %in% snps.not)]
TN <- length(true.negative)
## get false positives
false.positive <- test.positive[which(test.positive %in% snps.not)]
FP <- length(false.positive)
## get false negatives
false.negative <- test.negative[which(test.negative %in% snps.assoc)]
FN <- length(false.negative)


#################
## sensitivity ##
#################
## ie. How many truly ASSOCIATED SNPs did you manage to catch
## ~ Pr(Positive Test | SNP ASSOCIATED)
## --> Set 1: will be 0/0 = NaN
sensitivity <- (TP / (TP + FN))
sensitivity
#################
## specificity ##
#################
## ie. Of all the truly NOT associated SNPs, how many did you manage to rule out?
## ~ Pr(Negative Test | SNP NOT associated)
specificity <- (TN / (TN + FP)) ## = (1 - FPR)
specificity
#########
## PPV ##
#########
## ie. Of all the POSITIVE calls you made, how many were CORRECT/ identified truly ASSOCIATED SNPs
## ~ Pr(SNP ASSOCIATED | Positive Test)
## --> Set 1: will be 0 (UNLESS you made NO positive calls, then 0/0 = NaN)
PPV <- (TP / (TP + FP)) ## = (1 - FDR)
PPV


#################    #################    #################    #################

## COMPARE TO: ##

#################
## Run treeWAS ##
#################

## Second, we can try treeWAS with dist=dist
## (where dist comes from the .Rdata file loaded just before we ran coalescent.sim)
## (so it will use the true distribution to
## identify the number of substitutions per site to simulate)

treeWAS.output2 <- treeWAS(snps, phen, n.subs = dist,
                           tree = tree,
                           dist.dna.model = NULL, plot.tree = FALSE,
                           test = "score",
                           p.value = 0.001, p.value.correct = "bonf", p.value.by = "count",
                           sim.n.snps = 10000, n.reps = 1,
                           plot.null.dist = TRUE, plot.dist = FALSE)

str(treeWAS.output2)

##############
## EVALUATE ##
##############
test.positive <- treeWAS.output2$sig.snps$SNP.locus
test.negative <- snps.indices[-which(snps.indices %in% test.positive)]
## get true positives
snps.not <- snps.names[-which(snps.indices %in% snps.assoc)]
true.positive <- test.positive[which(test.positive %in% snps.assoc)]
TP <- length(true.positive)
## get true negatives
true.negative <- test.negative[which(test.negative %in% snps.not)]
TN <- length(true.negative)
## get false positives
false.positive <- test.positive[which(test.positive %in% snps.not)]
FP <- length(false.positive)
## get false negatives
false.negative <- test.negative[which(test.negative %in% snps.assoc)]
FN <- length(false.negative)

#################
## sensitivity ##
#################
## ie. How many truly ASSOCIATED SNPs did you manage to catch
## ~ Pr(Positive Test | SNP ASSOCIATED)
## --> Set 1: will be 0/0 = NaN
sensitivity <- (TP / (TP + FN))
sensitivity
#################
## specificity ##
#################
## ie. Of all the truly NOT associated SNPs, how many did you manage to rule out?
## ~ Pr(Negative Test | SNP NOT associated)
specificity <- (TN / (TN + FP)) ## = (1 - FPR)
specificity
#########
## PPV ##
#########
## ie. Of all the POSITIVE calls you made, how many were CORRECT/ identified truly ASSOCIATED SNPs
## ~ Pr(SNP ASSOCIATED | Positive Test)
## --> Set 1: will be 0 (UNLESS you made NO positive calls, then 0/0 = NaN)
PPV <- (TP / (TP + FP)) ## = (1 - FDR)
PPV

#################    #################    #################    #################

## COMPARE TO: ##

#################
## Run treeWAS ##
#################

## Third, we can try treeWAS with dist=NULL
## So we will use the Fitch parsimony functions from R pkg phangorn
## (reconfigured for our purposes in treeWAS function get.fitch.n.mts)
## to reconstruct the distribution of n.subs-per-site from the snps data and tree.

treeWAS.output3 <- treeWAS(snps, phen, n.subs = NULL,
                           tree = tree,
                           dist.dna.model = NULL, plot.tree = FALSE,
                           test = "score",
                           p.value = 0.001, p.value.correct = "bonf", p.value.by = "count",
                           sim.n.snps = 10000, n.reps = 1,
                           plot.null.dist = TRUE, plot.dist = FALSE)

str(treeWAS.output3)

##############
## EVALUATE ##
##############
test.positive <- treeWAS.output3$sig.snps$SNP.locus
test.negative <- snps.indices[-which(snps.indices %in% test.positive)]
## get true positives
snps.not <- snps.names[-which(snps.indices %in% snps.assoc)]
true.positive <- test.positive[which(test.positive %in% snps.assoc)]
TP <- length(true.positive)
## get true negatives
true.negative <- test.negative[which(test.negative %in% snps.not)]
TN <- length(true.negative)
## get false positives
false.positive <- test.positive[which(test.positive %in% snps.not)]
FP <- length(false.positive)
## get false negatives
false.negative <- test.negative[which(test.negative %in% snps.assoc)]
FN <- length(false.negative)

#################
## sensitivity ##
#################
## ie. How many truly ASSOCIATED SNPs did you manage to catch
## ~ Pr(Positive Test | SNP ASSOCIATED)
## --> Set 1: will be 0/0 = NaN
sensitivity <- (TP / (TP + FN))
sensitivity
#################
## specificity ##
#################
## ie. Of all the truly NOT associated SNPs, how many did you manage to rule out?
## ~ Pr(Negative Test | SNP NOT associated)
specificity <- (TN / (TN + FP)) ## = (1 - FPR)
specificity
#########
## PPV ##
#########
## ie. Of all the POSITIVE calls you made, how many were CORRECT/ identified truly ASSOCIATED SNPs
## ~ Pr(SNP ASSOCIATED | Positive Test)
## --> Set 1: will be 0 (UNLESS you made NO positive calls, then 0/0 = NaN)
PPV <- (TP / (TP + FP)) ## = (1 - FDR)
PPV
