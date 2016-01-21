
######################
## Checking treeWAS ##
######################

## Check diff btw coalescent.sim.R and coalescent.sim.SP.R
## --> Currently using coalescent.sim.R...

## Check SOURCE lines in treeWAS.R (for fitch.sim.R and tree.sim.R)
## Or, better yet, get package to compile properly...

## Check adegenet sequences.R source -- still needed for coalescent.sim?
## If SO, either push to adegenet and/or include in treeWAS (for now)

#############
## EXAMPLE ##
#############

## load sample distribution
## (currently using the ClonalFrame Saureus output
## just so we can see what happens when we
## simulate data based on this distribution
## AND then use it to inform treeWAS
## (as compared to treeWAS's performance with dist=NULL))

data(dist)

## NOTE: if data(dist) not working, you may load the file like this,
## but you may need to change the file path here for now.
#dist <- get(load("./data/example.output_distribution.Rdata"))

out <- coalescent.sim(n.ind=100, gen.size=10000, sim.by="locus",
                      theta=1*2, dist=dist,
                      theta_p=15, phen=NULL,
                      n.snp.assoc=20, assoc.option="all", assoc.prob=90,
                      haploid=TRUE, biallelic=TRUE, seed=1,
                      plot=TRUE, heatmap=FALSE, plot2=FALSE)

## side note: if you set plot2="UPGMA", we can see how the reconstructed tree
## can be quite different (at least in terms of branch lengths) from the
## real tree used to generate the data
## --> implications for treeWAS performance if based on reconstructed tree?

#############################
## isolate elements of out ##
#############################
snps <- out[[1]]
tree <- out[[2]]

## snps names:
snps.names <- colnames(snps)
## snps indices:
snps.indices <- c(1:ncol(snps))

##########################################
## isolate set-specific elements of foo ##
##########################################
phen <- out[[3]][[2]][[5]]
snps.assoc <- out[[3]][[2]][[6]]
snps.assoc

#################
## Run treeWAS ##
#################

## First, we'll try treeWAS with dist=NULL
## (so it will use the default Poisson with parameter 1 to
## get the number of substitutions per site to simulate)

foo <- treeWAS(x=snps, y=phen, tree=tree,
               mt.rate=NULL,
               p.value=0.001, mt.correct=FALSE, p.value.option="count",
               plot.null.dist=TRUE, plot.dist=FALSE, plot.Fitch=TRUE,
               plot.tree=FALSE, tree.method="ml",
               n.reps=1, sim.gen.size=10000, sim.by="locus", dist=NULL,
               method="Didelot", test="score", corr.dat.fisher=NULL)

str(foo)

##############
## EVALUATE ##
##############
test.positive <- foo[[4]][,1]
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



#################    #################    #################    #################

## COMPARE TO: ##

#################
## Run treeWAS ##
#################

## Second, we can try treeWAS with dist=dist
## (where dist comes from the .Rdata file loaded just before we ran coalescent.sim)
## (so it will use the true distribution to
## identify the number of substitutions per site to simulate)

foo2 <- treeWAS(x=snps, y=phen, tree=tree,
               mt.rate=NULL,
               p.value=0.001, mt.correct=FALSE, p.value.option="count",
               plot.null.dist=TRUE, plot.dist=FALSE, plot.Fitch=TRUE,
               plot.tree=FALSE, tree.method="ml",
               n.reps=1, sim.gen.size=10000, sim.by="locus", dist=dist,
               method="Didelot", test="score", corr.dat.fisher=NULL)

str(foo2)

##############
## EVALUATE ##
##############
test.positive <- foo2[[4]][,1]
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
