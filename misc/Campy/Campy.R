
################
## Campy data ##
################

require(ape)
require(adegenet)

###########
## ST 21 ##
###########

ST21.dna <- fasta2DNAbin("~/treeWAS/misc/Campy/ST21.fasta")
str(ST21.dna)

ST21.gen <- DNAbin2genind(ST21.dna)
str(ST21.gen@tab) # ncol = 94,711

###
ST21.phen <- read.table("~/treeWAS/misc/Campy/ST21-hosts.csv", header=FALSE, sep=",")
head(ST21.phen)
ST21.noms.phen <- ST21.phen$V1
ST21.phen <- ST21.phen$V2

prefix <- "~/treeWAS/misc/Campy/ST21.fasta.out"

ST21 <- read.CFML(prefix=prefix, plot=TRUE)

ST21.tree <- ST21$tree
tree <- ST21.tree
ST21.seqs <- ST21$seqs
seqs <- ST21.seqs
ST21.index <- ST21$index
index <- ST21.index
ST21.dist <- ST21$dist
dist <- ST21.dist

l <- length(seqs[1,])
l # l = 12,522

table(index)[1:10]

plot(tree)
title("Newick Tree (CFML)
      \n ST 21")

str(seqs)
tips <- which(rownames(seqs) %in% tree$tip.label)
snps <- DNAbin2genind(seqs[tips,], polyThres = 0) # SNPs already selected by CFML
snps <- snps@tab
str(snps) # ncol = 25,712
snps[1:10,1:10]

ST21.snps <- snps

## Save data for ST21
# save(ST21.snps, file="~/treeWAS/misc/Campy/ST21.snps.Rdata")
# save(ST21.phen, file="~/treeWAS/misc/Campy/ST21.phen.Rdata")
# save(ST21.tree, file="~/treeWAS/misc/Campy/ST21.tree.Rdata")
# save(ST21.seqs, file="~/treeWAS/misc/Campy/ST21.seqs.Rdata")
# save(ST21.index, file="~/treeWAS/misc/Campy/ST21.index.Rdata")
# save(ST21.dist, file="~/treeWAS/misc/Campy/ST21.dist.Rdata")
# save(ST21, file="~/treeWAS/misc/Campy/ST21.read.CFML.Rdata")

## load data for ST21
snps <- get(load("~/treeWAS/misc/Campy/ST21.snps.Rdata"))
phen <- get(load("~/treeWAS/misc/Campy/ST21.phen.Rdata"))
tree <- get(load("~/treeWAS/misc/Campy/ST21.tree.Rdata"))
index <- get(load("~/treeWAS/misc/Campy/ST21.index.Rdata"))
dist <- get(load("~/treeWAS/misc/Campy/ST21.dist.Rdata"))


## Check for non-binary SNPs loci:
coln <- colnames(snps)
coln <- removeLastN(coln, 2) # multiple columns for multiple bases...
coln <- as.numeric(coln)
length(unique(coln)) # 12,552
length(which(table(coln) > 2)) # 657
length(which(table(coln) > 3)) # 11
length(which(table(coln) > 4)) # 0

## WHAT TO DO ABOUT NON-BINARY SNPS??

## removing for now...
toRemove.4 <- which(table(coln) > 3)
# snps.ori <- snps
snps <- snps[,-toRemove.4]
# coln.ori <- coln
coln <- coln[-toRemove.4]
toRemove.3 <- which(table(coln) > 3)
# snps.ori <- snps
snps <- snps[,-toRemove.4]


## REMEMBER to consider the removed columns later
## when mapping the original sequences back using the index!


###########
## ST 45 ##
###########

ST45.dna <- fasta2DNAbin("~/treeWAS/misc/Campy/ST45.fasta")
str(ST45.dna)

ST45.gen <- DNAbin2genind(ST45.dna)
str(ST45.gen@tab) # ncol = 94,711

###
ST45.phen <- read.table("~/treeWAS/misc/Campy/ST45-hosts.csv", header=FALSE, sep=",")
head(ST45.phen)
ST45.noms.phen <- ST45.phen$V1
ST45.phen <- ST45.phen$V2

prefix <- "~/treeWAS/misc/Campy/ST45.fasta.out"

ST45 <- read.CFML(prefix=prefix, plot=TRUE)

ST45.tree <- ST45$tree
tree <- ST45.tree
ST45.seqs <- ST45$seqs
seqs <- ST45.seqs
ST45.index <- ST45$index
index <- ST45.index
ST45.dist <- ST45$dist
dist <- ST45.dist

l <- length(seqs[1,])
l # l = 9,138

table(index)[1:10]

plot(tree)
title("Newick Tree (CFML)
      \n ST 45")

str(seqs)
tips <- which(rownames(seqs) %in% tree$tip.label)
snps <- DNAbin2genind(seqs[tips,], polyThres = 0) # SNPs already selected by CFML
snps <- snps@tab
str(snps) # ncol = 19,029
snps[1:10,1:10]

ST45.snps <- snps

## Save data for ST45
# save(ST45.snps, file="~/treeWAS/misc/Campy/ST45.snps.Rdata")
# save(ST45.phen, file="~/treeWAS/misc/Campy/ST45.phen.Rdata")
# save(ST45.tree, file="~/treeWAS/misc/Campy/ST45.tree.Rdata")
# save(ST45.seqs, file="~/treeWAS/misc/Campy/ST45.seqs.Rdata")
# save(ST45.index, file="~/treeWAS/misc/Campy/ST45.index.Rdata")
# save(ST45.dist, file="~/treeWAS/misc/Campy/ST45.dist.Rdata")
# save(ST45, file="~/treeWAS/misc/Campy/ST45.read.CFML.Rdata")

## Load data for ST45
snps <- get(load("~/treeWAS/misc/Campy/ST45.snps.Rdata"))
phen <- get(load("~/treeWAS/misc/Campy/ST45.phen.Rdata"))
tree <- get(load("~/treeWAS/misc/Campy/ST45.tree.Rdata"))
index <- get(load("~/treeWAS/misc/Campy/ST45.index.Rdata"))
dist <- get(load("~/treeWAS/misc/Campy/ST45.dist.Rdata"))


## Check for non-binary SNPs loci:
coln <- colnames(snps)
coln <- removeLastN(coln, 2) # multiple columns for multiple bases...
coln <- as.numeric(coln)
length(unique(coln)) # 9,138
length(which(table(coln) > 2)) # 724
length(which(table(coln) > 3)) # 29
length(which(table(coln) > 4)) # 0

## WHAT TO DO ABOUT NON-BINARY SNPS??




