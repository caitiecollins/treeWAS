
#################
## ACE example ##
#################

setwd("/media/caitiecollins/Seagate Backup Plus Drive/SimBac")

prefix <- ("CFML_R_0")

tree<-read.tree(sprintf('%s.labelled_tree.newick',prefix))
seqs<-read.dna(sprintf('%s.ML_sequence.fasta',prefix),format='fasta')
mapping<-scan(sprintf('%s.position_cross_reference.txt',prefix),sep=',',quiet=T)
l<-length(seqs[1,])

## check & convert tree!
## (violations cause problems, eg. in phen.sim)
if(!is.binary.tree(tree)) tree <- multi2di(tree)
if(!identical(tree, reorder.phylo(tree,"postorder"))) tree <- reorder.phylo(tree,"postorder")
## seems like our tree$edge inds 1:n.ind are in plot order (ie. 1:100 botttom:top instead of corresponding to actual labs)
## --> added postorder check (as in fitch.testing) to readCFML code...

## examine seqs
str(seqs) # seqs for all inds (terminal and internal) & all unique site patterns
n.ind <- tree$Nnode + 1
n.nodes <- tree$Nnode + n.ind

############################################
## simulate example phenotype for testing ##
############################################
n.phen.subs <- 20 # 15
set.seed(1)
phen.list <- phen.sim(tree, n.subs = n.phen.subs)
## get phenotype for terminal nodes only
phen <- phen.list$phen
## get phenotype for all nodes,
## terminal and internal
phen.nodes <- phen.list$phen.nodes
## get the indices of phen.subs (ie. branches)
phen.loci <- phen.list$phen.loci

###############
## PLOT TREE ##
###############
## phylogram
phen.plot.col <- plot.phen(tree = tree,
                           phen.nodes = phen.nodes,
                           plot = TRUE)
## cladogram
phen.plot.col <- plot.phen(tree = tree,
                           phen.nodes = phen.nodes,
                           plot = TRUE,
                           type="c", use.edge.length=FALSE)



##################
## update names ##
##################

phen.ini <- phen

tree$tip.label <- as.character((as.numeric(tree$tip.label)+1))
names(phen) <- tree$tip.label

################
## CONTINUOUS ##
################

## Trying ACE after converting phen to continuous [0,100]
## (phen 0 --> 0, phen 1 --> 100)

## HYP: We could treat the inferred states of internal nodes as
## the probabilities of the true binary states??

phen.ori <- phen

phen <- as.character(phen)
phen <- replace(phen, which(phen=="A"), 0)
phen <- replace(phen, which(phen=="B"), 100)
phen <- as.numeric(phen)
phen <- as.factor(phen)
names(phen) <- names(phen.ori)
phen.f <- phen

phen <- as.numeric(as.character(phen))
names(phen) <- names(phen.f)

## Try ACE as continuous var...

#########
## ACE ##
#########
phen.ace <- ace(phen, tree, type="continuous")
str(phen.ace)
## get probs
## (ie. get the inferred continuous var value,
## which we will interpret as the probability
## of being phen state B).
probs <- phen.ace$ace
probs


###########################
## INDEPENDENT CONTRASTS ##
###########################
phen.ace2 <- ace(phen, tree, type="continuous", method="pic", scaled=TRUE)
str(phen.ace2)
probs2 <- phen.ace2$ace
identical(probs, probs2) # FALSE

phen.ace3 <- ace(phen, tree, type="continuous", method="pic", scaled=FALSE)
str(phen.ace3)
probs3 <- phen.ace3$ace
identical(probs2, probs3) # TRUE
## ?? Why doesn't scaling have any impact on the IC ace??


###############
## PLOT TREE ##
###############

## convert probs to show on tree pies:
probs.ori <- probs
df <- data.frame(100-probs, probs)
names(df) <- c("prob.A", "prob.B")
head(df)
## cladogram
phen.plot.col <- plot.phen(tree = tree,
                           phen.nodes = phen.nodes,
                           plot = TRUE,
                           type="c", use.edge.length=FALSE,
                           label.offset=1)
## try plotting ~ example: ##
#### Showing the probs and inverse probs on each node:
myCol <- c("blue", "red")
## convert phen to 1s and 2s (to indicate myCol 1 or 2)
tipCol.pattern <- as.numeric(as.factor(as.character(phen)))
## plot terminal squares
tiplabels(pch = 22, bg = myCol[tipCol.pattern], cex = 0.7, adj = 1)
## plot internal pies
nodelabels(thermo = df/100, piecol = transp(myCol, 0.8), cex = 0.25)




#####    #####    #####    #####    #####    #####    #####    #####    #####    #####    #####    #####

#####################################################
## MANUALLY RE-WRITE PHEN FOR PRESENTATION EXAMPLE ##
#####################################################

phen.old <- phen

phen

## want:
## 10 clusters of ~ 10 inds each
## 5 pure phen == R
## 5 pure phen == S, except 1 ind == R
### clusters should alternate to avoid meta-common ancestors...
phen[which(names(phen) %in%
             as.character(
               c(52, 100, 91, 38, 97, 75, 36, 76, 87)))] <- 100

phen[which(names(phen) %in%
             as.character(
               c(28, 90, 45, 30, 93, 3, 14, 49, 51, 95)))] <- 0

phen[which(names(phen) %in%
             as.character(
               c(18, 67, 19, 65, 66)))] <- 100

phen[which(names(phen) %in%
             as.character(57))] <- 100

phen[which(names(phen) %in%
             as.character(
               c(34, 39, 84, 2, 64, 70, 98, 26, 69, 40, 99, 59, 24, 46)))] <- 0

phen[which(names(phen) %in%
             as.character(
               c(44, 81, 63, 25, 56, 85, 80, 62, 20, 92, 60, 15)))] <- 100

phen[which(names(phen) %in%
             as.character(
               c(27, 6, 21, 17)))] <- 100
##
# ind <- which(names(phen) == "4")
# edgeCol[which(tree$edge[,2]==ind)] <- "green"
# phen[ind] <- "100"
#
# ind <- which(names(phen) == "14")
# edgeCol[which(tree$edge[,2]==ind)] <- "green"
# phen[ind] <- "100"
#
# ind <- which(names(phen) == "2")
# edgeCol[which(tree$edge[,2]==ind)] <- "green"
# phen[ind] <- "100"
#
# ind <- which(names(phen) == "40")
# edgeCol[which(tree$edge[,2]==ind)] <- "green"
# phen[ind] <- "100"
#
# ind <- which(names(phen) == "9")
# edgeCol[which(tree$edge[,2]==ind)] <- "green"
# phen[ind] <- "100"
#
# ind <- which(names(phen) == "52")
# edgeCol[which(tree$edge[,2]==ind)] <- "green"
# phen[ind] <- 0
#
# ind <- which(names(phen) == "73")
# edgeCol[which(tree$edge[,2]==ind)] <- "green"
# phen[ind] <- 100
#
# ind <- which(names(phen) == "28")
# edgeCol[which(tree$edge[,2]==ind)] <- "green"
# phen[ind] <- 100


## plot again ##
# phen.nodes.ori <- phen.nodes
phen.nodes <- as.character(phen)
phen.nodes <- replace(phen.nodes, which(phen.nodes=="100"), "B")
phen.nodes <- replace(phen.nodes, which(phen.nodes=="0"), "A")
phen.nodes <- as.factor(phen.nodes)

leafCol <- as.character(phen.nodes)
leafCol <- replace(leafCol, which(leafCol=="B"), "red")
leafCol <- replace(leafCol, which(leafCol=="A"), "blue")

edgeCol <- "black"
edgeCol <- rep("black", nrow(tree$edge))
edgeCol[c(1:22)] <- "red"
edgeCol[c(198)] <- "red"
edgeCol[c(23:48)] <- "blue"
edgeCol[c(49:57)] <- "red"
edgeCol[c(58:74)] <- "blue"
edgeCol[c(75:102)] <- "red"
edgeCol[c(103:122, 180:182)] <- "blue"
edgeCol[c(123:144)] <- "red"
edgeCol[c(145:159)] <- "blue"
edgeCol[c(160:177)] <- "red"
edgeCol[c(178:193)] <- "blue"
edgeCol[c(198, 177, 160, 194, 57)] <- "green"
edgeCol[c(195:197)] <- "blue"
edgeCol[c(111)] <- "green"# "blue"

# save(phen.nodes, file="./MRC_phen.nodes.Rdata")
# save(leafCol, file="./MRC_leafCol.Rdata")
# save(edgeCol, file="./MRC_edgeCol.Rdata")
# save(phen, file="./MRC_phen.Rdata")
# save(tree, file="./MRC_tree.Rdata")
# save(phen.cont, file="./MRC_phen.cont.Rdata")

### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###

###############
## load data ##
###############
setwd("/media/caitiecollins/Seagate Backup Plus Drive/MRC Seminar")

tree <- get(load("./MRC_tree.Rdata"))
phen <- get(load("./MRC_phen.Rdata"))
phen.nodes <- get(load("./MRC_phen.nodes.Rdata"))
leafCol <- get(load("./MRC_leafCol.Rdata"))
edgeCol <- get(load("./MRC_edgeCol.Rdata"))
snps <- get(load("./MRC_snps.Rdata"))
# phen.cont <- get(load("./MRC_phen.cont.Rdata"))
# phen.ori <- phen.cont


library(phangorn)
library(adegenet)

##########
## PLOT ##
##########
par(mar=c(1,1,4,2)+0.1)

lab <- as.character(tree$tip.label)
temp <- sapply(c(1:length(lab)),
               function(e)
                 paste(
                   paste(rep(c("-"), 3-nchar(lab[e])), collapse=""),
                   "-", lab[e], sep=""))
tiptext <- as.vector(unlist(temp))


plot(tree, show.tip=FALSE, edge.width=2,
     edge.color=edgeCol, type="c", use.edge.length=FALSE)
title("Coalescent tree w/ phenotypic changes")
tiplabels(text=tiptext, cex=0.6, adj=c(-.4, 0), col=leafCol, frame="none")

## Try ACE as continuous var...

#########
## ACE ##
#########

## CONTINUOUS

phen.ori <- phen
phen <- as.numeric(phen)
names(phen) <- names(phen.ori)
phen.ace <- ace(phen, tree, type="continuous")
str(phen.ace)
## get probs
## (ie. get the inferred continuous var value,
## which we will interpret as the probability
## of being phen state B).
probs <- phen.ace$ace
#probs


## DISCRETE

phen <- phen.ori
phen.ace.d <- ace(phen, tree, type="discrete")
str(phen.ace.d)
round(phen.ace.d$lik.anc, 3)
probs.d <- phen.ace.d$lik.anc[,2]*100

###########################
## INDEPENDENT CONTRASTS ##
###########################
phen.ace2 <- ace(phen, tree, type="continuous", method="pic", scaled=TRUE)
str(phen.ace2)
probs2 <- phen.ace2$ace
identical(probs, probs2) # FALSE

# phen.ace3 <- ace(phen, tree, type="continuous", method="pic", scaled=FALSE)
# str(phen.ace3)
# probs3 <- phen.ace3$ace
# identical(probs2, probs3) # TRUE
## ?? Why doesn't scaling have any impact on the IC ace??


###############
## PLOT TREE ##
###############

# ## cladogram
# phen.plot.col <- plot.phen(tree = tree,
#                            phen.nodes = phen.nodes,
#                            plot = TRUE,
#                            type="c", use.edge.length=FALSE,
#                            label.offset=1)

# par(mar=c(1,1,4,2)+0.1)
# lab <- as.character(tree$tip.label)
# temp <- sapply(c(1:length(lab)),
#                function(e)
#                  paste(
#                    paste(rep(c("-"), 3-nchar(lab[e])), collapse=""),
#                    "-", lab[e], sep=""))
# tiptext <- as.vector(unlist(temp))

## convert probs to show on tree pies:
probs.ori <- probs.d
df <- data.frame(100-probs.ori, probs.ori)
names(df) <- c("prob.A", "prob.B")
head(df)
# save(df, file="./MRC_probs_ACE.Rdata")

par(mar=c(1,1,1,2)-0.5)
plot(tree, show.tip=FALSE, edge.width=2,
     edge.color=edgeCol, type="c", use.edge.length=FALSE)
# title("Coalescent tree w/ phenotypic changes")
tiplabels(text=tiptext, cex=0.6, adj=c(-.4, 0), col=leafCol, frame="none")


## try plotting ~ example: ##
#### Showing the probs and inverse probs on each node:
myCol <- c("blue", "red")
## convert phen to 1s and 2s (to indicate myCol 1 or 2)
tipCol.pattern <- as.numeric(as.factor(as.character(phen)))
## plot terminal squares
tiplabels(pch = 22, bg = myCol[tipCol.pattern], cex = 0.7, adj = 1)
## plot internal pies
nodelabels(thermo = df/100, piecol = transp(myCol, 0.8), cex = 0.25)



##############
## get SNPs ##
##############

snps <- snp.sim(10000, n.subs=1, n.snps.assoc=0, tree=tree)
snps <- snps$snps

## get associated SNPS...
# phen <- phen.ori

## perfect corr w phen
snp <- phen
snp <- replace(snp, which(snp==100), 1)
snp1 <- snp

## corr only for lone phen R inds
snp <- rep(0, length(snp1))
snp <- replace(snp, which(names(phen) %in% as.character(c(9, 4, 14, 2, 40))), 1)
snp2 <- snp

## corr for all but lone phen R inds
snp <- snp1
snp <- replace(snp, which(names(phen) %in% as.character(c(9, 4, 14, 2, 40))), 0)
snp3 <- snp

## corr for bottom half only
snp <- snp1
snp <- replace(snp, c(which(names(phen) == "73"):100), 0)
snp4 <- snp

## corr for top half only
snp <- snp1
snp <- replace(snp, c(1:(which(names(phen) == "73")-1)), 0)
snp5 <- snp

snps[,1] <- snp1
snps[,2] <- snp2
snps[,3] <- snp3
snps[,4] <- snp4
snps[,5] <- snp5

#names(phen)
row.names(snps) <- names(phen)
## NOTE: use tree labels in plot order for snps names (or reorder snps rows)

# save(snps, file="./MRC_snps.Rdata")


# setwd("~/treeWAS/pkg")
# snps.ace <- snps.ori
# phen.ace <- phen
# tree.ace <- tree
# library(devtools)
# use_data(snps.ace)
# use_data(phen.ace)
# use_data(tree.ace)


### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###

### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###


###############
## LOAD DATA ##
###############

data("snps.ace")
data("phen.ace")
data("tree.ace")
snps <- snps.ace
phen <- phen.ace
tree <- tree.ace

################################################################################################################################################################

#########
## ACE ##
#########

## DISCRETE
var <- snps[,1]

system.time(
  snps.ace.d <- ace(var, tree, type="discrete")
  )

# x <- ace(var, tree, type="discrete")
# str(x)

snps.ace.d <- ace(var, tree, type="discrete")
# str(snps.ace.d)
# round(snps.ace.d$lik.anc, 3)
probs.d <- snps.ace.d$lik.anc[,2]*100
# probs.d <- ace.snp$lik.anc[,2]*100
##########
## PLOT ##
##########
## convert probs to show on tree pies:
probs.ori <- probs.d
df <- data.frame(100-probs.ori, probs.ori)
names(df) <- c("prob.A", "prob.B")
head(df)
# save(df, file="./MRC_probs_ACE.Rdata")

lab <- as.character(tree$tip.label)
temp <- sapply(c(1:length(lab)),
               function(e)
                 paste(
                   paste(rep(c("-"), 3-nchar(lab[e])), collapse=""),
                   "-", lab[e], sep=""))
tiptext <- as.vector(unlist(temp))

par(mar=c(1,1,1,2)-0.5)
edgeCol <- "black"
plot(tree, show.tip=FALSE, edge.width=2,
     edge.color=edgeCol, type="c", use.edge.length=FALSE)
## PHEN-coloured tip labels:
myCol <- c("blue", "red")
## convert phen to 1s and 2s (to indicate myCol 1 or 2)
tipCol.pattern <- as.numeric(as.factor(as.character(phen)))
tiplabels(text=tiptext, cex=0.6, adj=c(-.4, 0), col=myCol[tipCol.pattern], frame="none")


## try plotting ~ example: ##
#### Showing the probs and inverse probs on each node:
myCol <- c("blue", "red")
## convert phen to 1s and 2s (to indicate myCol 1 or 2)
tipCol.pattern <- as.numeric(as.factor(as.character(var[1:(tree$Nnode+1)])))
## plot terminal squares
tiplabels(pch = 22, bg = myCol[tipCol.pattern], cex = 0.7, adj = 1)
## plot internal pies
require(adegenet)
nodelabels(thermo = df/100, piecol = transp(myCol, 0.5), cex = 0.55)

#nodelabels(text=rev(unique(tree$edge[,1])), cex=0.5, bg=transp("blue", 0.3), adj=2)

################################################################################################################################################################

##########
## TEST ##
##########


#########
## ACE ##
#########

## DISCRETE ##

## get tree edge list
edges <- tree$edge

## get variable
var <- snps[,1]
var.terminal <- var
snps.ace.d <- ace(var, tree, type="discrete")
var.internal <- snps.ace.d$lik.anc[,2]
var <- c(var.terminal, var.internal)

#########################
## get diffs for SNPS: ##
#########################
diffs <- list()
for(i in 1:10){
  ## get variable
  var <- snps[,i]
  var.terminal <- var
  snps.ace.d <- ace(var, tree, type="discrete")
  var.internal <- snps.ace.d$lik.anc[,2]
  var <- c(var.terminal, var.internal)

  ## get differences for this variable's ace likelihoods
  diffs[[i]] <- get.ace.diffs(var, edges)
}

#########################
## get diffs for PHEN: ##
#########################
## check phen is numeric
# phen <- as.numeric(phen.cont)
# phen <- replace(phen, which(phen==100), 1)
# names(phen) <- names(phen.cont)
phen.terminal <- phen
phen.ace.d <- ace(phen, tree, type="discrete")
phen.internal <- phen.ace.d$lik.anc[,2]
var <- c(phen.terminal, phen.internal)

phen.diffs <-  get.ace.diffs(var, edges)

###########################################
## compare diffs for each SNPs vs. PHEN: ##
###########################################
## check w snps[,3] # all corr except lone phen.R inds:
snps.diffs <- diffs[[10]]

sp.diffs <- signs <- list()
for(i in 1:length(snps.diffs)){
  sp.diffs[[i]] <- snps.diffs[i] * phen.diffs[i]
  signs[[i]] <- ((snps.diffs[i] < 0 & phen.diffs[i] < 0) | (snps.diffs[i] > 0 & phen.diffs[i] > 0))
}

length(which(as.vector(unlist(signs)) == FALSE))
round(sum(as.vector(unlist(sp.diffs))), 4)

# hist(snps.diffs)
# hist(phen.diffs)

# par(mar=c(5,4,4,1)+0.1)


### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###

### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###

####################################  MOVED FNS and EXAMPLE TO R/ace.R  ####################################

### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###

### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###   ### ###



##############
## example: ##
##############

# load data
data("snps.ace")
data("phen.ace")
data("tree.ace")
snps <- snps.ace
phen <- phen.ace
tree <- tree.ace

## run ace.test on subset of data (otherwise SLOW!)
out <- ace.test(snps[,1:10], phen, tree, method="discrete")

## examine output
round(out, 4)



###########################################################################

###############
## PARSIMONY ##
###############

##############
## phangorn ##   ###   ###   ###   ###   ###   ###   ###   ###
##############
require(phangorn)

lev <- unique(as.vector(snps))
dna <- as.phyDat(snps, type="USER", levels=lev)
D <- dist.gene(snps)
## NOTE--dist.gene does not have the model arg that
## dist.dna has (allowing user to specify, eg. "raw", "JC69")
## Why not? Does it matter? If so, how to coerce binary snps into DNAbin??

####

snps.nt <- replace(snps, which(snps==0), "A")
snps.nt <- replace(snps.nt, which(snps==1), "G")

dna <- as.phyDat(snps.nt)
dna.bin <- as.DNAbin(dna)
D <- dist.dna(dna.bin, model="raw")

####

tree.ini <- nj(D)

plot(tree.ini, cex=0.5)
title("NJ tree")

plot(midpoint(tree.ini), cex=0.5)
title("NJ tree, midpoint")

## get parsimony (ie. overall parsimony cost)
cost <- parsimony(tree.ini, dna)
cost
#15004 # w snps 0/1
#5668 # w snps.nt

## get site-wise parsimony scores
cost <- parsimony(tree.ini, dna, site="site")
cost

#######################################
## ALL ABOVE STEPS CAN BE DONE WITH: ##
#######################################
cost <- get.fitch.n.mts(snps, tree)
## possibly replacing tree w nj...
cost <- get.fitch.n.mts(snps, tree.ini)

## OR w sankoff:
sank <- sankoff(tree, dna, site="site")
str(sank)

##############################################
## get optimised parsimony
tree.pars <- optim.parsimony(tree.ini, dna)
tree.pars
#Final p-score 15004 after  0 nni operations
## WHY DID IT NOT UNDERGO ANY NNIs???

plot(tree.pars, cex=0.5)
## AND YET WHY DOES THE OUTPUT TREE LOOK SO VERY DIFFERENT???

## phangorn doc claims the parsimony ratchet
## is the preferred way to search for the best tree
tree.prat <- pratchet(dna)
str(tree.prat)

########################################
## try parsimony on individual snps ? ##
########################################

## NB: branch lengths corrsond to n.subs
# foo.ori <- foo
# states.ori <- states
# mpr.ori <- mpr
# dna.ori <- dna

snp <- snps[,1:10]

lev <- unique(as.vector(snp))
dna <- as.phyDat(snp, type="USER", levels=lev)

D <- dist.gene(snps)
## NOTE--dist.gene does not have the model arg that
## dist.dna has (allowing user to specify, eg. "raw", "JC69")
## Why not? Does it matter? If so, how to coerce binary snps into DNAbin??
tree.ini <- nj(D)
plot(tree.ini, cex=0.5)
plot(midpoint(tree.ini), cex=0.5)

## get parsimony (ie. overall parsimony cost)
cost <- parsimony(tree.ini, dna)
cost
#39

## get optimised parsimony
tree.pars <- optim.parsimony(tree.ini, dna)
str(tree.pars)
#Final p-score 19 after  1 nni operations
## WHY DID IT UNDERGO 1 NNI NOW W ONE SNP BUT NOT UNDERGO ANY BEFORE????

plot(tree.pars, cex=0.5)
## AND YET WHY DOES THE OUTPUT TREE LOOK SO VERY DIFFERENT???

## phangorn doc claims the parsimony ratchet
## is the preferred way to search for the best tree
tree.prat <- pratchet(dna)
str(tree.prat)

#######################################
## ALL ABOVE STEPS CAN BE DONE WITH: ##
#######################################
cost <- get.fitch.n.mts(snps, tree)
## possibly replacing tree w nj...
cost1 <- get.fitch.n.mts(snps, tree.ini)

## OR w sankoff:
sank <- sankoff(tree, dna, site="site")
str(sank)

#############################################
## get SNP states of all internal nodes ?? ##
#############################################

#########
## MPR ##
#########

## pace == ancestral.pars
pa.MPR <- pace(tree, dna, type="MPR")
str(pa.MPR)

#############
## ACCTRAN ##
#############

## pace == ancestral.pars
pa.ACCTRAN <- pace(tree, dna, type="ACCTRAN")
str(pa.ACCTRAN)

## pace  --> diff resuls w MPR vs. ACCTRAN
diffs <- sapply(c(1:length(pa.ACCTRAN)), function(e) identical(pa.MPR[[e]], pa.ACCTRAN[[e]]))
which(diffs==FALSE) # 101:199

###########################################
## convert reconstruction back to snps.. ##
###########################################
## each of the n.ind elements of pa is a matrix w n.snps rows and 4 columns, each for the 4 nts possible (acgt)

# rec <- pa.MPR
rec <- pa.ACCTRAN

# str(rec)
snps.rec <- list()
for(i in 1:length(rec)){
  snps.rec[[i]] <- rec[[i]][,3]
}
snps.rec <- do.call("rbind", snps.rec)
rownames(snps.rec) <- c(rownames(snps), c((nrow(snps)+1):((nrow(snps)*2)-1)))
colnames(snps.rec) <- colnames(snps)
# identical(snps.rec[1:nrow(snps),], snps)


##############################################
## get LOCATIONS (branches) of snps subs ?? ##
##############################################
subs.edges <- rep(list(NULL), ncol(snps))
for(i in 1:ncol(snps)){
  snp <- snps.rec[, i]
  subs.logical <- sapply(c(1:nrow(edges)),
                       function(e)
                         snp[edges[e,1]]
                       ==
                         snp[edges[e,2]])
  ## get indices of all edges containing a substitution
  subs.total <- which(subs.logical == FALSE)
  ## get df of states of ancestor and descendants nodes on these edges
  df <- data.frame(snp[edges[subs.total,1]], snp[edges[subs.total,2]])
  names(df) <- c("anc", "dec")
  ## get indices of all edges w a positive sub (0 --> 1)
  subs.pos <- subs.total[which(df$anc==0)]
  ## get indices of all edges w a negative sub (1 --> 0)
  subs.neg <- subs.total[which(df$anc==1)]

  ## get output list
  subs.edges[[i]] <- rep(list(NULL), 3)
  names(subs.edges[[i]]) <- c("total", "pos", "neg")
  if(length(subs.total) > 0) subs.edges[[i]][["total"]] <- subs.total
  if(length(subs.pos) > 0) subs.edges[[i]][["pos"]] <- subs.pos
  if(length(subs.neg) > 0) subs.edges[[i]][["neg"]] <- subs.neg
}

cost3 <- sapply(c(1:length(subs.edges)), function(e) length(subs.edges[[e]]))
head(cost3)

## check:
## for SNP1, does it identify the correct/reasonable branches?
## see plot w edgeCol2--reasonable yes, but not the real simulated branches...
# temp <- rep(0, nrow(edges))
# temp <- replace(temp, subs.edges[[1]], 1)

#######################################
## test for association w phen (ACE) ##
#######################################
ace.score <- list()
## NB: length(subs.edges) == ncol(snps.unique)
for(i in 1:length(subs.edges)){
  ## get the absolute value of
  ## the sum of all the positive (0-->1) and negative (1-->0) subs
  ace.score[[i]] <- abs(sum(phen.diffs[subs.edges[[i]][["pos"]]])) +
    abs(sum(phen.diffs[subs.edges[[i]][["neg"]]]))
}
ace.score <- as.vector(unlist(ace.score))
max(ace.score) # 5.96
min(ace.score) # 0

#par(mar=c(5,4,4,1)+0.1)
hist(ace.score)

table(round(ace.score, 2))

head(round(ace.score, 2), 20)

###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###



#########################
## MAXIMUM LIKELIHOOD? ##
#########################
snps.nt <- replace(snps, which(snps==0), "A")
snps.nt <- replace(snps.nt, which(snps==1), "G")

dna.bin <- as.DNAbin(dna)
D <- dist.dna(dna.bin, model="raw")

dna4 <- as.phyDat(dna.bin)
tre.ini <- nj(D)
fit.ini <- pml(tre.ini, dna4, k=n.ind)
fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE,
                 optQ = TRUE, optGamma = TRUE)

# anova(fit.ini, fit)
# AIC(fit.ini)
# AIC(fit)

tree.ml <- fit$tree
#tree.ml <- midpoint(ladderize(tree.ml))
tree.ml <- midpoint(tree.ml)

## plot
plot(tree.ml, show.tip=TRUE, edge.width=2)
title("Maximum-likelihood tree")
axisPhylo()

str(tree.ml)

#fit$lv # what is this?
#



#


###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

#########
## ape ##   ###   ###   ###   ###   ###   ###   ###   ###   ###
#########
require(ape)

mpr <- MPR(phen, unroot(tree), outgroup=1)
str(mpr)

#############
## get.mpr ##
#############
get.mpr <- function(snps){

  ###############################
  ## get unique SNPs patterns: ##
  ###############################
  tab.out <- table.matrix(t(snps))
  snps.unique <- t(as.matrix(tab.out$unique.data))
  index <- tab.out$index
  #   if(ncol(snps.unique) == ncol(snps)){
  #     all.unique <- TRUE
  #   }else{
  #     all.unique <- FALSE
  #   }

  ## work w only unique snps:
  snps.ori <- snps
  snps <- snps.unique

  n.ind <- nrow(snps)

  mpr <- list()
  # system.time # 35 elapsed
  # w 2500 unique snps
  for(i in 1:ncol(snps)){
    mpr[[i]] <- MPR(snps[,i], unroot(tree), outgroup=n.ind)
    ## MPR returns a df w states for all internal nodes
    ## the df has 2 columns for the lower and upper bound of these states
    ## bc it treats the snp 0/1 as an integer value.
    ## NOTE--we hope that all lower/upper pairs will be the same
    ## (ie. that MPR is confident in the SNP selected... ).
    ## If not ??????
  }

  return(mpr)
} # end get.mpr

foo <- mpr[[1]]
foo

edges <- tree$edge
tail(edges)

## check if all upper and lower MPR estimates are the same:
temp <- sapply(c(1:nrow(foo)),
               function(e) foo[e,1]==foo[e,2])
## identify which are different...
diffs <- which(temp==FALSE)
## deal with these later!
foo[diffs,]

## get a vector of SNPi states for all internal nodes
## (use MPR lower estimate for now)
states.int <- foo[,1]
###########
## NOTE--THIS IS SPECIFIC TO YOUR DATASET ###################
## 1 = SNP state for ind n.ind (ie 100) thus c(x, 0)
states.int <- c(foo[,1], 1) ## NEED to FIX OUtGROUP PROBLEM!
## reorder??
noms <- names(states.int)
## remove "NODE_" from noms
## ie. NODE_112 --> 112
noms <- as.numeric(removeFirstN(noms, 5))

## get a vector of SNPi states for all terminal nodes
x <- 1
states.term <- snps[,x]
## concatenate to get vector of SNPi states for all nodes
states <- c(states.term, states.int)

## go through edges from root down to identify  branches where change occurs.
subs <- list()
for(i in 1:nrow(edges)){
  state.anc <- states[edges[i,1]]
  state.dec <- states[edges[i,2]]
  subs[[i]] <- state.anc==state.dec
}
subs <- as.vector(unlist(subs))
subs.edges <- which(subs==FALSE)
subs.edges
length(subs.edges) # 40...
## should only be 10 edges w a sub for snp 1...

################
## ape MPR eg ##
################
## the example in Narushima and Hanazawa (1997):
tr <- read.tree(text = "(((i,j)c,(k,l)b)a,(h,g)e,f)d;")
x <- c(1, 3, 0, 6, 5, 2, 4)
names(x) <- letters[6:12]
o <- MPR(x, tr, "f")
plot(tr)
nodelabels(paste("[", o[, 1], ",", o[, 2], "]", sep = ""))
tiplabels(x[tr$tip.label], adj = -2)

## some random data:
x <- rpois(30, 1)
tr <- rtree(30, rooted = FALSE)
p <- MPR(x, tr, "t1")
plot(tr, type="c", use.edge.length=FALSE)
nodelabels(paste("[", p[, 1], ",", p[, 2], "]", sep = ""), cex=0.7)
tiplabels(x[1:length(tr$tip.label)], adj = -5, cex=0.7)


#


###############
## plot tree ##
###############

## ACE ##
## DISCRETE
var <- snps[,1]

snps.ace.d <- ace(var, tree, type="discrete")
probs.d <- snps.ace.d$lik.anc[,2]*100
## convert probs to show on tree pies:
probs.ori <- probs.d
df <- data.frame(100-probs.ori, probs.ori)
names(df) <- c("prob.A", "prob.B")
lab <- as.character(tree$tip.label)
temp <- sapply(c(1:length(lab)),
               function(e)
                 paste(
                   paste(rep(c("-"), 3-nchar(lab[e])), collapse=""),
                   "-", lab[e], sep=""))
tiptext <- as.vector(unlist(temp))
par(mar=c(1,1,1,2)-0.5)
plot(tree, show.tip=FALSE, edge.width=2,
     edge.color=edgeCol, type="c", use.edge.length=FALSE)
tiplabels(text=tiptext, cex=0.6, adj=c(-.4, 0), col=leafCol, frame="none")
## try plotting ~ example: ##
#### Showing the probs and inverse probs on each node:
myCol <- c("blue", "red")
## convert phen to 1s and 2s (to indicate myCol 1 or 2)
tipCol.pattern <- as.numeric(as.factor(as.character(var[1:100])))
## plot terminal squares
tiplabels(pch = 22, bg = myCol[tipCol.pattern], cex = 0.7, adj = 1)
## plot internal pies
nodelabels(thermo = df/100, piecol = transp(myCol, 0.8), cex = 0.25, adj=1.5)
##########
# nodelabels(paste("[", foo[, 1], ",", foo[, 2], "]", sep = ""), cex=0.7)
# tiplabels(states.term[1:length(tree$tip.label)], adj = -5, cex=0.7)
##########
nodelabels(states.int, cex=0.5)
tiplabels(states.term, adj = -5, cex=0.7, frame="none")
##########
edgeCol2 <- replace(temp, which(temp==0), "black")
edgeCol2 <- replace(edgeCol2, which(temp==1), "yellow")
edgelabels(text=temp,
           cex=0.7, font=2, bg=transp(edgeCol2, 0.3))
##########
plot(tree, cex=0.5)
edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."),
          cex=0.5, font=2, bg=transp("yellow", 0.3), adj=c(1,1))
nodelabels(text=rev(unique(tree$edge[,1])), cex=0.5, bg=transp("blue", 0.3), adj=2)
## should be numbered s.t. the root node is n.term+1
## RECALL: terminal nodes are numbered 1:n.ind from bottom to top of plot of tree;
## edges are numbered 1:nrow(edges) by following the lowest trace on the plot??
## (starting from the root down to the lowermost tips);
## thus, internal nodes are numbered (n.ind+1):(n.ind+(n.ind-1)),
## from root to the top-most internal node to be connected (ie. the highest in the plot)



# ## branch and bound? ##
# ## WAYYY TO SlOW!
# # dna = phyDat
# tree.bab <- bab(dna) # trace=0
# str(tree.bab)
#
# cost.bab <- parsimony(tree.bab, dna)
#



#



###########################################################################




################
## QUESTIONS: ##
################

## Does directionality matter or can we work w abs val??
## --> it will matter for comparison btw snps and phen ace.diffs...

## How to ensure we are working w the right column of var.internal??

## Work w var indices or names?
## NB: tree$edge works w indices



##








#







#
