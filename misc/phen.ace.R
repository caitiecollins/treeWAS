

# setwd("/media/caitiecollins/Seagate Backup Plus Drive/SimBac")

prefix <- ("CFML_R_0.1")

require(ape)

tree<-read.tree(sprintf('%s.labelled_tree.newick',prefix))
seqs<-read.dna(sprintf('%s.ML_sequence.fasta',prefix),format='fasta')
mapping<-scan(sprintf('%s.position_cross_reference.txt',prefix),sep=',',quiet=T)
l<-length(seqs[1,])

## check & convert tree!
## (violations cause problems, eg. in phen.sim)
if(!is.binary.tree(tree)) tree <- multi2di(tree)
if(!identical(tree, reorder.phylo(tree,"postorder"))) tree <- reorder.phylo(tree,"postorder")

## plot tree
plot(tree, cex=0.6)
edgeLabCol <- "purple"
edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."),
           cex=0.7, font=2,
           frame="none", adj=c(1,1),
           col="purple") # bg=transp(edgeLabCol, 0.3), col="black"
nodelabels(text=c(1:nrow(tree$edge)),
           cex=0.7, font=2,
           frame="none", adj=c(1,1),
           col="black")
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
edgeLabCol <- "purple"
edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."),
           cex=0.7, font=2,
           frame="none", adj=c(1,1),
           col="purple") # bg=transp(edgeLabCol, 0.3), col="black"
nodelabels(text=c(1:tree$Nnode),
           cex=0.5, font=1,
           frame="rect", adj=c(1,1),
           col="black", bg=transp("yellow", 0.3))
## match names to tip labels
# tree$tip.label
names(phen) <- c(0:(length(phen)-1))
# table(phen)

###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

#########
## ACE ##
#########
phen.ace <- ace(phen, tree, type="discrete")
str(phen.ace)
lik1 <- phen.ace$lik.anc

phen.ace2 <- ace(phen, tree, type="discrete", marginal=TRUE)
str(phen.ace2)
lik2 <- phen.ace2$lik.anc


## change names to start at 1 not 0...
phen.ori <- phen
tree.ori <- tree # note: this is AFTER postordering the initial tree
names(phen)
as.numeric(tree$tip.label)

names(phen) <- as.numeric(names(phen))+1
phen.ori1 <- phen
tree$tip.label <- as.numeric(tree$tip.label)+1
tree$tip.label <- as.character(tree$tip.label)
str(tree)

## try reordering phen to match tree tips
## (by order of labels, not as plotted)
phen <- phen[as.numeric(tree$tip.label)]

## retry ACE w new phen order...
lik1.ori <- lik1
lik2.ori <- lik2
phen.ace <- ace(phen, tree, type="discrete")
str(phen.ace)
lik1 <- phen.ace$lik.anc
lik1

phen.ace2 <- ace(phen, tree, type="discrete", marginal=TRUE)
str(phen.ace2)
lik2 <- phen.ace2$lik.anc
lik2

## compare...
identical(lik1.ori, lik1) #TRUE ## Must be running ACE using names.
identical(lik2.ori, lik2)


## how to interpret ace lik output??
lik1
diffs <- sapply(c(1:nrow(lik1)),
                function(e)
                  abs(lik1[e,1] - lik1[e,2]))
which.max(diffs) # 52

## which node is row 52 identifying?
tree$edge
tree$tip.label[c(48, 49)] # 42 29
phen[c(42, 29)] # both A!
phen.ori1[c(42, 29)] # one A, one B

## Hyp: ace is confused and is reading tip.labels as phen indices?
## try renaming the reordered phen with 1:100
phen.ori2 <- phen
names(phen) <- c(1:100)
phen

## retry ACE w new phen order but old names...
lik1.ori1 <- lik1
lik2.ori1 <- lik2
phen.ace <- ace(phen, tree, type="discrete")
str(phen.ace)
lik1 <- phen.ace$lik.anc
lik1
## Oops... all 0.5 now...

phen.ace2 <- ace(phen, tree, type="discrete", marginal=TRUE)
str(phen.ace2)
lik2 <- phen.ace2$lik.anc
lik2

###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

##################
## ACE example: ##
##################
data("bird.orders") # a tree
### For discrete characters:
x <- factor(c(rep(0, 5), rep(1, 18)))
ans <- ace(x, bird.orders, type = "d")
ans$lik.anc

#### Showing the likelihoods on each node:
plot(bird.orders, type = "c", FALSE, label.offset = 1)
co <- c("blue", "yellow")
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 2, adj = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.75)

## interpret?
lik.eg <- ans$lik.anc
diffs <- sapply(c(1:nrow(lik.eg)),
                function(e)
                  abs(lik.eg[e,1] - lik.eg[e,2]))
which.max(diffs) # 16
diffs[16] #1

## BUT--on visual inspection, row 1 looks to contain the biggest
## (and only real sig diff) of the set!!
## SOLN--row 16 is identified bc the lik.anc
## (ie. the scaled log likelihood) of state 1 is SO insignificantly SMALL
## that it is interpreted as diff btw 0 and 1 --> abs(0-1) --> 1!!
## Whereas, the REAL biggest diff is on row 1 btw 0.67 and 0.33 --> 0.33!!!

str(bird.orders) # tree w 23 tips, 22 internal nodes
lik.eg # df w 2 cols, 22 rows (one for each internal node, presumably)

###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

## Trying again!
## working w phen w new order but original names (in their new location)
phen <- phen.ori2

## retry ACE
phen.ace <- ace(phen, tree, type="discrete")
str(phen.ace)
lik1 <- phen.ace$lik.anc
lik1

phen.ace2 <- ace(phen, tree, type="discrete", marginal=TRUE)
str(phen.ace2)
lik2 <- phen.ace2$lik.anc
lik2

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
## try plotting ~ example: ##
#### Showing the likelihoods on each node:
# plot(tree, type = "c", FALSE, label.offset = 1, cex=0.55)
myCol <- c("blue", "red")
tiplabels(pch = 22, bg = myCol[as.numeric(phen)], cex = 0.7, adj = 1)
nodelabels(thermo = lik1, piecol = transp(myCol, 0.8), cex = 0.25)
nodelabels(thermo = lik2, piecol = transp(myCol, 0.8), cex = 0.25)

###########

# ## retry ACE specifying transition rates...
# lik1.ori <- lik1
# phen.ace <- ace(phen, tree, type="discrete", model=matrix(c(0, 1, 1, 0), 2))
# str(phen.ace)
# lik1 <- phen.ace$lik.anc
# lik1
# identical(lik1, lik1.ori) # TRUE


###########

################
## CONTINUOUS ##
################

## Trying ACE after converting phen to continuous [0,100]
## (phen 0 --> 0, phen 1 --> 100)

## HYP: We could treat the inferred states of internal nodes as
## the probabilities of the true binary states??

phen.ini <- phen.ori
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
phen.ace <- ace(phen, tree, type="continuous")
str(phen.ace)
## get probs
## (ie. get the inferred continuous var value,
## which we will interpret as the probability
## of being phen state B).
probs <- phen.ace$ace
probs

# phen.ace2 <- ace(phen, tree, type="continuous", marginal=TRUE)
# str(phen.ace2)
# probs2 <- phen.ace2$ace
# identical(probs, probs2) # TRUE

# phen.ace2 <- ace(phen, tree, type="continuous", scaled=FALSE)
# str(phen.ace2)
# probs2 <- phen.ace2$ace
# identical(probs, probs2) # TRUE

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
tipCol.pattern.ori <- as.numeric(as.factor(as.character(phen)))
## reorder
for(i in 1:length(tipCol.pattern.ori)){
  tipCol.pattern[i] <- tipCol.pattern.ori[which(names(phen) == i)]
}
names(tipCol.pattern) <- c(1:length(tipCol.pattern))
## plot terminal squares
tiplabels(pch = 22, bg = myCol[tipCol.pattern], cex = 0.7, adj = 1)
## plot internal pies
nodelabels(thermo = df/100, piecol = transp(myCol, 0.8), cex = 0.25)


##########################################################
## UHOH--pies make it look like ACE cont is confused... ##
##########################################################
## Hint: pies make sense/match the previously-incorrect tip square colours.
## Soln: Need to rearrange phen (was in right order, but its names
## needed to be replaced with tree$tip.label) ...

phen.temp <- phen

phen.new.ori <- phen
#phen.new <- phen.new.ori
## reorder
for(i in 1:length(phen.new.ori)){
  phen.new[i] <- phen.new.ori[which(names(phen) == i)]
}
names(phen.new) <- tree$tip.label
phen.new

phen <- phen.new

tipCol.pattern ## represents phen in order plotted on tips of tree (bottom to top)
## in our case, this happens to be the order of the tree$tip.label

# tree$tip.label
# get.tip.order(tree)

##################################
## Try ACE as continuous var... ##
##################################
phen.ace <- ace(phen, tree, type="continuous")
str(phen.ace)
## get probs
probs <- phen.ace$ace
probs

###########################
## INDEPENDENT CONTRASTS ##
###########################
phen.ace2 <- ace(phen, tree, type="continuous", method="pic", scaled=TRUE)
str(phen.ace2)
probs2 <- phen.ace2$ace
identical(probs, probs2) # FALSE # (though not very different...)

phen.ace3 <- ace(phen, tree, type="continuous", method="pic", scaled=FALSE)
str(phen.ace3)
probs3 <- phen.ace3$ace
identical(probs2, probs3) # TRUE
## ?? Why doesn't scaling have any impact on the IC ace??

###############
## PLOT TREE ##
###############

## convert probs to show on tree pies:
probs.temp <- probs2
df <- data.frame(100-probs.temp, probs.temp)
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
## plot terminal squares
tiplabels(pch = 22, bg = myCol[tipCol.pattern], cex = 0.7, adj = 1)
## plot internal pies
nodelabels(thermo = df/100, piecol = transp(myCol, 0.8), cex = 0.25)


#################
## CONCLUSION: ##
#################

## Finally got ACE working and interpreting it correctly,
## BUT--as most of the scaled likelihoods are 0.5/0.5,
## it seems like this may not be the best method for our purposes
## (ie. for determining whether a phenotypic change has been likely
## across the same branch as a SNP substitution).

## I thought we were hoping for something a little more deterministic/
## concrete in its output (ie. something that could give us the inferred
## phen state at all internal nodes, even if we were only planning on
## working with the probabilities of those states). It seems like ACE
## can not do that (in the code implemented above, at least...).

## Continuing search for something that can do ASR differently...

## Side note: ACE (fn ace) can do Idependent Contrasts!!!
## Just set method = "pic" (and scaled = TRUE)
#

