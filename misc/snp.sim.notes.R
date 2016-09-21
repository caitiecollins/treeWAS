
####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####

##################################################
######## toy eg - reversability justification: ##
##################################################

# phen <- c(0,0,0,1,1,1,0,0,1,0,1,1,1,1,0,0,1,1,1,1) # n.0 = 8 # n.1 = 12
#
# snp1 <- c(0,0,0,1,1,1,0,0,1,0,1,1,0,0,1,1,1,1,1,1) # n.0 = 8 # n.1 = 12
# snp2 <- c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,1,1,1,1) # n.0 = 12 # n.1 = 8
# snp3 <- c(0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1) # n.0 = 4 # n.1 = 16
#
# cor(phen, snp1) # 0.58333
# cor(phen, snp2) # 0.66667
# cor(phen, snp3) # 0.61237

####################################################################################################################################




####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####

############    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####
## sim 2 ###    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####
############    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####

####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####

## COMPLEMENTARY ASSOCIATIONS:
## Each SNP-phen assoc only accounts for ~ 1/2 the assoc's on the tree (by major clade).

## get clades --> numbers of descendants in successive major clades:

## phytools (?) fn from Revell blog:
## http://blog.phytools.org/2012/01/function-to-get-descendant-node-numbers.html
getDescendants<-function(tree,node,curr=NULL){

  if(is.null(curr)) curr <- vector()

  daughters <- tree$edge[which(tree$edge[,1]==node),2]
  curr <- c(curr,daughters)

  w <- which(daughters>=length(tree$tip))
  if(length(w) > 0) for(i in 1:length(w))
    curr <- getDescendants(tree, daughters[w[i]], curr)

  return(curr)
} # end getDescendants

# getDescendants(tree, node=104)
getDescendants(tree, node=102)

tree <- tree.ori
tree$edge
node <- 199




## get/load tree: ##
# tree <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_31_tree.Rdata"))
tree <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_29_tree.Rdata"))

edges <- tree$edge

# dec <- getDescendants(tree, node=edges[nrow(edges), 1])
# length(dec)


## OR -- use cutree to get numbers at every level!

# require(graphics)
# hc <- hclust(dist(USArrests), "ave")
# plot(hc, hang = -1)
# str(hc)

# ?cutree

## Get tree as hclust tree:
tree.hc <- as.hclust.phylo(tree)
plot(tree.hc)
str(tree.hc)

## Want to divide tree into 2 sets of clades btw 1/3:2/3 and 1/2:1/2
clades <- tab <- grp.options <- clade1 <- clade2 <- sets.complete <- list()

min.size <- ceiling((tree$Nnode+1)*(1/3))
max.size <- floor((tree$Nnode+1)*(2/3))
grp1 <- tree$Nnode+1

i <- 2
counter <- 0

while(grp1 < min.size | grp1 > max.size){
## FOR LOOP to get size of clades:
# for(i in 2:(nrow(edges)/2)){
  clades[[i]] <- cutree(tree.hc, k=i)
  tab[[i]] <- table(clades[[i]])
  grp.opts <- grp.options[[i]] <- sapply(c(1:(i-1)), function(e) sum(tab[[i]][1:e]))
  ## make grp1 first clade in grp.options:
  group1 <- grp.opts[1]
  ## remove first clade from options:
  grp.opts <- grp.opts[-1]
  ## and record n.grps:
  n.grp <- 1
  ## try to identify a (set of) clade(s) that's big enough (but not too big):
  while(group1 < min.size){
    group1 <- sum(group1, grp.opts[1])
    grp.opts <- grp.opts[-1]
    n.grp <- n.grp+1
  }
  clade1[[i]] <- clades[[i]][which(clades[[i]] %in% (1:n.grp))]
  clade2[[i]] <- clades[[i]][which(!clades[[i]] %in% (1:n.grp))]
  sets.complete[[i]] <- replace(clades[[i]], which(clades[[i]] %in% (1:n.grp)), 1)
  sets.complete[[i]] <- replace(sets.complete[[i]], which(!clades[[i]] %in% (1:n.grp)), 2)
  grp1 <- sum(grp.options[[i]][1:n.grp])
  k <- i
  i <- i+1
  counter <- counter+1
# } # end for loop
} # end while loop
###########


# clades[[length(clades)]]
# clade1

sets <- sets.complete[[length(sets.complete)]]

set1 <- names(sets)[which(sets == 1)]
set2 <- names(sets)[which(sets == 2)]

## PLOT PHEN on tREE + CLADES ALONG TIPS: ##
############################################
## load/get phen.rec:
# res <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_31_res.Rdata"))
res <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_29_res.Rdata"))
phen.rec <- res$dat$phen.rec
# phen.rec <- out$res[[1]]$dat$phen.rec

## convert phen.rec to A/B:
phen.rec <- as.character(phen.rec)
phen.rec <- replace(phen.rec, which(phen.rec == 0), "A")
phen.rec <- replace(phen.rec, which(phen.rec == 1), "B")

## plot phen on tree:
plot.phen(tree, phen.nodes=phen.rec)


## Get CLADES:
set1 <- names(sets)[which(sets == 1)]
set2 <- names(sets)[which(sets == 2)]
cladeCol <- rep(NA, length(tree$tip.label))
cladeCol <- replace(cladeCol, which(tree$tip.label %in% set1), "darkgreen")
cladeCol <- replace(cladeCol, which(tree$tip.label %in% set2), "lightgreen")
tiplabels(text=NULL, cex=0.6, adj=c(0.65, 0.75), col=cladeCol, pch=15)
#


# perf <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_31_performance.Rdata"))
perf <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_29_performance.Rdata"))
snps.assoc <- perf$snps.assoc
#

snps <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_29_snps.Rdata"))




#





###############

####################################################################################################################################




####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####

############    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####
## sim 3 ###    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####
############    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####

####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####

## SUBSEQUENT/COMPENSATORY SUBS:
## Variable rates of simultaneous subs/maintained assoc's/subsequent subs contributing to associations...



library(phytools)

## adopt something like sim.history but with probs/rates not ML-based?
# http://blog.phytools.org/2014/12/r-function-for-pagels-1994-correlation.html

sim.history

?reorder.phylo
str(tree)
str(tree.hc)
tree.cw <- reorder.phylo(tree, order="cladewise")
str(tree.cw)

plot(tree.cw, cex=0.6)
edgelabels(text=c(1:nrow(tree.cw$edge)), frame="none", cex=0.5, col="blue")

## UNcorrelated transition rates btw SNP and phen states:
Q.u <- matrix(c(-1,1,1,-1),2,2)
Q.u
## correlated transition rates btw SNP and phen states:
Q.c <- matrix(c(0,0.5,0.5,0,2,0,0,2,2,0,0,2,0,0.5,0.5,0),4,4,byrow=TRUE)
rownames(Q.c)<-colnames(Q.c)<-c("aa","ab","ba","bb")
Q.c

#


##############################################################################
## Can we identify a formula for Exp(n.subs | s, af, n.ind, tree$edge...)?? ##
##############################################################################

## Can we work within the framework of our process as a Markov chain?

## EG.

Q <- matrix(c(0.9, 0.075, 0.025,
            0.15, 0.8, 0.05,
            0.25, 0.25, 0.5),
            3, 3, byrow=T)

Q%*%Q%*%Q

Q.mat <- Q
diag(Q.mat) <- sapply(c(1:nrow(Q.mat)), function(e) -sum(Q.mat[e, c(1:ncol(Q.mat))[-e]]))

matexpo(Q.mat*3)
?matexpo

temp <- matexpo(Q.mat*0.1)
rownames(temp) <- colnames(temp) <- c("0|0", "0|1", "1|0", "1|1")
temp







#






#





###############

####################################################################################################################################
