
#####################
## get.fitch.n.mts ##
#####################
## phangorn-based fitch fn

get.fitch.n.mts <- function(snps, tree){
  require(phangorn)
  ## returns only unique patterns...
  snps.phyDat <- as.phyDat(as.matrix(snps),
                           type="USER", levels=c(0,1))
  ## get index of all original snps columns to map to unique pattern
  phyDat.index <- attr(snps.phyDat, "index")

  fitch.phangorn <- phangorn::fitch
  ## get parsimony score for all unique patterns in snps
  fitch.unique <- fitch.phangorn(tree, snps.phyDat, site="site")

  ## get score for all original sites
  ## eg.
  ## There are 10000 origial snps columns,
  ## but only 6543 unique snps patterns.
  ## 174/6543 are not 1 (so they are repeated columns).
  ## If index 5 is repeated 96 times, we need to repeat
  ## the Fitch score for this pattern 96 times.
  fitch.complete <- rep(NA, ncol(snps))
  for(i in 1:length(fitch.unique)){
    fitch.complete[which(phyDat.index == i)] <- fitch.unique[i]
  }
  return(fitch.complete)
} # end get.fitch.n.mts



# str(fitch.complete)
# table(fitch.complete)







####################
## Xav's fitch fn ##
####################

# ## example:
# library(ape)
# set.seed(1)
# ## get tree
# # tree has to be binary and have postorder ordering. If the former is
# # not true, use tree=multi2di(tree). If the latter is not true use
# # tree=reorder.phylo(tree,'postorder')
# tree <- rtree(1000)
# ## get data
# # data is vector of 0s and 1s
# snps <- sample(c(0,1),1000,T)
# plot(tree)
# str(snps)
# ## run fitch fn
# # gets the total minimum cost of the root
# out <- fitch(tree,snps)
# out


fitch<-function(tree,data){
  #   tree has to be binary and have postorder ordering. If the former is
  #   not true, use tree=multi2di(tree). If the latter is not true use
  #   tree=reorder.phylo(tree,'postorder')
  #   data is vector of 0s and 1s
  ## CHECKS
  if(any(data > 1 | data < 0)) stop("Your data is not binary.")
  if(!is.binary.tree(tree)) tree <- multi2di(tree)
  if(!identical(tree, reorder.phylo(tree,"postorder"))) tree <- reorder.phylo(tree,"postorder")

  ## CORE FITCH FN
  ## (for a single vector of 0s and 1s)
  f1 <- function(tree, data){
    n<-tree$Nnode+1
    br<-tree$edge
    children=matrix(NA,n*2-1,2)#Matrix of children in tree
    for (i in 1:(nrow(br))) if (is.na(children[br[i,1],1]))
      children[br[i,1],1]=br[i,2] else children[br[i,1],2]=br[i,2]
    out<-matrix(Inf,n*2-1,2)#Store cost for a given node and type
    for (i in 1:n) out[i,data[i]+1]=0
    for (i in (2*n-1):(n+1)) {
      out[i,1]=min(out[children[i,1],]+c(0,1))+min(out[children[i,2],]+c(0,1))
      out[i,2]=min(out[children[i,1],]+c(1,0))+min(out[children[i,2],]+c(1,0))
    }
    return(min(out[n+1,]))#Final cost is minimum cost of root
  } # end f1

  ## RUN FN
  ## On single vector data:
  if(is.vector(data)){
    f1(tree, data)
  }else{
    if(is.matrix(data)){
      ## On matrix data
      sapply(c(1:ncol(data)), function(e) f1(tree, data[,e]))
    }
  }
} # end fitch


## get unique patterns from phyDat to matrix
lala <- do.call(rbind, snps.phyDat[1:100])
toto <- lala-1

t1 <- system.time(fitch.phangorn(tree, snps.phyDat, site="site"))
t2 <- system.time(fitch(tree, toto))
# > t1
# user  system elapsed
# 0.01    0.00    0.02
# > t2
# user  system elapsed
# 8.48    0.00    8.50

## Even running Xavier's fitch algorithm on just unique patterns,
## the phangorn fitch fn is much faster...
## Both algorithms return the same result.
## THOUGH-- we note that this result is much more parsimonious than
## the distribution used as initial input to coalescent.sim
## Again--this may reflect the fact that this dist was drawn from
## the (MUCH LARGER) S.aureus CFml example data...


hist(fitch.complete)
hist(x.fitch)

barplot(dist.prop)
hist(n.mts, breaks=20)
#########
## Xavier's method:
x.fitch <- sapply(c(1:ncol(snps)), function(e) fitch(tree, snps[,e]))
str(x.fitch)
table(x.fitch)

all(table(fitch.complete) == table(x.fitch)) # TRUE

#########
## OLD METHOD (in treeWAS, incomplete)
# mt.rate <- parsimony(tree.ori, snps.phyDat, method="fitch")
