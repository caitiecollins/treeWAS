


##################################################################################################################


##################################
## TRoUBLESHOOTING rtree ISSUES ##
##################################

## Run simTest up to treeWAS output:
## Running set 1 (seed 41) w assoc.prrob = 100, n.snps=1000, n.snps.assoc=10000...



## check terminal results:
nrow(res$terminal$pval.0.01.bonf.count.10.x.n.snps$sig.snps)


corr.sim <- out$vals$terminal$corr.sim
corr.dat <- out$vals$terminal$corr.dat

thresh <- res$terminal$pval.0.01.bonf.count.10.x.n.snps$sig.thresh

length(which(corr.dat > thresh))

## check all treeWAS scores:

hist(out$vals$terminal$corr.sim, breaks=10, xlim=c(0,1), col=transp("red", 0.5),  freq=T, add=T)
hist(out$vals$terminal$corr.dat, breaks=40, xlim=c(0,1), col=transp("blue", 0.5), freq=T) #  add=T


hist(out$vals$simultaneous$corr.sim[1:1000], breaks=40,
     xlim=c(0, max(max(out$vals$simultaneous$corr.sim), max(out$vals$simultaneous$corr.dat))),
     col=transp("red", 0.5),  freq=T)
hist(out$vals$simultaneous$corr.dat, breaks=40,
     xlim=c(0, max(max(out$vals$simultaneous$corr.sim), max(out$vals$simultaneous$corr.dat))),
     col=transp("blue", 0.5), add=T, freq=T)

hist(out$vals$subsequent$corr.sim[1:1000], breaks=40,
     xlim=c(0, max(max(out$vals$subsequent$corr.sim), max(out$vals$subsequent$corr.dat))),
     col=transp("red", 0.5),  freq=T)
hist(out$vals$subsequent$corr.dat, breaks=40,
     xlim=c(0, max(max(out$vals$subsequent$corr.sim), max(out$vals$subsequent$corr.dat))),
     col=transp("blue", 0.5), add=T, freq=T)


## Compare input n.subs (= dist_0) to output fitch n.subs: ##

## input n.subs dist:
barplot(dist_0, names.arg=c(1:length(dist_0)), col=transp("royalblue", 0.5))

## output n.subs dist:
fitch.complete <- get.fitch.n.mts(snps, tree)
hist(fitch.complete, col=transp("red", 0.5))

## Uhoh...


## Check fitch score for snps.assoc vs n.phen.subs:
fitch.complete[snps.assoc] # all = 23

n.phen.subs <- length(which(phen.plot.col$edges == "grey"))
n.phen.subs # = 11...

## Shit, if anything fitch should UNDER-estimate n.phen.subs!

fitch.complete.41 <- fitch.complete

###########
## Double check that fitch is working for coalescent tree first...

tree.ori <- tree

c.sim <- coalescent.sim(n.ind=100,
                        n.snps=10000,
                        n.subs=dist_0,
                        n.snps.assoc=10,
                        assoc.prob=100,
                        n.phen.subs=15,
                        phen=NULL,
                        plot=TRUE,
                        heatmap=FALSE,
                        reconstruct=FALSE,
                        dist.dna.model="JC69",
                        grp.min = 0.25,
                        row.names=NULL,
                        set=1,
                        coaltree = TRUE,
                        s = 1,
                        af = 2,
                        filename = NULL,
                        seed=41)

tree.c <- tree <- c.sim$tree
snps.c <- snps <- c.sim$snps
phen.c <- phen <- c.sim$phen
length(which(c.sim$phen.plot.col$edges == "grey")) # 11
snps.assoc.c <- snps.assoc <- c.sim$snps.assoc


## Compare input n.subs (= dist_0) to output fitch n.subs: ##

## input n.subs dist:
barplot(dist_0, names.arg=c(1:length(dist_0)), col=transp("royalblue", 0.5))

## output n.subs dist:
fitch.complete <- get.fitch.n.mts(snps, tree)
hist(fitch.complete, col=transp("red", 0.5)) ## Good! Approx same shape as input dist!

fitch.complete[snps.assoc] # 9
## Also good! Just under-estimates the real number (11)...

## QED something is definitely wrong w rtree that is fine for coaltree...

tree$tip.label
plot(tree, cex=0.5)

get.tip.order(tree)

tree.r <- tree.41
# tree <- tree.c
# tree.ord <- reorder.phylo(tree, order="postorder")
# str(tree.ord)
## None of the reorder order options --> tree$tip.label being in order presented!

str(tree.r)
tree.r2 <- reorder.phylo(tree.r, order="cladewise")
str(tree.r2)

tree.r.ori <- tree.r
str(snps.41)
tree <- tree.r
tree$tip.label <- c(1:100) ## REPLACE tree$tip.lab w simple c(1:100) before running fitch!

## try fitch again:
fitch.complete <- get.fitch.n.mts(snps.41, tree)
hist(fitch.complete) ## SOLVED! -- looks like the hist for the coalescent tree ( and thus ~ like the input n.subs dist_0)


############################### ############################################################################################
## BUT -- PROBLEM REMAINS !! ##
###############################
## WHAT IF USER'S REAL DATA CONTAINS rownames(snps) as a set of character labels in a given order
## and tree$tip.label cotains the same set of character labels in a DIFFerent order!

## You will NOT be able to solve their problem simply by replacing tree$tip.label w c(1:N)
## and WORSE, you WILL HAVE NO IDEA THERE IS A PROBLEM!

############################################################################################################################


# ## try fitch again w modified snps rownames?
# rownames(snps) <- rev(c(1:100))
# fitch.complete <- get.fitch.n.mts(snps, tree)
# hist(fitch.complete) ## still works fine :(

# ## try making tree tip.labs numeric?
# # tree.r.new <- tree
# tree <- tree.r
# ## if we can remove "t" to get numeric, do so:
# prefix <- keepFirstN(tree$tip.label, 1)
# if(all(tolower(prefix) == "t")){
#   temp <- removeFirstN(tree$tip.label, 1)
#   if(all.is.numeric(temp))  tree$tip.label <- as.numeric(temp)
# }
#
# snps <- snps.41
#
# ## try fitch again w numeric original tree$tip.labs?
# fitch.complete <- get.fitch.n.mts(snps, tree)
# hist(fitch.complete)
# ## NOPE--fails again!
#
# ## Hyp-- fitch interprets the indices in tree$edge as the rownames(snps) to map to.
# ## Changing the cells of tree$edge to be the numeric tree$tip.label they map to would probably work? ## NO THIS FAILS TOO!
# ## Or possibly changing the rownames(snsp) to match the numeric tree$tip.labs? (NO, fails)
# ## And doing BOTH of the above together ALSO FAILS!
#
# # rownames(snps) <- tree$tip.label
#
# tree$edge <- replace(tree$edge, which(tree$edge %in% c(1:length(tree$tip.label))), tree$tip.label)
# snps <- snps.41
# rownames(snps) <- tree$tip.label
#


#




#


##
##################################################################################################################

## Compare snps & snps.sim:
tree$tip.label <- c(1:100)
fitch.complete <- get.fitch.n.mts(snps, tree)
hist(fitch.complete)

snps.sim <- out$dat$snps.sim
fitch.complete <- get.fitch.n.mts(snps.sim, tree)
hist(fitch.complete)

names(out$dat)



##################################################################################################################
##################################################################################################################
##################################################################################################################

###################################################################
## Try running pipeline on a coalescent tree generated w rcoal?! ##
###################################################################

## (i) Check that rtree and rcoal give same tree label ordering
## (ii) Run w rcoal tree & coal=FALSE
######  (OR better yet, write in code so you can distinguish appropriate action given tree structure AND label order!)

set.seed(1)
tree.r <- rtree(n=10, rooted=TRUE)
set.seed(1)
tree.c <- rcoal(n=10, rooted=TRUE)

plot(tree.r)
plot(tree.c)


str(tree.r)
str(tree.c)

tree.ori <- tree
phen.ori <- phen
snps.ori <- snps

head(tree$edge)

###################################################################

# labs <- rownames(snps)
labs <- c(1:100)

if(all.is.numeric(tree$tip.label)){
  ## if we can convert to numeric, do so:
  if(!is.null(tree$tip.label)) tree$tip.label <- as.numeric(tree$tip.label)
}else{
  ## if we can remove "t" to get numeric, do so:
  prefix <- keepFirstN(tree$tip.label, 1)
  if(all(tolower(prefix) == "t")){
    temp <- removeFirstN(tree$tip.label, 1)
    if(all.is.numeric(temp)){
      tree$tip.label <- as.numeric(temp)
    }else{
      ## else, replace with numeric indices:
      # tree$tip.label <- c(1:length(tree$tip.label))
      warning("Site-wise parsimony scores (phangorn's
              fitch parsimony function) may not be calculated correctly
              when tip.labels are not numeric.
              Please change tree$tip.label to numeric values.")
    }
  }
}

## NODE labels ##
if(all.is.numeric(tree$node.label)){
  if(!is.null(tree$node.label)) tree$node.label <- as.numeric(tree$node.label)
}else{
  ## if we can remove "NODE_" to get numeric, do so:
  prefix <- keepFirstN(tree$node.label, 4)
  if(all(tolower(prefix) == "node")){
    temp <- removeFirstN(tree$node.label, 5)
    if(all.is.numeric(temp)){
      tree$node.label <- as.numeric(temp)
    }else{
      ## else, replace with numeric indices:
      # tree$node.label <- c((n.ind+1):(n.ind+tree$Nnode))
      warning("Site-wise parsimony scores (phangorn's
              fitch parsimony function) may not be calculated correctly
              when node.labels are not numeric.
              Please change tree$node.label to numeric values.")
    }
  }
}

treelabs <- c(tree$tip.label, tree$node.label)
edges <- tree$edge
# for (i in 1:nrow(edges)) for (j in 1:ncol(edges)) edges[i,j] <- which(labs==treelabs[edges[i,j]]) ## original (from readCFML)

## replace indices in cells of tree$edge with the corresponding index in rownames(snps)
for (i in 1:nrow(edges)){
  for (j in 1:ncol(edges)){
    if(!is.na(treelabs[edges[i,j]])){
      edges[i,j] <- which(labs==treelabs[edges[i,j]])
    }
  }
}



##################################################  Try it w rtree(10)     ##################################################

set.seed(1)
tree.r <- rtree(n=10, rooted=TRUE)

labs <- paste("t", c(1:10), sep="")

#


#




#



#




#






#







#






#
