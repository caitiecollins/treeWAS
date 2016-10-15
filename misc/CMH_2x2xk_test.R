
######################
## CMH TEST (2x2xk) ##
######################

?mantelhaen.test

## Use the same K as DAPC...


## sample sim data...

# snps <- get(load("E:/treeWAS_Sims/set1/n.subs_1/set1_1_snps.Rdata"))
# phen <- get(load("E:/treeWAS_Sims/set1/n.subs_1/set1_1_phen.Rdata"))
# snps.assoc <- get(load("E:/treeWAS_Sims/set1/n.subs_1/set1_1_performance.Rdata"))
# snps.assoc <- snps.assoc$snps.assoc
# tree <- get(load("E:/treeWAS_Sims/set1/n.subs_1/set1_1_tree.Rdata"))

snps <- get(load("/media/caitiecollins/Seagate Backup Plus Drive/treeWAS_Sims/set1/n.subs_1/set1_1_snps.Rdata"))
phen <- get(load("/media/caitiecollins/Seagate Backup Plus Drive/treeWAS_Sims/set1/n.subs_1/set1_1_phen.Rdata"))
snps.assoc <- get(load("/media/caitiecollins/Seagate Backup Plus Drive/treeWAS_Sims/set1/n.subs_1/set1_1_performance.Rdata"))
snps.assoc <- snps.assoc$snps.assoc
tree <- get(load("/media/caitiecollins/Seagate Backup Plus Drive/treeWAS_Sims/set1/n.subs_1/set1_1_tree.Rdata"))


snps.ori <- snps.ori.ori
phen.ori <- phen.ori.ori

## Get colnames(snps)
if(is.null(colnames(snps))) colnames(snps) <- c(1:ncol(snps))
snps.names <- colnames(snps)

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


## Get pop factor (same as in PCA/DAPC):

## Identify main pop clusters:
n.PCs <- 5
grp <- find.clusters(snps, n.pca=n.PCs, choose.n.clust=FALSE, max.n.clust=(n.PCs + 1)) # pca.select="percVar", perc.pca=60,
pop <- grp$grp # gives same result as cutree(clust, k=6)
n.grp <- length(levels(pop))

##############################################

snps.12 <- snps+1
phen.34 <- phen+3
mat <- t(matrix(as.numeric(paste(snps.12, ".", phen.34, sep="")), nrow=ncol(snps), byrow=T))

## get only unique columns of pasted mat:
mat.u <- get.unique.matrix(mat)
mat.unique <- mat.u$unique.data
index <- mat.u$index

## get all 2x2 combos of snps.12 and phen.34:
noms <- c("1.3", "1.4", "2.3", "2.4")

## get array from table, by pop:
arr.l <- list()
for(i in 1:ncol(mat.unique)){
  tab <- list()
  for(e in 1:length(levels(pop))){
    temp <- ftable(mat.unique[pop==e, i])
    tab[[e]] <- replace(rep(0, 4), which(noms %in% attr(temp, "col.vars")[[1]]), temp)
  } # end (e) loop
  arr.l[[i]] <- do.call(cbind, tab)
} # end for (i) loop
arr <- do.call(rbind, arr.l)


arr.complete <- arr
##############
## FOR LOOP ##
##############
## TO GET P-VALUES FROM CMH TEST for EACH SNPs COLUMN:
p.vals <- list()
for(i in 1:ncol(mat.unique)){
  ## get indices for this snp for all pops and all 4 2x2 combos:
  from <- seq(1, nrow(arr.complete), 4)[i]
  to <- from+3
  arr <- arr.complete[from:to,]
  dat <- array(arr,
               dim = c(2,2,ncol(arr)),
               dimnames = list(
                 phen = c("0", "1"),
                 SNP = c("0", "1"),
                 pop = levels(pop)
               ))
  ## Run CMH test on this unique snps column:
  CMH <- mantelhaen.test(dat)
  p.vals[[i]] <- CMH$p.value
} # end for loop
p.vals <- as.vector(unlist(p.vals))

## get full set of p-vals for non-unique columns:
p.vals.ori <- p.vals
p.vals <- p.vals[index]

## get Bonferroni threshold:
thresh <- 0.01/ncol(snps)
## identify sig snps:
sig.snps <- snps.names[which(p.vals < thresh)]

## STORE AND SAVE p.vals & sig.snps...



#




#
