

#############################
## treeWAS example dataset ##
#############################


foo <- coalescent.sim(n.ind = 100,
                      n.snps = 10000,
                      n.subs = dist_0  ,
                      n.snps.assoc = 10,
                      assoc.prob = 90,
                      n.phen.subs = 15,
                      phen = NULL,
                      plot = TRUE,
                      heatmap = FALSE,
                      reconstruct = FALSE,
                      dist.dna.model = "JC69",
                      grp.min = 0.4,
                      row.names = TRUE,
                      set = 1,
                      coaltree = TRUE,
                      s = 20,
                      af = 10,
                      filename.plot = NULL,
                      seed = 2)

FOO <- foo
foo <- FOO

snps.ori <- foo$snps
snps.assoc <- foo$snps.assoc
tree <- foo$tree
phen <- foo$phen
phen.plot.col <- foo$phen.plot.col

snps <- snps.ori


save(snps.ori, file="C:/Users/Caitlin/treeWAS/misc/snps.ori.eg.Rdata")
save(snps.assoc, file="C:/Users/Caitlin/treeWAS/misc/snps.assoc.eg.Rdata")
save(tree, file="C:/Users/Caitlin/treeWAS/misc/tree.eg.Rdata")
save(phen, file="C:/Users/Caitlin/treeWAS/misc/phen.eg.Rdata")
save(phen.plot.col, file="C:/Users/Caitlin/treeWAS/misc/phen.plot.col.eg.Rdata")

## Make snps more interesting... ###########################################################################
## (1) make a DNAbin(-type?) object (w nts not binary, and 2-columns for each locus)
## (2) add some NAs
## (3) add some tri-allelic and 4-allelic? loci
## (4) make row.names and tree labs more interesting
## (5) put snps rows out of order (NB - leave these as is, but use for label-checking eg)
############################################################################################################

## save snsp.assoc as NAMES:
snps.assoc <- colnames(snps)[snps.assoc]
save(snps.assoc, file="C:/Users/Caitlin/treeWAS/misc/snps.assoc.eg.Rdata")


############################################################################################################
## (1) make a DNAbin(-type?) object (w nts not binary, and 2-columns for each locus)

str(snps)

SNPS <- snps

snps.assoc
#  618 2017 2656 3721 5728 6287 6604 8980 9080 9442
#  "618"  "2017" "2656" "3721" "5728" "6287" "6604" "8980" "9080" "9442"

foo <- matrix(as.character(snps), nrow=nrow(snps), ncol=ncol(snps))
rownames(foo) <- rownames(snps)
colnames(foo) <- colnames(snps)

str(foo)

# select cols to change to nts:
loci <- sample(1:ncol(snps), replace=FALSE)

## ac
toReplace <- loci[1:1000]
foo[,toReplace] <- replace(foo[,toReplace], which(foo[,toReplace]=="0"), "a")
foo[,toReplace] <- replace(foo[,toReplace], which(foo[,toReplace]=="1"), "c")

## ag
toReplace <- loci[1001:3000]
foo[,toReplace] <- replace(foo[,toReplace], which(foo[,toReplace]=="0"), "a")
foo[,toReplace] <- replace(foo[,toReplace], which(foo[,toReplace]=="1"), "g")

## tc
toReplace <- loci[3001:4000]
foo[,toReplace] <- replace(foo[,toReplace], which(foo[,toReplace]=="0"), "t")
foo[,toReplace] <- replace(foo[,toReplace], which(foo[,toReplace]=="1"), "c")

## gt
toReplace <- loci[4001:7000]
foo[,toReplace] <- replace(foo[,toReplace], which(foo[,toReplace]=="0"), "g")
foo[,toReplace] <- replace(foo[,toReplace], which(foo[,toReplace]=="1"), "t")

## cg
toReplace <- loci[7001:9000]
foo[,toReplace] <- replace(foo[,toReplace], which(foo[,toReplace]=="0"), "c")
foo[,toReplace] <- replace(foo[,toReplace], which(foo[,toReplace]=="1"), "g")

## ga
toReplace <- loci[9001:10000]
foo[,toReplace] <- replace(foo[,toReplace], which(foo[,toReplace]=="0"), "g")
foo[,toReplace] <- replace(foo[,toReplace], which(foo[,toReplace]=="1"), "a")

unique(as.vector(unlist(foo)))
# "g" "c" "t" "a"

snps <- foo

bin <- as.DNAbin(snps)
gen <- DNAbin2genind(bin, polyThres=0)
gen.s <- gen@tab

str(gen.s)


snps.eg <- gen.s
save(snps.eg, file="C:/Users/Caitlin/treeWAS/misc/snps.eg.Rdata")

snps <- snps.eg

############################################################################################################
## (2) add some NAs

snps <- get(load("C:/Users/Caitlin/treeWAS/misc/snps.eg.Rdata"))


is.even <- function(x) x %% 2 == 0
is.odd <- function(x) x %% 2 != 0


toReplace <- sample(c(1:ncol(snps))[-which(removeLastN(colnames(snps), 2) %in% snps.assoc)],
                    3000, replace=F)
toReplace <- as.list(toReplace)

for(i in 1:length(toReplace)){
  x <- toReplace[[i]]
  if(is.odd(x)){
    x <- c(x, x+1)
  }else{
    x <- c(x-1, x)
  }
  toNA <- sample(c(1:nrow(snps)), sample(c(1:nrow(snps)), 1), replace=F)
  snps[, x[1]] <- replace(snps[, x[1]], toNA, NA)
  snps[, x[2]] <- replace(snps[, x[2]], toNA, NA)

} # end for loop



SNPS <- SNPS.ORI <- snps
############################################################################################################
## (3) add some tri-allelic and 4-allelic? loci

# toReplace <- sample(c(1:ncol(snps))[-snps.assoc], 200, replace=F)
# toReplace <- as.list(toReplace)

toReplace <- as.list(c(17, 1417, 5099)) ## 9 708 2549

# for(i in 1:length(toReplace)){

SNPS.NEW <- snps

for(i in 1:length(toReplace)){

  x <- toReplace[[i]]
  a <- removeLastN(colnames(snps)[x-1], 2)
  b <- removeLastN(colnames(snps)[x], 2)
  c <- removeLastN(colnames(snps)[x+1], 2)
  names(a) <- names(b) <- names(c) <- NULL
  if(identical(b,c)){
    x <- c(x, x+1)
  }else{
    x <- c(x-1, x)
  }
  to3 <- sample(c(1:nrow(snps)), sample(c(1:(nrow(snps)*0.25)), 1), replace=F)
  snps[, x[1]] <- replace(snps[, x[1]], to3, 0)
  snps[, x[2]] <- replace(snps[, x[2]], to3, 0)
  SNPS <- snps
  s <- replace(rep(0, nrow(snps)), to3, 1)
  temp <- cbind(snps[,(1:(x[2]))], s)
  toto <- cbind(temp, SNPS[, (x[2]+1):ncol(SNPS)])
  snps <- toto
  nom <- removeLastN(colnames(snps)[x], 2)
  suff <- keepLastN(colnames(snps)[x], 1)
  nts <- c("a", "c", "g", "t")
  nt <- sample(nts[-which(nts %in% suff)], 1)
  colnames(snps)[(x[2]+1)] <- paste(nom[1], nt, sep=".")

} # end for loop

save(snps, file="C:/Users/Caitlin/treeWAS/misc/snps.eg.Rdata")



############################################################################################################
## (4) make row.names and tree labs more interesting

head(row.names(snps))
head(tree$tip.label)

identical(as.character(row.names(snps)), as.character(tree$tip.label)) # TRUE


n1 <- as.character(phen)
n2 <- sample(c(1051:1150), nrow(snps), replace=F)

nom <- paste(n1, n2, sep="_")

row.names(snps) <- nom
names(phen) <- nom
tree$tip.label <- nom

save(snps, file="C:/Users/Caitlin/treeWAS/misc/snps.eg.Rdata")

############################################################################################################
## (5) put snps rows out of order (NB - leave these as is, but use for label-checking eg)

ord <- sample(c(1:nrow(snps)), nrow(snps), replace=F)

snps <- snps[ord,]
str(snps)

save(snps, file="C:/Users/Caitlin/treeWAS/misc/snps.eg.Rdata")

############################################################################################################

# snps[, which(removeLastN(colnames(snps), 2) %in% snps.assoc)]

## save as data objects:
library(devtools)

use_data(snps, overwrite=T)
use_data(phen, overwrite=T)
use_data(tree, overwrite=T)
use_data(snps.assoc, overwrite=T)
use_data(phen.plot.col, overwrite=T)

dat <- list("tree" = tree,
            "snps" = snps,
            "snps.assoc" = snps.assoc,
            "phen" = phen,
            "phen.plot.col" = phen.plot.col)

use_data(dat, overwrite=T)

## add snps.assoc as attribute?
# toto <- snps
# attr(toto, which="snps.assoc") <- snps.assoc


# snps[, which(removeLastN(colnames(snps), 2) %in% snps.assoc)]

############################################################################################################


SNPS.OUT <- snps

snps <- get.binary.snps(snps)

out <- treeWAS(snps = snps,
               phen = phen,
               tree = tree,
               n.subs = NULL,
               sim.n.snps = ncol(snps)*10,
               test = c("terminal", "simultaneous", "subsequent"),
               snps.reconstruction = "parsimony",
               snps.sim.reconstruction = "parsimony",
               phen.reconstruction = "parsimony",
               na.rm = TRUE,
               p.value = 0.01,
               p.value.correct = "bonf",
               p.value.by = "count",
               dist.dna.model = "JC69",
               plot.tree = FALSE,
               plot.manhattan = TRUE,
               plot.null.dist = TRUE,
               plot.dist = FALSE,
               snps.assoc = NULL,
               filename.plot = NULL,
               seed = 1)

treeWAS.out <- out
save(treeWAS.out, file="C:/Users/Caitlin/treeWAS/misc/treeWAS.out.eg.Rdata")
use_data(treeWAS.out, overwrite=T)

############################################################################################################
#



#




#


#


