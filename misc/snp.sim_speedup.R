


####################
## SPEEDING UP... ##
####################


###############
## snp.sim.R ##
###############

setwd("/home/caitiecollins/ClonalFrameML/src/pubMLST/Neisseria/B/UK/WG/phylip/")

snps <- get(load("./B_UK_WG.fas.out.snps_clean.Rdata"))
phen <- get(load("./B_UK_WG.fas.out.phen_clean.Rdata"))
tree <- get(load("./B_UK_WG.fas.out.tree_clean.Rdata"))
dat <- get(load("./B_UK_WG.fas.out.read.CFML_dat.Rdata"))

n.subs <- dat$n.subs


set.args(list(n.snps = 10000,
              n.subs = n.subs,
              snp.root = NULL,
              n.snps.assoc = 0,
              assoc.prob = 100,
              tree = tree,
              phen.loci = NULL,
              heatmap = FALSE,
              reconstruct = FALSE,
              dist.dna.model = "JC69",
              row.names = row.names(snps),
              set = NULL,
              seed = 1))

# rm(snps)
# rm(phen)
# rm(dat)

object.size(tree) # 24200 bytes





## snp.sim ##
ss <- snp.sim(n.snps = 10000,
              n.subs = n.subs,
              snp.root = NULL,
              n.snps.assoc = 0,
              assoc.prob = 100,
              tree = tree,
              phen.loci = NULL,
              heatmap = FALSE,
              reconstruct = FALSE,
              dist.dna.model = "JC69",
              row.names = row.names(snps),
              set = NULL,
              seed = 1)

ss <- ss$snps

nm <- get.fitch.n.mts(ss, tree)
table(nm)

str(n.subs)
barplot(n.subs)
barplot(table(nm))



## snp.sim (OLD) ##
# ss2 <- snp.sim(n.snps = 10000,
#               n.subs = n.subs,
#               snp.root = NULL,
#               n.snps.assoc = 0,
#               assoc.prob = 100,
#               tree = tree,
#               phen.loci = NULL,
#               heatmap = FALSE,
#               reconstruct = FALSE,
#               dist.dna.model = "JC69",
#               row.names = row.names(snps),
#               set = NULL,
#               seed = 1)
#
# ss2 <- ss2$snps
#
# nm2 <- get.fitch.n.mts(ss2, tree)
# table(nm2)
#
# str(n.subs)
# barplot(n.subs)
# barplot(table(nm2))

## All fine until WHILE LOOP replacing non-polymorphic loci (l 463)



#########################
## STORAGE EFFICIENCY? ##
#########################

## CHARACTER: ##
foo <- sample(c("a", "t"), 10000, replace=T)

l.c <- as.list(rep(list(foo), 100))
str(l.c)

m.c <- matrix(rep(foo, 100), nrow=100)
str(m.c)


## NUMERIC: ##
foo <- sample(c(0, 1), 10000, replace=T)

l.n <- as.list(rep(list(foo), 100))
str(l.n)

m.n <- matrix(rep(foo, 100), nrow=100)
str(m.n)


## LOGICAL?: ##
foo <- sample(c(FALSE, TRUE), 10000, replace=T)

l.l <- as.list(rep(list(foo), 100))
str(l.l)

m.l <- matrix(rep(foo, 100), nrow=100)
str(m.l)


#





#






#





#


#
