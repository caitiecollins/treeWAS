

######################
## Troubleshooting! ##
######################

## sample sim data for troubleshooting...
working.dir <- "E:/treeWAS_Sims/set"
set.number <- 3

snps <- get(load("E:/treeWAS_Sims/set3/n.subs_1/set3_31_snps.Rdata"))
phen <- get(load("E:/treeWAS_Sims/set3/n.subs_1/set3_31_phen.Rdata"))
snps.assoc <- get(load("E:/treeWAS_Sims/set3/n.subs_1/set3_31_performance.Rdata"))
snps.assoc <- snps.assoc$snps.assoc
tree <- get(load("E:/treeWAS_Sims/set3/n.subs_1/set3_31_tree.Rdata"))

score3 <- get(load("E:/treeWAS_Sims/set3/n.subs_1/set3_31_score3.Rdata"))
# SCORE3.ori <- SCORE3 <- score3


snps.ori <- snps
phen.ori <- phen
names(snps.assoc) <- as.character(snps.assoc)
snps.assoc.ori <- snps.assoc


corr.dat <- SCORE3$corr.dat
corr.sim <- SCORE3$corr.sim
corr.dat.ori <- corr.dat
corr.sim.ori <- corr.sim

## thresh method = pval 0.0001, bonf, count, 1x.n.snps
snps.names <- colnames(snps.ori)


# length(corr.dat)
int.dat <- corr.dat[[9]]


## Q without simultaneous subs??
Q <- corr.dat[[14]]

score3.p <- corr.dat[[7]]
score3.n <- corr.dat[[8]]

####################################


## FROM SCRATCH: ####
res <- get(load("E:/treeWAS_Sims/set3/n.subs_1/set3_31_res.Rdata"))
names(res)
snps.rec <- res$dat$snps.rec
phen.rec <- res$dat$phen.rec
snps.sim.rec <- res$dat$snps.sim.rec

##############
## FOR LOOP ##
##############
edges <- tree$edge
SCORE3.integral <- SCORE3.integral.NoEdgeL <- list()
for(i in 1:ncol(snps.rec)){
  ## HANDLE SNPS.REC:
  snp.rec <- snps.rec[,i]
  # snp.rec <- replace(snp.rec, which(! snp.rec %in% c(0, 1)), NA)

  score3.integral <- score3.integral.NoEdgeL <- list()
  for(e in 1:nrow(edges)){

    Pa <- phen.rec[edges[e,1]]
    Pd <- phen.rec[edges[e,2]]
    Sa <- snp.rec[edges[e,1]]
    Sd <- snp.rec[edges[e,2]]
    l <- tree$edge.length[e]

    ## Original integral-based score:
    score3.integral[[e]] <- l*(((4/3)*Pa*Sa) +
                                 ((2/3)*Pa*Sd) +
                                 ((2/3)*Pd*Sa) +
                                 ((4/3)*Pd*Sd) -
                                 Pa -
                                 Pd -
                                 Sa -
                                 Sd +
                                 1)

    ## integral-based score WITHOUT accounting for edge.length??
    score3.integral.NoEdgeL[[e]] <- (((4/3)*Pa*Sa) +
                                       ((2/3)*Pa*Sd) +
                                       ((2/3)*Pd*Sa) +
                                       ((4/3)*Pd*Sd) -
                                       Pa -
                                       Pd -
                                       Sa -
                                       Sd +
                                       1)

  } # end (e) for loop

  ## integral-based score
  score3.integral <- as.vector(unlist(score3.integral))
  ## save BOTH the sum and the max-only score:
  ## sum:
  SCORE3.integral[[i]] <- abs(sum(score3.integral[!is.na(score3.integral)]))

  score3.integral.NoEdgeL <- as.vector(unlist(score3.integral.NoEdgeL))
  SCORE3.integral.NoEdgeL[[i]] <- abs(sum(score3.integral.NoEdgeL[!is.na(score3.integral.NoEdgeL)]))
} # end (i) for loop


int <- as.vector(unlist(SCORE3.integral))
int.NoL <- as.vector(unlist(SCORE3.integral.NoEdgeL))

## With NAs:
# int.ori <- int
# int.NoL.ori <- int.NoL

## WithOUT NAs in sum:
# int.noNAs <- int
# int.NoL.noNAs <- int.NoL

## WithOUT converting non-binary to NA:
int.ini <- int
int.NoL.ini <- int.NoL

hist(int)
hist(int.NoL)

snps.assoc

# int[snps.assoc]
sapply(c(1:length(snps.assoc)), function(e) length(which(int > int[snps.assoc[e]])))
# [1]  2  0  0  0  0  1 19  0  0  3
# [1] 2540 3114  172    4  516  180 5241    2    3  424
# [1] 2906 3246  172    4 2766 2228 5246    2    3 2746

sapply(c(1:length(snps.assoc)), function(e) length(which(int.NoL > int.NoL[snps.assoc[e]])))
# [1]    6    0    0    0    3    0 2128    0    0    2
# [1]  6  2  0  4  3  5 38  7  1  8
# [1]  6  2  0  4  3  5 39  6  1  8

i <- snps.assoc[2]

#

#


