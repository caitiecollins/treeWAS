

######################
## Troubleshooting! ##
######################

## sample sim data for troubleshooting...
working.dir <- "E:/treeWAS_Sims/set"
set.number <- 3

# args <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_args.Rdata", sep="")))
# args
#
# snps <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_snps.Rdata", sep="")))
# phen <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_phen.Rdata", sep="")))
# snps.assoc <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_performance.Rdata", sep="")))
# snps.assoc <- snps.assoc$snps.assoc
# tree <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_tree.Rdata", sep="")))
#
# score3 <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_score3.Rdata", sep="")))
########

args <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_args.Rdata", sep="")))
args

snps <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_snps.Rdata", sep="")))
phen <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_phen.Rdata", sep="")))
snps.assoc <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_performance.Rdata", sep="")))
snps.assoc <- snps.assoc$snps.assoc
tree <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_tree.Rdata", sep="")))

score3 <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_score3.Rdata", sep="")))



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
# res <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_res.Rdata", sep="")))
# # names(res)
# snps.rec <- res$dat$snps.rec
# phen.rec <- res$dat$phen.rec
# snps.sim.rec <- res$dat$snps.sim.rec
# tree <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_tree.Rdata", sep="")))
# snps.assoc <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_performance.Rdata", sep="")))
# snps.assoc <- snps.assoc$snps.assoc
# # score3 <- get(load(paste("E:/treeWAS_Sims/set3/n.subs_1/set3_", n, "_score3.Rdata", sep="")))

## "C:/Cait 2016/Work/Xavier/Sims/set3/set3_"
n <- 32
res <- get(load(paste("C:/Cait 2016/Work/Xavier/Sims/set3/set3_", n, "_res.Rdata", sep="")))
# names(res)
snps.rec <- res$dat$snps.rec
phen.rec <- res$dat$phen.rec
snps.sim.rec <- res$dat$snps.sim.rec
tree <- get(load(paste("C:/Cait 2016/Work/Xavier/Sims/set3/set3_", n, "_tree.Rdata", sep="")))
snps.assoc <- get(load(paste("C:/Cait 2016/Work/Xavier/Sims/set3/set3_", n, "_performance.Rdata", sep="")))
snps.assoc <- snps.assoc$snps.assoc
# score3 <- get(load(paste("C:/Cait 2016/Work/Xavier/Sims/set3/set3_", n, "_score3.Rdata", sep="")))
args <- get(load(paste("C:/Cait 2016/Work/Xavier/Sims/set3/set3_", n, "_args.Rdata", sep="")))
args

##########
n <- 38
res <- get(load(paste("C:/Cait 2016/Work/Xavier/Sims/set3/temp/old_score3/set3_", n, "_res.Rdata", sep="")))
# names(res)
snps.rec <- res$dat$snps.rec
phen.rec <- res$dat$phen.rec
snps.sim.rec <- res$dat$snps.sim.rec
tree <- get(load(paste("C:/Cait 2016/Work/Xavier/Sims/set3/temp/old_score3/set3_", n, "_tree.Rdata", sep="")))
snps.assoc <- get(load(paste("C:/Cait 2016/Work/Xavier/Sims/set3/temp/old_score3/set3_", n, "_performance.Rdata", sep="")))
snps.assoc <- snps.assoc$snps.assoc
# score3 <- get(load(paste("C:/Cait 2016/Work/Xavier/Sims/set3/temp/old_score3/set3_", n, "_score3.Rdata", sep="")))
args <- get(load(paste("C:/Cait 2016/Work/Xavier/Sims/set3/temp/old_score3/set3_", n, "_args.Rdata", sep="")))
args
########

snps.rec.ori <- snps.rec
snps.sim.rec.ori <- snps.sim.rec
phen.rec.ori <- phen.rec
snps.assoc.ori <- snps.assoc
snps.assoc

# identical(snps.rec.ori, snps.rec)
# identical(snps.sim.rec.ori, snps.sim.rec)
# identical(phen.rec.ori, phen.rec)
# identical(tree.ori, tree)
# identical(args, args.41) # FALSE -- 41 was simulated w only 10,000 SIM.N.SNPS!
# args.41 <- args
# n <- 40


# snps.rec <- snps.rec.ori

#######################
## Calculate SCORE 3 ##    ############################################################################
#######################
edges <- tree$edge

Pa <- phen.rec[edges[,1]]
Pd <- phen.rec[edges[,2]]
Sa <- snps.rec[edges[,1], ]
Sd <- snps.rec[edges[,2], ]
l <- tree$edge.length

## Original integral-based score:
score3.integral <- l*(((4/3)*Pa*Sa) +
                        ((2/3)*Pa*Sd) +
                        ((2/3)*Pd*Sa) +
                        ((4/3)*Pd*Sd) -
                        Pa -
                        Pd -
                        Sa -
                        Sd +
                        1)

## integral-based score WITHOUT accounting for edge.length??
score3.integral.NoEdgeL <- (((4/3)*Pa*Sa) +
                              ((2/3)*Pa*Sd) +
                              ((2/3)*Pd*Sa) +
                              ((4/3)*Pd*Sd) -
                              Pa -
                              Pd -
                              Sa -
                              Sd +
                              1)

SCORE3.integral <- abs(colSums(score3.integral, na.rm=TRUE))
SCORE3.integral.NoEdgeL <- abs(colSums(score3.integral.NoEdgeL, na.rm=TRUE))

######## end ###########################################################################################

int <- SCORE3.integral
int.NoL <- SCORE3.integral.NoEdgeL

corr.dat.L <- int
corr.dat.NoL <- int.NoL
###

#########################################    ###################################    ##############################################
snps.rec <- snps.sim.rec

#######################
## Calculate SCORE 3 ##    ############################################################################
#######################
edges <- tree$edge

Pa <- phen.rec[edges[,1]]
Pd <- phen.rec[edges[,2]]
Sa <- snps.rec[edges[,1], ]
Sd <- snps.rec[edges[,2], ]
l <- tree$edge.length

## Original integral-based score:
score3.integral <- l*(((4/3)*Pa*Sa) +
                        ((2/3)*Pa*Sd) +
                        ((2/3)*Pd*Sa) +
                        ((4/3)*Pd*Sd) -
                        Pa -
                        Pd -
                        Sa -
                        Sd +
                        1)

## integral-based score WITHOUT accounting for edge.length??
score3.integral.NoEdgeL <- (((4/3)*Pa*Sa) +
                              ((2/3)*Pa*Sd) +
                              ((2/3)*Pd*Sa) +
                              ((4/3)*Pd*Sd) -
                              Pa -
                              Pd -
                              Sa -
                              Sd +
                              1)

SCORE3.integral <- abs(colSums(score3.integral, na.rm=TRUE))
SCORE3.integral.NoEdgeL <- abs(colSums(score3.integral.NoEdgeL, na.rm=TRUE))

######## end ###########################################################################################

int <- SCORE3.integral
int.NoL <- SCORE3.integral.NoEdgeL

corr.sim.L <- int
corr.sim.NoL <- int.NoL

#########################################    ###################################    ##############################################

# hist(corr.dat.L)
# hist(corr.sim.L)
hist(corr.dat.NoL)
hist(corr.sim.NoL)

snps.assoc
# int[snps.assoc]
# int.NoL[snps.assoc]

#################
### snps only ###
#################
# ## WITH edge-length
# sapply(c(1:length(snps.assoc)), function(e) length(which(corr.dat.L > corr.dat.L[snps.assoc[e]])))
# ## (n.subs=1) ##
# # [1]  2  0  0  0  0  1 19  0  0  3
# # [1] 2540 3114  172    4  516  180 5241    2    3  424
# # [1] 2906 3246  172    4 2766 2228 5246    2    3 2746
#
# ## set3_32
# # [1] 5670 1816    3    1  265   30 5318    2   33    0 # ini
#
# ## WITHOUT edge-length
# sapply(c(1:length(snps.assoc)), function(e) length(which(corr.sim.NoL > corr.dat.NoL[snps.assoc[e]])))
# ## (n.subs=1) ##
# ## set3_31
# # [1]    6    0    0    0    3    0 2128    0    0    2
# # [1]  6  2  0  4  3  5 38  7  1  8
# # [1]  6  2  0  4  3  5 39  6  1  8
#
# ## set3_32
# # [1] 12 11  3  1  5  7  0  4  6  2 # ini

##################
### w snps.sim ###
##################
## WITH edge-length
sapply(c(1:length(snps.assoc)), function(e) length(which(corr.sim.L > corr.dat.L[snps.assoc[e]])))
## (n.subs=1) ##
## set3_32
# [1] 57609 18742     5     0  2503   243 54571     2   286     0

## (n.subs = dist_R0.05)
## set3_38
# [1]     5  8956 95356  5668  8971 14745 12445  7275 11325  9806

## set3_39
# [1]  4204 21208    55  1152 16303  1746  5492     2  3704  3607

## set3_40
# [1]  1859  1104 98275     0  3560  1767 27611  4318 23684 99871

## (n.subs = dist_0.05) NEW:
# set3_32
# [1]    44  9106 11953 37336  6469 64655  7005 12719 37378 37426

## WITHOUT edge-length
sapply(c(1:length(snps.assoc)), function(e) length(which(corr.sim.NoL > corr.dat.NoL[snps.assoc[e]])))
## (n.subs=1) ##
## set3_32
# [1] 32 29  0  0  0  4  0  0  1  0
## set3_41
#  [1]  291   35   15    8  175   11 5359 7686 5498  253

## (n.subs = dist_R0.05)
## set3_38
# [1]    1  442 7118  118 1251  941  600  346  485  464

## set3_39
# [1]    1  342    0    1 1336   19   22    0    8  522

## set3_40
# [1] 0 0 0 0 0 0 0 0 0 3
## set3_41
# [1]  253   30    1   13  217   35  570 4172  524 1101

## (n.subs = dist_0.05) NEW:
# set3_32
# [1] 44  698  675   17  134 1557 1293 5867  291  517

## Properly:
corr.dat <- corr.dat.NoL
corr.sim <- corr.sim.NoL
p.value <- 0.01/length(corr.sim)
thresh <- quantile(corr.sim, probs=1-p.value)
hist(corr.sim)
corr.dat[snps.assoc]

## Identify (real) SNPs w correlations > thresh:
sig.snps <- which(corr.dat > thresh)

print("sig.snps"); print(sig.snps)
print("sig.snps.assoc"); print(snps.assoc[which(snps.assoc %in% sig.snps)])
print("N TRUE pos"); print(length(which(snps.assoc %in% sig.snps)))
# 6
# 9 # 3_40 OLD n.subs=dist_0.05
# 0 # 3_32 NEW
## set3_39
# 2
## set3_38
# 0

TP <- which(sig.snps %in% snps.assoc)
FP <- length(sig.snps)
if(!.is.integer0(TP)) FP <- sig.snps[-TP]
print("N FALSE pos"); print(length(FP))
# 0 #
# 0 # 3_40 OLD n.subs=dist_0.05
# 1 # 3_32 NEW
## set3_39
# 0
## set3_38
# 0

corr.dat.NoL[snps.assoc]
head(round(sort(corr.dat.NoL, decreasing=T), 2), 20)

#


snps.ori <- snps
phen.ori <- phen

identical(snps.ori, snps)  ## NO! (new) snps.ori has the REDUCED n.snps columns! ## ONLY for set 3 snps!
identical(phen.ori, phen) ## NO! (new) phen.ori has phen as character...
