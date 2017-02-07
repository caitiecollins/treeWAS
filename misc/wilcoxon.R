
###################
## WILCOXON TEST ##
###################

## (1) PARSIMONY vs. ACE/ML
## (2) SCORE 3 -- No Edge Length vs. with Edge Length


###############################
## (1) PARSIMONY vs. ACE/ML: ##
###############################

#######################################################################################################################

##############
## GET DATA ##    #######    #######    #######    #######    #######    #######    #######    #######    #######    #######
##############

# setwd("C:/Cait 2016/Work/Xavier/Sims/set1/")

set.n <- "set2"

dirname <- paste("E:/Cait 2016/Work/Xavier/Sims/", set.n, sep="") #  "/coaltree",
# dirname <- paste("C:/Cait 2016/Work/Xavier/Sims/", set.n, sep="") #  "/coaltree",
setwd(dirname)

#################################################################################################


evalStats.pars <- get(load("./coaltree/parsimony/set2_coal_1_80_treeWAS_ALL_all_tests_best_thresh_evalStats.Rdata"))

evalStats.ACE <- get(load("./coaltree/ACE/set2_ACE_coal_1_54_56_60_treeWAS_ALL_all_tests_best_thresh_evalStats.Rdata"))

## Needed to keep only subset of df.pars bc not done in df.ACE:
# df <- evalStats.pars
# toKeep <- c(1:54, 56:60)
# # toKeep <- rep(toKeep, 10)
# TOKEEP <- list()
# for(i in 1:10){
#   TOKEEP[[i]] <- toKeep+(80*(i-1))
# }
# TOKEEP <- as.vector(unlist(TOKEEP))
#
# df <- df[TOKEEP, ]
# nrow(df)
# nrow(evalStats.ACE)
# evalStats.pars <- df
#
# save(evalStats.pars, file="./coaltree/parsimony/set2_coal_1_54_56_60_treeWAS_ALL_all_tests_best_thresh_evalStats.Rdata")

evalStats.pars.ori.ori <- evalStats.pars.ori <- evalStats.pars
evalStats.ACE.ori.ori <- evalStats.ACE.ori <- evalStats.ACE

df.pars <- evalStats.pars
df.ACE <- evalStats.ACE

# ?wilcox.test

head(df.pars)
head(df.ACE)

treeWAS.tests <- c("terminal", "simultaneous", "subsequent", "treeWAS.all")

df <- df.pars
df <- df[df$assoc.test %in% treeWAS.tests, ]
df.pars <- df

df <- df.ACE
df <- df[df$assoc.test %in% treeWAS.tests, ]
df.ACE <- df

nrow(df.pars) # 236 # 320
nrow(df.ACE) # 236 # 320

df <- df.pars
table(df$assoc.test)

###################
## WILCOXON TEST ##
###################
## Run test for all 3 treeWAS tests + treeWAS.all
## And for each of the 4 main performance statistics
treeWAS.tests <- c("terminal", "simultaneous", "subsequent", "treeWAS.all")
perfStats <- c("F1.score", "PPV", "sensitivity", "FPR")

#################
## STORE in DF ##
#################
treeWAS.test <- as.vector(unlist(sapply(c(1:length(treeWAS.tests)),
                                        function(e)
                                          rep(treeWAS.tests[e], 4))))

stat <- rep(perfStats, 4)
p.value <- median.diff <- CI.L <- CI.U <- rep(NA, length(stat))

wilcox <- data.frame(treeWAS.test, stat, p.value, median.diff, CI.L, CI.U)

###############
## Get list: ##
###############
out <- list()
counter <- 0
for(i in 1:length(treeWAS.tests)){
  df <- df.pars
  x.tab <- df[df$assoc.test == treeWAS.tests[i], ]
  df <- df.ACE
  y.tab <- df[df$assoc.test == treeWAS.tests[i], ]

  out[[i]] <- list()
  for(j in 1:length(perfStats)){
    x <- x.tab[, perfStats[j]]
    y <- y.tab[, perfStats[j]]
    toRemove <- which(is.na(x))
    toRemove <- unique(c(toRemove, which(is.na(y))))
    if(length(toRemove) == length(x)){
      out[[i]][[j]] <- NULL
    }else{
      if(length(toRemove) > 0){
        x <- x[-toRemove]
        y <- y[-toRemove]
      }
      ## Check if we can get a conf.int:
      if(suppressWarnings(class(try(wilcox.test(x=x, y=y, alternative="two.sided", paired=TRUE, conf.int=TRUE), silent=TRUE)) == "try-error")){
        out[[i]][[j]] <- suppressWarnings(wilcox.test(x=x, y=y, alternative="two.sided", paired=TRUE, conf.int=FALSE))
      }else{
        out[[i]][[j]] <- suppressWarnings(wilcox.test(x=x, y=y, alternative="two.sided", paired=TRUE, conf.int=TRUE))
      }
    }

    counter <- counter+1

    vals <- out[[i]][[j]]
    p.value <- median.diff <- CI.L <- CI.U <- NA
    if(!is.null(vals$p.value)) p.value <- vals$p.value
    if(!is.null(vals$estimate)) median.diff <- vals$estimate
    if(!is.null(vals$conf.int)){
      CI.L <- vals$conf.int[1]
      CI.U <- vals$conf.int[2]
    }

    wilcox[counter, 3:6] <- c(p.value, median.diff, CI.L, CI.U)

  } # end for (j) loop
  names(out[[i]]) <- perfStats
} # end for (i) loop
names(out) <- treeWAS.tests


## Check out wilcox df:
wilcox
data.frame(wilcox[,1:2], round(wilcox[,3:6], 3))

##################
## SAVE wilcox: ##
##################
save(wilcox, file="./coaltree/set2_parsimony_vs_ACE_wilcox.Rdata")
#



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


#####################################################
## (1) SCORE 3 No Edge Length vs With Edge Length: ##
#####################################################

#######################################################################################################################

##############
## GET DATA ##    #######    #######    #######    #######    #######    #######    #######    #######    #######    #######
##############

# setwd("C:/Cait 2016/Work/Xavier/Sims/set1/")

set.n <- "set3"

# dirname <- paste("C:/Cait 2016/Work/Xavier/Sims/", set.n, sep="") #  "/coaltree",
dirname <- paste("E:/treeWAS_Sims/", set.n, "/coaltree/parsimony/", sep="") #
setwd(dirname)

#################################################################################################

#################
## get score3: ##
#################

foo <- dir(dirname)

## get all performance Rdata names
toKeep <- grep("_score3", foo) ##??
foo <- foo[toKeep]
foo

## load performance data
dat <- list()
setwd(dirname)
system.time(
  for(i in 1:length(foo)){ # length(foo)
    print(i)
    dat[[i]] <- get(load(paste("./", foo[i], sep="")))
    gc()
  }
)

## REORDER dat s.t. order is numeric not character...
# dat <- snps.assoc
inds.new <- c(1:length(dat))
inds.ori <- sort(as.character(inds.new))
ord <- sapply(c(1:length(inds.new)), function(e) which(as.numeric(inds.ori) == inds.new[e]))

dat.ori <- dat
dat <- dat[ord]

score3.ori <- score3 <- dat

save(score3, file="./set3_1_80_ALL_SCORE3.Rdata")

#################################################################################################

######################
## get performance: ##
######################

foo <- dir(dirname)

## get all performance Rdata names
toKeep <- grep("performance", foo) ##??
foo <- foo[toKeep]
foo <- foo[1:80]
foo

## load performance data
dat <- list()
setwd(dirname)
system.time(
  for(i in 1:length(foo)){ # length(foo)
    print(i)
    dat[[i]] <- get(load(paste("./", foo[i], sep="")))
    gc()
  }
)

## REORDER dat s.t. order is numeric not character...
# dat <- snps.assoc
inds.new <- c(1:length(dat))
inds.ori <- sort(as.character(inds.new))
ord <- sapply(c(1:length(inds.new)), function(e) which(as.numeric(inds.ori) == inds.new[e]))

dat.ori <- dat
dat <- dat[ord]

perf.ori <- perf <- dat

#################################################################################################

#####################
## get snps.assoc: ##
#####################
# dat = performance
snps.assoc <- sapply(c(1:length(perf)), function(e) perf[[e]][[1]])

snps.assoc.ori <- snps.assoc <- as.list(as.data.frame(snps.assoc))

## SAVE:
# save(snps.assoc, file="./set3_1.100_ALL_snps.assoc.Rdata")

#################################################################################################

######################################################
## Get performance for WITHOUT and WITH edge.length ##
######################################################

snps.names <- c(1:10000)

performance <- list()
for(i in 1:length(score3)){

  print(i)

  # corr.dat <- score3[[i]]$corr.dat$score3.NoL
  # corr.sim <- score3[[i]]$corr.sim$score3.NoL

  corr.dat <- score3[[i]]$corr.dat$score3.L
  corr.sim <- score3[[i]]$corr.sim$score3.L

  base.p.value <- 0.01

  ## Get threshold:
  p.value <- base.p.value/(length(corr.dat)*3)
  thresh <- quantile(corr.sim, probs=1-p.value)

  ## Identify (real) SNPs w correlations > thresh:
  sig.snps <- which(corr.dat > thresh)
  #   p.vals <- .get.p.vals(corr.sim = corr.sim,
  #                         corr.dat = corr.dat,
  #                         fisher.test = FALSE)
  #   sig.p.vals <- p.vals[sig.snps]



  ## get list of those correlation values
  sig.corrs <- corr.dat[sig.snps]
  ## get the list of those SNPs (ie. their locus names)
  # sig.snps <- dimnames(snps)[[2]][sig.snps]
  sig.snps.names <- snps.names[sig.snps]

  ## re-order list of sig.snps and sig.corrs by value of sig.corr
  NWO <- order(sig.corrs, decreasing=TRUE)
  sig.snps <- sig.snps[NWO]
  sig.corrs <- sig.corrs[NWO]
  # sig.p.vals <- sig.p.vals[NWO]
  sig.snps.names <- sig.snps.names[NWO]



  #####################
  ## get PERFORMANCE ##
  #####################

  test.positive <- sig.snps.names
  # test.positive <- SNP.loci[[i]]$terminal
  # test.positive <- SNP.loci[[i]]$simultaneous
  # test.positive <- SNP.loci[[i]]$subsequent

  if(length(test.positive) > 0){
    test.negative <- snps.names[-test.positive]
  }else{
    test.negative <- snps.names
  }

  n.test.positive <- length(test.positive)
  n.test.negative <- length(test.negative)
  n.tests <- 10000

  snps.associated <- snps.assoc[[i]]

  ## get true positives
  true.positive <- test.positive[which(test.positive %in% snps.associated)]
  TP <- length(true.positive)

  ## get true negatives
  snps.not <- snps.names[-which(snps.names %in% snps.associated)]
  true.negative <- test.negative[which(test.negative %in% snps.not)]
  TN <- length(true.negative)

  ## get false positives
  false.positive <- test.positive[which(test.positive %in% snps.not)]
  FP <- length(false.positive)

  ## get false negatives
  false.negative <- test.negative[which(test.negative %in% snps.associated)]
  FN <- length(false.negative)

  #####################################
  ## CALCULATE METRICS OF EVALUATION ##
  #####################################
  ## accuracy ##
  accuracy <- ((TP + TN) / (TP + TN + FP + FN))

  ## specificity ##
  specificity <- (TN / (TN + FP)) ## = (1 - FPR)

  ## FPR
  FPR <- (FP / (FP + TN)) ## = (1 - specificity)

  ## FNR
  FNR <- (FN / (FN + TP))

  ## sensitivity ##
  sensitivity <- (TP / (TP + FN))

  ## PPV
  PPV <- (TP / (TP + FP)) ## = (1 - FDR)

  ## FDR ##
  FDR <- (FP / (FP + TP)) ## = (1 - PPV)

  ## F1.score ##
  ## Balanced accuracy-like score considering both sensitivity and PPV:
  F1.score <- 2*((sensitivity*PPV) / (sensitivity+PPV))

  ## combine eval metrics into df ##
  performance[[i]] <- data.frame(accuracy, specificity, FPR, FNR, sensitivity, PPV, FDR, F1.score)

} # end for loop


performance <- do.call(rbind, performance)

# performance.NoL <- performance
performance.L <- performance

#


#################################################################################################
#################################################################################################
#################################################################################################


save(performance.NoL, file="./set3_coal_1_80_score3_performance_NoL.Rdata")

save(performance.L, file="./set3_coal_1_80_score3_performance_L.Rdata")

performance.NoL.ori.ori <- performance.NoL.ori <- performance.NoL
performance.L.ori.ori <- performance.L.ori <- performance.L

df.NoL <- performance.NoL
df.L <- performance.L

# ?wilcox.test

nrow(df.NoL) # 80
nrow(df.L) # 80

df <- df.NoL
table(df$assoc.test)

###################
## WILCOXON TEST ##
###################
## Run test for all 3 treeWAS tests + treeWAS.all
## And for each of the 4 main performance statistics
# treeWAS.tests <- c("terminal", "simultaneous", "subsequent", "treeWAS.all")
perfStats <- c("F1.score", "PPV", "sensitivity", "FPR")

#################
## STORE in DF ##
#################
# treeWAS.test <- as.vector(unlist(sapply(c(1:length(treeWAS.tests)),
#                                         function(e)
#                                           rep(treeWAS.tests[e], 4))))

# stat <- rep(perfStats, 4)
stat <- perfStats
p.value <- median.diff <- CI.L <- CI.U <- rep(NA, length(stat))

wilcox <- data.frame(stat, p.value, median.diff, CI.L, CI.U)

###############
## Get list: ##
###############
out <- list()
counter <- 0

x.tab <- df.NoL
y.tab <- df.L

for(j in 1:length(perfStats)){
  x <- x.tab[, perfStats[j]]
  y <- y.tab[, perfStats[j]]
  toRemove <- which(is.na(x))
  toRemove <- unique(c(toRemove, which(is.na(y))))
  if(length(toRemove) == length(x)){
    out[[j]] <- NULL
  }else{
    if(length(toRemove) > 0){
      x <- x[-toRemove]
      y <- y[-toRemove]
    }
    ## Check if we can get a conf.int:
    if(suppressWarnings(class(try(wilcox.test(x=x, y=y, alternative="two.sided", paired=TRUE, conf.int=TRUE), silent=TRUE)) == "try-error")){
      out[[j]] <- suppressWarnings(wilcox.test(x=x, y=y, alternative="two.sided", paired=TRUE, conf.int=FALSE))
    }else{
      out[[j]] <- suppressWarnings(wilcox.test(x=x, y=y, alternative="two.sided", paired=TRUE, conf.int=TRUE))
    }
  }

  counter <- counter+1

  vals <- out[[j]]
  p.value <- median.diff <- CI.L <- CI.U <- NA
  if(!is.null(vals$p.value)) p.value <- vals$p.value
  if(!is.null(vals$estimate)) median.diff <- vals$estimate
  if(!is.null(vals$conf.int)){
    CI.L <- vals$conf.int[1]
    CI.U <- vals$conf.int[2]
  }

  wilcox[counter, 2:5] <- c(p.value, median.diff, CI.L, CI.U)

} # end for (j) loop
names(out) <- perfStats


## Check out wilcox df:
wilcox
data.frame(stat=wilcox[,1], round(wilcox[,2:5], 4))

##################
## SAVE wilcox: ##
##################
save(wilcox, file="./set3_score3_NoL_vs_L_wilcox.Rdata")

#





#





#




#
