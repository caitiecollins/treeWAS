

## see "testing_treeWAS.Rnw" for original boxplots, ggplot2 etc. code...

#########################################
## fn getting order of tips as plotted ##
#########################################
get.tip.order <- function(tree){
  tree2 <- read.tree(text=write.tree(tree))
  out <- as.numeric(tree2$tip.label)
  out <- rev(out)
  return(out)
} # end get.tip.order

##################
## tree (set 1) ##
##################

tree <- get(load("./set1_12_tree.Rdata"))
plot(tree)

phen <- get(load("./set1_12_phen.Rdata"))

ord <- get.tip.order(tree)
#phen.ori <- phen
phen <- phen[ord]
myCol <- as.character(phen)
myCol <- replace(myCol, which(myCol == "B"), "blue")
myCol <- replace(myCol, which(myCol == "A"), "red")
leafCol <- myCol
myCol
edgeCol <- "black"

plot(tree, show.tip=FALSE, edge.width=2, edge.color=edgeCol) # edgeCol
axisPhylo()
tiplabels(text=tree$tip.label, cex=0.6, adj=c(-0.5, 0), col=leafCol, frame="none", font=1)
#

#######################################################################################################################

## Replace boxpots with:
## box < violin < bean plot < beeswarm plot...

####################
## BEESWARM PLOTS ##
####################

install.packages("beeswarm", dep=TRUE)
library(beeswarm)

?beeswarm


###################
## BEESWARM PLOT ##   #########################   #########################   #########################
###################
#
# #################################
# ## EG from beeswarm package... ##
# #################################
# data(breast)
#
# #########################
# ## reg. beeswarm plot: ##
# #########################
# breast2 <- breast[order(breast$event_survival, breast$ER),]
#
# beeswarm(time_survival ~ event_survival, data = breast2, pch = 16,
#          pwcol = as.numeric(ER), xlab = '',
#          ylab = 'Follow-up time (months)',
#          labels = c('Censored', 'Metastasis'))
# legend('topright', legend = levels(breast$ER), title = 'ER',
#        pch = 16, col = 1:2)
#
# ##########################################
# ## Reg. beeswarm PLUS boxplot overlaid: ##
# ##########################################
# beeswarm(time_survival ~ event_survival, data = breast2, pch = 16,
#          pwcol = as.numeric(ER), xlab = '',
#          ylab = 'Follow-up time (months)',
#          labels = c('Censored', 'Metastasis'))
# boxplot(time_survival ~ event_survival, data = breast2, add = T,
#         names = c("",""), col="#0000ff22")
# legend('topright', legend = levels(breast$ER), title = 'ER',
#        pch = 16, col = 1:2)
#
# #################################################
# ## GGPLOT2 beeswarm plot PLUS boxplot overlaid ##
# #################################################
#
# ## The trick is to use the beeswarm call
# ## to get the x and y position.
# ## Beeswarm creates a dataframe
# ## from which we can get the necessary positionings...
#
# beeswarm <- beeswarm(time_survival ~ event_survival,
#                      data = breast, method = 'swarm',
#                      pwcol = ER)[, c(1, 2, 4, 6)]
# colnames(beeswarm) <- c("x", "y", "ER", "event_survival")
#
# library(ggplot2)
# library(plyr)
#
# ## Do not forget to REMOVE the OUTLIERS
# ## from your boxplot or they will superimpose
# ## with the points created by geom_point.
#
# beeswarm.plot <- ggplot(beeswarm, aes(x, y)) +
#   xlab("") +
#   scale_y_continuous(expression("Follow-up time (months)"))
#
# beeswarm.plot2 <- beeswarm.plot +
#   geom_boxplot(aes(x, y, group = round_any(x, 1, round)), outlier.shape = NA)
#
# beeswarm.plot3 <- beeswarm.plot2 + geom_point(aes(colour = ER)) +
#   scale_colour_manual(values = c("black", "red")) +
#   scale_x_continuous(breaks = c(1:2),
#                      labels = c("Censored", "Metastasis"), expand = c(0, 0.5))
#
# plot(beeswarm.plot3)



#


#######################################################################################################################


###############
#### SET X ####
###############

#######################################################################################################################

# ## p.value ##
# if(i %in% 1:8) p.value <- 0.05
# if(i %in% 9:16) p.value <- 0.01
# if(i %in% 17:24) p.value <- 0.001
# if(i %in% 25:32) p.value <- 0.0001
#
# ## p.value.correct ##
# if(i %in% c(1:4, 9:12, 17:20, 25:28)){
#   p.value.correct <- "bonf"
# }else{
#   p.value.correct <- "fdr"
# }
#
# ## p.value.by ##
# if(i %in% c(1,2,9,10,17,18,25,26, ## bonf
#             5,6,13,14,21,22,29,30)){ ## fdr
#   p.value.by <- "count"
# }else{
#   p.value.by <- "density"
# }
#
# ## n.snps.sim ##
# if(i %in% seq(1, 32, 2)){
#   corr.sim <- corr.sim.ori
# }else{
#   corr.sim <- corr.sim.ori[1:length(corr.dat)]
# }

##############
## GET DATA ##    #######    #######    #######    #######    #######    #######    #######    #######    #######    #######
##############

setwd("C:/Cait 2016/Work/Xavier/Sims/set2/")

set.n <- "set2"

dirname <- paste("C:/Cait 2016/Work/Xavier/Sims/", set.n, sep="")
# dirname <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/", set.n, sep="")

foo <- dir(dirname)
foo
## get all performance Rdata names
# toKeep <- grep("performance", foo) ##??
# toKeep <- grep("phen.plot.col", foo) ##??
toKeep <- grep("_res", foo) ##??
# toKeep <- grep("_args", foo) ##??
foo <- foo[toKeep]
foo

## keep only subset w same args:
# toFind <- paste(set.n, "_", c(51:81), "_", sep="")
# toFind <- paste(set.n, "_", c(21:30), "_", sep="")
# toKeep <- sapply(c(1:length(toFind)), function(e) grep(toFind[e], foo)) ##??
# foo <- foo[toKeep]
# foo

## load performance data
dat <- list()
setwd(dirname)
system.time(
for(i in 1:length(foo)){ # length(foo)
  print(i)
  # dat[[i]] <- get(load(paste("./", foo[i], sep="")))
  temp <- get(load(paste("./", foo[i], sep="")))
  dat[[i]] <- temp$vals$terminal$corr.dat[snps.assoc[, i]]
  gc()
}
)

## REORDER dat s.t. order is numeric not character...
inds.new <- c(1:length(dat))
inds.ori <- sort(as.character(inds.new))
ord <- sapply(c(1:length(inds.new)), function(e) which(as.numeric(inds.ori) == inds.new[e]))

dat.ori <- dat
dat <- dat[ord]

#####################
## get snps.assoc: ##
#####################
# dat = performance
snps.assoc <- sapply(c(1:length(dat)), function(e) dat[[e]][[1]])

snps.assoc <- as.list(as.data.frame(snps.assoc))
# snps.assoc[[1]]
snps.assoc <- snps.assoc[ord]
# res <- dat

######################
## get score1.mean: ##
######################
# dat = res
score1.mean <- sapply(c(1:length(dat)), function(e) mean(dat[[e]]))
score1 <- dat

## SAVE ##
save(score1, file="./set2_1.42_ALL_terminal.score_snps.assoc.Rdata")
save(snps.assoc, file="./set2_1.42_ALL_snps.assoc.Rdata")

## bind score1.mean to args in evalStats:
evalStats <- cbind(evalStats[,1:10], score1.mean, evalStats[,11:18])
# table(round(evalStats$score1.mean, 1), evalStats$s, evalStats$af, evalStats$tree.type)

## SAVE: ##
# filename <- "./set3_1.120_ALL_s.ALL_af.ALL_treeWAS_evalStats.df.Rdata"
# save(evalStats, file=filename)

######################
## get n.phen.subs: ##
######################
# dat = phen.plot.col
n.phen.subs <- sapply(c(1:length(dat)), function(e) length(which(dat[[e]]$edges == "grey")))

head(evalStats)

## bind n.phen.subs to args in evalStats:
evalStats <- cbind(evalStats[,1:9], n.phen.subs, evalStats[,10:17])
# table(evalStats$n.phen.subs, evalStats$s, evalStats$tree.type)

## SAVE: ##
# filename <- "./set3_1.120_ALL_s.ALL_af.ALL_treeWAS_evalStats.df.Rdata"
# save(evalStats, file=filename)

# evalStats <- get(load("./set3_1.120_ALL_s.ALL_af.ALL_treeWAS_evalStats.df.Rdata"))

###############
## get args: ##
###############
# dat = args
args <- dat




##############
## CAREFUL: ##
##############
## (1) accuracy may need to be multiplied by 2 (if ncol(snps) was doubled for PLINK!)
## (2) names of performance etc. may need to be changed (if you want to use them...)
## ### (eg. if all treeWAS == "terminal" --> 2nd 1/3 = "simultaneous", 3rd 1/3 = "subsequent")

names(dat[[1]])

# ## CORRECT NAMES: ##
## treeWAS = 4:99
(99-3)/3 # 32
inds.terminal <- c(4:(4+31))
inds.simultaneous <- c((max(inds.terminal)+1):((max(inds.terminal)+1)+31))
inds.subsequent <- c((max(inds.simultaneous)+1):((max(inds.simultaneous)+1)+31))
noms.ori <- names(dat[[1]])
noms <- noms.ori

temp <- noms[inds.simultaneous]
temp2 <- paste("treeWAS.simultaneous", removeFirstN(temp, nchar("treeWAS.terminal")), sep="")

temp <- noms[inds.subsequent]
temp3 <- paste("treeWAS.subsequent", removeFirstN(temp, nchar("treeWAS.terminal")), sep="")

noms <- replace(noms, inds.simultaneous, temp2)
noms <- replace(noms, inds.subsequent, temp3)
noms

for(i in 1:length(dat)){
  names(dat[[i]]) <- noms
}

treeWAS <- vector("list", length=3)
names(treeWAS) <- c("terminal", "simultaneous", "subsequent")
treeWAS$terminal <- treeWAS$simultaneous <- treeWAS$subsequent <- list()


fisher.bonf <- fisher.fdr <-
  plink.bonf <- plink.fdr <- plink.gc.bonf <- plink.gc.fdr <- list()

for(i in 1:length(dat)){
  treeWAS$terminal[[i]] <- treeWAS$simultaneous[[i]] <- treeWAS$subsequent[[i]] <- list()

  for(e in 2:length(names(dat[[1]]))){
    ## Fisher:
    if(e == 2){
      fisher.bonf[[i]] <- dat[[i]][[e]]
    }
    if(e == 3){
      fisher.fdr[[i]] <- dat[[i]][[e]]
    }
    ## treeWAS:
    if(e %in% 4:35){
      treeWAS$terminal[[i]][[(length(treeWAS$terminal[[i]])+1)]] <- dat[[i]][[e]]
    }
    if(e %in% 36:67){
      treeWAS$simultaneous[[i]][[(length(treeWAS$simultaneous[[i]])+1)]] <- dat[[i]][[e]]
    }
    if(e %in% 68:99){
      treeWAS$subsequent[[i]][[(length(treeWAS$subsequent[[i]])+1)]] <- dat[[i]][[e]]
    }
    ## PLINK:
    if(e == 100){
      plink.bonf[[i]] <- dat[[i]][[e]]
    }
    if(e == 101){
      plink.fdr[[i]] <- dat[[i]][[e]]
    }
    if(e == 102){
      plink.gc.bonf[[i]] <- dat[[i]][[e]]
    }
    if(e == 103){
      plink.gc.fdr[[i]] <- dat[[i]][[e]]
    }
  } # end e for loop
  names(treeWAS$terminal[[i]]) <- removeFirstN(noms[4:35], nchar("treeWAS.terminal."))
  names(treeWAS$simultaneous[[i]]) <- removeFirstN(noms[36:67], nchar("treeWAS.simultaneous."))
  names(treeWAS$subsequent[[i]]) <- removeFirstN(noms[68:99], nchar("treeWAS.subsequent."))
} # end i for loop





## combine
# treeWAS <- do.call("rbind", treeWAS)
fisher.bonf <- do.call("rbind", fisher.bonf)
fisher.fdr <- do.call("rbind", fisher.fdr)
plink.bonf <- do.call("rbind", plink.bonf)
plink.fdr <- do.call("rbind", plink.fdr)
plink.gc.bonf <- do.call("rbind", plink.gc.bonf)
plink.gc.fdr <- do.call("rbind", plink.gc.fdr)

## Get dfs of nrow=length(dat),
## for each of the 32 thresh methods and each of the 3 treeWAS tests:
df <- list()
for(t in 1:length(treeWAS)){
  df[[t]] <- list()
  for(e in 1:32){
    df[[t]][[e]] <- list()
    for(i in 1:length(treeWAS[[t]])){
      df[[t]][[e]][[i]] <- treeWAS[[t]][[i]][[e]]
      # df[[t]][[i]] <- do.call("rbind", sapply(c(1:length(treeWAS[[t]])), function(e) treeWAS[[t]][[e]][[i]]))
    }
    df[[t]][[e]] <- do.call("rbind", df[[t]][[e]])
    ################
    ## TEMP -- FOR SOME DATASETS, NEED TO MULTIPLY ACCURACY BY TWO!!!!!!!!!!!!!!
    ## (though arguably accuracy is not a useful metric for genetic data.. )
    # df[[t]][[e]]$accuracy <- df[[t]][[e]]$accuracy*2

    #     acc <- df[[t]][[e]]$accuracy
    #     for(a in toChange){
    #       acc[a] <- acc[a]*2
    #     }
    #     df[[t]][[e]]$accuracy <- acc
    ################
  }
  names(df[[t]]) <- names(treeWAS$terminal[[1]])
}
names(df) <- c("terminal", "simultaneous", "subsequent")

treeWAS.df <- df
str(treeWAS.df[[1]])

str(df[[1]])
length(df[[1]])
nrow(df[[1]][[1]])

## CHECK -- accuracy:
# DF <- df[[1]][[1]]
# acc.ori <- DF$accuracy
# acc.new <- ((DF$sensitivity*10) + (DF$specificity*9990))/10000
# toChange <- which(acc.ori <= 0.5)



## CORRECT ACCURACY:
# fisher.bonf$accuracy <- fisher.bonf$accuracy*2
# fisher.fdr$accuracy <- fisher.fdr$accuracy*2
#
# plink.bonf$accuracy <- plink.bonf$accuracy*2
# plink.fdr$accuracy <- plink.fdr$accuracy*2
# plink.gc.bonf$accuracy <- plink.gc.bonf$accuracy*2
# plink.gc.fdr$accuracy <- plink.gc.fdr$accuracy*2

## OR -- correct a SUBSET of accuracy...
# fisher.bonf$accuracy[toChange] <- fisher.bonf$accuracy[toChange]*2
# fisher.fdr$accuracy[toChange] <- fisher.fdr$accuracy[toChange]*2
#
# plink.bonf$accuracy[toChange] <- plink.bonf$accuracy[toChange]*2
# plink.fdr$accuracy[toChange] <- plink.fdr$accuracy[toChange]*2
# plink.gc.bonf$accuracy[toChange] <- plink.gc.bonf$accuracy[toChange]*2
# plink.gc.fdr$accuracy[toChange] <- plink.gc.fdr$accuracy[toChange]*2

##########
## save ##    #######    #######    #######    #######    #######    #######
##########
tree.type <- "ALL"
dir.n <- "1.40"
# getwd()
filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_treeWAS_evalStats.Rdata", sep="")
save(treeWAS.df, file=filename)
filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_fisher.bonf_evalStats.Rdata", sep="")
save(fisher.bonf, file=filename)
filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_fisher.fdr_evalStats.Rdata", sep="")
save(fisher.fdr, file=filename)
filename <- paste("./", set.n, "_", dir.n, "_", tree.type,  "_plink.bonf_evalStats.Rdata", sep="")
save(plink.bonf, file=filename)
filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_plink.fdr_evalStats.Rdata", sep="")
save(plink.fdr, file=filename)
filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_plink.gc.bonf_evalStats.Rdata", sep="")
save(plink.gc.bonf, file=filename)
filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_plink.gc.fdr_evalStats.Rdata", sep="")
save(plink.gc.fdr, file=filename)

## summarise
summary(fisher.bonf)
summary(fisher.fdr)
summary(plink.bonf)
summary(plink.fdr)
summary(plink.gc.bonf)
summary(plink.gc.fdr)

## summarise treeWAS by test & thresh...
treeWAS.sum <- list()
for(t in 1:length(treeWAS.df)){
  treeWAS.sum[[t]] <- list()
  for(e in 1:length(treeWAS.df[[t]])){
    treeWAS.sum[[t]][[e]] <- summary(treeWAS.df[[t]][[e]])
  }
  names(treeWAS.sum[[t]]) <- names(treeWAS.df[[t]])
}
names(treeWAS.sum) <- names(treeWAS.df)

treeWAS.sum



##############################################################################################################
##########
## plot ##    #######    #######    #######    #######    #######    #######    #######    #######    #######
##########
##############################################################################################################

###################################
## BY THRESHOLD-SELECTION METHOD ##
###################################

###################
## GET evalStats ##
###################
## (for each sel method, by test)
evalStats <- vector("list", length=3)
names(evalStats) <- c("terminal", "simultaneous", "subsequent")
evalStats$terminal <- evalStats$simultaneous <- evalStats$subsequent <- list()

## Get repeated thresh-sel names column: ##
nom.ori <- names(treeWAS.df[[1]])
nom <- list()
for(i in 1:length(nom.ori)){
  nom[[i]] <- rep(nom.ori[[i]], nrow(treeWAS.df[[1]][[1]]))
}
thresh.sel <- as.vector(unlist(nom))
test <- thresh.sel # using test instead of thresh.sel for convenience w plot code below...

## PLUS: break test into its composite parts:
################

temp <- strsplit(test, "[.]")
temp <- do.call("rbind", temp)

# pval <- temp[,3] # 05
pval <- paste(temp[,2], temp[,3], sep=".") # 0.05
pval.correct <- temp[,4] # bonf # fdr
pval.by  <- temp[,5] # count # density
n.snps.sim <- paste(temp[,6], "x", sep="") # 1x # 10x.n.snps

thresh.sel <- data.frame(test, pval, pval.correct, pval.by, n.snps.sim)

tree.type <- c(rep("coal", 20), rep("rtree", 20))

# args <- data.frame(tree.type, s, af)

params <- cbind(thresh.sel, tree.type)

################


## Get evalStats and make df: ##
for(t in 1:length(treeWAS.df)){
  df <- treeWAS.df[[t]]

  accuracy <- specificity <- FPR <- FNR <- sensitivity <- PPV <- FDR <- list()

  for(i in 1:length(df)){
    accuracy[[i]] <- df[[i]]$accuracy
    specificity[[i]] <- df[[i]]$specificity
    FPR[[i]] <- df[[i]]$FPR
    FNR[[i]] <- df[[i]]$FNR
    sensitivity[[i]] <- df[[i]]$sensitivity
    PPV[[i]] <- df[[i]]$PPV
    FDR[[i]] <- df[[i]]$FDR
  } # end (i) for loop

  accuracy <- as.vector(unlist(accuracy))
  specificity <- as.vector(unlist(specificity))
  FPR <- as.vector(unlist(FPR))
  FNR <- as.vector(unlist(FNR))
  sensitivity <- as.vector(unlist(sensitivity))
  PPV <- as.vector(unlist(PPV))
  FDR <- as.vector(unlist(FDR))

  F1.score <- 2*((PPV*sensitivity) / (PPV+sensitivity))

  evalStats[[t]] <- data.frame(params, accuracy, specificity, FPR, FNR, sensitivity, PPV, FDR, F1.score)

} # end (t) for loop

##########
## save ##    #######    #######    #######    #######    #######    #######
##########
tree.type <- "ALL"
dir.n <- "1.40"
# getwd()
filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_treeWAS_evalStats.df.Rdata", sep="")
save(evalStats, file=filename)

#######    #######    #######    #######    #######    #######    #######

## append treeWAS.test to treeWAS evalStats and bind:
treeWAS.test <- c(rep("terminal", nrow(evalStats[[1]])),
                  rep("simultaneous", nrow(evalStats[[2]])),
                  rep("subsequent", nrow(evalStats[[3]])))
temp <- do.call("rbind", evalStats)
temp <- data.frame(treeWAS.test, temp)
rownames(temp) <- NULL
evalStats <- temp

## append args to other tests:
fisher.bonf <- data.frame(tree.type, fisher.bonf)
fisher.fdr <- data.frame(tree.type, fisher.fdr)
plink.bonf <- data.frame(tree.type, plink.bonf)
plink.fdr <- data.frame(tree.type, plink.fdr)
plink.gc.bonf <- data.frame(tree.type, plink.gc.bonf)
plink.gc.fdr <- data.frame(tree.type, plink.gc.fdr)


##########
## save ##    #######    #######    #######    #######    #######    #######
##########
tree.type <- "ALL"
dir.n <- "1.40"

filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_treeWAS_evalStats.df.all.tests.Rdata", sep="")
save(evalStats, file=filename)

filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_fisher.bonf_evalStats.df.Rdata", sep="")
save(fisher.bonf, file=filename)
filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_fisher.fdr_evalStats.df.Rdata", sep="")
save(fisher.fdr, file=filename)
filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_plink.bonf_evalStats.df.Rdata", sep="")
save(plink.bonf, file=filename)
filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_plink.fdr_evalStats.df.Rdata", sep="")
save(plink.fdr, file=filename)
filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_plink.gc.bonf_evalStats.df.Rdata", sep="")
save(plink.gc.bonf, file=filename)
filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_plink.gc.fdr_evalStats.df.Rdata", sep="")
save(plink.gc.fdr, file=filename)

###################
## BEESWARM PLOT ##   #########################   #########################   #########################
###################
## uses BOTH df and beeswarm dataframes...

# library(beeswarm)
# library(plyr)
# library(ggplot2)

#############
## MY DATA ##
#############

evalStats <- get(load("C:/Cait 2016/Work/Xavier/Sims/set3/set3_1.120_ALL_s.ALL_af.ALL_treeWAS_evalStats.df.Rdata"))
df <- evalStats

df$test <- mapvalues(df$test, from = levels(df$test), to = c(1:32))
evalStats.ori <- evalStats
evalStats <- df

# t <- 3
# df <- evalStats[[t]]


treeWAS.tests <- c("terminal", "simultaneous", "subsequent")
Y <- c("F1.score", "PPV", "sensitivity", "FPR")

for(y in 1:length(Y)){ #
  for(t in 1:3){
    df <- evalStats
    df <- df[df$treeWAS.test==treeWAS.tests[t] & df$tree.type=="rtree" & df$s==0.5 & df$af==5, ] #
    beeswarm.plot(y=Y[y], x="test", df, #y.lab="Sensitivity",
                  pt.size=4, x.text=TRUE)
  }
}
#

## SAVE ##

## TO DO: (RE)SAVE ALL PLOTS WITH (a) NO y-label and (b) REPLACE x-axis text with numbers 1:n.tests
## so we can make a 4- (or 2?) panel figure....
## (later, if needed, could put labels corresponding to numbers and colours alongside in a legend.. )

# set3_101_110_rtree_s_0_5_af_5_treeWAS_terminal_F1score
# set3_101_110_rtree_s_0_5_af_5_treeWAS_simultaneous_F1score
# set3_101_110_rtree_s_0_5_af_5_treeWAS_subsequent_F1score
#
# set3_101_110_rtree_s_0_5_af_5_treeWAS_terminal_PPV
# set3_101_110_rtree_s_0_5_af_5_treeWAS_simultaneous_PPV
# set3_101_110_rtree_s_0_5_af_5_treeWAS_subsequent_PPV
#
# set3_101_110_rtree_s_0_5_af_5_treeWAS_terminal_sensitivity
# set3_101_110_rtree_s_0_5_af_5_treeWAS_simultaneous_sensitivity
# set3_101_110_rtree_s_0_5_af_5_treeWAS_subsequent_sensitivity
#
# set3_101_110_rtree_s_0_5_af_5_treeWAS_terminal_FPR
# set3_101_110_rtree_s_0_5_af_5_treeWAS_simultaneous_FPR
# set3_101_110_rtree_s_0_5_af_5_treeWAS_subsequent_FPR

#################

# set3_111_120_rtree_s_0_1_af_5_treeWAS_terminal_F1score
# set3_111_120_rtree_s_0_1_af_5_treeWAS_simultaneous_F1score
# set3_111_120_rtree_s_0_1_af_5_treeWAS_subsequent_F1score
#
# set3_111_120_rtree_s_0_1_af_5_treeWAS_terminal_PPV
# set3_111_120_rtree_s_0_1_af_5_treeWAS_simultaneous_PPV
# set3_111_120_rtree_s_0_1_af_5_treeWAS_subsequent_PPV
#
# set3_111_120_rtree_s_0_1_af_5_treeWAS_terminal_sensitivity
# set3_111_120_rtree_s_0_1_af_5_treeWAS_simultaneous_sensitivity
# set3_111_120_rtree_s_0_1_af_5_treeWAS_subsequent_sensitivity
#
# set3_111_120_rtree_s_0_1_af_5_treeWAS_terminal_FPR
# set3_111_120_rtree_s_0_1_af_5_treeWAS_simultaneous_FPR
# set3_111_120_rtree_s_0_1_af_5_treeWAS_subsequent_FPR

####
# set3_1_120_ALL_treeWAS_terminal_F1score
# set3_1_120_ALL_treeWAS_simultaneous_F1score
# set3_1_120_ALL_treeWAS_subsequent_F1score
#
# set3_1_120_ALL_treeWAS_terminal_PPV
# set3_1_120_ALL_treeWAS_simultaneous_PPV
# set3_1_120_ALL_treeWAS_subsequent_PPV
#
# set3_1_120_ALL_treeWAS_terminal_sensitivity
# set3_1_120_ALL_treeWAS_simultaneous_sensitivity
# set3_1_120_ALL_treeWAS_subsequent_sensitivity
#
# set3_1_120_ALL_treeWAS_terminal_FPR
# set3_1_120_ALL_treeWAS_simultaneous_FPR
# set3_1_120_ALL_treeWAS_subsequent_FPR

####
# set1_21_40_rtree_treeWAS_terminal_F1score_leg
# set1_21_40_rtree_treeWAS_simultaneous_F1score_leg
# set1_21_40_rtree_treeWAS_subsequent_F1score_leg # NA

# set1_21_40_rtree_treeWAS_terminal_PPV_leg
# set1_21_40_rtree_treeWAS_simultaneous_PPV_leg
# set1_21_40_rtree_treeWAS_subsequent_PPV_leg
#
# set1_21_40_rtree_treeWAS_terminal_sensitivity_leg
# set1_21_40_rtree_treeWAS_simultaneous_sensitivity_leg
# set1_21_40_rtree_treeWAS_subsequent_sensitivity_leg
#
# set1_21_40_rtree_treeWAS_terminal_FPR_leg
# set1_21_40_rtree_treeWAS_simultaneous_FPR_leg
# set1_21_40_rtree_treeWAS_subsequent_FPR_leg

#########



####
# set3_1_120_ALL_treeWAS_terminal_F1score_leg
# set3_1_120_ALL_treeWAS_simultaneous_F1score_leg
# set3_1_120_ALL_treeWAS_subsequent_F1score_leg
#
# set3_1_120_ALL_treeWAS_terminal_PPV_leg
# set3_1_120_ALL_treeWAS_simultaneous_PPV_leg
# set3_1_120_ALL_treeWAS_subsequent_PPV_leg
#
# set3_1_120_ALL_treeWAS_terminal_sensitivity_leg
# set3_1_120_ALL_treeWAS_simultaneous_sensitivity_leg
# set3_1_120_ALL_treeWAS_subsequent_sensitivity_leg
#
# set3_1_120_ALL_treeWAS_terminal_FPR_leg
# set3_1_120_ALL_treeWAS_simultaneous_FPR_leg
# set3_1_120_ALL_treeWAS_subsequent_FPR_leg


######


##




###################
## beeswarm.plot ##
###################
beeswarm.plot <- function(y="sensitivity", x="test", df, y.lab=NULL,
                          pt.size=4, x.text=FALSE){

  require(beeswarm)

  if(is.null(y.lab)) y.lab <- y

  ## Y ~ X ??
  fm <- as.formula(paste(y, x, sep=" ~ "))

  beeswarm <- beeswarm(fm,
                       data = df,
                       #method="swarm", # swarm square hex center
                       #priority="descending", ## ONLY for SWARM method...
                       method="center", # swarm square hex center
                       #priority="descending", ## ONLY for SWARM method...
                       pwcol = eval(parse(text=x)),
                       #col = myCol, ## to set w funky colours (INSTEAD of pwcol = test)
                       ylim = c(-0.001,1), # otherwise ggplot can't plot ZERO values --> NAs
                       las=2,
                       cex=0.8,
                       corral = "omit",
                       do.plot = FALSE) # none gutter wrap omit
  # head(beeswarm)

  ######################################################
  ## Find and Replace OUTLIERS(' symbols in plot...): ##
  ######################################################
  outliers <- outlier.vals <- PCH <- list()

  if(!all(as.character(beeswarm$col) %in% levels(df[,x]))){
    foo <- beeswarm$col
    foo <- levels(df[,x])[foo]
    beeswarm$col <- factor(foo, levels=levels(df[,x]))
  }else{
    beeswarm$col <- factor(beeswarm$col)
  }

  noms <- as.character(levels(beeswarm$col))

  ## FOR LOOP ##
  for(i in 1:length(noms)){
    #i <- 1
    # get vals for variable (and boxplot)
    val <- beeswarm$y[which(beeswarm$col==noms[i])]
    #boxplot(val, ylim=c(-0.001, 1))
    if(length(val) == 0){
      PCH[[i]] <- NULL
      outliers[[i]] <- NULL
    }else{
      PCH[[i]] <- rep(16, length(val)) # standard filled circle...

      ## get median
      M <- as.numeric(quantile(val, 0.5))
      # get lower 25 of box
      Q25 <- as.numeric(quantile(val, 0.25))
      # get upper 75 of box
      Q75 <- as.numeric(quantile(val, 0.75))
      # get box length
      box <- Q75-Q25

      if(box == 0) box <- 0.0000001

      # with a coef of 1.5 (the default for boxplots), identify outlying values
      outliers[[i]] <- c(which(val < Q25-(1.5*box)), which(val > Q75+(1.5*box)))
      # get values of outliers
      if(length(outliers[[i]]) > 0){
        outlier.vals[[i]] <- val[outliers[[i]]]
        PCH[[i]] <- replace(PCH[[i]], outliers[[i]], 17) # replace with triangle...
      }else{
        outlier.vals[[i]] <- NULL
      }
    }
  } # end for loop

  #outliers
  PCH <- as.vector(unlist(PCH))
  # PCH


  #########################################################################################################
  ######################
  ## plots, layers... ##
  ######################

  if(x.text == FALSE){

    ################
    ## NO X-TEXT: ##
    ################

    beeswarm.plot <- ggplot(beeswarm, aes(x, y)) +
      xlab("") +
      guides(fill=FALSE) +
      scale_x_discrete(drop=FALSE) +
      scale_y_continuous(y.lab, limits=c(0,1))  # expression("char")

    beeswarm.plot2 <- beeswarm.plot +
      guides(fill=FALSE) +
      geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.25,
                   outlier.shape = 17,
                   outlier.size=0) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), axis.title.x=element_blank(),
          axis.text.y = element_text(size=13),
          # axis.title.y=element_text(size=18),
          axis.title.y=element_blank(),
          legend.position="none")

    beeswarm.plot3 <- beeswarm.plot2 +
      geom_point(data=beeswarm, aes(colour = col), pch = PCH, size=pt.size, na.rm=TRUE, alpha=0.6) +
      guides(fill=FALSE) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
          axis.text.y = element_text(size=13),
          # axis.title.y=element_text(size=18),
          axis.title.y=element_blank(),
          legend.position="none")

    beeswarm.plot4 <- beeswarm.plot3 +
      guides(fill=FALSE) +
      geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.0, fatten=3,
                   outlier.shape = 17,
                   outlier.size=0) +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
            axis.text.y = element_text(size=13),
            # axis.title.y=element_text(size=18),
            axis.title.y=element_blank(),
            legend.position="none")
  }else{

    ##################
    ## WITH X-TEXT: ##
    ##################

    beeswarm.plot <- ggplot(beeswarm, aes(x, y)) +
      xlab("") +
      guides(fill=FALSE) +
      scale_x_discrete(drop=FALSE) +
      scale_y_continuous(y.lab, limits=c(0,1)) # expression("char")

    beeswarm.plot2 <- beeswarm.plot +
      guides(fill=FALSE) +
      geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.25,
                   outlier.shape = 17,
                   outlier.size=0) +
      theme(axis.text.x = element_text(size=13),
            # axis.text.x = element_text(angle=65, hjust=1, vjust=0.95, size=10),
            axis.text.y = element_text(size=13),
            axis.title.y = element_text(size=18),
            axis.title.x = element_blank(),
            legend.position="none")

    beeswarm.plot3 <- beeswarm.plot2 +
      geom_point(data=beeswarm, aes(colour = col), pch = PCH, size=pt.size, na.rm=TRUE, alpha=0.6) +
      guides(fill=FALSE) +
      theme(axis.text.x = element_text(size=13),
            # axis.text.x = element_text(angle=65, hjust=1, vjust=0.95, size=10),
            axis.text.y = element_text(size=13),
            axis.title.y = element_text(size=18),
            axis.title.x = element_blank(),
            legend.position="none")

    beeswarm.plot4 <- beeswarm.plot3 +
      guides(fill=FALSE) +
      geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.0, fatten=3,
                   outlier.shape = 17,
                   outlier.size=0) +
      theme(axis.text.x = element_text(size=13),
            # axis.text.x = element_text(angle=65, hjust=1, vjust=0.95, size=10),
            axis.text.y = element_text(size=13),
            axis.title.y = element_text(size=18),
            axis.title.x = element_blank(),
            legend.position="none")
  }


  ## PRINT PLOT ##
  plot(beeswarm.plot4)

} # end beeswarm.plot



## CAREFUL: CHECK that your outliers are really in (ALMOST) the right place by plotting beeswarm2 w oultier.cex=2, and not outlier.size=0!
beeswarm.plot.terminal <- beeswarm.plot4
beeswarm.plot.simultaneous <- beeswarm.plot4
beeswarm.plot.subsequent <- beeswarm.plot4


# set3_1.120_ALL_treeWAS.subsequent_FPR

##########################################################################################################













# ##############
# ## GET DATA ##    #######    #######    #######    #######    #######    #######    #######    #######    #######    #######
# ##############
#
# setwd("C:/Cait 2016/Work/Xavier/Sims/set3/")
#
# set.n <- "set3"
#
# dirname <- paste("C:/Cait 2016/Work/Xavier/Sims/", set.n, sep="")
# # dirname <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/", set.n, sep="")
#
# foo <- dir(dirname)
# foo
# ## get all performance Rdata names
# # toKeep <- grep("performance", foo) ##??
# # toKeep <- grep("phen.plot.col", foo) ##??
# toKeep <- grep("_res", foo) ##??
# foo <- foo[toKeep]
# foo
#
# ## keep only subset w same args:
# # toFind <- paste(set.n, "_", c(51:81), "_", sep="")
# # toFind <- paste(set.n, "_", c(21:30), "_", sep="")
# # toKeep <- sapply(c(1:length(toFind)), function(e) grep(toFind[e], foo)) ##??
# # foo <- foo[toKeep]
# # foo
#
# ## load performance data
# dat <- list()
# setwd(dirname)
# system.time(
#   for(i in 1:length(foo)){ # length(foo)
#     print(i)
#     # dat[[i]] <- get(load(paste("./", foo[i], sep="")))
#     temp <- get(load(paste("./", foo[i], sep="")))
#     dat[[i]] <- temp$vals$terminal$corr.dat[snps.assoc[, i]]
#     gc()
#   }
# )
#
# ## REORDER dat s.t. order is numeric not character...
# inds.new <- c(1:length(dat))
# inds.ori <- sort(as.character(inds.new))
# ord <- sapply(c(1:length(inds.new)), function(e) which(as.numeric(inds.ori) == inds.new[e]))
#
# dat.ori <- dat
# dat <- dat[ord]
#
# #####################
# ## get snps.assoc: ##
# #####################
# # dat = performance
# snps.assoc <- sapply(c(1:length(dat)), function(e) dat[[e]][[1]])
#
# snps.assoc <- as.list(as.data.frame(snps.assoc))
# # snps.assoc[[1]]
# snps.assoc <- snps.assoc[ord]
# # res <- dat
#
# ######################
# ## get score1.mean: ##
# ######################
# # dat = res
# score1.mean <- sapply(c(1:length(dat)), function(e) mean(dat[[e]]))
# score1 <- dat
#
# ## SAVE ##
# save(score1, file="./set3_1.120_ALL_terminal.score_snps.assoc.Rdata")
# save(snps.assoc, file="./set3_1.120_ALL_snps.assoc.Rdata")
#
# ## bind score1.mean to args in evalStats:
# evalStats <- cbind(evalStats[,1:10], score1.mean, evalStats[,11:18])
# # table(round(evalStats$score1.mean, 1), evalStats$s, evalStats$af, evalStats$tree.type)
#
# ## SAVE: ##
# # filename <- "./set3_1.120_ALL_s.ALL_af.ALL_treeWAS_evalStats.df.Rdata"
# # save(evalStats, file=filename)
#
# ######################
# ## get n.phen.subs: ##
# ######################
# n.phen.subs <- sapply(c(1:length(dat)), function(e) length(which(dat[[e]]$edges == "grey")))
#
# head(evalStats)
#
# ## bind n.phen.subs to args in evalStats:
# evalStats <- cbind(evalStats[,1:9], n.phen.subs, evalStats[,10:17])
# # table(evalStats$n.phen.subs, evalStats$s, evalStats$tree.type)
#
# ## SAVE: ##
# # filename <- "./set3_1.120_ALL_s.ALL_af.ALL_treeWAS_evalStats.df.Rdata"
# # save(evalStats, file=filename)
#
# evalStats <- get(load("./set3_1.120_ALL_s.ALL_af.ALL_treeWAS_evalStats.df.Rdata"))
#
# ######################
# ## get avg. score1: ##
# ######################
#
#
#
#
# ##############
# ## CAREFUL: ##
# ##############
# ## (1) accuracy may need to be multiplied by 2 (if ncol(snps) was doubled for PLINK!)
# ## (2) names of performance etc. may need to be changed (if you want to use them...)
# ## ### (eg. if all treeWAS == "terminal" --> 2nd 1/3 = "simultaneous", 3rd 1/3 = "subsequent")
#
# names(dat[[1]])
#
# # ## CORRECT NAMES: ##
# ## treeWAS = 4:99
# (99-3)/3 # 32
# inds.terminal <- c(4:(4+31))
# inds.simultaneous <- c((max(inds.terminal)+1):((max(inds.terminal)+1)+31))
# inds.subsequent <- c((max(inds.simultaneous)+1):((max(inds.simultaneous)+1)+31))
# noms.ori <- names(dat[[1]])
# noms <- noms.ori
#
# temp <- noms[inds.simultaneous]
# temp2 <- paste("treeWAS.simultaneous", removeFirstN(temp, nchar("treeWAS.terminal")), sep="")
#
# temp <- noms[inds.subsequent]
# temp3 <- paste("treeWAS.subsequent", removeFirstN(temp, nchar("treeWAS.terminal")), sep="")
#
# noms <- replace(noms, inds.simultaneous, temp2)
# noms <- replace(noms, inds.subsequent, temp3)
# noms
#
# for(i in 1:length(dat)){
#   names(dat[[i]]) <- noms
# }
#
# treeWAS <- vector("list", length=3)
# names(treeWAS) <- c("terminal", "simultaneous", "subsequent")
# treeWAS$terminal <- treeWAS$simultaneous <- treeWAS$subsequent <- list()
#
#
# fisher.bonf <- fisher.fdr <-
#   plink.bonf <- plink.fdr <- plink.gc.bonf <- plink.gc.fdr <- list()
#
# for(i in 1:length(dat)){
#   treeWAS$terminal[[i]] <- treeWAS$simultaneous[[i]] <- treeWAS$subsequent[[i]] <- list()
#
#   for(e in 2:length(names(dat[[1]]))){
#     ## Fisher:
#     if(e == 2){
#       fisher.bonf[[i]] <- dat[[i]][[e]]
#     }
#     if(e == 3){
#       fisher.fdr[[i]] <- dat[[i]][[e]]
#     }
#     ## treeWAS:
#     if(e %in% 4:35){
#       treeWAS$terminal[[i]][[(length(treeWAS$terminal[[i]])+1)]] <- dat[[i]][[e]]
#     }
#     if(e %in% 36:67){
#       treeWAS$simultaneous[[i]][[(length(treeWAS$simultaneous[[i]])+1)]] <- dat[[i]][[e]]
#     }
#     if(e %in% 68:99){
#       treeWAS$subsequent[[i]][[(length(treeWAS$subsequent[[i]])+1)]] <- dat[[i]][[e]]
#     }
#     ## PLINK:
#     if(e == 100){
#       plink.bonf[[i]] <- dat[[i]][[e]]
#     }
#     if(e == 101){
#       plink.fdr[[i]] <- dat[[i]][[e]]
#     }
#     if(e == 102){
#       plink.gc.bonf[[i]] <- dat[[i]][[e]]
#     }
#     if(e == 103){
#       plink.gc.fdr[[i]] <- dat[[i]][[e]]
#     }
#   } # end e for loop
#   names(treeWAS$terminal[[i]]) <- removeFirstN(noms[4:35], nchar("treeWAS.terminal."))
#   names(treeWAS$simultaneous[[i]]) <- removeFirstN(noms[36:67], nchar("treeWAS.simultaneous."))
#   names(treeWAS$subsequent[[i]]) <- removeFirstN(noms[68:99], nchar("treeWAS.subsequent."))
# } # end i for loop
#
#
#
#
#
# ## combine
# # treeWAS <- do.call("rbind", treeWAS)
# fisher.bonf <- do.call("rbind", fisher.bonf)
# fisher.fdr <- do.call("rbind", fisher.fdr)
# plink.bonf <- do.call("rbind", plink.bonf)
# plink.fdr <- do.call("rbind", plink.fdr)
# plink.gc.bonf <- do.call("rbind", plink.gc.bonf)
# plink.gc.fdr <- do.call("rbind", plink.gc.fdr)
#
# ## Get dfs of nrow=length(dat),
# ## for each of the 32 thresh methods and each of the 3 treeWAS tests:
# df <- list()
# for(t in 1:length(treeWAS)){
#   df[[t]] <- list()
#   for(e in 1:32){
#     df[[t]][[e]] <- list()
#     for(i in 1:length(treeWAS[[t]])){
#       df[[t]][[e]][[i]] <- treeWAS[[t]][[i]][[e]]
#       # df[[t]][[i]] <- do.call("rbind", sapply(c(1:length(treeWAS[[t]])), function(e) treeWAS[[t]][[e]][[i]]))
#     }
#     df[[t]][[e]] <- do.call("rbind", df[[t]][[e]])
#     ################
#     ## TEMP -- FOR SOME DATASETS, NEED TO MULTIPLY ACCURACY BY TWO!!!!!!!!!!!!!!
#     ## (though arguably accuracy is not a useful metric for genetic data.. )
#     # df[[t]][[e]]$accuracy <- df[[t]][[e]]$accuracy*2
#
#     acc <- df[[t]][[e]]$accuracy
#     for(a in toChange){
#       acc[a] <- acc[a]*2
#     }
#     df[[t]][[e]]$accuracy <- acc
#     ################
#   }
#   names(df[[t]]) <- names(treeWAS$terminal[[1]])
# }
# names(df) <- c("terminal", "simultaneous", "subsequent")
#
# treeWAS.df <- df
# str(treeWAS.df[[1]])
#
# str(df[[1]])
# length(df[[1]])
# nrow(df[[1]][[1]])
#
# ## CHECK -- accuracy:
# # DF <- df[[1]][[1]]
# # acc.ori <- DF$accuracy
# # acc.new <- ((DF$sensitivity*10) + (DF$specificity*9990))/10000
# # toChange <- which(acc.ori <= 0.5)
#
#
#
# ## CORRECT ACCURACY:
# # fisher.bonf$accuracy <- fisher.bonf$accuracy*2
# # fisher.fdr$accuracy <- fisher.fdr$accuracy*2
# #
# # plink.bonf$accuracy <- plink.bonf$accuracy*2
# # plink.fdr$accuracy <- plink.fdr$accuracy*2
# # plink.gc.bonf$accuracy <- plink.gc.bonf$accuracy*2
# # plink.gc.fdr$accuracy <- plink.gc.fdr$accuracy*2
#
# ## OR -- correct a SUBSET of accuracy...
# # fisher.bonf$accuracy[toChange] <- fisher.bonf$accuracy[toChange]*2
# # fisher.fdr$accuracy[toChange] <- fisher.fdr$accuracy[toChange]*2
# #
# # plink.bonf$accuracy[toChange] <- plink.bonf$accuracy[toChange]*2
# # plink.fdr$accuracy[toChange] <- plink.fdr$accuracy[toChange]*2
# # plink.gc.bonf$accuracy[toChange] <- plink.gc.bonf$accuracy[toChange]*2
# # plink.gc.fdr$accuracy[toChange] <- plink.gc.fdr$accuracy[toChange]*2
#
# ##########
# ## save ##    #######    #######    #######    #######    #######    #######
# ##########
# tree.type <- "ALL"
# s <- "ALL"
# af <- "ALL"
# dir.n <- "1.120"
# # getwd()
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_treeWAS_evalStats.Rdata", sep="")
# save(treeWAS.df, file=filename)
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_fisher.bonf_evalStats.Rdata", sep="")
# save(fisher.bonf, file=filename)
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_fisher.fdr_evalStats.Rdata", sep="")
# save(fisher.fdr, file=filename)
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_plink.bonf_evalStats.Rdata", sep="")
# save(plink.bonf, file=filename)
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_plink.fdr_evalStats.Rdata", sep="")
# save(plink.fdr, file=filename)
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_plink.gc.bonf_evalStats.Rdata", sep="")
# save(plink.gc.bonf, file=filename)
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_plink.gc.fdr_evalStats.Rdata", sep="")
# save(plink.gc.fdr, file=filename)
#
# ## summarise
# summary(fisher.bonf)
# summary(fisher.fdr)
# summary(plink.bonf)
# summary(plink.fdr)
# summary(plink.gc.bonf)
# summary(plink.gc.fdr)
#
# ## summarise treeWAS by test & thresh...
# treeWAS.sum <- list()
# for(t in 1:length(treeWAS.df)){
#   treeWAS.sum[[t]] <- list()
#   for(e in 1:length(treeWAS.df[[t]])){
#     treeWAS.sum[[t]][[e]] <- summary(treeWAS.df[[t]][[e]])
#   }
#   names(treeWAS.sum[[t]]) <- names(treeWAS.df[[t]])
# }
# names(treeWAS.sum) <- names(treeWAS.df)
#
# treeWAS.sum
#
#
#
# ##############################################################################################################
# ##########
# ## plot ##    #######    #######    #######    #######    #######    #######    #######    #######    #######
# ##########
# ##############################################################################################################
#
# ###################################
# ## BY THRESHOLD-SELECTION METHOD ##
# ###################################
#
# ###################
# ## GET evalStats ##
# ###################
# ## (for each sel method, by test)
# evalStats <- vector("list", length=3)
# names(evalStats) <- c("terminal", "simultaneous", "subsequent")
# evalStats$terminal <- evalStats$simultaneous <- evalStats$subsequent <- list()
#
# ## Get repeated thresh-sel names column: ##
# nom.ori <- names(treeWAS.df[[1]])
# nom <- list()
# for(i in 1:length(nom.ori)){
#   nom[[i]] <- rep(nom.ori[[i]], nrow(treeWAS.df[[1]][[1]]))
# }
# thresh.sel <- as.vector(unlist(nom))
# test <- thresh.sel # using test instead of thresh.sel for convenience w plot code below...
#
# ## PLUS: break test into its composite parts:
# ################
#
# temp <- strsplit(test, "[.]")
# temp <- do.call("rbind", temp)
#
# # pval <- temp[,3] # 05
# pval <- paste(temp[,2], temp[,3], sep=".") # 0.05
# pval.correct <- temp[,4] # bonf # fdr
# pval.by  <- temp[,5] # count # density
# n.snps.sim <- paste(temp[,6], "x", sep="") # 1x # 10x.n.snps
#
# thresh.sel <- data.frame(test, pval, pval.correct, pval.by, n.snps.sim)
#
# tree.type <- c(rep("coal", 40), rep("rtree", 41), rep("coal", 19), rep("rtree", 20))
# s <- c(rep(1, 10), rep(10, 10), rep(1, 10), rep(10, 10), rep(2, 10), rep(1, 31), rep(0.5, 10), rep(0.1, 9), rep(0.5, 10), rep(0.1, 10))
# af <- c(rep(2, 20), rep(5, 100))
#
# args <- data.frame(tree.type, s, af)
#
# params <- cbind(thresh.sel, args)
#
# ################
#
#
# ## Get evalStats and make df: ##
# for(t in 1:length(treeWAS.df)){
#   df <- treeWAS.df[[t]]
#
#   accuracy <- specificity <- FPR <- FNR <- sensitivity <- PPV <- FDR <- list()
#
#   for(i in 1:length(df)){
#     accuracy[[i]] <- df[[i]]$accuracy
#     specificity[[i]] <- df[[i]]$specificity
#     FPR[[i]] <- df[[i]]$FPR
#     FNR[[i]] <- df[[i]]$FNR
#     sensitivity[[i]] <- df[[i]]$sensitivity
#     PPV[[i]] <- df[[i]]$PPV
#     FDR[[i]] <- df[[i]]$FDR
#   } # end (i) for loop
#
#   accuracy <- as.vector(unlist(accuracy))
#   specificity <- as.vector(unlist(specificity))
#   FPR <- as.vector(unlist(FPR))
#   FNR <- as.vector(unlist(FNR))
#   sensitivity <- as.vector(unlist(sensitivity))
#   PPV <- as.vector(unlist(PPV))
#   FDR <- as.vector(unlist(FDR))
#
#   F1.score <- 2*((PPV*sensitivity) / (PPV+sensitivity))
#
#   evalStats[[t]] <- data.frame(params, accuracy, specificity, FPR, FNR, sensitivity, PPV, FDR, F1.score)
#
# } # end (t) for loop
#
# ##########
# ## save ##    #######    #######    #######    #######    #######    #######
# ##########
# tree.type <- "ALL"
# s <- "ALL"
# af <- "ALL"
# dir.n <- "1.120"
# # getwd()
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_treeWAS_evalStats.df.Rdata", sep="")
# save(evalStats, file=filename)
#
# #######    #######    #######    #######    #######    #######    #######
#
# ## append treeWAS.test to treeWAS evalStats and bind:
# treeWAS.test <- c(rep("terminal", nrow(evalStats[[1]])),
#                   rep("simultaneous", nrow(evalStats[[2]])),
#                   rep("subsequent", nrow(evalStats[[3]])))
# temp <- do.call("rbind", evalStats)
# temp <- data.frame(treeWAS.test, temp)
# rownames(temp) <- NULL
# evalStats <- temp
#
# ## append args to other tests:
# fisher.bonf <- data.frame(args, fisher.bonf)
# fisher.fdr <- data.frame(args, fisher.fdr)
# plink.bonf <- data.frame(args, plink.bonf)
# plink.fdr <- data.frame(args, plink.fdr)
# plink.gc.bonf <- data.frame(args, plink.gc.bonf)
# plink.gc.fdr <- data.frame(args, plink.gc.fdr)
#
#
# ##########
# ## save ##    #######    #######    #######    #######    #######    #######
# ##########
# tree.type <- "ALL"
# s <- "ALL"
# af <- "ALL"
# dir.n <- "1.120"
#
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_treeWAS_evalStats.df.all.tests.Rdata", sep="")
# save(evalStats, file=filename)
#
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_fisher.bonf_evalStats.df.Rdata", sep="")
# save(fisher.bonf, file=filename)
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_fisher.fdr_evalStats.df.Rdata", sep="")
# save(fisher.fdr, file=filename)
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_plink.bonf_evalStats.df.Rdata", sep="")
# save(plink.bonf, file=filename)
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_plink.fdr_evalStats.df.Rdata", sep="")
# save(plink.fdr, file=filename)
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_plink.gc.bonf_evalStats.df.Rdata", sep="")
# save(plink.gc.bonf, file=filename)
# filename <- paste("./", set.n, "_", dir.n, "_", tree.type, "_s.", s, "_af.", af, "_plink.gc.fdr_evalStats.df.Rdata", sep="")
# save(plink.gc.fdr, file=filename)
#
# ###################
# ## BEESWARM PLOT ##   #########################   #########################   #########################
# ###################
# ## uses BOTH df and beeswarm dataframes...
#
# # library(beeswarm)
# # library(plyr)
# # library(ggplot2)
#
# #############
# ## MY DATA ##
# #############
#
# # t <- 3
# # df <- evalStats[[t]]
#
# treeWAS.tests <- c("terminal", "simultaneous", "subsequent")
#
# # t <- 1
# for(t in 1:3){
#   df <- evalStats
#   df <- df[df$treeWAS.test==treeWAS.tests[t] & df$tree.type=="coal" & df$s==0.1 & df$af==5, ]
#   beeswarm.plot(y="F1.score", x="test", df, #y.lab="Sensitivity",
#                 pt.size=4, legend=TRUE)
# }
# #
#
# ## SAVE
#
# # set3_92.100_rtree_s.0.1_af.5_treeWAS.terminal_PPV_leg
# # set3_92.100_rtree_s.0.1_af.5_treeWAS.simultaneous_PPV_leg
# # set3_92.100_rtree_s.0.1_af.5_treeWAS.subsequent_PPV_leg
# #
# # set3_92.100_rtree_s.0.1_af.5_treeWAS.terminal_sensitivity_leg
# # set3_92.100_rtree_s.0.1_af.5_treeWAS.simultaneous_sensitivity_leg
# # set3_92.100_rtree_s.0.1_af.5_treeWAS.subsequent_sensitivity_leg
# #
# # set3_92.100_rtree_s.0.1_af.5_treeWAS.terminal_FPR_leg
# # set3_92.100_rtree_s.0.1_af.5_treeWAS.simultaneous_FPR_leg
# # set3_92.100_rtree_s.0.1_af.5_treeWAS.subsequent_FPR_leg
# #
# # set3_92.100_rtree_s.0.1_af.5_treeWAS.terminal_F1.score_leg
# # set3_92.100_rtree_s.0.1_af.5_treeWAS.simultaneous_F1.score_leg
# # set3_92.100_rtree_s.0.1_af.5_treeWAS.subsequent_F1.score_leg
#
# ## REDO ALL PPV
# # set3_1.120_ALL_treeWAS.terminal_F1.score_leg
# # set3_1.120_ALL_treeWAS.simultaneous_F1.score_leg
# # set3_1.120_ALL_treeWAS.subsequent_F1.score_leg
#
# # set3_1.120_ALL_treeWAS.terminal_PPV_leg
# # set3_1.120_ALL_treeWAS.simultaneous_PPV_leg
# # set3_1.120_ALL_treeWAS.subsequent_PPV_leg
#
# # set3_1.120_ALL_treeWAS.terminal_sensitivity_leg
# # set3_1.120_ALL_treeWAS.simultaneous_sensitivity_leg
# # set3_1.120_ALL_treeWAS.subsequent_sensitivity_leg
#
# # set3_1.120_ALL_treeWAS.terminal_FPR_leg
# # set3_1.120_ALL_treeWAS.simultaneous_FPR_leg
# # set3_1.120_ALL_treeWAS.subsequent_FPR_leg
#
# ###################
# ## beeswarm.plot ##
# ###################
# beeswarm.plot <- function(y="sensitivity", x="test", df, y.lab=NULL,
#                           pt.size=4, legend=FALSE){
#
#   if(is.null(y.lab)) y.lab <- y
#
#   ## Y ~ X ??
#   fm <- as.formula(paste(y, x, sep=" ~ "))
#
#   beeswarm <- beeswarm(fm,
#                        data = df,
#                        #method="swarm", # swarm square hex center
#                        #priority="descending", ## ONLY for SWARM method...
#                        method="center", # swarm square hex center
#                        #priority="descending", ## ONLY for SWARM method...
#                        pwcol = eval(parse(text=x)),
#                        #col = myCol, ## to set w funky colours (INSTEAD of pwcol = test)
#                        ylim = c(-0.001,1), # otherwise ggplot can't plot ZERO values --> NAs
#                        las=2,
#                        cex=0.8,
#                        corral = "omit",
#                        do.plot = FALSE) # none gutter wrap omit
#   # head(beeswarm)
#
#   ######################################################
#   ## Find and Replace OUTLIERS(' symbols in plot...): ##
#   ######################################################
#   outliers <- outlier.vals <- PCH <- list()
#
#   if(!all(beeswarm$col %in% levels(df[,x]))){
#     foo <- beeswarm$col
#     foo <- levels(df[,x])[foo]
#     beeswarm$col <- factor(foo, levels=levels(df[,x]))
#   }
#
#   noms <- as.character(levels(beeswarm$col))
#
#   ## FOR LOOP ##
#   for(i in 1:length(noms)){
#     #i <- 1
#     # get vals for variable (and boxplot)
#     val <- beeswarm$y[which(beeswarm$col==noms[i])]
#     #boxplot(val, ylim=c(-0.001, 1))
#     if(length(val) == 0){
#       PCH[[i]] <- NULL
#       outliers[[i]] <- NULL
#     }else{
#       PCH[[i]] <- rep(16, length(val)) # standard filled circle...
#
#       ## get median
#       M <- as.numeric(quantile(val, 0.5))
#       # get lower 25 of box
#       Q25 <- as.numeric(quantile(val, 0.25))
#       # get upper 75 of box
#       Q75 <- as.numeric(quantile(val, 0.75))
#       # get box length
#       box <- Q75-Q25
#
#       if(box == 0) box <- 0.0000001
#
#       # with a coef of 1.5 (the default for boxplots), identify outlying values
#       outliers[[i]] <- c(which(val < Q25-(1.5*box)), which(val > Q75+(1.5*box)))
#       # get values of outliers
#       if(length(outliers[[i]]) > 0){
#         outlier.vals[[i]] <- val[outliers[[i]]]
#         PCH[[i]] <- replace(PCH[[i]], outliers[[i]], 17) # replace with triangle...
#       }else{
#         outlier.vals[[i]] <- NULL
#       }
#     }
#   } # end for loop
#
#   #outliers
#   PCH <- as.vector(unlist(PCH))
#   # PCH
#
#
#   #########################################################################################################
#   ######################
#   ## plots, layers... ##
#   ######################
#
#   if(legend == FALSE){
#
#     ################
#     ## NO LEGEND: ##
#     ################
#
#     beeswarm.plot <- ggplot(beeswarm, aes(x, y)) +
#       xlab("") +
#       guides(fill=FALSE) +
#       scale_x_discrete(drop=FALSE) +
#       scale_y_continuous(y.lab, limits=c(0,1))  # expression("char")
#
#     beeswarm.plot2 <- beeswarm.plot +
#       guides(fill=FALSE) +
#       geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.25,
#                    outlier.shape = 17,
#                    outlier.size=0) +
#       theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
#             axis.text.y = element_text(size=13),
#             axis.title.y=element_text(size=18),
#             legend.position="none")
#
#     beeswarm.plot3 <- beeswarm.plot2 +
#       geom_point(data=beeswarm, aes(colour = col), pch = PCH, size=pt.size, na.rm=TRUE, alpha=0.6) +
#       guides(fill=FALSE) +
#       theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
#             axis.text.y = element_text(size=13),
#             axis.title.y=element_text(size=18),
#             legend.position="none")
#
#     beeswarm.plot4 <- beeswarm.plot3 +
#       guides(fill=FALSE) +
#       geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.0, fatten=3,
#                    outlier.shape = 17,
#                    outlier.size=0) +
#       theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
#             axis.text.y = element_text(size=13),
#             axis.title.y=element_text(size=18),
#             legend.position="none")
#   }else{
#
#     ##################
#     ## WITH LEGEND: ##
#     ##################
#
#     beeswarm.plot <- ggplot(beeswarm, aes(x, y)) +
#       xlab("") +
#       guides(fill=FALSE) +
#       scale_x_discrete(drop=FALSE) +
#       scale_y_continuous(y.lab, limits=c(0,1)) # expression("char")
#
#     beeswarm.plot2 <- beeswarm.plot +
#       guides(fill=FALSE) +
#       geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.25,
#                    outlier.shape = 17,
#                    outlier.size=0) +
#       theme(axis.text.x = element_text(angle=65, hjust=1, vjust=0.95, size=10),
#             axis.text.y = element_text(size=13),
#             axis.title.y = element_text(size=18),
#             axis.title.x = element_blank(),
#             legend.position="none")
#
#     beeswarm.plot3 <- beeswarm.plot2 +
#       geom_point(data=beeswarm, aes(colour = col), pch = PCH, size=pt.size, na.rm=TRUE, alpha=0.6) +
#       guides(fill=FALSE) +
#       theme(axis.text.x = element_text(angle=65, hjust=1, vjust=0.95, size=10),
#             axis.text.y = element_text(size=13),
#             axis.title.y = element_text(size=18),
#             axis.title.x = element_blank(),
#             legend.position="none")
#
#     beeswarm.plot4 <- beeswarm.plot3 +
#       guides(fill=FALSE) +
#       geom_boxplot(data=df, aes(x=eval(parse(text=x)), y=eval(parse(text=y)), fill=eval(parse(text=x))), alpha=0.0, fatten=3,
#                    outlier.shape = 17,
#                    outlier.size=0) +
#       theme(axis.text.x = element_text(angle=65, hjust=1, vjust=0.95, size=10),
#             axis.text.y = element_text(size=13),
#             axis.title.y = element_text(size=18),
#             axis.title.x = element_blank(),
#             legend.position="none")
#   }
#
#
#   ## PRINT PLOT ##
#   plot(beeswarm.plot4)
#
# } # end beeswarm.plot
#
#
#
# ## CAREFUL: CHECK that your outliers are really in (ALMOST) the right place by plotting beeswarm2 w oultier.cex=2, and not outlier.size=0!
# beeswarm.plot.terminal <- beeswarm.plot4
# beeswarm.plot.simultaneous <- beeswarm.plot4
# beeswarm.plot.subsequent <- beeswarm.plot4
#
#
# # set3_1.120_ALL_treeWAS.subsequent_FPR
#
# ##########################################################################################################





















###############################
## All 3 tests together now! ##
###############################

## save as png, then...
# install.packages("png", dep=T)
library(png)

filenames <- c("./set3_21.30_coal_s.1_af.5_treeWAS.terminal_FPR.png",
               "./set3_21.30_coal_s.1_af.5_treeWAS.simultaneous_FPR.png",
               "./set3_21.30_coal_s.1_af.5_treeWAS.subsequent_FPR.png")
foo<-list()
for(j in 1:3) foo[[j]] <- readPNG(filenames[j])

layout(matrix(1:3,nr=3,byr=T))
for (j in 1:3) plot(foo[[j]])


################

#
# install.packages("gridExtra", dep=T)
# library("gridExtra")
#
# arrangeGrob(beeswarm.plot.terminal,
#           beeswarm.plot.simultaneous,
#           beeswarm.plot.subsequent,
#           labels=c("terminal", "simultaneous", "subsequent"), ncol = 3, nrow = 1)


# par(mfrow=c(3,1))
# beeswarm.plot.terminal
# beeswarm.plot.simultaneous
# beeswarm.plot.subsequent
# par(mfrow=c(1,1))




##########################################################
## BY TEST (w BEST treeWAS thresh method only x3 tests) ##
##########################################################

dat.ori <- dat

test <- c(rep("treeWAS", length(treeWAS$accuracy)),
          rep("fisher.bonf", length(fisher.bonf$accuracy)),
          rep("fisher.fdr", length(fisher.fdr$accuracy)),
          rep("plink.bonf", length(plink.bonf$accuracy)),
          rep("plink.fdr", length(plink.bonf$accuracy)),
          rep("plink.gc.bonf", length(plink.bonf$accuracy)),
          rep("plink.gc.fdr", length(plink.bonf$accuracy))
          )

# accuracy <- c(treeWAS$accuracy, fisher.bonf$accuracy, fisher.fdr$accuracy,
#               plink.bonf$accuracy, plink.fdr$accuracy, plink.gc.bonf$accuracy, plink.gc.fdr$accuracy)


sensitivity <- c(treeWAS$sensitivity, fisher.bonf$sensitivity, fisher.fdr$sensitivity,
                 plink.bonf$sensitivity, plink.fdr$sensitivity, plink.gc.bonf$sensitivity, plink.gc.fdr$sensitivity)


specificity <- c(treeWAS$specificity, fisher.bonf$specificity, fisher.fdr$specificity,
                 plink.bonf$specificity, plink.fdr$specificity, plink.gc.bonf$specificity, plink.gc.fdr$specificity)

FPR <- c(treeWAS$FPR, fisher.bonf$FPR, fisher.fdr$FPR,
         plink.bonf$FPR, plink.fdr$FPR, plink.gc.bonf$FPR, plink.gc.fdr$FPR)


## WITH NAs...
df <- data.frame(test, specificity, FPR, sensitivity)
#filename <- paste("./", set.n, "theta_p_50_combined_df.Rdata", sep="")
# filename <- paste("./", set.n, "_combined_df.Rdata", sep="")
# save(df, file=filename)

length(which(df$sensitivity[which(df$test == "treeWAS")] == 0))


##########################################################################################################

#####################################################

##########################################################################################################


## ACCURACY ##

# ## DENSITY + HISTOGRAM
# library(ggplot2)
#
# myCol <- funky(7)
#
# ## box plots:
# # A basic box with the conditions colored
# bp <- ggplot(df, aes(x=test, y=accuracy, fill=test)) + geom_boxplot()
# bp + ggtitle("Accuracy")+
#   theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
#         plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))
#
#
# # The above adds a redundant legend. With the legend removed:
# ggplot(dat, aes(x=cond, y=rating, fill=cond)) + geom_boxplot() +
#   guides(fill=FALSE)

##########################################################################################################


#####################################################

##########################################################################################################



## FPR ##

## DENSITY + HISTOGRAM??
#
# library(adegenet)
# library(ggplot2)
#
# myCol <- funky(7)

###########################
## box plots: ## ## FPR #########################      #########################
###########################

# A basic box with the conditions colored
bp <- ggplot(df, aes(x=test, y=FPR, fill=test)) + geom_boxplot()+scale_y_continuous(limits=c(0,1))
bp + ggtitle("FPR")+
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))


# The above adds a redundant legend. With the legend removed:
ggplot(df, aes(x=test, y=FPR, fill=test, color=factor(test))) + geom_boxplot() + scale_y_continuous(limits=c(0,1)) +
#   ggtitle("FPR") +
  guides(fill=FALSE)+
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        plot.title = element_text(size = 20, face = "bold", vjust =1),
        legend.position="none") # colour = "black",

## set 1 FPR:
## plink.fdr has one FPR dot at ~ 0.3
## fisher.fdr, plink.fdr, and plink.gc.fdr all have one FPR dot at ~ 0.02 ... which one?

which.max(df$FPR[which(df$test=="fisher.fdr")]) # 82
which.max(df$FPR[which(df$test=="plink.fdr")]) # 69
which(df$FPR[which(df$test=="plink.fdr")] == sort(df$FPR[which(df$test=="plink.fdr")], decreasing=TRUE)[2]) # 82
which.max(df$FPR[which(df$test=="plink.gc.fdr")]) # 82

args82 <- get(load("./set1_82_args.Rdata"))
args82


###################
## BEESWARM PLOT ##   #########################   #########################   #########################
###################
## uses BOTH df and beeswarm dataframes...

# library(beeswarm)
# library(plyr)

#############
## MY DATA ##
#############
beeswarm <- beeswarm(FPR ~ test,
                     data = df,
                     #method="swarm", # swarm square hex center
                     #priority="descending", ## ONLY for SWARM method...
                     method="center", # swarm square hex center
                     #priority="descending", ## ONLY for SWARM method...
                     pwcol = test,
                     #col = myCol, ## to set w funky colours (INSTEAD of pwcol = test)
                     ylim = c(-0.001,1), # otherwise ggplot can't plot ZERO values --> NAs
                     las=2,
                     cex=0.8,
                     corral = "omit") # none gutter wrap omit
head(beeswarm)

######################################################
## Find and Replace OUTLIERS(' symbols in plot...): ##
######################################################
outliers <- outlier.vals <- PCH <- list()
noms <- as.character(levels(beeswarm$col))

## FOR LOOP ##
for(i in 1:length(noms)){
  #i <- 1
  # get vals for variable (and boxplot)
  val <- beeswarm$y[which(beeswarm$col==noms[i])]
  #boxplot(val, ylim=c(-0.001, 1))
  PCH[[i]] <- rep(16, length(val)) # standard filled circle...

  ## get median
  M <- as.numeric(quantile(val, 0.5))
  # get lower 25 of box
  Q25 <- as.numeric(quantile(val, 0.25))
  # get upper 75 of box
  Q75 <- as.numeric(quantile(val, 0.75))
  # get box length
  box <- Q75-Q25

  if(box == 0) box <- 0.0000001

  # with a coef of 1.5 (the default for boxplots), identify outlying values
  outliers[[i]] <- c(which(val < Q25-(1.5*box)), which(val > Q75+(1.5*box)))
  # get values of outliers
  if(length(outliers[[i]]) > 0){
    outlier.vals[[i]] <- val[outliers[[i]]]
    PCH[[i]] <- replace(PCH[[i]], outliers[[i]], 17) # replace with triangle...
  }else{
    outlier.vals[[i]] <- NULL
  }
} # end for loop

#outliers
PCH <- as.vector(unlist(PCH))
PCH


#########################################################################################################
######################
## plots, layers... ##
######################
beeswarm.plot <- ggplot(beeswarm, aes(x, y)) +
  xlab("") +
  guides(fill=FALSE) +
  scale_y_continuous(expression("False Positive Rate (FPR)"), limits=c(0,1))

beeswarm.plot2 <- beeswarm.plot +
  guides(fill=FALSE) +
  geom_boxplot(data=df, aes(x=test, y=FPR, fill=test), alpha=0.25,
               outlier.shape = 17,
               outlier.size=0,
               #outlier.cex=2,
               ylim=c(-0.001, 1)) +# , varwidth=FALSE, outlier.shape = 17
  theme(axis.text.x = element_text(angle=45, hjust=0.9, vjust=0.85, size=18),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=18),
        legend.position="none")

beeswarm.plot3 <- beeswarm.plot2 +
  geom_point(data=beeswarm, aes(colour = col), pch = PCH, size=5, na.rm=TRUE, alpha=0.6) +
  #scale_colour_manual(values = myCol) +
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle=45, hjust=0.9, vjust=0.85, size=18),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=18),
        legend.position="none")

beeswarm.plot4 <- beeswarm.plot3 +
  guides(fill=FALSE) +
  geom_boxplot(data=df, aes(x=test, y=FPR, fill=test), alpha=0.0, fatten=3,
               outlier.shape = 17,
               outlier.size=0,
               #outlier.cex=2,
               ylim=c(-0.001, 1)) +# , varwidth=FALSE, outlier.shape = 17
  theme(axis.text.x = element_text(angle=45, hjust=0.9, vjust=0.85, size=18),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=18),
        legend.position="none")

plot(beeswarm.plot4)
## CAREFUL: CHECK that your outliers are really in (ALMOST) the right place by plotting beeswarm2 w oultier.cex=2, and not outlier.size=0!



##########################################################################################################






#####################################
## box plots: ## ## SENSITIVITY    #########################      #########################
#####################################

# A basic box with the conditions colored
bp <- ggplot(df, aes(x=test, y=sensitivity, fill=test)) + geom_boxplot()+scale_y_continuous(limits=c(0,1))
bp + ggtitle("Sensitivity")+
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))


# The above adds a redundant legend. With the legend removed:
ggplot(df, aes(x=test, y=sensitivity, fill=test)) + geom_boxplot() + scale_y_continuous(limits=c(0,1)) +
  #   ggtitle("Sensitivity") +
  guides(fill=FALSE)+
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))

# which.min(df$sensitivity[which(df$test == "treeWAS")]) # 13
# length(which(df$sensitivity[which(df$test == "treeWAS")] == 0)) # 19
# which(df$sensitivity[which(df$test == "treeWAS")] == 0)
## 13 18 23 24 28 32 33 38 40 43 44 56 60 76 79 87 91 93 94


###################
## BEESWARM PLOT ##   #########################   #########################   #########################
###################
## uses BOTH df and beeswarm dataframes...

# library(beeswarm)
# library(plyr)

#############
## MY DATA ##
#############
beeswarm <- beeswarm(sensitivity ~ test,
                     data = df,
                     #method="swarm", # swarm square hex center
                     #priority="descending", ## ONLY for SWARM method...
                     method="center", # swarm square hex center
                     #priority="descending", ## ONLY for SWARM method...
                     pwcol = test,
                     #col = myCol, ## to set w funky colours (INSTEAD of pwcol = test)
                     ylim = c(-0.001, 1), # otherwise ggplot can't plot ZERO values --> NAs
                     las=2,
                     cex=0.8,
                     corral = "omit") # none gutter wrap omit
head(beeswarm)

######################################################
## Find and Replace OUTLIERS(' symbols in plot...): ##
######################################################
outliers <- outlier.vals <- PCH <- list()
noms <- as.character(levels(beeswarm$col))

## FOR LOOP ##
for(i in 1:length(noms)){
  #i <- 1
  # get vals for variable (and boxplot)
  val <- beeswarm$y[which(beeswarm$col==noms[i])]
  #boxplot(val, ylim=c(-0.001, 1))
  PCH[[i]] <- rep(16, length(val)) # standard filled circle...

  ## get median
  M <- as.numeric(quantile(val, 0.5))
  # get lower 25 of box
  Q25 <- as.numeric(quantile(val, 0.25))
  # get upper 75 of box
  Q75 <- as.numeric(quantile(val, 0.75))
  # get box length
  box <- Q75-Q25

  if(box == 0) box <- 0.0000001

  # with a coef of 1.5 (the default for boxplots), identify outlying values
  outliers[[i]] <- c(which(val < Q25-(1.5*box)), which(val > Q75+(1.5*box)))
  # get values of outliers
  if(length(outliers[[i]]) > 0){
    outlier.vals[[i]] <- val[outliers[[i]]]
    PCH[[i]] <- replace(PCH[[i]], outliers[[i]], 17) # replace with triangle...
  }else{
    outlier.vals[[i]] <- NULL
  }
} # end for loop

#outliers
PCH <- as.vector(unlist(PCH))
PCH
#PCH <- replace(PCH, which(PCH==1), 16)

######################
## plots, layers... ##
######################
beeswarm.plot <- ggplot(beeswarm, aes(x, y)) +
  xlab("") +
  guides(fill=FALSE) +
  scale_y_continuous(expression("Sensitivity"), limits=c(0,1.001))

beeswarm.plot2 <- beeswarm.plot +
  guides(fill=FALSE) +
  geom_boxplot(data=df, aes(x=test, y=sensitivity, fill=test), alpha=0.25,
               outlier.shape = 17,
               outlier.size=0,
               #outlier.cex=2,
               ylim=c(-0.001, 1)) +# , varwidth=FALSE, outlier.shape = 17
  theme(axis.text.x = element_text(angle=45, hjust=0.9, vjust=0.85, size=18),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=18),
        legend.position="none")

beeswarm.plot3 <- beeswarm.plot2 +
  geom_point(data=beeswarm, aes(colour = col), pch = PCH, size=5, na.rm=TRUE, alpha=0.6) +
  #scale_colour_manual(values = myCol) +
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle=45, hjust=0.9, vjust=0.85, size=18),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=18),
        legend.position="none")

beeswarm.plot4 <- beeswarm.plot3 +
  guides(fill=FALSE) +
  geom_boxplot(data=df, aes(x=test, y=sensitivity, fill=test), alpha=0.0, fatten=3,
               outlier.shape = 17,
               outlier.size=0,
               #outlier.cex=2,
               ylim=c(-0.001, 1)) +# , varwidth=FALSE, outlier.shape = 17
  theme(axis.text.x = element_text(angle=45, hjust=0.9, vjust=0.85, size=18),
        axis.text.y = element_text(size=13),
        axis.title=element_text(size=18),
        legend.position="none")

plot(beeswarm.plot4)


##########################################################################################################
# beeswarm.plot <- ggplot(beeswarm, aes(x, y)) +
#   xlab("") +
#   guides(fill=FALSE) +
#   scale_y_continuous(expression("Sensitivity"), limits=c(-0.001,1))
# # +
# #   scale_x_continuous(labels=as.character(unique(beeswarm$x.orig)),
# #                      breaks=(c(1:length(unique(beeswarm$x.orig)))))
#
#
# beeswarm.plot2 <- beeswarm.plot +
#   guides(fill=FALSE) +
#   geom_boxplot(data=df, aes(x=test, y=sensitivity, fill=test), alpha=0.25,
#                outlier.shape = 17,
#                outlier.size=0,
#                #outlier.cex=2,
#                ylim=c(-0.001, 1)) +# , varwidth=FALSE, outlier.shape = 17
#   theme(axis.text.x = element_text(angle=45, hjust=0.9, vjust=0.85, size=13),
#         axis.text.y = element_text(size=11),
#         axis.title=element_text(size=14),
#         legend.position="none")
# #ylim = c(-0.001,1)
#
# beeswarm.plot3 <- beeswarm.plot2 +
#   geom_point(data=beeswarm, aes(colour = col), pch = PCH, cex=5, alpha=0.7, na.rm=TRUE) +
#   #scale_colour_manual(values = myCol) +
#   guides(fill=FALSE) +
#   theme(axis.text.x = element_text(angle=45, hjust=0.9, vjust=0.85, size=18),
#         axis.text.y = element_text(size=13),
#         axis.title=element_text(size=18),
#         legend.position="none")
#
# plot(beeswarm.plot3)
# ## CAREFUL: CHECK that your outliers are really in (ALMOST) the right place by plotting beeswarm2 w oultier.cex=2, and not outlier.size=0!
#


##########################################################################################################

#
#
# ######################
# ## plots, layers... ##
# ######################
# beeswarm.plot <- ggplot(beeswarm, aes(x, y)) +
#   xlab("") +
#   guides(fill=FALSE) +
#   scale_y_continuous(expression("False Positive Rate (FPR)"), limits=c(0,1)) +
#   scale_x_continuous(labels=as.character(unique(beeswarm$x.orig)),
#                      breaks=(c(1:7)))
#
#
# beeswarm.plot2 <- beeswarm.plot +
#   guides(fill=FALSE) +
#   geom_boxplot(aes(x, y, group = round_any(x, 1, round))) # , varwidth=FALSE, outlier.shape = 17
#
# beeswarm.plot3 <- beeswarm.plot2 +
#   geom_point(aes(colour = col, cex = 2), pch = PCH, na.rm=TRUE) +
#   #scale_colour_manual(values = myCol) +
#   guides(fill=FALSE) +
#   theme(axis.text.x = element_text(angle=45, hjust=0, vjust=0, size=12), legend.position="none",
#         plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))
# #   theme(axis.text.x = element_text(angle=90, hjust=0, vjust=1, size=12), legend.position="none",
# #       plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))
#
# plot(beeswarm.plot3)
# #######################################################
# bp <- ggplot(df, aes(x=test, y=FPR)) + # , fill=test
#   geom_boxplot() + scale_y_continuous(limits=c(0,1)) +
#   #   ggtitle("False Positive Rate (FPR)") +
#   guides(fill=FALSE)+
#   theme(axis.text.x = element_text(angle=45, hjust=0, vjust=0, size=12), legend.position="none",
#         plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))
#
# plot(bp)
# points(beeswarm$x, beeswarm$y, type="p", pch=PCH, add=TRUE)
#
# #####################################################
# # The above adds a redundant legend. With the legend removed:
# p1 <- ggplot(df, aes(x=test, y=FPR, fill=test)) + geom_boxplot() + scale_y_continuous(limits=c(0,1)) +
#         #   ggtitle("False Positive Rate (FPR)") +
#         guides(fill=FALSE)+
#         theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
#         plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))
#
# plot(p1)
# str(p1)

##########################################################################################################
## RESET YOUR OPTIONS AS:
#options(warn=warn.ori, error=error.ori)

# ##############
# beeswarm.plot3 <- beeswarm.plot2 + geom_point(aes(colour = col)) +
#   scale_colour_manual(values = c("black", "red")) +
#   scale_x_continuous(breaks = c(1:2),
#                      labels = c("Censored", "Metastasis"), expand = c(0, 0.5))
# plot(beeswarm.plot3)
#####################################################
#####################################################
## plots, layers...
beeswarm.plot <- ggplot(beeswarm, aes(x, y)) +
  xlab("") +
  guides(fill=FALSE) +
  scale_y_continuous(expression("False Positive Rate (FPR)"), limits=c(0,1)) +
  scale_x_discrete(#breaks= seq(0, 1, length.out=9)[2:8],  limits=c(0,1), #c(1:7),
    labels = as.character(unique(beeswarm$col)), las=2) #, expand=c(0,0.5))

beeswarm.plot2 <- beeswarm.plot +
  guides(fill=FALSE) +
  geom_boxplot(aes(x, y, group = round_any(x, 1, round)), outlier.shape = NA) #  #??

beeswarm.plot3 <- beeswarm.plot2 +
  geom_point(aes(colour = col, cex=2)) +
  #scale_colour_manual(values = myCol) +
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle=90, hjust=0, vjust=1),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1),
        legend.position="none")

plot(beeswarm.plot3)

#####################################################

##########################################################################################################



## SPECIFICITY ##

## DENSITY + HISTOGRAM
library(ggplot2)

myCol <- funky(7)

## box plots:
# A basic box with the conditions colored
bp <- ggplot(df, aes(x=test, y=specificity, fill=test)) + geom_boxplot()
bp + ggtitle("Specificity")+
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
  plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))


# The above adds a redundant legend. With the legend removed:
ggplot(df, aes(x=test, y=specificity, fill=test)) + geom_boxplot() +
  ggtitle("Specificity") +
  guides(fill=FALSE)+
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))

###########################
## OR with individual data frames (NAs removed --> diff nrow per test problem... )
bp <- ggplot(df.specificity, aes(x=test, y=dat.specificity, fill=test)) + geom_boxplot()
bp + ggtitle("Specificity")+
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))


##########################################################################################################

#####################################################

##########################################################################################################



## FPR ##

## DENSITY + HISTOGRAM
library(ggplot2)
library(adegenet)

myCol <- funky(7)

## box plots:
# A basic box with the conditions colored
bp <- ggplot(df, aes(x=test, y=FPR, fill=test)) + geom_boxplot()
bp + ggtitle("False Positive Rate (FPR)")+
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))


# The above adds a redundant legend. With the legend removed:
ggplot(df, aes(x=test, y=FPR, fill=test)) + geom_boxplot() +  scale_y_continuous(limits=c(0,1)) +
#   ggtitle("False Positive Rate (FPR)") +
  guides(fill=FALSE)+
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))


## all together now:

## can't get all 3 on same plot...
#ggplot(df, aes(x=FPR, fill=test)) + geom_density(alpha=0.3)
## works for fishers only:
#ggplot(df[length(treeWAS$FPR):nrow(df),], aes(x=FPR, fill=test)) + geom_density(alpha=0.3)


##########################################################################################################

#####################################################

##########################################################################################################
#' @











# <<echo=FALSE>>=

  ###############
#### SET 2 ####
###############

#######################################################################################################################

#################
## theta_p = 5 ##
#################

foo <- dir("C:/Cait 2012/Work/Xavier/Sims/set2/theta_p = 25")
foo
## get all performance Rdata names
toKeep <- grep("performance", foo) ##??
foo <- foo[toKeep]

## load performance data
dat <- list()
setwd("C:/Cait 2012/Work/Xavier/Sims/set2/theta_p = 25")
for(i in 1:length(foo)){
  dat[[i]] <- get(load(paste("./", foo[i], sep="")))
}

treeWAS <- list()
for(i in 1:length(dat)){
  treeWAS[[i]] <- dat[[i]][[1]]
}

fisher.bonf <- list()
for(i in 1:length(dat)){
  fisher.bonf[[i]] <- dat[[i]][[2]]
}

fisher.fdr <- list()
for(i in 1:length(dat)){
  fisher.fdr[[i]] <- dat[[i]][[3]]
}

## combine
treeWAS <- do.call("rbind", treeWAS)
fisher.bonf <- do.call("rbind", fisher.bonf)
fisher.fdr <- do.call("rbind", fisher.fdr)

## save
save(treeWAS, file="./set2_theta_p_25_treeWAS_performance.Rdata")
save(fisher.bonf, file="./set2_theta_p_25_fisher.bonf_performance.Rdata")
save(fisher.fdr, file="./set2_theta_p_25_fisher.fdr_performance.Rdata")

## summarise
summary(treeWAS)
summary(fisher.bonf)
summary(fisher.fdr)

##########
## plot ##
##########

## GET DATA
dat.ori <- dat
test <- c(rep("treeWAS", length(treeWAS$accuracy)),
          rep("fisher.bonf", length(fisher.bonf$accuracy)),
          rep("fisher.fdr", length(fisher.fdr$accuracy)))
accuracy <- c(treeWAS$accuracy, fisher.bonf$accuracy, fisher.fdr$accuracy)
specificity <- c(treeWAS$specificity, fisher.bonf$specificity, fisher.fdr$specificity)
FPR <- c(treeWAS$FPR, fisher.bonf$FPR, fisher.fdr$FPR)

df <- data.frame(test, accuracy, specificity, FPR)
save(df, file="./set2_theta_p_25_combined_df.Rdata")



##########################################################################################################

#####################################################

##########################################################################################################


## ACCURACY ##

## DENSITY + HISTOGRAM
library(ggplot2)

myCol <- funky(5)

# Histogram overlaid with kernel density curve
p <- ggplot(df[which(test=="treeWAS"),], aes(x=accuracy)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.00015,
                 colour="black", fill="white", xlim=c(0.9995, 1.001)) +
  geom_density(alpha=.5, fill=myCol[1], xlim=c(0.9995, 1.001))  # Overlay with transparent density plot

p + ggtitle(expression(atop("Accuracy: treeWAS", atop(italic("(theta_p = 25)"), "")))) +
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        #plot.margin = unit(c(1.5, 1, 1, 1), "cm"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust = -1))


# Histogram overlaid with kernel density curve
p2 <- ggplot(df[which(test=="fisher.bonf"),], aes(x=accuracy)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.05,
                 colour="black", fill="white") +
  geom_density(alpha=.5, fill=myCol[3])  # Overlay with transparent density plot

p2 + ggtitle(expression(atop("Accuracy: Fisher's Exact Test with Bonferonni Correction"), atop(italic("(theta_p = 5)"), ""))) +
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        #plot.margin = unit(c(1.5, 1, 1, 1), "cm"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust = -1))


p3 <- ggplot(df[which(test=="fisher.fdr"),], aes(x=accuracy)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.05,
                 colour="black", fill="white") +
  geom_density(alpha=.5, fill=myCol[2])  # Overlay with transparent density plot

p3 + ggtitle(expression(atop("Accuracy: Fisher's Exact Test with FDR Correction"), atop(italic("(theta_p = 5)"), ""))) +
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        #plot.margin = unit(c(1.5, 1, 1, 1), "cm"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust = -1))


## box plots:
# A basic box with the conditions colored
bp <- ggplot(df, aes(x=test, y=accuracy, fill=test)) + geom_boxplot()
bp + ggtitle("Accuracy")+
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))


# The above adds a redundant legend. With the legend removed:
ggplot(dat, aes(x=cond, y=rating, fill=cond)) + geom_boxplot() +
  guides(fill=FALSE)


## all together now:

## can't get all 3 on same plot...
#ggplot(df, aes(x=accuracy, fill=test)) + geom_density(alpha=0.3)
## works for fishers only:
#ggplot(df[length(treeWAS$accuracy):nrow(df),], aes(x=accuracy, fill=test)) + geom_density(alpha=0.3)


##########################################################################################################

#####################################################

##########################################################################################################



# ## SPECIFICITY ##
#
# ## DENSITY + HISTOGRAM
# library(ggplot2)
#
# myCol <- funky(5)
#
# # Histogram overlaid with kernel density curve
# p <- ggplot(df[which(test=="treeWAS"),], aes(x=accuracy)) +
#     geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
#                    binwidth=.00015,
#                    colour="black", fill="white", xlim=c(0.9995, 1.001)) +
#     geom_density(alpha=.5, fill=myCol[1], xlim=c(0.9995, 1.001))  # Overlay with transparent density plot
#
# p + ggtitle(expression(atop("Accuracy: treeWAS", atop(italic("(theta_p = 5)"), "")))) +
#    theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
#    #plot.margin = unit(c(1.5, 1, 1, 1), "cm"),
#    plot.title = element_text(size = 20, face = "bold", colour = "black", vjust = -1))
#
#
# # Histogram overlaid with kernel density curve
# p2 <- ggplot(df[which(test=="fisher.bonf"),], aes(x=accuracy)) +
#     geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
#                    binwidth=.05,
#                    colour="black", fill="white") +
#     geom_density(alpha=.5, fill=myCol[3])  # Overlay with transparent density plot
#
# p2 + ggtitle(expression(atop("Accuracy: Fisher's Exact Test with Bonferonni Correction"), atop(italic("(theta_p = 5)"), ""))) +
#    theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
#    #plot.margin = unit(c(1.5, 1, 1, 1), "cm"),
#    plot.title = element_text(size = 20, face = "bold", colour = "black", vjust = -1))
#
#
# p3 <- ggplot(df[which(test=="fisher.fdr"),], aes(x=accuracy)) +
#     geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
#                    binwidth=.05,
#                    colour="black", fill="white") +
#     geom_density(alpha=.5, fill=myCol[2])  # Overlay with transparent density plot
#
# p3 + ggtitle(expression(atop("Accuracy: Fisher's Exact Test with FDR Correction"), atop(italic("(theta_p = 5)"), ""))) +
#    theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
#    #plot.margin = unit(c(1.5, 1, 1, 1), "cm"),
#    plot.title = element_text(size = 20, face = "bold", colour = "black", vjust = -1))
#
#
# ## box plots:
# # A basic box with the conditions colored
# bp <- ggplot(df, aes(x=test, y=accuracy, fill=test)) + geom_boxplot()
# bp + ggtitle("Accuracy")+
#   theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
#   plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))
#
#
# # The above adds a redundant legend. With the legend removed:
# ggplot(dat, aes(x=cond, y=rating, fill=cond)) + geom_boxplot() +
#     guides(fill=FALSE)
#
#
# ## all together now:
#
# ## can't get all 3 on same plot...
# #ggplot(df, aes(x=accuracy, fill=test)) + geom_density(alpha=0.3)
# ## works for fishers only:
# #ggplot(df[length(treeWAS$accuracy):nrow(df),], aes(x=accuracy, fill=test)) + geom_density(alpha=0.3)


##########################################################################################################

#####################################################

##########################################################################################################



## FPR ##

## DENSITY + HISTOGRAM
library(ggplot2)

myCol <- funky(5)

# Histogram overlaid with kernel density curve
p <- ggplot(df[which(test=="treeWAS"),], aes(x=FPR)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.00015,
                 colour="black", fill="white", xlim=c(0.9995, 1.001)) +
  geom_density(alpha=.5, fill=myCol[1], xlim=c(0.9995, 1.001))  # Overlay with transparent density plot

p + ggtitle(expression(atop("FPR: treeWAS", atop(italic("(theta_p = 25)"), "")))) +
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        #plot.margin = unit(c(1.5, 1, 1, 1), "cm"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust = -1))


# Histogram overlaid with kernel density curve
p2 <- ggplot(df[which(test=="fisher.bonf"),], aes(x=FPR)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.05,
                 colour="black", fill="white") +
  geom_density(alpha=.5, fill=myCol[3])  # Overlay with transparent density plot

p2 + ggtitle(expression(atop("FPR: Fisher's Exact Test with Bonferonni Correction"), atop(italic("(theta_p = 5)"), ""))) +
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        #plot.margin = unit(c(1.5, 1, 1, 1), "cm"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust = -1))


p3 <- ggplot(df[which(test=="fisher.fdr"),], aes(x=FPR)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.05,
                 colour="black", fill="white") +
  geom_density(alpha=.5, fill=myCol[2])  # Overlay with transparent density plot

p3 + ggtitle(expression(atop("FPR: Fisher's Exact Test with FDR Correction"), atop(italic("(theta_p = 5)"), ""))) +
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        #plot.margin = unit(c(1.5, 1, 1, 1), "cm"),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust = -1))


## box plots:
# A basic box with the conditions colored
bp <- ggplot(df, aes(x=test, y=FPR, fill=test)) + geom_boxplot()
bp + ggtitle("False Positive Rate (FPR)")+
  theme(axis.text.x = element_text(angle=0, hjust=0, vjust=1),
        plot.title = element_text(size = 20, face = "bold", colour = "black", vjust =1))


# The above adds a redundant legend. With the legend removed:
ggplot(dat, aes(x=cond, y=rating, fill=cond)) + geom_boxplot() +
  guides(fill=FALSE)


## all together now:

## can't get all 3 on same plot...
#ggplot(df, aes(x=FPR, fill=test)) + geom_density(alpha=0.3)
## works for fishers only:
#ggplot(df[length(treeWAS$FPR):nrow(df),], aes(x=FPR, fill=test)) + geom_density(alpha=0.3)


##########################################################################################################

#####################################################

##########################################################################################################

## histogram
hist(treeWAS$accuracy, xlab="Accuracy", main="ACCURACY\n(treeWAS)", sub="theta_p = 5", col.sub="red", font.sub=2)
hist(treeWAS$accuracy, xlim=c(0,1), xlab="Accuracy", main="ACCURACY\n(treeWAS)", sub="theta_p = 5", col.sub="red", font.sub=2)
hist(fisher.bonf$accuracy, xlim=c(0,1), xlab="Accuracy",
     main="ACCURACY\n(Fisher's Exact Test w/ Bonferonni Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)
hist(fisher.fdr$accuracy, xlim=c(0,1), xlab="Accuracy",
     main="ACCURACY\n(Fisher's Exact Test w/ FDR Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)

## density dist
plot(density(treeWAS$accuracy), xlab="Accuracy", main="ACCURACY\n(treeWAS)", sub="theta_p = 5", col.sub="red", font.sub=2)
plot(density(fisher.bonf$accuracy), xlab="Accuracy",
     main="ACCURACY\n(Fisher's Exact Test w/ Bonferonni Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)
plot(density(fisher.fdr$accuracy), xlab="Accuracy",
     main="ACCURACY\n(Fisher's Exact Test w/ FDR Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)


# ## DENSITY + HISTOGRAM
# myCol <- funky(5)
#
# hist(treeWAS$accuracy, xlab="Accuracy",
#      col=transp(myCol[1], 0.05),
#      main="ACCURACY\n(treeWAS)", sub="theta_p = 5", col.sub="red", font.sub=2)
# polygon(density(treeWAS$accuracy), col=transp(myCol[1], .3), border=myCol[1])
#
# hist(fisher.bonf$accuracy, xlim=c(0,1), xlab="Accuracy",
#      col=transp(myCol[2], 0.05),
#      main="ACCURACY\n(Fisher's Exact Test w/ Bonferonni Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)
# polygon(density(fisher.bonf$accuracy), col=transp(myCol[2], .3), border=myCol[2])
#
# hist(fisher.fdr$accuracy, xlim=c(0,1), xlab="Accuracy",
#      main="ACCURACY\n(Fisher's Exact Test w/ FDR Correction)", sub="theta_p = 5", col.sub=myCol[3], font.sub=2)
#
# plot(density(fisher.fdr$accuracy), xlab="Accuracy",
#      main="ACCURACY\n(Fisher's Exact Test w/ FDR Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)


## SPECIFICITY ##

## histogram
hist(treeWAS$specificity, xlab="specificity", main="SPECIFICITY\n(treeWAS)", sub="theta_p = 5", col.sub="red", font.sub=2)
hist(treeWAS$specificity, xlim=c(0,1), xlab="specificity", main="SPECIFICITY\n(treeWAS)", sub="theta_p = 5", col.sub="red", font.sub=2)
hist(fisher.bonf$specificity, xlim=c(0,1), xlab="specificity",
     main="SPECIFICITY\n(Fisher's Exact Test w/ Bonferonni Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)
hist(fisher.fdr$specificity, xlim=c(0,1), xlab="specificity",
     main="SPECIFICITY\n(Fisher's Exact Test w/ FDR Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)

## density dist
plot(density(treeWAS$specificity), xlab="specificity", main="SPECIFICITY\n(treeWAS)", sub="theta_p = 5", col.sub="red", font.sub=2)
plot(density(fisher.bonf$specificity), xlab="specificity",
     main="SPECIFICITY\n(Fisher's Exact Test w/ Bonferonni Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)
plot(density(fisher.fdr$specificity), xlab="specificity",
     main="SPECIFICITY\n(Fisher's Exact Test w/ FDR Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)


## FPR ##

## histogram
hist(treeWAS$FPR, xlab="FPR", main="FALSE POSITIVE RATE\n(treeWAS)", sub="theta_p = 5", col.sub="red", font.sub=2)
hist(treeWAS$FPR, xlim=c(0,1), xlab="FPR", main="FALSE POSITIVE RATE\n(treeWAS)", sub="theta_p = 5", col.sub="red", font.sub=2)
hist(fisher.bonf$FPR, xlim=c(0,1), xlab="FPR",
     main="FALSE POSITIVE RATE\n(Fisher's Exact Test w/ Bonferonni Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)
hist(fisher.fdr$FPR, xlim=c(0,1), xlab="FPR",
     main="FALSE POSITIVE RATE\n(Fisher's Exact Test w/ FDR Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)

## density dist
plot(density(treeWAS$FPR), xlab="FPR", main="FALSE POSITIVE RATE\n(treeWAS)", sub="theta_p = 5", col.sub="red", font.sub=2)
plot(density(fisher.bonf$FPR), xlab="FPR",
     main="FALSE POSITIVE RATE\n(Fisher's Exact Test w/ Bonferonni Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)
plot(density(fisher.fdr$FPR), xlab="FPR",
     main="FALSE POSITIVE RATE\n(Fisher's Exact Test w/ FDR Correction)", sub="theta_p = 5", col.sub="red", font.sub=2)







#######################################################################################################################

#################
## theta_p = 25 ##
#################

foo <- dir("C:/Cait 2012/Work/Xavier/Sims/set2/theta_p = 25")
foo
## get all performance Rdata names
toKeep <- grep("performance", foo) ##??
foo <- foo[toKeep]

## load performance data
dat <- list()
setwd("C:/Cait 2012/Work/Xavier/Sims/set2/theta_p = 25")
for(i in 1:length(foo)){
  dat[[i]] <- get(load(paste("./", foo[i], sep="")))
}

treeWAS <- list()
for(i in 1:length(dat)){
  treeWAS[[i]] <- dat[[i]][[1]]
}

fisher.bonf <- list()
for(i in 1:length(dat)){
  fisher.bonf[[i]] <- dat[[i]][[2]]
}

fisher.fdr <- list()
for(i in 1:length(dat)){
  fisher.fdr[[i]] <- dat[[i]][[3]]
}

## combine
treeWAS <- do.call("rbind", treeWAS)
fisher.bonf <- do.call("rbind", fisher.bonf)
fisher.fdr <- do.call("rbind", fisher.fdr)

## save
save(treeWAS, file="./set2_theta_p_25_treeWAS_performance.Rdata")
save(fisher.bonf, file="./set2_theta_p_25_fisher.bonf_performance.Rdata")
save(fisher.fdr, file="./set2_theta_p_25_fisher.fdr_performance.Rdata")

## summarise
summary(treeWAS)
summary(fisher.bonf)
summary(fisher.fdr)

##########
## plot ##
##########

## ACCURACY ##

## histogram
hist(treeWAS$accuracy)
hist(treeWAS$accuracy, xlim=c(0,1))
hist(fisher.bonf$accuracy, xlim=c(0,1))
hist(fisher.fdr$accuracy, xlim=c(0,1))

## density dist
plot(density(treeWAS$accuracy), main="ACCURACY\n(treeWAS)")
plot(density(fisher.bonf$accuracy), main="ACCURACY\n(Fisher's Exact Test w/ Bonferonni Correction)")
plot(density(fisher.fdr$accuracy), main="ACCURACY\n(Fisher's Exact Test w/ FDR Correction)")


## SPECIFICITY ##

## histogram
hist(treeWAS$specificity)
hist(treeWAS$specificity, xlim=c(0,1))
hist(fisher.bonf$specificity, xlim=c(0,1))
hist(fisher.fdr$specificity, xlim=c(0,1))

## density dist
plot(density(treeWAS$specificity), main="SPECIFICITY\n(treeWAS)")
plot(density(fisher.bonf$specificity), main="SPECIFICITY\n(Fisher's Exact Test w/ Bonferonni Correction)")
plot(density(fisher.fdr$specificity), main="SPECIFICITY\n(Fisher's Exact Test w/ FDR Correction)")


## FPR ##

## histogram
hist(treeWAS$FPR)
hist(treeWAS$FPR, xlim=c(0,1))
hist(fisher.bonf$FPR, xlim=c(0,1))
hist(fisher.fdr$FPR, xlim=c(0,1))

## density dist
plot(density(treeWAS$FPR), main="FALSE POSITIVE RATE\n(treeWAS)")
plot(density(fisher.bonf$FPR), main="FALSE POSITIVE RATE\n(Fisher's Exact Test w/ Bonferonni Correction)")
plot(density(fisher.fdr$FPR), main="FALSE POSITIVE RATE\n(Fisher's Exact Test w/ FDR Correction)")







#######################################################################################################################

#################
## theta_p = 5 ##
#################

foo <- dir("C:/Cait 2012/Work/Xavier/Sims/set2/theta_p = 5")
foo
## get all performance Rdata names
toKeep <- grep("performance", foo) ##??
foo <- foo[toKeep]

## load performance data
dat <- list()
setwd("C:/Cait 2012/Work/Xavier/Sims/set2/theta_p = 5")
for(i in 1:length(foo)){
  dat[[i]] <- get(load(paste("./", foo[i], sep="")))
}

treeWAS <- list()
for(i in 1:length(dat)){
  treeWAS[[i]] <- dat[[i]][[1]]
}

fisher.bonf <- list()
for(i in 1:length(dat)){
  fisher.bonf[[i]] <- dat[[i]][[2]]
}

fisher.fdr <- list()
for(i in 1:length(dat)){
  fisher.fdr[[i]] <- dat[[i]][[3]]
}

## combine
treeWAS <- do.call("rbind", treeWAS)
fisher.bonf <- do.call("rbind", fisher.bonf)
fisher.fdr <- do.call("rbind", fisher.fdr)

## save
save(treeWAS, file="./set2_theta_p_5_treeWAS_performance.Rdata")
save(fisher.bonf, file="./set2_theta_p_5_fisher.bonf_performance.Rdata")
save(fisher.fdr, file="./set2_theta_p_5_fisher.fdr_performance.Rdata")

## summarise
summary(treeWAS)
summary(fisher.bonf)
summary(fisher.fdr)

##########
## plot ##
##########

## ACCURACY ##

## histogram
hist(treeWAS$accuracy)
hist(treeWAS$accuracy, xlim=c(0,1))
hist(fisher.bonf$accuracy, xlim=c(0,1))
hist(fisher.fdr$accuracy, xlim=c(0,1))

## density dist
plot(density(treeWAS$accuracy), main="ACCURACY\n(treeWAS)")
plot(density(fisher.bonf$accuracy), main="ACCURACY\n(Fisher's Exact Test w/ Bonferonni Correction)")
plot(density(fisher.fdr$accuracy), main="ACCURACY\n(Fisher's Exact Test w/ FDR Correction)")


## SPECIFICITY ##

## histogram
hist(treeWAS$specificity)
hist(treeWAS$specificity, xlim=c(0,1))
hist(fisher.bonf$specificity, xlim=c(0,1))
hist(fisher.fdr$specificity, xlim=c(0,1))

## density dist
plot(density(treeWAS$specificity), main="SPECIFICITY\n(treeWAS)")
plot(density(fisher.bonf$specificity), main="SPECIFICITY\n(Fisher's Exact Test w/ Bonferonni Correction)")
plot(density(fisher.fdr$specificity), main="SPECIFICITY\n(Fisher's Exact Test w/ FDR Correction)")


## FPR ##

## histogram
hist(treeWAS$FPR)
hist(treeWAS$FPR, xlim=c(0,1))
hist(fisher.bonf$FPR, xlim=c(0,1))
hist(fisher.fdr$FPR, xlim=c(0,1))

## density dist
plot(density(treeWAS$FPR), main="FALSE POSITIVE RATE\n(treeWAS)")
plot(density(fisher.bonf$FPR), main="FALSE POSITIVE RATE\n(Fisher's Exact Test w/ Bonferonni Correction)")
plot(density(fisher.fdr$FPR), main="FALSE POSITIVE RATE\n(Fisher's Exact Test w/ FDR Correction)")
