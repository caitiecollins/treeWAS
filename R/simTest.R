


## GENERIC SIM TESTING FUNCTION ##


#############
## simTest ##
#############

#
# ## CHECKING SNP.SIM.Q vs. SNP.SIM: ##
# #####################################
#
# ## load data:
# data(tree)
# data(phen)
# data(dist_0)
# data(dist_0.01)
# data(dist_0.05)
# data(dist_0.1)
#
# ## get phen loci:
# n.phen.subs <- 15
# grp.min <- 0.25
# seed <- 1
# phen.list <- phen.sim(tree, n.subs = n.phen.subs, grp.min = grp.min, seed = seed)
#
# ## set args:
# n.snps <- 10000
# n.subs <- dist_0.01
# n.snps.assoc <- 10
# assoc.prob <- 90
# tree <- tree
# phen.loci <- phen.list$phen.loci
# heatmap <- FALSE
# reconstruct <- FALSE
# dist.dna.model <- "JC69"
# row.names <- NULL
# set <- 1
# # seed <- 1
#
# #############
# ## SNP.SIM ##
# #############
# system.time(
#   snps.list <- snp.sim(n.snps=n.snps,
#                        n.subs=n.subs,
#                        n.snps.assoc=n.snps.assoc,
#                        assoc.prob=assoc.prob,
#                        tree=tree,
#                        phen.loci=phen.loci,
#                        heatmap=heatmap,
#                        reconstruct=reconstruct,
#                        dist.dna.model=dist.dna.model,
#                        row.names = NULL,
#                        set=set,
#                        seed=seed)
# ) # 0.764 (w n.snps.assoc = 10), or 0.256 (w NO snps.assoc)
#
# str(snps.list)
# snps.list$snps[1:10,snps.list$snps.assoc]
#
# ###############
# ## SNP.SIM.Q ##
# ###############
# ## args:
# # Q.mat <- matrix(c(2, 0.75, 0.75, 1, 3, 0.5,
# #                   0.25, 3, 3, 0.25, 0.5, 3, 1, 0.75, 0.75, 2), nrow = 4,
# #                 byrow = T, dimnames = rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
# s <- 20 # n.subs
# af <- 10 # association factor
# s <- s/sum(tree$edge.length)
# Q.mat <- matrix(c(NA, 1*s, 1*s, 0,
#                   1*af*s, NA, 0, 1*af*s,
#                   1*af*s, 0, NA, 1*af*s,
#                   0, 1*s, 1*s, NA),
#                 nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
# diag(Q.mat) <- sapply(c(1:nrow(Q.mat)), function(e) -sum(Q.mat[e, c(1:ncol(Q.mat))[-e]]))
#
#
# n.phen.subs <- 15
# grp.min <- 0.25
# set <- 3
#
#
# ## run snp.sim.Q:
# system.time(
#   snps.list.Q <- snp.sim.Q(n.snps=n.snps,
#                            n.subs=n.subs,
#                            snp.root = NULL,
#                            n.snps.assoc=n.snps.assoc,
#                            assoc.prob=assoc.prob,
#                            Q = Q.mat,
#                            tree=tree,
#                            n.phen.subs = n.phen.subs,
#                            phen.loci=NULL,
#                            heatmap=heatmap,
#                            reconstruct=reconstruct,
#                            dist.dna.model=dist.dna.model,
#                            grp.min = grp.min,
#                            row.names = NULL,
#                            set=set,
#                            seed=seed)
# ) # 1.132 (w n.snps.assoc = 10), or 0.256 (w NO snps.assoc)
#
# str(snps.list.Q)
# snps.list.Q$snps[1:10,snps.list.Q$snps.assoc]



###################
## R PROFILING?? ##
###################
## profiling of time (and memory??)
## COALESCENT TREE (w n.snps.sim = 10000)
# Rprof("E:/treeWAS_Sims/Rprof_simTest", memory.profiling=T)
# foo <- simTest()
# Rprof(NULL)
# summaryRprof("E:/treeWAS_Sims/Rprof_simTest") # , memory="both"

## R-TREE (w n.snps.sim = 10000)
# Rprof("E:/treeWAS_Sims/Rprof_simTest_rtree", memory.profiling=T)
# # foo <- simTest()
# Rprof(NULL)
# summaryRprof("E:/treeWAS_Sims/Rprof_simTest_rtree") # , memory="both"

###################

## PCA/DAPC error?
# Error in weights * y : non-numeric argument to binary operator

## set 3 error:
# Error in text.default(x = (max(h.null$breaks) * 3/4), y = (max(h.null$counts) *  :
#    object 'myCol' not found
# In addition: Warning message:
#   In matrix(as.numeric(snps), nrow = nrow(snps.ori), ncol = ncol(snps.ori)) :
#   data length [47700] is not a sub-multiple or multiple of the number of columns [10000]



# ## rtree error:
# Error in if (snp.node[[edges[e, 1]]] == "0" & phen.node[[edges[e, 1]]] ==  :
#              missing value where TRUE/FALSE needed

# HYP -- problem in snp.sim.Q snps.assc sim -- edges not going in ordered pairs
### SOLN? --> need to add which(edge == e) as in snps.sim instead of just using snp.node[[edges[e,1]]] to get index...




#################
## PARALLELIZE ##
#################

# ## Load libraries:
# library(foreach)
# library(doParallel)
#
# ## Calculate the number of cores:
# no_cores <- detectCores() - 1 # 3
#
# ## Initiate implicit cluster:
# registerDoParallel(no_cores)

################################################################################################################################
################################################################################################################################



################################################################################################################################
################################################################################################################################


#######################
## SET 2 (ACE / ML!) ##
#######################

# out <- foreach(n.reps=rep(1, 20), file.n=c(1:20), .packages="treeWAS") %dopar%
#   simTest(
#     ## simTest args:
#     set.number = 2,
#     n.reps = n.reps,
#     set.seed.as = "file.number",
#     working.dir = "/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims",
#
#     ## data from file args:
#     from.file = FALSE,
#     file.n = file.n,
#     Windows=TRUE,
#
#     ## coalescent.sim args:
#     n.ind = 100,
#     n.snps = 10000, # gen.size
#     # sim.by = "locus",
#     n.subs = dist_0, # 15, # theta (*2)
#     n.phen.subs = 15, # theta_p = NULL # 15
#     n.snps.assoc = 10, #
#     # assoc.option = "all",
#     assoc.prob = NULL, #100, #  90, #
#     grp.min = 0.25,
#     s = 20,
#     af = 10,
#     coaltree = TRUE,
#
#     ## treeWAS args:
#     ## RUNNING ALL OF THESE OPTIONS (FOR NOW):
#     p.value = 0.01, # REQUIRED FOR FISHER TEST
#     #   p.value.correct = c("bonf", "fdr", FALSE), #mt.correct = FALSE
#     #   p.value.by = c("count", "density"),
#     sim.n.snps = 100000, # 10*n.snps #sim.gen.size = NULL ###################### CAREFUL (10,000) !!
#     treeWAS.test = c("terminal", "simultaneous", "subsequent"), # "score"
#     snps.reconstruction = "ace",
#     phen.reconstruction = "ace"
#   )
# ###############################################################################################################################
#
#
#
# out <- foreach(n.reps=rep(1, 20), file.n=c(21:40), .packages="treeWAS") %dopar%
#   simTest(
#     ## simTest args:
#     set.number = 2,
#     n.reps = n.reps,
#     set.seed.as = "file.number",
#     working.dir = "/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims",
#
#     ## data from file args:
#     from.file = FALSE,
#     file.n = file.n,
#     Windows=TRUE,
#
#     ## coalescent.sim args:
#     n.ind = 100,
#     n.snps = 10000, # gen.size
#     # sim.by = "locus",
#     n.subs = dist_0.01, # 15, # theta (*2)
#     n.phen.subs = 15, # theta_p = NULL # 15
#     n.snps.assoc = 10, #
#     # assoc.option = "all",
#     assoc.prob = NULL, #100, #  90, #
#     grp.min = 0.25,
#     s = 20,
#     af = 10,
#     coaltree = TRUE,
#
#     ## treeWAS args:
#     ## RUNNING ALL OF THESE OPTIONS (FOR NOW):
#     p.value = 0.01, # REQUIRED FOR FISHER TEST
#     #   p.value.correct = c("bonf", "fdr", FALSE), #mt.correct = FALSE
#     #   p.value.by = c("count", "density"),
#     sim.n.snps = 100000, # 10*n.snps #sim.gen.size = NULL ###################### CAREFUL (10,000) !!
#     treeWAS.test = c("terminal", "simultaneous", "subsequent"), # "score"
#     snps.reconstruction = "ace",
#     phen.reconstruction = "ace"
#   )
# ###############################################################################################################################
#
#
# out <- foreach(n.reps=rep(1, 20), file.n=c(41:60), .packages="treeWAS") %dopar%
#   simTest(
#     ## simTest args:
#     set.number = 2,
#     n.reps = n.reps,
#     set.seed.as = "file.number",
#     working.dir = "/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims",
#
#     ## data from file args:
#     from.file = FALSE,
#     file.n = file.n,
#     Windows=TRUE,
#
#     ## coalescent.sim args:
#     n.ind = 100,
#     n.snps = 10000, # gen.size
#     # sim.by = "locus",
#     n.subs = dist_0.05, # 15, # theta (*2)
#     n.phen.subs = 15, # theta_p = NULL # 15
#     n.snps.assoc = 10, #
#     # assoc.option = "all",
#     assoc.prob = NULL, #100, #  90, #
#     grp.min = 0.25,
#     s = 20,
#     af = 10,
#     coaltree = TRUE,
#
#     ## treeWAS args:
#     ## RUNNING ALL OF THESE OPTIONS (FOR NOW):
#     p.value = 0.01, # REQUIRED FOR FISHER TEST
#     #   p.value.correct = c("bonf", "fdr", FALSE), #mt.correct = FALSE
#     #   p.value.by = c("count", "density"),
#     sim.n.snps = 100000, # 10*n.snps #sim.gen.size = NULL ###################### CAREFUL (10,000) !!
#     treeWAS.test = c("terminal", "simultaneous", "subsequent"), # "score"
#     snps.reconstruction = "ace",
#     phen.reconstruction = "ace"
#   )
# ###############################################################################################################################
#
#
# out <- foreach(n.reps=rep(1, 20), file.n=c(61:80), .packages="treeWAS") %dopar%
#   simTest(
#     ## simTest args:
#     set.number = 2,
#     n.reps = n.reps,
#     set.seed.as = "file.number",
#     working.dir = "/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims",
#
#     ## data from file args:
#     from.file = FALSE,
#     file.n = file.n,
#     Windows=TRUE,
#
#     ## coalescent.sim args:
#     n.ind = 100,
#     n.snps = 10000, # gen.size
#     # sim.by = "locus",
#     n.subs = dist_0.1, # 15, # theta (*2)
#     n.phen.subs = 15, # theta_p = NULL # 15
#     n.snps.assoc = 10, #
#     # assoc.option = "all",
#     assoc.prob = NULL, #100, #  90, #
#     grp.min = 0.25,
#     s = 20,
#     af = 10,
#     coaltree = TRUE,
#
#     ## treeWAS args:
#     ## RUNNING ALL OF THESE OPTIONS (FOR NOW):
#     p.value = 0.01, # REQUIRED FOR FISHER TEST
#     #   p.value.correct = c("bonf", "fdr", FALSE), #mt.correct = FALSE
#     #   p.value.by = c("count", "density"),
#     sim.n.snps = 100000, # 10*n.snps #sim.gen.size = NULL ###################### CAREFUL (10,000) !!
#     treeWAS.test = c("terminal", "simultaneous", "subsequent"), # "score"
#     snps.reconstruction = "ace",
#     phen.reconstruction = "ace"
#   )
# ###############################################################################################################################


################################################################################################################################
################################################################################################################################


## AT THE END, RUN:
# stopImplicitCluster()




## PROBLEM -- WHY WERE THERE SO FEW UNIQUE COLUMNS BEFORE AND NOW SO MANY? (just bc dist_0 -> dist_0.01???)














###############################################################################################################################################

###############################################################################################################################################



# out <- simTest(
#   ## simTest args:
#   set.number = 1,
#   n.reps = 1,
#   set.seed.as = "file.number",
#   working.dir = "/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims",
#
#   ## data from file args:
#   from.file = FALSE,
#   file.n = 21,
#   Windows=TRUE,
#
#   ## coalescent.sim args:
#   n.ind = 100,
#   n.snps = 10000, # gen.size
#   # sim.by = "locus",
#   n.subs = dist_0.01, # 15, # theta (*2)
#   n.phen.subs = 15, # theta_p = NULL # 15
#   n.snps.assoc = 10, #
#   # assoc.option = "all",
#   assoc.prob = NULL, #100, #  90, #
#   grp.min = 0.25,
#   s = 20,
#   af = 10,
#   coaltree = TRUE,
#
#   ## treeWAS args:
#   ## RUNNING ALL OF THESE OPTIONS (FOR NOW):
#   p.value = 0.01, # REQUIRED FOR FISHER TEST
#   #   p.value.correct = c("bonf", "fdr", FALSE), #mt.correct = FALSE
#   #   p.value.by = c("count", "density"),
#   sim.n.snps = 100000, # 10*n.snps #sim.gen.size = NULL ###################### CAREFUL (10,000) !!
#   treeWAS.test = c("terminal", "simultaneous", "subsequent"), # "score"
#   snps.reconstruction = "ace",
#   phen.reconstruction = "ace"
# )


###############################################################################################################################################

#
# # simTest args:
# set.number = 3
# n.reps = 80
# set.seed.as = "file.number"
# working.dir = "/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims"
#
# ## data from file args:
# from.file = FALSE
# file.n = c(1:80) # Need to update this if not starting from 1!
# Windows = TRUE
# cluster = TRUE
#
# ## coalescent.sim args:
# # n.ind = sample(round(seq(50, 200, length.out = 80), 0), 80, replace = FALSE) # 100
# # n.snps = sample(round(seq(10000, 100000, length.out = 80), 0), 80, replace = FALSE) # 10000
# # n.subs =  dist_0.01 # 1 #theta (*2)
# # n.phen.subs = 15 # 5 #theta_p = NULL ###
# # n.snps.assoc = round((n.snps/1000), 0) # = 0
# n.ind = 100
# n.snps = 2000
# n.subs <- c(rep(list(dist_0), 20), rep(list(dist_0.01), 20), rep(list(dist_0.05), 20), rep(list(dist_0.1), 20))
# n.phen.subs = 15 # 5 #theta_p = NULL ###
# n.snps.assoc = 10
# assoc.prob = 100
# grp.min = 0.25
# s = 20
# af = 10
# coaltree = TRUE
#
# ## treeWAS args:
# p.value = 0.01
# p.value.correct = "bonf"
# p.value.by = "count"
# sim.n.snps = n.snps*10 # 100000 # 10*n.snps #sim.gen.size = NULL
# treeWAS.test = c("terminal", "simultaneous", "subsequent") # "score"
# snps.reconstruction = "parsimony" # "ace" #
# phen.reconstruction = "parsimony" # "ace" #




########################################################################

###################
## DOCUMENTATION ##
###################

#' Simulation Testing.
#'
#' Generic simulation-testing function used to validate treeWAS performance on simulated datasets. Not designed for public use!
#'
#' @param test A character string or vector containing one or more of the following available tests of association:
#' "terminal", "simultaneous", "subsequent", "cor", "fisher". By default, the first three tests are run.
#' See details for more information on what these tests do and when they may be appropriate.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#  adegenet ape
#' @importFrom Hmisc all.is.numeric
#'
#' @export

########################################################################


simTest <- function(

  ## simTest args:
  set.number = 3,
  n.reps = 1,
  set.seed.as = "file.number",
  working.dir = "/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims",

  ## data from file args:
  from.file = FALSE,
  file.n = NULL,
  Windows = FALSE,
  cluster = FALSE,

  ## coalescent.sim args:
  n.ind = 100,
  n.snps = 10000, # gen.size
  n.subs = 1, # theta (*2)
  n.phen.subs = 15, # theta_p = NULL
  n.snps.assoc = 10, # = 0
  assoc.prob = 90, # 100 (set2)
  grp.min = 0.25,
  s = 20,
  af = 10,
  coaltree = TRUE,

  ## treeWAS args:
  ## RUNNING ALL OF THESE OPTIONS (FOR NOW):
  p.value = 0.01, # REQUIRED FOR FISHER TEST
  p.value.correct = "bonf", # c("bonf", "fdr", FALSE), #mt.correct = FALSE
  p.value.by = "count", #  c("count", "density"),
  sim.n.snps = 100000, # 10*n.snps #sim.gen.size = NULL
  treeWAS.test = c("terminal", "simultaneous", "subsequent"), # "score"
  snps.reconstruction = "parsimony",
  phen.reconstruction = "parsimony"

){


  ###############################################
  ## make lists in which to store all ###########
  ## data and output from each of n.reps runs: ##
  ###############################################
  SNPS <- PHEN <- PHEN.PLOT.COL <-  TREE <- OUT <- RES <-
    FISHER.RESULTS <- PLINK.RESULTS <- PCA.RESULTS <- DAPC.RESULTS <- CMH.RESULTS <-
    ARGS <- PERFORMANCE <- SCORE3 <- list()
  ## and make lists for saving filenames
  filename.snps <- filename.phen <- filename.phen.plot.col <- filename.tree <-
    filename.out <- filename.res <- filename.fisher.results <-
    filename.plink.results <- filename.pca <- filename.dapc <- filename.cmh <-
    filename.args <- filename.performance <- filename.score3 <-
    filename.plot <- filename.tree.plot <- list()


  ####################################################################################################################################
  ################################### *** DEFINE ARGUMENTS | SET *** #################################################################
  ####################################################################################################################################

  ## coalescent.sim args:
  if(missing(n.reps)) n.reps <- 1
  if(missing(set.seed.as)) set.seed.as <- "file.number"
  if(missing(n.ind)) n.ind <- 100
  if(missing(n.snps)) n.snps <- 10000
  if(missing(n.subs)) n.subs <- 1
  if(missing(n.phen.subs)) n.phen.subs <- NULL
  if(missing(n.snps.assoc)) n.snps.assoc <- 10
  if(missing(assoc.prob)) assoc.prob <- 90
  if(missing(grp.min)) grp.min <- 0.25

  ## treeWAS args:
  #   if(missing(p.value)) p.value <- 0.001
  #   if(missing(p.value.correct)) p.value.correct <- c("bonf", "fdr", FALSE)
  #   if(missing(p.value.by)) p.value.by <- c("count", "density")
  if(missing(sim.n.snps)) sim.n.snps <- 100000
  if(missing(treeWAS.test)) treeWAS.test <- c("terminal", "simultaneous", "subsequent")
  if(missing(snps.reconstruction)) snps.reconstruction <- "parsimony"
  if(missing(phen.reconstruction)) phen.reconstruction <- "parsimony"

  ## ensure phen is NULL (will be generated by sim)
  phen <- NULL

  ###########
  ## SET 1 ##
  ###########
  if(set.number == 1){
    if(!is.null(phen)) phen <- NULL
    ### ensure n.phen.subs is NOT null
    if(is.null(n.phen.subs)) n.phen.subs <- 15
    if(is.null(assoc.prob)) assoc.prob <- 90
  } # end set 1

  ###########
  ## SET 2 ##
  ###########
  if(set.number == 2){
    if(!is.null(phen)) phen <- NULL
    ## ensure n.phen.subs is NOT null
    if(is.null(n.phen.subs)) n.phen.subs <- 15
    assoc.prob <- 100
  } # end set 2

  ###########
  ## SET 3 ##
  ###########
  if(set.number == 3){
    if(!is.null(phen)) phen <- NULL
    ## n.phen.subs & assoc.prob not actually being used in snp.sim.Q
    if(is.null(n.phen.subs)) n.phen.subs <- 15
    if(is.null(assoc.prob)) assoc.prob <- 90
    ## ensure 10 associated SNPs
    if(is.null(n.snps.assoc)){
      if(length(n.snps) == 1 & n.snps[1] == 2000){
        n.snps.assoc <- 10
      }else{
        n.snps.assoc <- round((n.snps/1000), 0)
      }
    }
    ## set sim.n.snps:
    if(is.null(sim.n.snps)){
      sim.n.snps <- n.snps*10
    }

  } # end set 3

  ###########
  ## SET 4 ##
  ###########
  if(set.number == 4){
    if(!is.null(phen)) phen <- NULL
    ## ensure n.phen.subs is NOT null
    if(is.null(n.phen.subs)) n.phen.subs <- 15
    ## ensure MULTIPLE associated SNPs
    if(is.null(n.snps.assoc)) n.snps.assoc <- 10
    if(is.null(assoc.prob)) assoc.prob <- 90
  } # end set 4

  ## If cluster, no from.file!
  if(cluster == TRUE){
    from.file <- FALSE
    working.dir <- "" # just in case?
  }


  args <- snps.assoc <- NULL
  sim.n.snps.ori <- sim.n.snps

  ## change working dir if WINDOWS:
  if(Windows == TRUE){
    working.dir <- c("C:/Cait 2016/Work/Xavier/Sims")
  }

  ##############################################################################################################################
  ########################################### *** DATA FROM FILE *** ###########################################################
  ##############################################################################################################################

  if(from.file==TRUE){

    ###################################################
    ## READING IN DATA (OPTION for POST HOC TESTING) ##
    ###################################################

    ## Get from file: snps, phen, tree, performance, snps.assoc

    if(is.null(file.n)) stop("Use argument file.n
                             to specify which files to read in.")

    ## get n.reps for main for loop
    n.reps <- length(file.n)

    ## set working directory for the set specified
    wd <- paste(working.dir, "/", "set", set.number, sep="")
    setwd(wd)

    ## for loop:
    for(i in 1:n.reps){

      ## get filename prefix for this round of this set:
      filename.prefix <- paste("set", set.number, "_", file.n[i], "_", sep="")

      ## get snps
      filename <- paste("./", filename.prefix, "snps.Rdata", sep="")
      snps <- get(load(filename))

      ## get phen
      filename <- paste("./", filename.prefix, "phen.Rdata", sep="")
      phen <- get(load(filename))

      ## get tree
      filename <- paste("./", filename.prefix, "tree.Rdata", sep="")
      tree <- get(load(filename))

      ## get performance
      filename <- paste("./", filename.prefix, "performance.Rdata", sep="")
      performance <- get(load(filename))

      ## get snps.assoc (if any)
      snps.assoc <- NULL
      if(!is.null(performance$snps.assoc)) snps.assoc <- performance$snps.assoc

    } # end for loop


  } # end from.file == TRUE


  ##############################################################################################################################

  ##################################
  ## CONVERT VECTOR ARGS --> LIST ##
  ##################################
  ## For args of length n.reps, make lists:
  if(length(n.ind) == 1){
    N.IND <- as.list(rep(n.ind[1], n.reps))
    cat("n.ind was of length 1. Using n.ind[1] for all reps.\n")
  }else{
    N.IND <- as.list(n.ind)
    if(length(n.ind) != n.reps) warning("n.ind was not of length 1 or n.reps.\n")
  }

  if(length(n.snps) == 1){
    N.SNPS <- as.list(rep(n.snps[1], n.reps))
    cat("n.snps was of length 1. Using n.snps[1] for all reps.\n")
  }else{
    N.SNPS <- as.list(n.snps)
    if(length(n.snps) != n.reps) warning("n.snps was not of length 1 or n.reps.\n")
  }

  if(length(n.snps.assoc) == 1){
    N.SNPS.ASSOC <- as.list(rep(n.snps.assoc[1], n.reps))
    cat("n.snps.assoc was of length 1. Using n.snps.assoc[1] for all reps.\n")
  }else{
    N.SNPS.ASSOC <- as.list(n.snps.assoc)
    if(length(n.snps.assoc) != n.reps) warning("n.snps.assoc was not of length 1 or n.reps.\n")
  }

  if(length(sim.n.snps) == 1){
    SIM.N.SNPS <- as.list(rep(sim.n.snps[1], n.reps))
    cat("sim.n.snps was of length 1. Using sim.n.snps[1] for all reps.\n")
  }else{
    SIM.N.SNPS <- as.list(sim.n.snps)
    if(length(sim.n.snps) != n.reps) warning("sim.n.snps was not of length 1 or n.reps.\n")
  }

  if(length(n.subs) == 1 | !is.list(n.subs)){
    N.SUBS <- as.list(rep(list(n.subs), n.reps))
    cat("n.subs was of length 1. Using same n.subs for all reps.\n")
  }else{
    N.SUBS <- as.list(n.subs)
    if(length(N.SUBS) != n.reps){
      N.SUBS <- as.list(rep(n.subs, n.reps))[1:n.reps]
      warning("n.subs was not of length 1 or n.reps: repeating sequence until n.subs met.\n")
    }
  }
  ##############################################################################################################################

  ##############
  ## FOR LOOP ##
  ##############

  for(i in 1:n.reps){

    ## get list args i:
    n.ind <- N.IND[[i]]
    n.snps <- N.SNPS[[i]]
    n.snps.assoc <- N.SNPS.ASSOC[[i]]
    sim.n.snps <- SIM.N.SNPS[[i]]
    n.subs <- N.SUBS[[i]]

    wd <- paste(working.dir, "/", "set", set.number, sep="")
    if(cluster == FALSE) setwd(wd)

    ## get file number:
    if(from.file==FALSE){
      if(cluster == FALSE){
        ## get number for group | number of set3_snps already in file:
        number <- (length(grep("_snps", dir("./")))+1)
      }

      ## Unless file.n is provided as argument!
      if(!is.null(file.n)){
        if(!is.null(file.n[i])){
          number <- file.n[i]
        }
      }

    }else{
      number <- file.n[i]
    }


    ################
    ## dummy plot ##
    ################
    ## to give user indication of what round of simTest.set3 we are on:
    par(mfrow=c(1,1))
    round.marker <- paste("\n\n\n\n\n\n\n\n\n\n\n\nROUND", i,
                          "\nset", set.number,
                          "\n(file number", number, ")",
                          sep=" ")
    plot.new()
    title(round.marker, adj=0.5)

    ###############
    ## set seed? ##
    ###############
    if(!is.null(set.seed.as)){
      if(set.seed.as == "file.number"){
        seed <- number
        set.seed(seed)
      }else{
        if(length(set.seed.as) == n.reps){
          seed <- set.seed.as[i]
          set.seed(seed)
        }else{
          warning("seed is not of length n.reps; seed will not be set.")
          seed <- NULL
        }
      }
    }
    # end set.seed.as

    ##############################################################################################################################
    ################################################ *** COALESCENT.SIM *** ######################################################
    ##############################################################################################################################

    if(from.file==FALSE){

      ########################
      ## get PHEN for SET 1 ##
      ########################
      #       if(set.number == 1){
      #         ## simulate phen first, by random sampling ##
      #         ## enforce even split of cases and controls??
      #         #phen <- sample(c(rep("A", floor(n.ind/2)), rep("B", ceiling(n.ind/2)), replace=FALSE))
      #         ## or just draw phen by purely random sampling...?
      #         phen <- sample(c("A", "B"), n.ind, replace=TRUE)
      #         phen <- as.factor(phen)
      #       }else{
      #         phen <- NULL
      #       }

      ############################
      ## simulate data and tree ##
      ############################

      ## TESTING -- used to be phen = phen (but Hyp: multiple rounds causing problems--if so, could rename to phen.prior, eg.)
      ## CHECK:
      # print("NUMBER"); print(number)
      # print("PHEN BEFORE"); print(phen)

      # gc()

      filename.tree.plot[[i]] <- paste("./set", set.number, "_", number, "_tree_plot", ".pdf", sep="")

      filename.panel.plot <- paste("./set", set.number, "_", number, "_panel_plot", ".pdf", sep="")

      foo <- coalescent.sim(n.ind = n.ind,
                            n.snps = n.snps,
                            n.subs = n.subs,
                            n.snps.assoc = n.snps.assoc,
                            assoc.prob = assoc.prob,
                            n.phen.subs = n.phen.subs,
                            phen = NULL,
                            plot = TRUE,
                            heatmap = FALSE,
                            reconstruct = FALSE,
                            dist.dna.model = "JC69",
                            row.names = NULL,
                            grp.min = grp.min,
                            coaltree = coaltree,
                            set=set.number,
                            s = s,
                            af = af,
                            filename=list(filename.tree.plot[[i]],
                                          filename.panel.plot),
                            seed = seed)


      print("coalescent done")
      # gc()
      ####################################
      ## isolate common elements of foo ##
      ####################################
      snps <- snps.ori <- snps.ori.ori <- foo$snps
      if(!is.null(n.snps.assoc)) if(n.snps.assoc > 0){
        snps.assoc <- snps.assoc.ori <- snps.assoc.loci <- foo$snps.assoc
      }else{
        snps.assoc <- NULL
      }
      phen <- phen.ori <- phen.ori.ori <- foo$phen
      tree <- tree.ori <- foo$tree
      phen.plot.col <- foo$phen.plot.col

      ## snps names:
      if(is.null(colnames(snps))) colnames(snps) <- c(1:ncol(snps))
      snps.names <- colnames(snps)

      ##########################################
      ## isolate set-specific elements of foo ##
      ##########################################
      # if(is.null(phen)) phen <- foo$phen
      ## MAKE SURE PHEN IS IN CORRECT ORDER OF INDS NOT IN LEAF ORDER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ######################
      ## save plot as pdf ##
      ######################
      ## CHECK-- Not sure, but may be safer to write this using pdf() instead of dev.copy ~ treeWAS plots...?
      # filename.tree.plot[[i]] <- paste("./set", set.number, "_", number, "_tree_plot", ".pdf", sep="")
      # dev.copy(pdf, file=filename.tree.plot[[i]], width=7, height=11) # , pointsize=1
      # dev.off()


    }# end from.file == FALSE



    ##############################################################################################################################
    ############################################### *** treeWAS *** ##############################################################
    ##############################################################################################################################

    sim.n.snps.ori <- sim.n.snps
    if(is.null(sim.n.snps)) sim.n.snps <- ncol(snps)*10

    #######################
    ## save treeWAS plot ##
    #######################
    ## NB: plot.png will not be viewable until fn has finished running...
    # filename.plot[[i]] <- list()
    # for(t in 1:length(treeWAS.test)){
    #
    #   ## Save both Manhattan and Hist per test:
    #   filename.plot[[i]][[t]] <- c(## manhattan:
    #     paste("./set",
    #           set.number,
    #           "_", number,
    #           "_plot_manhattan_",
    #           treeWAS.test[t],
    #           ".pdf", sep=""),
    #
    #     ## null.dist:
    #     paste("./set",
    #           set.number,
    #           "_", number,
    #           "_plot_",
    #           treeWAS.test[t],
    #           ".pdf", sep="")
    #   )
    # }

    ## Save both Manhattan and Hist per test:
    filename.plot[[i]] <-
                              ## null.dist:
                              paste("./set",
                                    set.number,
                                    "_", number,
                                    "_plot",
                                    ".pdf", sep="")

    #################
    ## RUN treeWAS ##
    #################
    print("treeWAS started")
    # gc()
    set.seed(seed)

    syst.time <- system.time( # 341
      out <- treeWAS(snps = snps,
                     phen = phen,
                     tree = tree,
                     n.subs = NULL,
                     n.snps.sim = sim.n.snps,
                     chunk.size = ncol(snps),
                     test = treeWAS.test, # c("terminal", "simultaneous", "subsequent")
                     snps.reconstruction =  snps.reconstruction, # "parsimony",
                     snps.sim.reconstruction = "parsimony",
                     phen.reconstruction = phen.reconstruction,
                     phen.type = NULL,
                     na.rm = TRUE,
                     p.value = p.value, # 0.01
                     p.value.correct = p.value.correct, # "bonf
                     p.value.by = p.value.by, # "count"
                     dist.dna.model = "JC69",
                     plot.tree = FALSE,
                     plot.manhattan = TRUE,
                     plot.null.dist = TRUE,
                     plot.dist = FALSE,
                     snps.assoc = snps.assoc,
                     filename.plot = filename.plot[[i]],
                     seed = 1)
    )

    print("treeWAS done")
    # gc()

    # i <- 1 #####

    #####

    # dev.copy(pdf, file=filename.plot[[i]], width=7, height=11) # , pointsize=12
    # dev.off()

    ##################################
    ## isolate df of Sig SNPs found ##
    ##################################
    res.complete <- out
    res <- out$res

    treeWAS.all <- out$treeWAS.combined ## also stored in res.complete, which gets saved later.


    ##############
    ## get call ##
    ##############
    ## get arguments inputed to simTest for this round
    #call <- match.call()

    ## TO DO: update arguments to contain actual values used
    ## (ie. return actual value of seed instead of set.seed.as)
    ## (eg. if n.phen.subs was NULL and got set to 25 due to set.number,
    ## report 25)

    args <- NA
    # args <- mget(names(formals()), sys.frame(sys.nframe()))


    ###########################################################################################################################
    ############################################# *** FISHER TEST *** #########################################################
    ###########################################################################################################################

    #############################
    ## RUN FISHER'S EXACT TEST ##
    #############################

    ## NOTE--we will run and save all 3 (uncorrected, Bonferroni-corrected, and FDR-corrected),
    ######### BUT we will only calculate performance metrics for the latter two (or only FDR??)...

    ## RUN TEST
    pval.fisher <- sapply(c(1:ncol(snps)),
                          function(e) fisher.test(snps[,e], y=phen,
                                                  alternative="two.sided")$p.value)
    ## two.sided bc we want to know if
    ## inds w the phen have EITHER more 1s or 0s

    p.thresh <- p.value

    ## WITHOUT CORRECTION, identify sig.snps
    fisher.snps.uncorr <- colnames(snps)[which(pval.fisher < p.thresh)]
    n.fisher.snps.uncorr <- length(fisher.snps.uncorr)

    ## w BONFERRONI CORRECTION, identify sig.snps
    pval.bonf <- p.adjust(pval.fisher, method="bonferroni", n=length(pval.fisher))
    fisher.snps.bonf <- colnames(snps)[which(pval.bonf < p.thresh)]
    n.fisher.snps.bonf <- length(fisher.snps.bonf)

    ## w FDR CORRECTION, identify sig.snps
    pval.fdr <- p.adjust(pval.fisher, method="fdr", n=length(pval.fisher))
    fisher.snps.fdr <- colnames(snps)[which(pval.fdr < p.thresh)]
    n.fisher.snps.fdr <- length(fisher.snps.fdr)

    ## CONVERT 0-LENGTH RESULTS TO NULL
    if(length(fisher.snps.uncorr) == 0) fisher.snps.uncorr <- NULL
    if(length(fisher.snps.bonf) == 0) fisher.snps.bonf <- NULL
    if(length(fisher.snps.fdr) == 0) fisher.snps.fdr <- NULL

    ## STORE FISHER TEST RESULTS ##
    fisher.results <- list(pval.fisher, fisher.snps.uncorr, fisher.snps.bonf, fisher.snps.fdr)
    names(fisher.results) <- c("pval.fisher", "fisher.snps.uncorr", "fisher.snps.bonf", "fisher.snps.fdr")

    ## end fisher tests


    ###########################################################################################################################
    ################################################ *** PLINK  *** ###########################################################
    ###########################################################################################################################


    ###################
    ## RUNNING PLINK ##
    ###################

    # ## Set working directory to run plink program
    # if(Windows == FALSE){
    #   setwd("/media/caitiecollins/88CC9BCECC9BB4C2/Program Files/plink-1.07-dos")
    # }else{
    #   setwd("C:/Program Files/plink-1.07-dos")
    # }
    # plink.wd <- getwd()
    # # setwd(plink.wd)
    #
    # #################
    # ## Handle DATA ##
    # #################
    #
    # ## STORE ORIGINAL DATA FOR LATER ##
    # snps.ori <- snps.ori.ori
    # phen.ori <- phen.ori.ori
    #
    # #     #######################
    # #     ## from.file = FALSE ##
    # #     #######################
    # #     if(from.file == FALSE){
    #
    # ## If data was created during THIS round of simTest...
    #
    # ###########################
    # ## CONVERT AND SAVE DATA ##
    # ###########################
    #
    # ## set working directory for saving
    # dat.wd <- paste(working.dir, "/set", set.number, "/", sep="")
    # setwd(dat.wd)
    #
    # ##################################
    # ## convert snps --> ped format: ##
    # ##################################
    #
    # ## PED FORMAT: ##
    # ## 6 "MANDATORY" COLUMNS first:
    # ## FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype
    # #### ... EXCLUDE OTHER "mandatory" columns w PLINK commands
    # #### --> (IndividualID, Phenotype)
    # ## 7-to-p: GENOTYPE COLUMNS
    # ## NO HEADER ROW (belongs in .MAP file)
    #
    # ## SNPs can NOT be 0 --> convert from 0/1 to 1/2:
    # snps <- snps+1
    # row.names(snps) <- paste("ind", row.names(snps), sep=".")
    # ## save row and column names for elsewhere:
    # individualID <- row.names(snps)
    # loci.names <- colnames(snps)
    #
    # ## convert to matrix:
    # snps <- matrix(snps, byrow=FALSE, ncol=ncol(snps))
    # ## Replace every other column with a copy of the column before it:
    # ## NOTE--THIS MAY CHANGE (ie. IF WORKING INSTEAD WITH DATA AS MT DNA....)
    # #snps[,seq(2, ncol(snps), 2)] <- snps[,seq(1, ncol(snps), 2)]
    # ## DUPLICATE columns:
    # snps.ori <- snps
    # snps.new <- matrix(NA, nrow=nrow(snps.ori), ncol=2*ncol(snps.ori))
    # snps.new[, seq(1, 2*ncol(snps.ori), 2)] <- snps.ori
    # snps.new[, seq(2, 2*ncol(snps.ori), 2)] <- snps.ori
    # snps <- snps.new
    #
    # ## recode phen: S/A as 0, R/B as 1:
    # phen <- as.character(phen)
    # # replace phen coding w 1/2:
    # phen <- replace(phen, which(phen %in% c("R", "B", "1")), 2)
    # phen <- replace(phen, which(phen %in% c("S", "A", "0")), 1)
    # # phen <- as.numeric(phen) #?
    # # names(phen) <- names(phen.ori) #?
    #
    # ## bind individualID, phen, snps
    # dat <- cbind(individualID, phen, snps)
    # colnames(dat) <- NULL
    #
    # ped <- dat
    #
    # ## get filename
    # uniqueID <- paste("set", set.number, "_", number, sep="")
    # filename <- paste("./", uniqueID, "_ped.Rdata", sep="")
    #
    # ## save dat.ped as Rdata
    # save(ped, file=filename)
    # #ped <- get(load("C:/Cait 2012/Work/Xavier/Sims/set3/set3_1_ped.Rdata"))
    #
    # ## save as text!
    # #ped <- dat
    # if(Windows == FALSE){
    #   filename <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".txt", sep="")
    # }else{
    #   filename <- paste("C:/PLINK/", uniqueID, ".txt", sep="")
    # }
    # write.table(ped, file=filename, quote=FALSE, row.names=FALSE, col.names=FALSE)
    #
    # ## Do NOT save as .PED
    # ## convert from text to PED!!
    #
    # ## UHOH -- STOPPED HERE -- SHELL COMMAND DOESNT SEEM TO BE WORKING FROM LINUX!!!!!!!!!!!!!!!!!!
    # ## SYSTEM COMMAND NOT READING FILE NAMES CORRECTLY....????????????????///
    #
    # #     filename.ori <- paste("\\media\\caitiecollins\\88CC9BCECC9BB4C2\\PLINK\\", uniqueID, ".txt", sep="")
    # #     filename.new <- paste("\\media\\caitiecollins\\88CC9BCECC9BB4C2\\PLINK\\", uniqueID, ".ped", sep="")
    #
    # if(Windows == FALSE){
    #   filename.ori <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".txt", sep="")
    #   filename.new <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".ped", sep="")
    #   command <- paste("mv", filename.ori, filename.new, sep=" ")
    # }else{
    #   filename.ori <- paste("C:\\PLINK\\", uniqueID, ".txt", sep="")
    #   filename.new <- paste("C:\\PLINK\\", uniqueID, ".ped", sep="")
    #   command <- paste("move", filename.ori, filename.new, sep=" ")
    # }
    #
    # ## run command
    # if(Windows == FALSE){
    #   system(command)
    # }else{
    #   shell(command)
    # }
    #
    #
    # ####################################
    # ## convert snp meta-data --> .map ##
    # ####################################
    #
    # ## MAP format: ##
    # ## EACH LINE of a .map file describes a SINGLE marker.
    # ## Contains 4 "MANDATORY" COLUMNS:
    # ## Chromosome, rs#/SNP identifier, (Genetic distance), Base-pair position.
    #
    #
    # ## loci.names:
    # ## SHOULD BE ONE NAME PER SITE (ie. PER TWO LOCI): ie. L001, NOT L001.1, L001.2 !!!!
    # ## remove last TWO characters (ie. decimal and trailing digit):
    # #       loci.names <- substr(loci.names, 1, nchar(loci.names)-2)
    # #       ## keep only every other:
    # #       loci.names <- loci.names[seq(1, length(loci.names), 2)]
    # ## OR: # loci.names <- unique(loci.names)
    #
    # ## make dummy variables for irrelevant fields:
    # chromosome <- rep(26, length(loci.names)) # 26 = human mitochondrial (haploid)
    # gen.dist <- rep(0, length(loci.names))
    # ## get base-pair posi (ie. loci name - L):
    # bp <- loci.names
    # bp <- as.numeric(gsub("L", "", bp))
    # dat <- data.frame(chromosome, loci.names, gen.dist, bp)
    #
    # ## as matrix, no header:
    # dat <- as.matrix(dat, byrow=FALSE, ncol=ncol(dat))
    # colnames(dat) <- NULL
    #
    # map <- dat
    #
    # ## get filename
    # uniqueID <- paste("set", set.number, "_", number, sep="")
    # filename <- paste("./", uniqueID, "_map.Rdata", sep="")
    #
    # ## save dat.map as Rdata
    # save(map, file=filename)
    # #map <- get(load("C:/Cait 2012/Work/Xavier/Sims/set3/set3_1_map.Rdata"))
    #
    # ## save as text!
    # if(Windows == FALSE){
    #   filename <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".txt", sep="")
    # }else{
    #   filename <- paste("C:/PLINK/", uniqueID, ".txt", sep="")
    # }
    # write.table(map, file=filename, quote=FALSE, row.names=FALSE, col.names=FALSE)
    #
    # ## Do NOT save as .MAP
    # ## convert from text to MAP!!
    #
    # #     filename.ori <- paste("C:\\PLINK\\", uniqueID, ".txt", sep="")
    # #     filename.new <- paste("C:\\PLINK\\", uniqueID, ".map", sep="")
    #
    # if(Windows == FALSE){
    #   filename.ori <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".txt", sep="")
    #   filename.new <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".map", sep="")
    #   command <- paste("mv", filename.ori, filename.new, sep=" ")
    # }else{
    #   filename.ori <- paste("C:\\PLINK\\", uniqueID, ".txt", sep="")
    #   filename.new <- paste("C:\\PLINK\\", uniqueID, ".map", sep="")
    #   command <- paste("move", filename.ori, filename.new, sep=" ")
    # }
    #
    # ## Run command
    # if(Windows == FALSE){
    #   system(command)
    # }else{
    #   shell(command)
    # }
    #
    #
    #
    # ###### ###### ###### ######
    #
    # #     }else{ # end from.file = FALSE
    # #
    # #       ######################
    # #       ## from.file = TRUE ##
    # #       ######################
    # #
    # #     } # end from.file = TRUE
    #
    #
    # ########## #################### ########## #################### #################### ########## ####################
    #
    #
    # #####################
    # ## GWAS with PLINK ##
    # #####################
    #
    # ## set wd for PLINK program
    # setwd(plink.wd)
    #
    # ## get filename
    # if(Windows == FALSE){
    #   filename <- paste("\\media\\caitiecollins\\88CC9BCECC9BB4C2\\PLINK\\", uniqueID, sep="")
    # }else{
    #   filename <- paste("C:\\PLINK\\", uniqueID, sep="")
    # }
    #
    #
    # ## inspect file?
    # command <- paste("plink --file ", filename, " --no-fid --no-parents --no-sex --allow-no-sex", sep="")
    # #     ## Run command
    # #     if(Windows == FALSE){
    # #       system(command)
    # #     }else{
    # #       shell(command)
    # #     }
    #
    # ## make a binary PED file ##
    # ## (provide the full path, not just the file name)
    # command <- paste("plink --file ", filename, " --make-bed --no-fid --no-parents --no-sex --allow-no-sex --out ", filename, sep="")
    #
    # ## Run command
    # if(Windows == FALSE){
    #   system(command)
    # }else{
    #   shell(command)
    # }
    # ########################################################################################################################################
    # ## SYSTEM(COMMAND) FAILS HERE WITH THE FOLLOWING ERROR #################################################    ####    ####    ####    ####
    #
    # # sh: 1: plink: not found
    #
    # ## PROBABLY NEED TO FIND SOLN AND APPLY IT TO ALL PLINK-RELATED CODE, STARTING WELL ABOVE HERE....................
    # ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####
    # ########################################################################################################################################
    #
    # ## use --bfile to work with the BINARY file
    # # (same as --file, but loads the binary one and prints summary stats)
    # command <- paste("plink --bfile ", filename, " --no-fid --no-parents --no-sex --allow-no-sex", sep="")
    #
    # ## Run command
    # if(Windows == FALSE){
    #   system(command)
    # }else{
    #   shell(command)
    # }
    #
    # #######################
    # ## basic association ##
    # #######################
    # ##--> 1df chi-square test
    #
    # ## check freq of SNPs...?
    # command <- paste("plink --file ", filename, " --no-fid --no-parents --no-sex --allow-no-sex --freq --out ", filename, sep="")
    #
    # ## Run command
    # if(Windows == FALSE){
    #   system(command)
    # }else{
    #   shell(command)
    # }
    #
    # ## yay!
    # if(Windows == FALSE){
    #   filename <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".frq", sep="")
    # }else{
    #   filename <- paste("C:/PLINK/", uniqueID, ".frq", sep="")
    # }
    # freq <- read.table(filename, header=TRUE)
    # #head(freq)
    #
    # ## perform a basic association analysis on the disease trait for all single SNPs
    # if(Windows == FALSE){
    #   filename <- paste("\\media\\caitiecollins\\88CC9BCECC9BB4C2\\PLINK\\", uniqueID, sep="")
    # }else{
    #   filename <- paste("C:\\PLINK\\", uniqueID, sep="")
    # }
    # command <- paste("plink --bfile ",  filename, " --assoc --counts --allow-no-sex --out ", filename, sep="")
    #
    # ## Run command
    # if(Windows == FALSE){
    #   system(command)
    # }else{
    #   shell(command)
    # }
    #
    # ## to view the file you created, just read it in with R:
    # if(Windows == FALSE){
    #   filename <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".assoc", sep="")
    # }else{
    #   filename <- paste("C:/PLINK/", uniqueID, ".assoc", sep="")
    # }
    # plink.res <- read.table(filename, header=TRUE)
    # # head(plink.res)
    #
    #
    # ## get p.vals
    # pval.plink.assoc <- plink.res$P
    #
    # ## get sig ##
    #
    # ## p.thresh:
    # p.thresh <- p.value # 0.05 # 0.01 # 0.001 # ??
    #
    #
    # ## Uncorrected ##
    # plink.assoc.snps.uncorr <- snps.names[which(pval.plink.assoc < p.thresh)]
    #
    # ## Bonferonni ##
    # p.vals.bonf <- p.adjust(pval.plink.assoc, "bonferroni")
    # p.bonf <- which(p.vals.bonf < p.thresh)
    # plink.assoc.snps.bonf <- snps.names[p.bonf]
    #
    # ## FDR ##
    # p.vals.fdr <- p.adjust(pval.plink.assoc, "fdr")
    # p.fdr <- which(p.vals.fdr < p.thresh)
    # plink.assoc.snps.fdr <- snps.names[p.fdr]
    #
    # ## CONVERT 0-LENGTH RESULTS TO NULL
    # if(length(plink.assoc.snps.uncorr) == 0) plink.assoc.snps.uncorr <- NULL
    # if(length(plink.assoc.snps.bonf) == 0) plink.assoc.snps.bonf <- NULL
    # if(length(plink.assoc.snps.fdr) == 0) plink.assoc.snps.fdr <- NULL
    #
    # ## STORE PLINK TEST RESULTS ##
    # plink.assoc.results <- list(pval.plink.assoc, plink.assoc.snps.uncorr, plink.assoc.snps.bonf, plink.assoc.snps.fdr)
    # names(plink.assoc.results) <- c("pval.plink.assoc", "plink.assoc.snps.uncorr", "plink.assoc.snps.bonf", "plink.assoc.snps.fdr")
    #
    # ############################################
    #
    # ########################################################
    # ## association w control for genomic inflation factor ##
    # ########################################################
    # ##--> 1df chi-square test
    #
    # ## perform a basic association analysis on the disease trait for all single SNPs
    # if(Windows == FALSE){
    #   filename <- paste("\\media\\caitiecollins\\88CC9BCECC9BB4C2\\PLINK\\", uniqueID, sep="")
    # }else{
    #   filename <- paste("C:\\PLINK\\", uniqueID, sep="")
    # }
    # command <- paste("plink --bfile ",  filename, " --assoc --adjust --gc --counts --allow-no-sex --out ", filename, sep="")
    #
    # ## Run command
    # if(Windows == FALSE){
    #   system(command)
    # }else{
    #   shell(command)
    # }
    #
    #
    # ## to view the file you created, just read it in with R:
    # if(Windows == FALSE){
    #   filename <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".assoc.adjusted", sep="")
    # }else{
    #   filename <- paste("C:/PLINK/", uniqueID, ".assoc.adjusted", sep="")
    # }
    # plink.res <- read.table(filename, header=TRUE)
    # # head(plink.res)
    # ## NOT SURE WHY, BUT THE "UNADJ" p-values and the "GC" p-values are the same
    # ## in this table (even though the "UNADJ" p-values are not actually the same
    # ## as those in the plink.res from the basic association test above,
    # ## AND, in this case, lambdaGC was 5.12 and the mean chi-squared was 4.63!!!!!)
    #
    #
    # ## get p.vals
    # pval.plink.assoc.gc <- plink.res$GC
    #
    # ## REORDER!!! ##
    # ## NOTE: plink.gc returns results in order of SIGNIFICANCE!
    # pval.plink.assoc.gc <- pval.plink.assoc.gc[order(plink.res$SNP)]
    #
    #
    # ## get sig ##
    #
    # ## p.thresh:
    # p.thresh <- p.value # 0.05 # 0.01 # 0.001 # ??
    #
    #
    # ## Uncorrected ##
    # plink.assoc.gc.snps.uncorr <- snps.names[which(pval.plink.assoc.gc < p.thresh)]
    #
    # ## Bonferonni ##
    # p.vals.bonf <- p.adjust(pval.plink.assoc.gc, "bonferroni")
    # p.bonf <- which(p.vals.bonf < p.thresh)
    # plink.assoc.gc.snps.bonf <- snps.names[p.bonf]
    #
    # ## FDR ##
    # p.vals.fdr <- p.adjust(pval.plink.assoc.gc, "fdr")
    # p.fdr <- which(p.vals.fdr < p.thresh)
    # plink.assoc.gc.snps.fdr <- snps.names[p.fdr]
    #
    # ## CONVERT 0-LENGTH RESULTS TO NULL
    # if(length(plink.assoc.snps.uncorr) == 0) plink.assoc.gc.snps.uncorr <- NULL
    # if(length(plink.assoc.snps.bonf) == 0) plink.assoc.gc.snps.bonf <- NULL
    # if(length(plink.assoc.snps.fdr) == 0) plink.assoc.gc.snps.fdr <- NULL
    #
    #
    # ## STORE PLINK TEST RESULTS ##
    # plink.assoc.gc.results <- list(pval.plink.assoc.gc, plink.assoc.gc.snps.uncorr, plink.assoc.gc.snps.bonf, plink.assoc.gc.snps.fdr)
    # names(plink.assoc.gc.results) <- c("pval.plink.assoc.gc", "plink.assoc.gc.snps.uncorr", "plink.assoc.gc.snps.bonf", "plink.assoc.gc.snps.fdr")
    #
    #
    # ############################################
    #
    # ## STORE COMBINED PLINK RESULTS ##
    # plink.results <- list(plink.assoc.results,
    #                       plink.assoc.gc.results)


    ###########################################################################################################################
    ################################################# ***   PCA   *** #########################################################
    ###########################################################################################################################

    #########
    ## PCA ##
    #########

    ## STEPS: ##
    ## (1) Run PCA
    ## (2) Select sig. number of PCs (--> HOW??)
    ## (3) Regress snps along sig residuals
    ## (4) Run assoc test on adjusted snps dataset.

    ## (Below: taken from Glasgow/practical/practical-GWAS_before_cuts.Rnw ~ practical-GWAS_day4.pdf)

    print("PCA started")

    ## get snps:
    snps <- snps.ori.ori
    phen <- phen.ori.ori
    if(is.null(colnames(snps))) colnames(snps) <- c(1:ncol(snps))

    ## run PCA: ##
    ## Keep only n.PCs significant axes:
    n.PCs <- 5
    # n.PCs <- 15
    pca1 <- dudi.pca(snps, scale=FALSE, scannf=FALSE, nf=n.PCs)

    ## Identify main pop clusters: ##
    set.seed(seed)
    grp <- find.clusters(snps, n.pca=n.PCs, choose.n.clust=FALSE, max.n.clust=(n.PCs + 1))
    # grp <- find.clusters(snps, choose.n.clust=FALSE, max.n.clust=(n.PCs + 1), pca.select="percVar", perc.pca=60)
    pop <- grp$grp # gives same result as cutree(clust, k=6)
    n.grp <- length(levels(pop))

    ## CHECK-- ensure no populations contain only 1 individual:
    ## (UNLESS highest tree division leads to a single-node clade?)
    # snps.mat <- matrix(as.character(snps), nrow=nrow(snps), ncol=ncol(snps), dimnames=list(rownames(snps), colnames(snps)))
    # snps.mat <- replace(snps.mat, which(snps.mat == "0"), "a")
    # snps.mat <- replace(snps.mat, which(snps.mat == "1"), "t")
    # snps.dna <- as.DNAbin(snps.mat)
    # clust <- hclust(dist.dna(snps.dna))
    # tab <- table(cutree(clust, k=2))
    # ## If 1 of 2 smallest pops contains only 1 ind, CMH will fail...
    # if(any(tab == 1)) cmh.fails <- TRUE

    cmh.fails <- FALSE
    ## (1) reduce max.n.clust by 1:
    if(any(table(pop) < 2)){
      max.k <- n.PCs
      set.seed(seed)
      grp <- find.clusters(snps, n.pca=n.PCs, choose.n.clust=FALSE, max.n.clust=max.k) # pca.select="percVar", perc.pca=60,
      pop <- grp$grp # gives same result as cutree(clust, k=6)
      n.grp <- length(levels(pop))

      ## (2) repeat by re-setting seed with same reduced max.n.clust
      seed.new <- seed+1
      counter <- 0
      while(any(table(pop) < 2)){
        if(counter > 5){
          max.k <- max.k - 1
          counter <- 0
        }
        set.seed(seed.new)
        grp <- find.clusters(snps, n.pca=n.PCs, choose.n.clust=FALSE, max.n.clust=max.k) # pca.select="percVar", perc.pca=60,
        pop <- grp$grp # gives same result as cutree(clust, k=6)
        n.grp <- length(levels(pop))
        seed.new <- seed.new+1
        counter <- counter+1
      } # end while
    } # end pop check
    if(n.grp == 1) cmh.fails <- TRUE


    ########################################################
    ## Correct for pop strat by regressing along sig PCs: ##
    ########################################################

    ## NOTE-- if any MISSING DATA contained in SNPs dataset, must REPLACE it here (eg. with the mean for that snps column).

    ## get formula (get string up to n.PCs): lm(e ~ pca1$li[,1] + pca1$li[,2] + pca1$li[,3] + pca1$li[,4] + pca1$li[,5])
    PC.string <- sapply(c(1:n.PCs), function(e) paste("pca1$li[, ", e, "]", sep=""))
    PC.string <- paste0(PC.string, collapse=" + ")
    var.string <- paste("e ~ ", PC.string)

    ## Correct SNPs w PCA!
    ## SLOW STEP..!
    snps.corrected <- apply(snps, 2, function(e) residuals(do.call(lm, list(as.formula(var.string))))) # may take a minute

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

    if(!is.numeric(phen)) phen <- as.numeric(phen)

    ## Get UNIQUE snps.corrected ##
    snps.corrected.ori <- snps.corrected
    temp <- get.unique.matrix(snps.corrected, MARGIN=2)
    temp.unique <- temp$unique.data
    index <- temp$index

    if(ncol(temp.unique) == ncol(snps.corrected.ori)){
      all.unique <- TRUE
    }else{
      all.unique <- FALSE
    }

    ## work w only unique snps:
    snps.corrected <- temp.unique

    ## SLOW STEP..
    pval2 <- numeric(0)
    # system.time( # 120.78
    for(n in 1:ncol(snps.corrected)){
      foo <- suppressWarnings(glm(phen ~ snps.corrected[,n], family="binomial"))
      ANOVA <- anova(foo, test="Chisq")
      pval2[n] <- ANOVA$"Pr(>Chi)"[2]
    } # end for loop
    # )


    ## Get all non-unique values:
    snps.corrected <- snps.corrected.ori
    if(all.unique == FALSE) pval2 <- pval2[index]

    ## Store pvals and snps.corrected from pca-corrected ANOVA association test:
    pval.pca <- pval2
    snps.corrected.pca <- snps.corrected

    ## Get results:
    p.thresh <- p.value # 0.01
    p.vals.bonf <- p.adjust(pval.pca, "bonferroni")
    p.bonf <- which(p.vals.bonf < p.thresh)
    pca.snps.bonf <- snps.names[p.bonf]

    ## Store results:
    pca.results <- list(snps.corrected.pca, pval.pca, pca.snps.bonf)
    names(pca.results) <- c("snps.corrected.pca", "pval.pca", "pca.snps.bonf")

    print("PCA done")

    ###########################################################################################################################
    ################################################ ***   DAPC   *** #########################################################
    ###########################################################################################################################

    print("DAPC started")

    ## Using our pop clusters as the group factor in DAPC,
    ## we can generate a new DAPC object ...
    ## after performing cross-validation to optimise the discrimination between these subpopulations:

    ## NEW: make sure we get as many PCs  as DAs?
    n.pca.min <- (n.grp)
    n.pca.max <- min(min(dim(snps))*0.5, pca1$rank)
    runs <- 10
    n.pca <- round(pretty(n.pca.min:n.pca.max, runs))

    ## NOTE -- If this is SLOW and always results in 5-10 PCs, we may want to skip it... ??
    xval.pop <- xvalDapc(snps, pop, n.pca=n.pca)
    # str(xval.pop)
    # xval.pop[2:6]



    ## store DAPC object:
    dapc.pop <- xval.pop$DAPC

    ###################################
    ## Correct for pop strat w DAPC: ##
    ###################################

    ## As we did when correcting with PCA, we regress along the axes of DAPC.
    ## When correcting with the DAPC approach, we do not need to determine
    ## how many axes are ``significant'': we will always correct with (k - 1) axes.
    ## NOTE that (k - 1) is not necessarily = n.PCs
    ## ... though we may want to limit the K we work with depending on how system.time scales with K... (?)

    ## get formula (get string up to n.PCs):
    ## lm(e ~ dapc.pop$ind.coord[,1] + dapc.pop$ind.coord[,2] + dapc.pop$ind.coord[,3] + dapc.pop$ind.coord[,4])
    # DAPC.string <- sapply(c(1:n.PCs), function(e) paste("dapc.pop$ind.coord[, ", e, "]", sep=""))
    DAPC.string <- sapply(c(1:(n.grp - 1)), function(e) paste("dapc.pop$ind.coord[, ", e, "]", sep=""))
    DAPC.string <- paste0(DAPC.string, collapse=" + ")
    var.string <- paste("e ~ ", DAPC.string)

    ## Correct SNPs w DAPC!
    snps.corrected <- apply(snps, 2, function(e) residuals(do.call(lm, list(as.formula(var.string))))) # may take a minute


    ## Get UNIQUE snps.corrected ##
    snps.corrected.ori <- snps.corrected
    temp <- get.unique.matrix(snps.corrected, MARGIN=2)
    temp.unique <- temp$unique.data
    index <- temp$index

    if(ncol(temp.unique) == ncol(snps.corrected.ori)){
      all.unique <- TRUE
    }else{
      all.unique <- FALSE
    }

    ## work w only unique snps:
    snps.corrected <- temp.unique


    ## Run association test on DAPC-corrected snps:
    pval3 <- numeric(0)
    system.time(
      for(n in 1:ncol(snps.corrected)){
        foo <- suppressWarnings(glm(phen ~ snps.corrected[,n], family="binomial"))
        ANOVA <- anova(foo, test="Chisq")
        pval3[n] <- ANOVA$"Pr(>Chi)"[2]
      } # end for loop
    )


    ## Get all non-unique values:
    snps.corrected <- snps.corrected.ori
    if(all.unique == FALSE) pval3 <- pval3[index]


    ## store pvals and snps.corrected for DAPC-corrected assoc test:
    pval.dapc <- pval3
    snps.corrected.dapc <- snps.corrected

    ## Get results:
    p.thresh <- p.value # 0.01
    p.vals.bonf <- p.adjust(pval.dapc, "bonferroni")
    p.bonf <- which(p.vals.bonf < p.thresh)
    dapc.snps.bonf <- snps.names[p.bonf]

    ## Store results:
    dapc.results <- list(snps.corrected.dapc, pval.dapc, dapc.snps.bonf)
    names(dapc.results) <- c("snps.corrected.dapc", "pval.dapc", "dapc.snps.bonf")

    print("DAPC done")

    ###########################################################################################################################
    ################################################ ***   CMH   *** ##########################################################
    ###########################################################################################################################

    ## NB: Only running CMH test if it is possible to get at least 2 pops with > 1 individual in them...
    ## If not --> set results = CMH find 0 sig snps...
    ## (If problem recurs too often, could make a genuine soln w drop.tip etc... )
    if(cmh.fails == FALSE){

      print("CMH started")


      snps <- snps.ori.ori
      phen <- phen.ori.ori

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
      } # end make phen binary..
      phen <- as.numeric(phen)
      names(phen) <- names(phen.ini)

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
      for(n in 1:ncol(mat.unique)){
        tab <- list()
        for(e in 1:length(levels(pop))){
          temp <- ftable(mat.unique[pop==e, n])
          tab[[e]] <- replace(rep(0, 4), which(noms %in% attr(temp, "col.vars")[[1]]), temp)
        } # end (e) loop
        arr.l[[n]] <- do.call(cbind, tab)
      } # end for (i) loop
      arr <- do.call(rbind, arr.l)


      arr.complete <- arr
      ##############
      ## FOR LOOP ##
      ##############
      ## TO GET P-VALUES FROM CMH TEST for EACH SNPs COLUMN:
      p.vals <- list()
      for(n in 1:ncol(mat.unique)){
        ## get indices for this snp for all pops and all 4 2x2 combos:
        from <- seq(1, nrow(arr.complete), 4)[n]
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
        p.vals[[n]] <- CMH$p.value
      } # end for loop
      p.vals <- as.vector(unlist(p.vals))

      ## get full set of p-vals for non-unique columns:
      pval.cmh <- p.vals[index]


      ## Get results:
      p.thresh <- p.value # 0.01
      p.vals.bonf <- p.adjust(pval.cmh, "bonferroni")
      p.bonf <- which(p.vals.bonf < p.thresh)
      cmh.snps.bonf <- snps.names[p.bonf]

      ## Store results:
      cmh.results <- list(pval.cmh, cmh.snps.bonf)
      names(cmh.results) <- c("pval.cmh", "cmh.snps.bonf")

      print("CMH done")

    }else{

      print("CMH SKIPPED (single-node major clade)")

      ## set results to null results?
      # pval.cmh <- rep(1, ncol(snps))
      # cmh.snps.bonf <- snps.names[which(1 < 0)] # character(0)

      ## set results to fisher results:
      pval.cmh <- pval.fisher
      cmh.snps.bonf <- fisher.snps.bonf

      ## Store results:
      cmh.results <- list(pval.cmh, cmh.snps.bonf)
      names(cmh.results) <- c("pval.cmh", "cmh.snps.bonf")

      print("CMH (null) done")
    }

    ###########################################################################################################################
    ############################################# *** PERFORMANCE *** #########################################################
    ###########################################################################################################################

    ##########################
    ## EVALUATE PERFORMANCE ##
    ##########################

    performance <- list()

    ####################
    ## common metrics ##
    ####################
    ## get n.tests
    snps <- snps.ori.ori
    phen <- phen.ori.ori
    n.tests <- dim(snps)[2]

    ###########################################
    ## performance[[1]] contains snps.assoc: ##
    ###########################################
    performance[[1]] <- snps.assoc


    ##############################
    ## FOR LOOP FOR ALL 3 TESTS ##
    ##############################
    for(j in 2:10){

      if(j == 2) test <- "treeWAS.combined"
      if(j == 3) test <- "treeWAS.terminal"
      if(j == 4) test <- "treeWAS.simultaneous"
      if(j == 5) test <- "treeWAS.subsequent"

      if(j == 6) test <- "fisher.bonf"
      if(j == 7) test <- "fisher.fdr"

      if(j == 8) test <- "pca"
      if(j == 9) test <- "dapc"
      if(j == 10) test <- "cmh"


      ######################
      ## treeWAS.combined ##
      ######################
      if(test == "treeWAS.combined"){
        test.positive <- treeWAS.all$treeWAS.combined
      }

      ######################
      ## treeWAS.terminal ##
      ######################
      if(test == "treeWAS.terminal"){
        test.positive <- treeWAS.all$treeWAS$terminal
      }

      ##########################
      ## treeWAS.simultaneous ##
      ##########################
      if(test == "treeWAS.simultaneous"){
        test.positive <- treeWAS.all$treeWAS$simultaneous
      }

      ########################
      ## treeWAS.subsequent ##
      ########################
      if(test == "treeWAS.subsequent"){
        test.positive <- treeWAS.all$treeWAS$subsequent
      }


      ########### ########### ########### ########### ########### ########### ###########

      #################
      ## fisher.bonf ##
      #################
      if(test == "fisher.bonf"){
        ## get test.positive
        test.positive <- fisher.snps.bonf
      } # end test = fisher.bonf
      ################
      ## fisher.fdr ##
      ################
      if(test == "fisher.fdr"){
        test.positive <- fisher.snps.fdr
      } # end test = fisher.fdr


      ################
      ## PCA & DAPC ## ########### ########### ########### ########### ########### ###########
      ################

      #########
      ## PCA ##
      #########
      if(test == "pca"){
        test.positive <- pca.snps.bonf
      }

      ##########
      ## DAPC ##
      ##########
      if(test == "dapc"){
        test.positive <- dapc.snps.bonf
      }

      #########
      ## CMH ##
      #########
      if(test == "cmh"){
        test.positive <- cmh.snps.bonf
      }

      ###########
      ## PLINK ## ########### ########### ########### ########### ########### ###########
      ###########

      # ## Basic association ##
      #
      # ######################
      # ## plink.assoc.bonf ##
      # ######################
      # if(test == "plink.assoc.bonf"){
      #   test.positive <- plink.assoc.snps.bonf
      # } # end test = plink.assoc.bonf
      #
      # #####################
      # ## plink.assoc.fdr ##
      # #####################
      # if(test == "plink.assoc.fdr"){
      #   test.positive <- plink.assoc.snps.fdr
      # } # end test = plink.assoc.fdr
      #
      #
      # ## Corrected w Genomic Control ##
      #
      # #########################
      # ## plink.assoc.gc.bonf ##
      # #########################
      # if(test == "plink.assoc.gc.bonf"){
      #   test.positive <- plink.assoc.gc.snps.bonf
      # } # end test = plink.assoc.gc.bonf
      #
      # ########################
      # ## plink.assoc.gc.fdr ##
      # ########################
      # if(test == "plink.assoc.gc.fdr"){
      #   test.positive <- plink.assoc.gc.snps.fdr
      # } # end test = plink.assoc.gc.fdr


      ########### ########### ########### ########### ########### ########### ###########



      if(is.null(names(snps.assoc))) names(snps.assoc) <- as.character(snps.assoc)
      if(is.null(snps.names)){
        if(is.null(colnames(snps))) colnames(snps) <- as.character(1:ncol(snps))
        snps.names <- colnames(snps)
      }

      #########################
      ## common calculations ##
      #########################
      ## get test.negative
      if(length(which(snps.names %in% test.positive)) != 0){
        test.negative <- snps.names[-which(snps.names %in% test.positive)]
      }else{
        test.negative <- snps.names
      }


      ## get n.test.positive
      n.test.positive <- length(test.positive)
      ## get n.test.negative
      n.test.negative <- length(test.negative) ## == n.tests - n.test.positive

      n.tests <- ncol(snps)

      ########################
      ## GET TP, TN, FP, FN ##
      ########################
      ## get true positives
      snps.associated <- names(snps.assoc)
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

      ## Do NOT need anything to find (no associated SNPs required) ######################################

      ##############
      ## accuracy ##
      ##############
      ## ie. SUMMARY STATISTIC--Of all the CALLS/tests you made, how many of them were CORRECT
      ## ~ Pr(Correct Call | Call)
      # accuracy <- ((TP + TN) / n.tests) ### WHY IS THIS GIVING ME 0.5 (when all other metrics seem to be working....) ??!?!
      accuracy <- ((TP + TN) / (TP + TN + FP + FN))
      # acc <- (sensitivity*length(snps.associated) + specificity*length(snps.not))/ncol(snps)

      #################
      ## specificity ##
      #################
      ## ie. Of all the truly NOT associated SNPs, how many did you manage to rule out?
      ## ~ Pr(Negative Test | SNP NOT associated)
      specificity <- (TN / (TN + FP)) ## = (1 - FPR)

      #########
      ## FPR ##
      #########
      ## ie. How many truly NOT associated SNPs did you accidentally call significant
      ## ~ Pr(Positive Test | SNP NOT associated)
      FPR <- (FP / (FP + TN)) ## = (1 - specificity)


      ## NEED something to FIND, else uninformative! (True ASSOCIATED SNPs ~ required) ###################

      #########
      ## FNR ##
      #########
      ## ie. How many truly ASSOCIATED SNPs did you accidentally miss
      ## ~ Pr(Negative Test | SNP ASSOCIATED)
      ## --> Set 1: will be 0/0 = NaN
      FNR <- (FN / (FN + TP))

      #################
      ## sensitivity ##
      #################
      ## ie. How many truly ASSOCIATED SNPs did you manage to catch
      ## ~ Pr(Positive Test | SNP ASSOCIATED)
      ## --> Set 1: will be 0/0 = NaN
      sensitivity <- (TP / (TP + FN))

      #########
      ## PPV ##
      #########
      ## ie. Of all the POSITIVE calls you made, how many were CORRECT/ identified truly ASSOCIATED SNPs
      ## ~ Pr(SNP ASSOCIATED | Positive Test)
      ## --> Set 1: will be 0 (UNLESS you made NO positive calls, then 0/0 = NaN)
      PPV <- (TP / (TP + FP)) ## = (1 - FDR)

      #########
      ## FDR ##
      #########
      ## ie. Of all the POSITIVE calls you made, how many were WRONG/ identified truly NOT associated SNPs
      ## ~ Pr(SNP NOT associated | Positive Test)
      ## --> Set 1: will be 1 (UNLESS you made NO positive calls, then 0/0 = NaN)
      FDR <- (FP / (FP + TP)) ## = (1 - PPV)


      ##############
      ## F1.score ##
      ##############
      ## Balanced accuracy-like score considering both sensitivity and PPV:
      F1.score <- 2*((sensitivity*PPV) / (sensitivity+PPV))

      ##################################
      ## combine eval metrics into df ##
      ##################################
      performance[[j]] <- data.frame(accuracy, specificity, FPR, FNR, sensitivity, PPV, FDR, F1.score)

    } # end for loop

    names(performance) <- c("snps.assoc",
                            "treeWAS.combined", "treeWAS.terminal", "treeWAS.simultaneous", "treeWAS.subsequent",
                            "fisher.bonf", "fisher.fdr",
                            # "plink.assoc.bonf", "plink.assoc.fdr",
                            # "plink.assoc.gc.bonf", "plink.assoc.gc.fdr",
                            "pca", "dapc", "cmh")

    ################################    ################################    ################################


    ###########################################################################################################################
    ######################################### *** SAVING & RETURNING *** ######################################################
    ###########################################################################################################################


    ########################
    ## SAVE DATA & OUTPUT ##
    ########################
    ## set wd
    if(cluster == FALSE) setwd(wd)
    ## get uniqueID
    uniqueID <- paste("set", set.number, "_", number, sep="")

    ## save snps, phen, tree, out, res, fisher.results, plink.assoc.results, performance
    snps <- snps.ori.ori
    phen <- phen.ori.ori

    ## save snps
    filename.snps[[i]] <- paste("./", uniqueID, "_snps", ".Rdata", sep="")
    save(snps, file=filename.snps[[i]])
    ## save phen
    filename.phen[[i]] <- paste("./", uniqueID, "_phen", ".Rdata", sep="")
    save(phen, file=filename.phen[[i]])
    ## save phen.plot.col
    filename.phen.plot.col[[i]] <- paste("./", uniqueID, "_phen.plot.col", ".Rdata", sep="")
    save(phen.plot.col, file=filename.phen.plot.col[[i]])
    ## save tree
    filename.tree[[i]] <- paste("./", uniqueID, "_tree", ".Rdata", sep="")
    save(tree, file=filename.tree[[i]])
    ## save out
    filename.out[[i]] <- paste("./", uniqueID, "_out", ".Rdata", sep="")
    save(out, file=filename.out[[i]])
    ## save res
    res <- res.complete # includes data from reconstructions, values from treeWAS tests
    filename.res[[i]] <- paste("./", uniqueID, "_res", ".Rdata", sep="")
    save(res, file=filename.res[[i]])
    ## save fisher.results
    filename.fisher.results[[i]] <- paste("./", uniqueID, "_fisher.results", ".Rdata", sep="")
    save(fisher.results, file=filename.fisher.results[[i]])
    ## save plink.assoc.results
    # filename.plink.results[[i]] <- paste("./", uniqueID, "_plink.results", ".Rdata", sep="")
    # save(plink.results, file=filename.plink.results[[i]])

    ## save pca
    filename.pca[[i]] <- paste("./", uniqueID, "_pca", ".Rdata", sep="")
    save(pca.results, file=filename.pca[[i]])
    ## save dapc
    filename.dapc[[i]] <- paste("./", uniqueID, "_dapc", ".Rdata", sep="")
    save(dapc.results, file=filename.dapc[[i]])

    ## save cmh
    filename.cmh[[i]] <- paste("./", uniqueID, "_cmh", ".Rdata", sep="")
    save(cmh.results, file=filename.cmh[[i]])


    ## save performance
    filename.args[[i]] <- paste("./", uniqueID, "_args", ".Rdata", sep="")
    save(args, file=filename.args[[i]])
    ## save performance
    filename.performance[[i]] <- paste("./", uniqueID, "_performance", ".Rdata", sep="")
    save(performance, file=filename.performance[[i]])





    #########################
    ## STORE DATA & OUTPUT ##
    #########################
    SNPS[[i]] <- snps
    names(SNPS)[[i]] <- uniqueID
    PHEN[[i]] <- phen
    names(PHEN)[[i]] <- uniqueID
    PHEN.PLOT.COL[[i]] <- phen.plot.col
    names(PHEN.PLOT.COL[[i]]) <- uniqueID
    TREE[[i]] <- tree
    names(TREE)[[i]] <- uniqueID
    OUT[[i]] <- out
    names(OUT)[[i]] <- uniqueID
    RES[[i]] <- res
    names(RES)[[i]] <- uniqueID
    FISHER.RESULTS[[i]] <- fisher.results
    names(FISHER.RESULTS)[[i]] <- uniqueID
    # PLINK.RESULTS[[i]] <- plink.results
    # names(PLINK.RESULTS)[[i]] <- uniqueID

    PCA.RESULTS[[i]] <- pca.results
    names(PCA.RESULTS)[[i]] <- uniqueID
    DAPC.RESULTS[[i]] <- dapc.results
    names(DAPC.RESULTS)[[i]] <- uniqueID

    CMH.RESULTS[[i]] <- cmh.results
    names(CMH.RESULTS)[[i]] <- uniqueID

    ARGS[[i]] <- args
    names(ARGS)[[i]] <- uniqueID
    PERFORMANCE[[i]] <- performance
    names(PERFORMANCE)[[i]] <- uniqueID

  } # end for loop




  ##########################
  ## RETURN DATA & OUTPUT ##
  ##########################

  toReturn <- list(SNPS, PHEN, PHEN.PLOT.COL, TREE, RES,
                   FISHER.RESULTS, PCA.RESULTS, DAPC.RESULTS, CMH.RESULTS,
                   ARGS, PERFORMANCE)
  names(toReturn) <- c("snps", "phen", "phen.plot.col", "tree", "res",
                       "fisher.results", "pca.results", "dapc.results", "cmh.results",
                       "arguments", "performance")

  ################################
  ## SAVE PERFORMANCE & OUTPUT: ##
  ################################
  ## get timestamp for filename:
  from <- number-length(PERFORMANCE)+1
  to <- number
  time <- Sys.time()
  t1 <- keepFirstN(time, 10)
  t2 <- strsplit(keepLastN(keepFirstN(time, 19), 8), split =":")
  names(t2) <- NULL
  t2 <- t2[[1]]
  t2 <- paste(t2, collapse="_", sep="")
  timestamp <- paste(t1, t2, sep="_")

  ## save PERFORMANCE:
  nom1 <- paste("./", "set", set.number, "_", from, "_",  to, "_PERFORMANCE", "_", timestamp, ".Rdata", sep="")
  save(PERFORMANCE, file=nom1)

  ## save all OUTPUT:
  output <- toReturn
  nom2 <- paste("./", "set", set.number, "_", from, "_",  to, "_OUTPUT", "_", timestamp, ".Rdata", sep="")
  save(output, file=nom2)

  ## return output:
  return(toReturn)

} # end simTest (generic) function
