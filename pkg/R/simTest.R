


## GENERIC SIM TESTING FUNCTION ##


#############
## simTest ##
#############




# out <- simTest(
#
#   ## simTest args:
#   set.number = 2,
#   n.reps = 1,
#   set.seed.as = "file.number",
#   working.dir = "/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims",
#
#   ## data from file args:
#   from.file = FALSE,
#   file.n =NULL,
#
#   ## coalescent.sim args:
#   n.ind = 100,
#   n.snps = 10000, # gen.size
#   # sim.by = "locus",
#   n.subs = 1, # theta (*2)
#   n.phen.subs = 15, # theta_p = NULL
#   n.snps.assoc = 10, # = 0
#   # assoc.option = "all",
#   assoc.prob = 100, # 90,
#   grp.min = 0.25,
#
#   ## treeWAS args:
#   ## RUNNING ALL OF THESE OPTIONS (FOR NOW):
#   p.value = 0.05, # REQUIRED FOR FISHER TEST
#   #   p.value.correct = c("bonf", "fdr", FALSE), #mt.correct = FALSE
#   #   p.value.by = c("count", "density"),
#   sim.n.snps = 100000, # 10*n.snps #sim.gen.size = NULL
#   treeWAS.test = c("terminal", "simultaneous", "subsequent"), # "score"
#   snps.reconstruction = "parsimony",
#   phen.reconstruction = "parsimony"
# )


## old simTest args:
# set.number = 1,
# n.reps = 1, set.seed.as = "file.number",
# p.value = 0.0001, n.phen.subs=15,
# n.snps.assoc=NULL, assoc.prob = 90



## new args:

# ## simTest args:
# set.number = 1
# n.reps = 1
# set.seed.as = "file.number"
# working.dir = "/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims"
#
# ## data from file args:
# from.file = FALSE
# file.n = NULL
#
# ## coalescent.sim args:
# n.ind = 100
# n.snps = 10000 # gen.size
# # sim.by = "locus"
# n.subs = 1 #theta (*2)
# n.phen.subs = 15 #theta_p = NULL
# n.snps.assoc = 10 # = 0
# # assoc.option = "all"
# assoc.prob = 90
# grp.min = 0.25
#
# ## treeWAS args:
# p.value = 0.05
# # p.value.correct = c("bonf", "fdr", FALSE) #mt.correct = FALSE
# # p.value.by = c("count", "density")
# sim.n.snps = 100000 # 10*n.snps #sim.gen.size = NULL
# treeWAS.test = c("terminal", "simultaneous", "subsequent") # "score"
# snps.reconstruction = "parsimony"
# phen.reconstruction = "parsimony"




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
#' @examples
#'
#' ## load data
#' data(dist)
#' str(dist)
#'
#' ## basic use of fn
#' fn(arg1, arg2)
#'
#' #' ## more elaborate use of fn
#' fn(arg1, arg2)
#'
#' @import adegenet ape phangorn
#' @importFrom Hmisc all.is.numeric
#'
#' @export

########################################################################


simTest <- function(

  ## simTest args:
  set.number = 1,
  n.reps = 1,
  set.seed.as = "file.number",
  working.dir = "/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims",

  ## data from file args:
  from.file = FALSE,
  file.n = NULL,

  ## coalescent.sim args:
  n.ind = 100,
  n.snps = 10000, # gen.size
  # sim.by = "locus",
  n.subs = 1, # theta (*2)
  n.phen.subs = 15, # theta_p = NULL
  n.snps.assoc = 10, # = 0
  # assoc.option = "all",
  assoc.prob = 90, # 100 (set2)
  grp.min = 0.25,

  ## treeWAS args:
  ## RUNNING ALL OF THESE OPTIONS (FOR NOW):
  p.value = 0.05, # REQUIRED FOR FISHER TEST
  #   p.value.correct = c("bonf", "fdr", FALSE), #mt.correct = FALSE
  #   p.value.by = c("count", "density"),
  sim.n.snps = 100000, # 10*n.snps #sim.gen.size = NULL
  treeWAS.test = c("terminal", "simultaneous", "subsequent"), # "score"
  snps.reconstruction = "parsimony",
  phen.reconstruction = "parsimony"

){

  ####################
  ## load packages: ##
  ####################
  # require(adegenet)
  # require(phangorn)
  # require(ape)
  # # require(ade4) #?

  #################################################
  ## set working directory & source required fns ##
  #################################################
  # setwd("C:/Users/Caitlin")
  # source("./adegenet/R/sequences.R") ## need for DNAbin2genind fn

  #   setwd("C:/Users/Caitlin/treeWAS")
  #   source("./misc/coalescent.sim.R")
  #   source("./pkg/R/tree.sim.R")
  #   source("./pkg/R/treeWAS.R")


  ###############################################
  ## make lists in which to store all ###########
  ## data and output from each of n.reps runs: ##
  ###############################################
  SNPS <- PHEN <- PHEN.PLOT.COL <-  TREE <- OUT <- RES <-
    FISHER.RESULTS <- PLINK.RESULTS <-
    ARGS <- PERFORMANCE <- list()
  ## and make lists for saving filenames
  filename.snps <- filename.phen <- filename.phen.plot.col <- filename.tree <-
    filename.out <- filename.res <- filename.fisher.results <-
    filename.plink.results <-
    filename.args <- filename.performance <-
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
    ### ensure n.phen.subs is NULL
    if(!is.null(n.phen.subs)) n.phen.subs <- 15
  } # end set 1

  ###########
  ## SET 2 ##
  ###########
  if(set.number == 2){
    ## ensure n.phen.subs is NOT null
    if(is.null(n.phen.subs)) n.phen.subs <- 15
    if(is.null(assoc.prob)) assoc.prob <- 100
  } # end set 2

  ###########
  ## SET 3 ##
  ###########
  if(set.number == 3){
    ## ensure n.phen.subs is NOT null
    if(is.null(n.phen.subs)) n.phen.subs <- 15
    if(!is.null(phen)) phen <- NULL
    ## ensure ONE associated SNP
    if(is.null(n.snps.assoc)) n.snps.assoc <- 10
    if(is.null(assoc.prob)) assoc.prob <- 90
  } # end set 3

  ###########
  ## SET 4 ##
  ###########
  if(set.number == 4){
    ## ensure n.phen.subs is NOT null
    if(is.null(n.phen.subs)) n.phen.subs <- 15
    if(!is.null(phen)) phen <- NULL
    ## ensure MULTIPLE associated SNPs
    if(is.null(n.snps.assoc)) n.snps.assoc <- 10
    if(is.null(assoc.prob)) assoc.prob <- 90
  } # end set 4



  args <- snps.assoc <- NULL
  sim.n.snps.ori <- sim.n.snps

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

  ##############
  ## FOR LOOP ##
  ##############

  for(i in 1:n.reps){

    wd <- paste(working.dir, "/", "set", set.number, sep="")
    setwd(wd)

    ## get file number:
    if(from.file==FALSE){
      ## get number for group | number of set3_snps already in file:
      number <- (length(grep("_snps", dir("./")))+1)
    }else{
      number <- file.n[i]
    }


    ################
    ## dummy plot ##
    ################
    ## to give user indication of what round of simTest.set3 we are on:
    round.marker <- paste("\n\n\n\n\n\n\n\n\n\n\n\nROUND", i, "\n(file number", number, ")", sep=" ")
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
                            set=set.number,
                            seed = seed)

      ####################################
      ## isolate common elements of foo ##
      ####################################
      snps <- snps.ori <- foo$snps
      if(!is.null(n.snps.assoc)) if(n.snps.assoc > 0){
        snps.assoc <- snps.assoc.ori <- snps.assoc.loci <- foo$snps.assoc
      }else{
        snps.assoc <- NULL
      }
      phen <- phen.ori <- foo$phen
      tree <- tree.ori <- foo$tree
      phen.plot.col <- foo$phen.plot.col

      ## snps names:
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
      filename.tree.plot[[i]] <- paste("./set", set.number, "_", number, "_tree_plot", ".pdf", sep="")
      dev.copy(pdf, file=filename.tree.plot[[i]], width=7, height=11) # , pointsize=1
      dev.off()


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
    filename.plot[[i]] <- list()
    for(t in 1:length(treeWAS.test)){
      ## Save only one plot per test:
      # filename.plot[[i]][[t]] <- paste("./set",
      #                             set.number,
      #                             "_", number,
      #                             "_plot_",
      #                             treeWAS.test[t],
      #                             ".pdf", sep="")

      ## Save both Manhattan and Hist per test:
      filename.plot[[i]][[t]] <- c(## manhattan:
                                  paste("./set",
                                       set.number,
                                       "_", number,
                                       "_plot_manhattan_",
                                       treeWAS.test[t],
                                       ".pdf", sep=""),

                                   ## null.dist:
                                   paste("./set",
                                         set.number,
                                         "_", number,
                                         "_plot_",
                                         treeWAS.test[t],
                                         ".pdf", sep="")
                                   )
    }

    #################
    ## RUN treeWAS ##
    #################
    set.seed(seed)

    syst.time <- system.time( # 341
    out <- treeWAS(snps = snps,
                  phen = phen,
                  n.subs = NULL,
                  tree = tree,
                  dist.dna.model = "JC69",
                  plot.tree = FALSE,
                  test = treeWAS.test,
                  p.value = p.value,
                  p.value.correct = p.value.correct,
                  p.value.by = p.value.by,
                  sim.n.snps = sim.n.snps,
                  n.reps = 1,
                  plot.manhattan = TRUE,
                  plot.null.dist = TRUE,
                  plot.dist = FALSE,
                  snps.assoc = snps.assoc,
                  snps.reconstruction = "parsimony",
                  phen.reconstruction = "parsimony",
                  filename.plot=filename.plot[[i]])
    )


    # dev.copy(pdf, file=filename.plot[[i]], width=7, height=11) # , pointsize=12
    # dev.off()

    ##################################
    ## isolate df of Sig SNPs found ##
    ##################################
    res.complete <- out
    res <- out$res


    ##############
    ## get call ##
    ##############
    ## get arguments inputed to simTest for this round
    #call <- match.call()

    ## TO DO: update arguments to contain actual values used
    ## (ie. return actual value of seed instead of set.seed.as)
    ## (eg. if n.phen.subs was NULL and got set to 25 due to set.number,
    ## report 25)
    args <- list(set.number,
                 seed,
                 n.ind,
                 n.snps,
                 n.subs,
                 n.phen.subs,
                 n.snps.assoc,
                 assoc.prob,
                 sim.n.snps,
                 treeWAS.test
    )

    names(args) <- c("set.number",
                     "seed",
                     "n.ind",
                     "n.snps",
                     "n.subs",
                     "n.phen.subs",
                     "n.snps.assoc",
                     "assoc.prob",
                     "sim.n.snps",
                     "treeWAS.test"
    )





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

    #
    #     ###################
    #     ## RUNNING PLINK ##
    #     ###################
    #
    #     ## Set working directory to run plink program
    #     setwd("/media/caitiecollins/88CC9BCECC9BB4C2/Program Files/plink-1.07-dos")
    #     plink.wd <- getwd()
    #     #setwd(plink.wd)
    #
    #     #################
    #     ## Handle DATA ##
    #     #################
    #
    #     ## STORE ORIGINAL DATA FOR LATER ##
    #     snps.ori <- snps
    #     phen.ori <- phen
    #
    #     #     #######################
    #     #     ## from.file = FALSE ##
    #     #     #######################
    #     #     if(from.file == FALSE){
    #
    #     ## If data was created during THIS round of simTest...
    #
    #     ###########################
    #     ## CONVERT AND SAVE DATA ##
    #     ###########################
    #
    #     ## set working directory for saving
    #     dat.wd <- paste(working.dir, "/set", set.number, "/", sep="")
    #     setwd(dat.wd)
    #
    #     ##################################
    #     ## convert snps --> ped format: ##
    #     ##################################
    #
    #     ## PED FORMAT: ##
    #     ## 6 "MANDATORY" COLUMNS first:
    #     ## FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype
    #     #### ... EXCLUDE OTHER "mandatory" columns w PLINK commands
    #     #### --> (IndividualID, Phenotype)
    #     ## 7-to-p: GENOTYPE COLUMNS
    #     ## NO HEADER ROW (belongs in .MAP file)
    #
    #     ## SNPs can NOT be 0 --> convert from 0/1 to 1/2:
    #     snps <- snps+1
    #     row.names(snps) <- paste("ind", row.names(snps), sep=".")
    #     ## save row and column names for elsewhere:
    #     individualID <- row.names(snps)
    #     loci.names <- colnames(snps)
    #
    #     ## convert to matrix:
    #     snps <- matrix(snps, byrow=FALSE, ncol=ncol(snps))
    #     ## Replace every other column with a copy of the column before it:
    #     ## NOTE--THIS MAY CHANGE (ie. IF WORKING INSTEAD WITH DATA AS MT DNA....)
    #     #snps[,seq(2, ncol(snps), 2)] <- snps[,seq(1, ncol(snps), 2)]
    #     ## DUPLICATE columns:
    #     snps.ori <- snps
    #     snps.new <- matrix(NA, nrow=nrow(snps.ori), ncol=2*ncol(snps.ori))
    #     snps.new[, seq(1, 2*ncol(snps.ori), 2)] <- snps.ori
    #     snps.new[, seq(2, 2*ncol(snps.ori), 2)] <- snps.ori
    #     snps <- snps.new
    #
    #     ## recode phen: S/A as 0, R/B as 1:
    #     phen <- as.character(phen)
    #     # set1 was coded S/R originally
    #     phen <- replace(phen, which(phen=="S"), 1)
    #     phen <- replace(phen, which(phen=="R"), 2)
    #     # most should be coded A/B
    #     phen <- replace(phen, which(phen=="A"), 1)
    #     phen <- replace(phen, which(phen=="B"), 2)
    #
    #     ## bind individualID, phen, snps
    #     dat <- cbind(individualID, phen, snps)
    #     colnames(dat) <- NULL
    #
    #     ped <- dat
    #
    #     ## get filename
    #     uniqueID <- paste("set", set.number, "_", number, sep="")
    #     filename <- paste("./", uniqueID, "_ped.Rdata", sep="")
    #
    #     ## save dat.ped as Rdata
    #     save(ped, file=filename)
    #     #ped <- get(load("C:/Cait 2012/Work/Xavier/Sims/set3/set3_1_ped.Rdata"))
    #
    #     ## save as text!
    #     #ped <- dat
    #     filename <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".txt", sep="")
    #     write.table(ped, file=filename, quote=FALSE, row.names=FALSE, col.names=FALSE)
    #
    #     ## Do NOT save as .PED
    #     ## convert from text to PED!!
    #
    #     ## UHOH -- STOPPED HERE -- SHELL COMMAND DOESNT SEEM TO BE WORKING FROM LINUX!!!!!!!!!!!!!!!!!!
    #     ## SYSTEM COMMAND NOT READING FILE NAMES CORRECTLY....????????????????///
    #
    #     #     filename.ori <- paste("\\media\\caitiecollins\\88CC9BCECC9BB4C2\\PLINK\\", uniqueID, ".txt", sep="")
    #     #     filename.new <- paste("\\media\\caitiecollins\\88CC9BCECC9BB4C2\\PLINK\\", uniqueID, ".ped", sep="")
    #
    #     filename.ori <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".txt", sep="")
    #     filename.new <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".ped", sep="")
    #
    #     # /media/caitiecollins/88CC9BCECC9BB4C2/PLINK
    #     command <- paste("mv", filename.ori, filename.new, sep=" ")
    #     # shell(command)
    #     system(command)
    #
    #
    #     ####################################
    #     ## convert snp meta-data --> .map ##
    #     ####################################
    #
    #     ## MAP format: ##
    #     ## EACH LINE of a .map file describes a SINGLE marker.
    #     ## Contains 4 "MANDATORY" COLUMNS:
    #     ## Chromosome, rs#/SNP identifier, (Genetic distance), Base-pair position.
    #
    #
    #     ## loci.names:
    #     ## SHOULD BE ONE NAME PER SITE (ie. PER TWO LOCI): ie. L001, NOT L001.1, L001.2 !!!!
    #     ## remove last TWO characters (ie. decimal and trailing digit):
    #     #       loci.names <- substr(loci.names, 1, nchar(loci.names)-2)
    #     #       ## keep only every other:
    #     #       loci.names <- loci.names[seq(1, length(loci.names), 2)]
    #     ## OR: # loci.names <- unique(loci.names)
    #
    #     ## make dummy variables for irrelevant fields:
    #     chromosome <- rep(26, length(loci.names)) # 26 = human mitochondrial (haploid)
    #     gen.dist <- rep(0, length(loci.names))
    #     ## get base-pair posi (ie. loci name - L):
    #     bp <- loci.names
    #     bp <- as.numeric(gsub("L", "", bp))
    #     dat <- data.frame(chromosome, loci.names, gen.dist, bp)
    #
    #     ## as matrix, no header:
    #     dat <- as.matrix(dat, byrow=FALSE, ncol=ncol(dat))
    #     colnames(dat) <- NULL
    #
    #     map <- dat
    #
    #     ## get filename
    #     uniqueID <- paste("set", set.number, "_", number, sep="")
    #     filename <- paste("./", uniqueID, "_map.Rdata", sep="")
    #
    #     ## save dat.map as Rdata
    #     save(map, file=filename)
    #     #map <- get(load("C:/Cait 2012/Work/Xavier/Sims/set3/set3_1_map.Rdata"))
    #
    #     ## save as text!
    #     filename <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".txt", sep="")
    #     write.table(map, file=filename, quote=FALSE, row.names=FALSE, col.names=FALSE)
    #
    #     ## Do NOT save as .MAP
    #     ## convert from text to MAP!!
    #
    #     #     filename.ori <- paste("C:\\PLINK\\", uniqueID, ".txt", sep="")
    #     #     filename.new <- paste("C:\\PLINK\\", uniqueID, ".map", sep="")
    #
    #     filename.ori <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".txt", sep="")
    #     filename.new <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, ".map", sep="")
    #     command <- paste("mv", filename.ori, filename.new, sep=" ")
    #     # shell(command)
    #     system(command)
    #
    #     ###### ###### ###### ######
    #
    #     #     }else{ # end from.file = FALSE
    #     #
    #     #       ######################
    #     #       ## from.file = TRUE ##
    #     #       ######################
    #     #
    #     #     } # end from.file = TRUE
    #
    #
    #     ########## #################### ########## #################### #################### ########## ####################
    #
    #
    #     #####################
    #     ## GWAS with PLINK ##
    #     #####################
    #
    #     ## set wd for PLINK program
    #     setwd(plink.wd)
    #
    #     ## get filename
    #     # filename <- paste("C:\\PLINK\\", uniqueID, sep="")
    #     filename <- paste("/media/caitiecollins/88CC9BCECC9BB4C2/PLINK/", uniqueID, sep="")
    #
    #
    #     ## inspect file?
    #     command <- paste("plink --file ", filename, " --no-fid --no-parents --no-sex --allow-no-sex", sep="")
    #     #shell(command)
    #
    #     ## make a binary PED file ##
    #     ## (provide the full path, not just the file name)
    #     command <- paste("plink --file ", filename, " --make-bed --no-fid --no-parents --no-sex --allow-no-sex --out ", filename, sep="")
    #     # shell(command)
    #     system(command)
    #     ## SYSTEM(COMMAND) FAILS HERE WITH THE FOLLOWING ERROR #################################################    ####    ####    ####    ####    ####    ####
    #
    #     # sh: 1: plink: not found
    #
    #     ## PROBABLY NEED TO FIND SOLN AND APPLY IT TO ALL PLINK-RELATED CODE, STARTING WELL ABOVE HERE....................................................................
    #     ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####    ####
    #     ## use --bfile to work with the BINARY file
    #     # (same as --file, but loads the binary one and prints summary stats)
    #     command <- paste("plink --bfile ", filename, " --no-fid --no-parents --no-sex --allow-no-sex", sep="")
    #     #shell(command)
    #
    #     #######################
    #     ## basic association ##
    #     #######################
    #     ##--> 1df chi-square test
    #
    #     ## check freq of SNPs...?
    #     command <- paste("plink --file ", filename, " --no-fid --no-parents --no-sex --allow-no-sex --freq --out ", filename, sep="")
    #     shell(command)
    #     ## yay!
    #     filename <- paste("C:/PLINK/", uniqueID, ".frq", sep="")
    #     freq <- read.table(filename, header=TRUE)
    #     #head(freq)
    #
    #     ## perform a basic association analysis on the disease trait for all single SNPs
    #     filename <- paste("C:\\PLINK\\", uniqueID, sep="")
    #     command <- paste("plink --bfile ",  filename, " --assoc --counts --allow-no-sex --out ", filename, sep="")
    #     shell(command)
    #
    #     ## to view the file you created, just read it in with R:
    #     filename <- paste("C:/PLINK/", uniqueID, ".assoc", sep="")
    #     plink.res <- read.table(filename, header=TRUE)
    #     #head(plink.res)
    #
    #
    #     ## get p.vals
    #     pval.plink.assoc <- plink.res$P
    #
    #     ## get sig ##
    #
    #     ## p.thresh:
    #     p.thresh <- p.value # 0.05 # 0.01 # 0.001 # ??
    #
    #
    #     ## Uncorrected ##
    #     plink.assoc.snps.uncorr <- snps.names[which(pval.plink.assoc < p.thresh)]
    #
    #     ## Bonferonni ##
    #     p.vals.bonf <- p.adjust(pval.plink.assoc, "bonferroni")
    #     p.bonf <- which(p.vals.bonf < p.thresh)
    #     plink.assoc.snps.bonf <- snps.names[p.bonf]
    #
    #     ## FDR ##
    #     p.vals.fdr <- p.adjust(pval.plink.assoc, "fdr")
    #     p.fdr <- which(p.vals.fdr < p.thresh)
    #     plink.assoc.snps.fdr <- snps.names[p.fdr]
    #
    #     ## CONVERT 0-LENGTH RESULTS TO NULL
    #     if(length(plink.assoc.snps.uncorr) == 0) plink.assoc.snps.uncorr <- NULL
    #     if(length(plink.assoc.snps.bonf) == 0) plink.assoc.snps.bonf <- NULL
    #     if(length(plink.assoc.snps.fdr) == 0) plink.assoc.snps.fdr <- NULL
    #
    #     ## STORE PLINK TEST RESULTS ##
    #     plink.assoc.results <- list(pval.plink.assoc, plink.assoc.snps.uncorr, plink.assoc.snps.bonf, plink.assoc.snps.fdr)
    #     names(plink.assoc.results) <- c("pval.plink.assoc", "plink.assoc.snps.uncorr", "plink.assoc.snps.bonf", "plink.assoc.snps.fdr")
    #
    #     ############################################
    #
    #     ########################################################
    #     ## association w control for genomic inflation factor ##
    #     ########################################################
    #     ##--> 1df chi-square test
    #
    #     ## perform a basic association analysis on the disease trait for all single SNPs
    #     filename <- paste("C:\\PLINK\\", uniqueID, sep="")
    #     command <- paste("plink --bfile ",  filename, " --assoc --adjust --gc --counts --allow-no-sex --out ", filename, sep="")
    #     shell(command)
    #
    #     ## to view the file you created, just read it in with R:
    #     filename <- paste("C:/PLINK/", uniqueID, ".assoc.adjusted", sep="")
    #     plink.res <- read.table(filename, header=TRUE)
    #     #head(plink.res)
    #     ## NOT SURE WHY, BUT THE "UNADJ" p-values and the "GC" p-values are the same
    #     ## in this table (even though the "UNADJ" p-values are not actually the same
    #     ## as those in the plink.res from the basic association test above,
    #     ## AND, in this case, lambdaGC was 5.12 and the mean chi-squared was 4.63!!!!!)
    #
    #
    #     ## get p.vals
    #     pval.plink.assoc.gc <- plink.res$P
    #
    #     ## get sig ##
    #
    #     ## p.thresh:
    #     p.thresh <- p.value # 0.05 # 0.01 # 0.001 # ??
    #
    #
    #     ## Uncorrected ##
    #     plink.assoc.gc.snps.uncorr <- snps.names[which(pval.plink.assoc.gc < p.thresh)]
    #
    #     ## Bonferonni ##
    #     p.vals.bonf <- p.adjust(pval.plink.assoc.gc, "bonferroni")
    #     p.bonf <- which(p.vals.bonf < p.thresh)
    #     plink.assoc.gc.snps.bonf <- snps.names[p.bonf]
    #
    #     ## FDR ##
    #     p.vals.fdr <- p.adjust(pval.plink.assoc.gc, "fdr")
    #     p.fdr <- which(p.vals.fdr < p.thresh)
    #     plink.assoc.gc.snps.fdr <- snps.names[p.fdr]
    #
    #     ## CONVERT 0-LENGTH RESULTS TO NULL
    #     if(length(plink.assoc.snps.uncorr) == 0) plink.assoc.gc.snps.uncorr <- NULL
    #     if(length(plink.assoc.snps.bonf) == 0) plink.assoc.gc.snps.bonf <- NULL
    #     if(length(plink.assoc.snps.fdr) == 0) plink.assoc.gc.snps.fdr <- NULL
    #
    #
    #     ## STORE PLINK TEST RESULTS ##
    #     plink.assoc.gc.results <- list(pval.plink.assoc.gc, plink.assoc.gc.snps.uncorr, plink.assoc.gc.snps.bonf, plink.assoc.gc.snps.fdr)
    #     names(plink.assoc.gc.results) <- c("pval.plink.assoc.gc", "plink.assoc.gc.snps.uncorr", "plink.assoc.gc.snps.bonf", "plink.assoc.gc.snps.fdr")
    #
    #
    #     ############################################
    #
    #     ## STORE COMBINED PLINK RESULTS ##
    #     plink.results <- list(plink.assoc.results,
    #                           plink.assoc.gc.results)

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
    n.tests <- dim(snps)[2]

    ###########################################
    ## performance[[1]] contains snps.assoc: ##
    ###########################################
    performance[[1]] <- snps.assoc


    ##############################
    ## FOR LOOP FOR ALL 3 TESTS ##
    ##############################
    for(j in 2:99){


      if(j==2) test <- "fisher.bonf"
      if(j==3) test <- "fisher.fdr"

      # if(j==5) test <- "plink.assoc.bonf"
      # if(j==6) test <- "plink.assoc.fdr"
      # if(j==7) test <- "plink.assoc.gc.bonf"
      # if(j==8) test <- "plink.assoc.gc.fdr"

      if(j %in% 4:99) test <- "treeWAS"

      ## get test run:
      if(j %in% 4:35)  t <- "terminal"
      if(j %in% 36:67) t <- "simultaneous"
      if(j %in% 68:99) t <- "subsequent"

      # ###########
      # ## PLINK ## ########### ########### ########### ########### ########### ###########
      # ###########
      #
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
      # if(test == "plink.assoc.bonf"){
      #   test.positive <- plink.assoc.gc.snps.bonf
      # } # end test = plink.assoc.bonf
      #
      # ########################
      # ## plink.assoc.gc.fdr ##
      # ########################
      # if(test == "plink.assoc.fdr"){
      #   test.positive <- plink.assoc.gc.snps.fdr
      # } # end test = plink.assoc.fdr
      #
      # ########### ########### ########### ########### ########### ########### ###########



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


      ########### ########### ########### ########### ########### ########### ###########

      #############
      ## treeWAS ##
      #############
      if(test == "treeWAS"){

        if(j %in% 4:35)  N <- 3
        if(j %in% 36:67) N <- 35
        if(j %in% 68:99) N <- 67

        ## get test.positive
        if(class(res[[t]][[(j-N)]]$sig.snps)=="data.frame"){
          # test.positive <- as.character(res$SNP.locus)
          test.positive <- as.character(res[[t]][[(j-N)]]$sig.snps$SNP.locus)
        }else{
          test.positive <- NULL
        }
      } # end test = treeWAS

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
      ## for Set 3 there is 1 associated SNP, so 1 true positives...
      # if(set.number < 3){
      #   snps.associated <- NULL
      # }else{
      #   #         snps.associated <- paste(sapply(c(1:length(snps.assoc)),
      #   #                                         function(e)
      #   #                                           rep(names(snps.assoc)[e], 2)), c(1, 2), sep=".")
        snps.associated <- names(snps.assoc)
      # }
      true.positive <- test.positive[which(test.positive %in% snps.associated)]
      TP <- length(true.positive)

      ## get true negatives
      ## for Set 3 all but ONE SNPs are NOT (intentionally) associated with the phenotype
      # if(set.number < 3){
      #   snps.not <- snps.names
      # }else{
        snps.not <- snps.names[-which(snps.names %in% snps.associated)]
      # }
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
      accuracy <- ((TP + TN) / n.tests)

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


      ##

      ##################################
      ## combine eval metrics into df ##
      ##################################
      performance[[j]] <- data.frame(accuracy, specificity, FPR, FNR, sensitivity, PPV, FDR)

    } # end for loop

    ## get names for treeWAS tests:
    treeWAS.names <- list()

    for(r in 1:length(res)){
      if(i == 1) t <- "terminal"
      if(i == 2) t <- "simultaneous"
      if(i == 3) t <- "subsequent"

      treeWAS.names[[r]] <- paste("treeWAS", t, names(res[[t]]), sep=".")
    }

    treeWAS.names <- as.vector(unlist(treeWAS.names))

    names(performance) <- c("snps.assoc",
                            "fisher.bonf", "fisher.fdr",
                            # "plink.assoc.bonf", "plink.assoc.fdr",
                            # "plink.assoc.gc.bonf", "plink.assoc.gc.fdr"
                            treeWAS.names)

    ################################    ################################    ################################


    ###########################################################################################################################
    ######################################### *** SAVING & RETURNING *** ######################################################
    ###########################################################################################################################


    ########################
    ## SAVE DATA & OUTPUT ##
    ########################
    ## set wd
    setwd(wd)
    ## get uniqueID
    uniqueID <- paste("set", set.number, "_", number, sep="")

    ## save snps, phen, tree, out, res, fisher.results, plink.assoc.results, performance

    ## save snps
    filename.snps[[i]] <- paste("./", uniqueID, "_snps", ".Rdata", sep="")
    save(snps, file=filename.snps[[i]])
    ## save phen
    filename.phen[[i]] <- paste("./", uniqueID, "_phen", ".Rdata", sep="")
    save(phen, file=filename.phen[[i]])
    ## save phen.plot.col
    filename.phen.plot.col[[i]] <- paste("./", uniqueID, "_phen", ".Rdata", sep="")
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
    ARGS[[i]] <- args
    names(ARGS)[[i]] <- uniqueID
    PERFORMANCE[[i]] <- performance
    names(PERFORMANCE)[[i]] <- uniqueID


  } # end for loop


  ##########################
  ## RETURN DATA & OUTPUT ##
  ##########################

  toReturn <- list(SNPS, PHEN, PHEN.PLOT.COL, TREE, RES, FISHER.RESULTS, ARGS, PERFORMANCE) #  PLINK.ASSOC.RESULTS
  names(toReturn) <- c("snps", "phen", "phen.plot.col", "tree", "res", "fisher.results", "arguments", "performance") # "plink.assoc.results"

  return(toReturn)

} # end simTest (generic) function
