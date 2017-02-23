
#############
## treeWAS ##
#############


## Re-implementation of GWAS method developed in Sheppard et al 2014,
## based on the generation of a phylogenetically-correct p-value
## derived by comparing correlation btw SNPs
## and a phenotype of interest to a null distribution
## (Poisson dist w parameter 1) OR
## (modification based on similar method developed in Farhat et al 2013)
## based on permutation of empirically-derived
## n.mts per site (except we determine this with the Fitch algorithm).

############
## TO DO: ##
############
## 1) Implement an internal protocol to get the n.subs distribution
## from the data either by
## (A) linking to/calling ClonalFrameML,
## (B) using the code used in ClonalFrameML without calling the program,
## (C) writing separate code using R fns (available or self-generated)
## 2)


########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
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

################################
## EXAMPLE (data/parameters): ##
################################

##############
## DATA #1: ##
##############
# data("snps.ace")
# data("phen.ace")
# data("tree.ace")
#
# snps <- snps.ori <- snps.ace
# phen <- phen.ori <- phen.ace
# tree <- tree.ori <- tree.ace


##############
## DATA #2: ##
##############
# foo <- coalescent.sim(n.ind = 100,
#                       n.snps = 10000,
#                       n.subs = 1,
#                       n.snps.assoc = 10,
#                       assoc.prob = 90,
#                       n.phen.subs = 15,
#                       phen = NULL,
#                       plot = TRUE,
#                       heatmap = FALSE,
#                       reconstruct = FALSE,
#                       dist.dna.model = "JC69",
#                       row.names = NULL,
#                       grp.min = 0.25,
#                       seed = 4)
#
# snps <- snps.ori <- foo$snps
# snps.assoc.loci <- snps.assoc <- foo$snps.assoc
# phen <- phen.ori <- foo$phen
# tree <- tree.ori <- foo$tree

#################
## PARAMETERS: ##
#################

# n.subs <- NULL
# dist.dna.model <- "JC69"
# plot.tree <- FALSE
# test <- c("terminal", "simultaneous", "subsequent")
# p.value <- 0.001
# p.value.correct <- "fdr"
# p.value.by <- "count"
# sim.n.snps <- ncol(snps)*10
# n.reps <- 1
# plot.manhattan <- TRUE
# plot.null.dist <- TRUE
# plot.dist <- FALSE
# snps.reconstruction <- "parsimony"
# phen.reconstruction <- "parsimony"

treeWAS <- function(snps,
                    phen,
                    n.subs = NULL,
                    tree = c("UPGMA", "nj", "ml"),
                    dist.dna.model = "JC69",
                    plot.tree = FALSE,
                    test = c("terminal", "simultaneous", "subsequent"),
                    p.value = 0.001,
                    p.value.correct = c("bonf", "fdr", FALSE), ## DO WE WANT TO ALLOW USERS TO RUN MANY DIFFERENT MULTIPLE TESTING CORRECTION METHODS FOR EACH TEST?????????
                    p.value.by = c("count", "density"),
                    sim.n.snps = ncol(snps),
                    n.reps = 1,
                    plot.manhattan = TRUE,
                    plot.null.dist = TRUE,
                    plot.dist = FALSE,
                    snps.assoc = NULL, # for (manhattan) plot
                    snps.reconstruction = "parsimony",
                    phen.reconstruction = "parsimony",
                    filename.plot = NULL){

  ###################
  ## LOAD PACKAGES ##
  ###################
  # require(adegenet)
  # require(phangorn)
  # require(ape)
  # # require(ade4) #?
  # require(Hmisc) # all.is.numeric

  #####################################################################
  ## 0) HANDLE INPUT DATA #############################################
  #####################################################################

  #####################
  ## HANDLE TEST ARG ##
  #####################

  ## Allow partial matching of argument names:
  test <- match.arg(arg = test,
                    choices = c("terminal", "simultaneous", "subsequent", "cor", "fisher"),
                    several.ok = TRUE)

  ########################
  ## HANDLE SNPS & PHEN ##
  ########################
  if(!is.matrix(snps)) snps <- as.matrix(snps)
  # x <- snps
  n.snps <- ncol(snps)

  ## convert phenotype to factor
  phen <- as.factor(phen)
  # y <- phen

  ## set n.ind:
  n.ind <- length(phen)
  inds <- c(1:n.ind)

  #################
  ## HANDLE TREE ##
  #################

  ## RECONSTRUCTED TREE ##

  if(class(tree) == "character"){

    tree <- tolower(tree)

    if(!any(c("upgma", "nj", "ml") %in% tree)){
      warning("If tree is not a phylo object,
              it should be one of 'UPGMA', 'NJ', 'ML',
              specifying which method is to be used to
              reconstruct the phylogenetic tree from the snps provided.
              Choosing 'UPGMA' by default.")
      tree <- "upgma"
    }
    tree <- tree.reconstruct(snps,
                             method = tree,
                             dist.dna.model = dist.dna.model,
                             plot = plot.tree)
  }else{

    ## USER-PROVIDED TREE ##

    ## If user has already submitted a tree as input:
    ## Work with a centered phylo tree for
    ## consistency and visualisation's sake:
    if(class(tree) != "phylo") tree <- as.phylo(tree)
    ## if the tree is not already rooted, root it:
    if(!is.rooted(tree)) tree <- midpoint(tree)

    ## HANDLE TREE: ##
    ## Always work with trees in "pruningwise" order:
    tree <- reorder.phylo(tree, order="pruningwise")


    if(plot.tree==TRUE){
      plot(tree)
      title("Phylogenetic tree (original)")
      axisPhylo()
    } # end plot.tree

  }# end tree...


  ####################################################################
  ## Check for COALESCENT or RTREE-TYPE ORDERING before SIMULATING: ##
  ####################################################################

  # if(coaltree == FALSE){
  ## Simulation should start from the lowest internal node index (ie n.terminal+1):
  if(unique(tree$edge[,1])[1] == (tree$Nnode+2)){
    ## Simulate from top:bottom?
    x <- 1:nrow(tree$edge)
  }else{
    ## Extra check:
    if(unique(tree$edge[,1])[length(unique(tree$edge[,1]))] == (tree$Nnode+2)){
      ## Simulate from bottom:top?
      x <- rev(c(1:nrow(tree$edge)))
    }else{
      stop("This simulation procedure expects to find the root node/first internal node
           (ie. n.terminal+1) in either the FIRST or LAST row of tree$edge[,1],
           once the tree has been reordered to be in 'pruningwise' configuration.
           This is NOT the case with your tree. Please check.")
    }
    }
  ####################################################################


  ##############################
  ## HANDLE TIP & NODE LABELS ##
  ##############################

  ## Use Fitch Parsimony to get homoplasy distribution (n.subs per site)
  ## unless n.subs distribution has been provided by the user.

  ## tree's tip.labels and node.labels must be numeric (for Fitch parsimony step)
  tree.ori <- tree

  ## NOTE--COERCING TO NUMERIC CAN CAUSE BIG PROBLEMS!!!!
  ## TO DO--FIX FITCH PARSIMONY FN S.T IT CAN WORK WITH
  ## TIP.LABELS & NODE.LABELS THAT ARE NOT NUMERIC!

  if(is.null(n.subs)){

    ## TIP labels ##

    #######################################################################################
    ## CAREFUL -- YOU MAY STILL HAVE PROBLEMS HERE!!! (see misc/rtree.troubleshooting.R) ##
    #######################################################################################


    ## RUN CHECKS TO ENSURE tree$tip.label and rownames(snps) CONTAIN SAME SET OF LABELS!
    ## Check snps vs. tree$tip.labs:
    if(is.null(tree$tip.label)) stop("Trees must have tip.labels corresponding to rownames(snps).")
    if(is.null(rownames(snps))) stop("SNPs must have rownames corresponding to tree$tip.label.")
    if(!all(tree$tip.label %in% rownames(snps))) stop("tree$tip.label and rownames(snps)
                                                      must contain the same set of labels
                                                      so that individuals can be correctly identified.")
    ## Check phen vs. tree$tip.labs:
    if(is.null(names(phen))) stop("Phen must have names corresponding to tree$tip.label.")
    if(!all(tree$tip.label %in% names(phen))) stop("tree$tip.label and names(phen)
                                                      must contain the same set of labels
                                                      so that individuals can be correctly identified.")

  ###################
  ## HANDLE N.SUBS ##
  ###################

  ## if n.subs is a vector (ie. distribution) ##
  ## we use this distribution directly (but in proportion with the number of sites)
  ## to specify the n.subs per site. (Handled within snp.sim fn.)

  ## if n.subs is NULL ##
  ## we compute the distribution of the n.subs-per-site
  ## using the Fitch parsimony score calculation fns from phangorn.


    ###########
    ## TO DO ##   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
    ###########

    ## if either test 2 or test 3 will be run with parsimonious/user-provided reconstruction,
    ## get n.subs from this to avoid duplication..?
    #     if(any(c("simultaneous", "subsequent") %in% test) & snps.reconstruction != "ace"){
    #
    #       ## run get.ancestral.pars
    #       snps.pars <- get.ancestral.pars(var=snps, tree=tree)
    #
    #       ## get elements of output
    #       snps.rec <- snps.pars$var.rec
    #       snps.subs.edges <- snps.pars$subs.edges
    #
    #       ## CHECK--Compare costs:
    #       cost1 <- get.fitch.n.mts(snps=snps, tree=tree)
    #       cost2 <- sapply(c(1:length(snps.subs.edges)), function(e) length(snps.subs.edges[[e]][["total"]]))
    #       table(cost1)
    #       table(cost2) ## longer tail...
    #
    #     }else{

    ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

    ## get parsimomy cost for each SNP locus using fitch:
    n.subs <- get.fitch.n.mts(snps=snps, tree=tree)
    n.subs <- table(n.subs)
    ## handle n.subs "levels" with 0 SNP loci at those levels:
    noms <- as.numeric(names(n.subs))
    temp <- rep(0, max(noms))
    for(i in 1:max(noms)){
      if(i %in% noms) temp[i] <- n.subs[which(noms==i)]
    }
    n.subs <- temp
    # }
  }



  #####################################################
  ## 1) Simulate multiple snps datasets to compare your
  ## real correlations w phen to  #####################
  #####################################################

  if(is.null(sim.n.snps)) sim.n.snps <- n.snps
  out <- genomes <- snps.mat <- list()

  for(i in 1:n.reps){
    ## SIMULATE A DATASET | your tree ##
    out[[i]] <- snp.sim(n.snps = sim.n.snps,
                        n.subs = n.subs,
                        n.snps.assoc = 0,
                        assoc.prob = 100,
                        tree = tree,
                        phen.loci = NULL,
                        heatmap = FALSE,
                        reconstruct = FALSE,
                        dist.dna.model = dist.dna.model,
                        row.names = rownames(snps),
                        seed = NULL)

    genomes[[i]] <- out[[i]][[1]]

    ## Modify genomes/snps matrices
    if(!is.null(genomes[[i]])){
      snps.mat[[i]] <- genomes[[i]]
    }else{
      snps.mat[[i]] <- NULL
    }

    gc()

  } # end for loop

  print("treeWAS snps sim done.")


  ################################################################
  ## 3) Get results:##############################################
  #### Determine the phylogenetially correct p-values for SNPs | #
  ##   null distributions of correlations from simulated data ####
  #### Synthesize results output: List of all significant SNPs, ##
  ##   their names/locations, their p-values for this phenotype ##
  ################################################################

  ##################################################
  ## RUN CHECKS ONCE BEFORE get.sig.snps FOR LOOP ##
  ##################################################

  ## NOTE: These checks are repeated within the get.sig.snps fn
  ## as an extra layer of safety/ in case users want to use it alone,
  ## but it is more economical to run them once outside of the for loop..

  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

  snps.unique <- snps.index <- snps.sim.unique <- snps.sim.index <- NULL

  #################
  ## Handle snps ##
  #################
  ## Check snps column names
  if(is.null(colnames(snps))) colnames(snps) <- c(1:ncol(snps))

  ################################
  ## Handle snps.sim --> matrix ##
  ################################
  snps.sim <- snps.mat

  ## Handle matrix/list input:
  if(class(snps.sim) == "list"){
    ## If list of length 1...
    if(length(snps.sim) == 1){
      ## keep matrix:
      snps.sim <- snps.sim[[1]]
    }else{
      ## If list of multiple matrices...
      ## merge all elements into one big matrix
      ## by pasting columns together:
      snps.sim <- do.call("cbind", snps.sim)
    }
  }


  #################
  ## Handle phen ##
  #################
  ## convert phenotype to numeric:
  ## NOTE--this is also necessary for returning results in step (5)!
  phen.ori <- phen
  if(!is.numeric(phen)) phen <- as.numeric(phen)
  ## for ease of interpretation,
  ## if phen has 2 levels, 1 and 2,
  ## make these 0 and 1:
  if(length(unique(phen))!=2){
    stop("This function is only designed for phenotypes with two levels.")
  }else{
    if(length(phen[-c(which(phen==1), which(phen==2))])==0){
      phen <- replace(phen, which(phen==1), 0)
      phen <- replace(phen, which(phen==2), 1)
    }
  }
  ## ensure ind names not lost
  names(phen) <- names(phen.ori)

  ##############################################################################################
  ## Reconstruct ancestral SNPs & phen by parsimony/ACE (for tests simultaneous & subsequent) ##
  ##############################################################################################

  ## Ensure we are only reconstructing ancestral states ONCE here, to be used in MULTIPLE tests later.
  snps.REC <- snps.sim.REC <- phen.REC <- NULL

  if(any(c("simultaneous", "subsequent") %in% test)){


    ############
    ## TO DO: ##   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
    ############

    ## ADD CODE TO HANDLE USER-INPUTTED SNPS/PHEN.RECONSTRUCTIONS:
    ## - Extract code to get subs.edges from within get.ancestral.pars fn in ace.R --> make a get.subs.edges fn.
    ## --> Determine if user input is from ACE or PARSIMONY (eg. Are all values 0/1/0.5 (= parsimony) or are any in between (=ace))
    ## --> NOTE--even if reconstruction provided, will still need to perform reconstruction on SNPS.SIM (& probably phen). Will use
    ##     inferred snps reconstruction method on snps.sim as well... (OR could add another argument to control this??)
    ## --> If input = from PARSIMONY, run get.subs.edges fn
    ## --> Store info as list of snps/phen.REC containing inputted data as var.rec (and get.subs.edges output as subs.edges if parsimony)
    ## --> Proceed by handling this output as you would if snps/phen were reconstructed within treeWAS as below...


    ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

    #######################
    ## Reconstruct SNPs: ##
    #######################

    ## By PARSIMONY or ACE: ##

    ## Reconstruct REAL SNPs: ##
    system.time( # 274
    snps.REC <- asr(var = snps, tree = tree, type = snps.reconstruction)
    )
    snps.rec <- snps.REC$var.rec

    ## Reconstruct SIMULATED SNPs: ##
    system.time(
    snps.sim.REC <- asr(var = snps.sim, tree = tree, type = snps.reconstruction)
    )
    snps.sim.rec <- snps.sim.REC$var.rec


    #######################
    ## Reconstruct phen: ##
    #######################

    ## By PARSIMONY or ACE: ##

    phen.REC <- asr(var = phen, tree = tree, type = phen.reconstruction)
    phen.rec <- phen.REC$var.rec


  } # end reconstruction for tests 2 & 3


  ## !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ###


  ###########################
  ## GET UNIQUE SNPS(.SIM) ##
  ###########################

  ## Get UNIQUE snps + index
  snps.complete <- snps
  temp <- get.unique.matrix(snps, MARGIN=2)
  snps.unique <- temp$unique.data
  snps.index <- temp$index

  ## Get UNIQUE snps.sim + index
  snps.sim.complete <- snps.sim
  temp <- get.unique.matrix(snps.sim, MARGIN=2)
  snps.sim.unique <- temp$unique.data
  snps.sim.index <- temp$index

  ## Get UNIQUE snps.reconstruction
  snps.rec.complete <- snps.rec
  temp <- get.unique.matrix(snps.rec, MARGIN=2)
  snps.rec <- temp$unique.data
  snps.rec.index <- temp$index
  if(!identical(snps.rec.index, snps.index)){
    warning("Careful-- snps and snps.rec should have the same index when reduced
              to their unique forms!") ## SHOULD THIS BE A "STOP" INSTEAD? OR IS THIS ERROR NOT FATAL OR NOT POSSIBLE????
  }

  ## Get UNIQUE snps.sim.reconstruction
  snps.sim.rec.complete <- snps.sim.rec
  temp <- get.unique.matrix(snps.sim.rec, MARGIN=2)
  snps.sim.rec <- temp$unique.data
  snps.sim.rec.index <- temp$index
  if(!identical(snps.sim.rec.index, snps.sim.index)){
    warning("Careful-- snps.sim and snps.sim.rec should have the same index when reduced
              to their unique forms!") ## SHOULD THIS BE A "STOP" INSTEAD? OR IS THIS ERROR NOT FATAL OR NOT POSSIBLE????
  }

  print("reconstructions done")

  #######################
  ## identify sig.snps ##
  #######################

  ## Note: UNIQUE snps & snps.sim are identified WITHIN the get.sig.snps fn
  ## to reduce computational time, but results are identified on the basis of all
  ## ORIGINAL snps & snps.sim columns inputted.

  sig.list <- list()

  # test <- c("terminal", "simultaneous", "subsequent")
  TEST <- as.list(test)

  ## Run get.sig.snps fn once for each association test:
  system.time( # 100 - 164 (why such a difference?)
  for(i in 1:length(TEST)){
    sig.list[[i]] <- get.sig.snps(snps = snps,
                                  snps.unique = snps.unique,
                                  snps.index = snps.index,
                                  snps.sim = snps.sim,
                                  snps.sim.unique = snps.sim.unique,
                                  snps.sim.index = snps.sim.index,
                                  phen = phen,
                                  tree = tree,
                                  test = TEST[[i]],
                                  n.tests = length(TEST),
                                  p.value = p.value,
                                  p.value.correct = p.value.correct,
                                  p.value.by = p.value.by,
                                  snps.reconstruction = snps.rec,
                                  snps.sim.reconstruction = snps.sim.rec,
                                  phen.reconstruction = phen.rec,
                                  rec = snps.reconstruction)
  }
  )

  SCORE3 <- sig.list[[3]]$SCORE3
  sig.list[[3]] <- sig.list[[3]]$res

  names(sig.list) <- test
  # str(sig.list)

  print("get sig snps done.")

  ## DOUBLE CHECKING ##
  #   str(sig.list[[i]])
  #   sig.list[[i]]$sig.snps
  #   sig.list[[i]]$sig.corrs
  #   ## plot
  # hist(sig.list[[i]][[1]]$corr.sim)
  # hist(sig.list[[i]][[1]]$corr.dat)

  # sig.list[[2]][[1]]$corr.dat[snps.assoc]
  # sig.list[[3]][[1]]$corr.dat[snps.assoc]

  ## BUG CHECKING ##
  ## get.sig.snps
  #     snps <-  snps
  #     snps.sim <- snps.sim
  #     phen <- phen
  #     tree <- tree
  #     test <- "simultaneous"
  #     p.value <- p.value
  #     p.value.correct <- p.value.correct
  #     p.value.by <- p.value.by
  #     snps.reconstruction <- snps.rec
  #     snps.sim.reconstruction <- snps.sim.rec
  #     phen.reconstruction <- phen.rec


  #################
  ## GET RESULTS ##
  #################

  RES <- thresholds <- VALS <- DAT <- list()

  ## get results for each test run:
  for(j in 1:length(sig.list)){

    RES[[j]] <- list()

    ## get corr.dat, corr.sim, p.vals x2:
    VALS[[j]] <- list()
    VALS[[j]][[1]] <- sig.list[[j]][[1]]$corr.dat
    VALS[[j]][[2]] <- sig.list[[j]][[1]]$corr.sim

    VALS[[j]][[3]] <- list()
    VALS[[j]][[3]][[1]] <- sig.list[[j]][[1]]$p.vals
    VALS[[j]][[3]][[2]] <- sig.list[[j]][[2]]$p.vals
    names(VALS[[j]][[3]]) <- c("10.x.n.snps", "1.x.n.snps")


    names(VALS[[j]]) <- c("corr.dat",
                         "corr.sim",
                         "p.vals")

    ## isolate thresholds for plot:
    THRESH <- list()
    for(n in 1:length(sig.list[[j]])){
      THRESH[[n]] <- sig.list[[j]][[n]]$sig.thresh
    }
    names(THRESH) <- names(sig.list[[j]])
    thresholds[[j]] <- THRESH

    ## Only plot the threshold that our sims show is most-consistently performing the best:
    thresh.best <- THRESH[["pval.0.01.bonf.count.10.x.n.snps"]]

    ##########################
    ## NEW: MANHATTAN PLOT! ##
    ##########################
    if(plot.manhattan == TRUE){

      ## save next plot:
      if(!is.null(filename.plot)){
        if(length(filename.plot) == length(sig.list)){
          if(class(filename.plot) != "list") filename.plot <- as.list(filename.plot)

          ## save whatever plots before dev.off:
          pdf(file=filename.plot[[j]][1], width=7, height=11)
        }
      }


      manhattan.plot(p.vals = sig.list[[j]][[1]]$corr.dat,
                     col = "wasp",
                     transp = 0.75,
                     sig.thresh = thresh.best,
                     thresh.col="red",
                     snps.assoc = snps.assoc,
                     snps.assoc.col = "red",
                     jitter.amount = 0.00001,
                     min.p = NULL,
                     log10=FALSE,
                     ylab=paste(TEST[[j]], "score", sep=" "))

      ## End saving:
      ## CHECK-- Do we need if statements??
      dev.off() ## Not sure what happens if you run this without having used pdf or dev.copy previously..


      ####
      ## NOTE-- if you want to see the plot, you need to plot it again (dev.copy not working consistently!)
      manhattan.plot(p.vals = sig.list[[j]][[1]]$corr.dat,
                     col = "wasp",
                     transp = 0.75,
                     sig.thresh = thresh.best,
                     thresh.col="red",
                     snps.assoc = snps.assoc,
                     snps.assoc.col = "red",
                     jitter.amount = 0.00001,
                     min.p = NULL,
                     log10=FALSE,
                     ylab=paste(TEST[[j]], "score", sep=" "))

      # ## save plot:
      # if(!is.null(filename.plot)){
      #   if(length(filename.plot) == length(sig.list)){
      #     if(class(filename.plot) != "list") filename.plot <- as.list(filename.plot)
      #     dev.copy(pdf, file=filename.plot[[j]][1], width=7, height=11) # , pointsize=12
      #     dev.off()
      #   }
      # } # end save pdf

    } # end plot.manhattan


    ########################################################
    ## Plot the null distribution w thresholds & findings ##
    ########################################################

    ## NOTE: For simplicity & clarity, only plotting truly significant SNPs,
    ## instead of plotting all findings from all tests/thresholds.


    ## save next plot:
    if(!is.null(filename.plot)){
      if(length(filename.plot) == length(sig.list)){
        if(class(filename.plot) != "list") filename.plot <- as.list(filename.plot)

        if(length(filename.plot[[j]]) > 1){
          pdf(file=filename.plot[[j]][2], width=7, height=11) # , pointsize=12
        }else{
          pdf(file=filename.plot[[j]][1], width=7, height=11) # , pointsize=12
        }
      }
    }

    ## Generate one histogram per test:
    plot.sig.snps(corr.dat = sig.list[[j]][[1]]$corr.dat,
                  corr.sim = sig.list[[j]][[1]]$corr.sim,
                  corr.sim.subset = sig.list[[j]][[1]]$corr.sim[1:10000],
                  sig.corrs = corr.dat[snps.assoc],
                  sig.snps = snps.assoc,
                  sig.thresh = thresh.best,
                  test = TEST[[j]],
                  sig.snps.col = "blue",
                  hist.col = rgb(0,0,1,0.5), # rgb(0,0,1,0.5) # blue ## OR ## rgb(0.1,0.1,0.1,0.5) # darkgrey
                  hist.subset.col = rgb(1,0,0,0.5), # rgb(1,0,0,0.5) # red ## OR ## rgb(0.8,0.8,0.8,0.5) # lightgrey
                  thresh.col = "red",
                  snps.assoc = snps.assoc,
                  snps.assoc.col = "black",
                  bg = "lightgrey",
                  grid = TRUE,
                  freq = FALSE,
                  plot.null.dist = TRUE,
                  plot.dist = FALSE)

    ## End saving:
    dev.off()


    #####
    ## Again-- to see this plot as simTest runs, need to plot again bc dev.copy failing/corrupting at random...
    plot.sig.snps(corr.dat = sig.list[[j]][[1]]$corr.dat,
                  corr.sim = sig.list[[j]][[1]]$corr.sim,
                  corr.sim.subset = sig.list[[j]][[1]]$corr.sim[1:10000],
                  sig.corrs = corr.dat[snps.assoc],
                  sig.snps = snps.assoc,
                  sig.thresh = thresh.best,
                  test = TEST[[j]],
                  sig.snps.col = "blue",
                  hist.col = rgb(0,0,1,0.5), # rgb(0,0,1,0.5) # blue ## OR ## rgb(0.1,0.1,0.1,0.5) # darkgrey
                  hist.subset.col = rgb(1,0,0,0.5), # rgb(1,0,0,0.5) # red ## OR ## rgb(0.8,0.8,0.8,0.5) # lightgrey
                  thresh.col = "red",
                  snps.assoc = snps.assoc,
                  snps.assoc.col = "black",
                  bg = "lightgrey",
                  grid = TRUE,
                  freq=FALSE,
                  plot.null.dist = TRUE,
                  plot.dist = FALSE)

    # ## save plot:
    # if(!is.null(filename.plot)){
    #   if(length(filename.plot) == length(sig.list)){
    #     if(class(filename.plot) != "list") filename.plot <- as.list(filename.plot)
    #     if(length(filename.plot[[j]]) > 1){
    #       dev.copy(pdf, file=filename.plot[[j]][2], width=7, height=11) # , pointsize=12
    #     }else{
    #       dev.copy(pdf, file=filename.plot[[j]][1], width=7, height=11) # , pointsize=12
    #     }
    #     dev.off()
    #   }
    # } # edn save pdf


    ## legend for thresholds? ##

  for(i in 1:length(sig.list[[j]])){

    #############################################
    ## isolate elements of get.sig.snps output ##
    #############################################

    corr.dat <- sig.list[[j]][[i]]$corr.dat
    corr.sim <- sig.list[[j]][[i]]$corr.sim
    p.vals <- sig.list[[j]][[i]]$p.vals
    sig.snps.names <- sig.list[[j]][[i]]$sig.snps.names
    sig.snps <- sig.list[[j]][[i]]$sig.snps
    sig.corrs <- sig.list[[j]][[i]]$sig.corrs
    sig.p.vals <- sig.list[[j]][[i]]$sig.p.vals
    min.p <- sig.list[[j]][[i]]$min.p
    sig.thresh <- sig.list[[j]][[i]]$sig.thresh

    ########################################


    ##################################
    ## 4) (A) Plot the distribution ##
    ##################################

    # plot.sig.snps(corr.dat, corr.sim, sig.corrs, sig.snps,
    #               sig.thresh=sig.thresh, test=TEST[[j]],
    #               plot.null.dist = plot.null.dist,
    #               plot.dist = plot.dist)


    ########################################
    ## 5) Return results list ##############
    ########################################

    if(length(sig.snps)==0) sig.snps <- sig.corrs <- NULL

    ###########
    ## make a data.frame containing all relevant output for sig.snps
    if(length(sig.snps) > 0){

      ## Get counts for n.sig.snps in each cell of the contingency table:
      #     toKeep <- sapply(c(1:length(sig.snps)),
      #                      function(e)
      #                        which(dimnames(snps)[[2]] == sig.snps))
      toKeep <- sig.snps
      snps.toKeep <- snps[,toKeep]

      ##
      if(length(toKeep) > 1){
        S1P1 <- sapply(c(1:ncol(snps.toKeep)),
                       function(e)
                         length(which(snps.toKeep[which(phen==1),e]==1)))
        S0P0 <- sapply(c(1:ncol(snps.toKeep)),
                       function(e)
                         length(which(snps.toKeep[which(phen==0),e]==0)))
        S1P0 <- sapply(c(1:ncol(snps.toKeep)),
                       function(e)
                         length(which(snps.toKeep[which(phen==0),e]==1)))
        S0P1 <- sapply(c(1:ncol(snps.toKeep)),
                       function(e)
                         length(which(snps.toKeep[which(phen==1),e]==0)))
      }else{
        ## if only ONE sig snp (haploid) identified:
        S1P1 <- length(which(snps.toKeep[which(phen==1)]==1))
        S0P0 <- length(which(snps.toKeep[which(phen==0)]==0))
        S1P0 <- length(which(snps.toKeep[which(phen==0)]==1))
        S0P1 <- length(which(snps.toKeep[which(phen==1)]==0))

      }
      df <- data.frame(sig.snps,
                       sig.p.vals,
                       sig.corrs,
                       S1P1, S0P0, S1P0, S0P1)
      names(df) <- c("SNP.locus",
                     "p.value",
                     "Test.statistic",
                     "S1P1", "S0P0", "S1P0", "S0P1")

      ## NOTE: Could return sig.snps.names somewhere here
      ## in addition to sig.snps loci ####    ####    ####    ####

    }else{
      df <- "No significant SNPs found."
    }

    ## 0 p.vals
    #   min.p <- paste("p-values listed as 0 are <",
    #                  1/length(corr.sim), sep=" ")
    min.p <- 1/length(corr.sim)
    names(min.p) <- c("p-values listed as 0 are less than:")

    ## TO DO:
    ## ADD MANHATTAN PLOT

    results <- list()
    results[[1]] <- sig.thresh
    results[[2]] <- df
    results[[3]] <- min.p

    names(results) <- c("sig.thresh",
                        "sig.snps",
                        "min.p.value")

    RES[[j]][[i]] <- results
  } # end for loop (i)

    names(RES[[j]]) <- names(sig.list[[j]])

  } # end for loop (j)

  ## assign test names to main list components:
  names(RES) <- names(VALS) <- test


  ################################
  ## NEW! Get COMBINED results: ##
  ################################
  ## get uniques(sig.snps) for terminal, simultaneous, subsequent results combined
  ## for best thresh only (ie. pval.0.01.bonf.count.10.x.n.snps)
  SNP.loci <- vector("list", length=length(test))
  names(SNP.loci) <- test
  for(t in 1:length(test)){
    temp <- RES[[t]]$pval.0.01.bonf.count.10.x.n.snps$sig.snps
    if(is.vector(temp)){
      SNP.loci[[t]] <- NULL
    }else{
      SNP.loci[[t]] <- temp$SNP.locus
    }
  }
  ## store 3 tests in separate list:
  SNP.loci.ori <- SNP.loci

  ## make list of length 2 (all, separately):
  SNP.loci <- vector("list", length=2)
  names(SNP.loci) <- c("treeWAS.combined", "treeWAS")
  if(length(as.vector(unlist(SNP.loci.ori))) > 0){
    SNP.loci[[1]] <- sort(unique(as.vector(unlist(SNP.loci.ori))), decreasing=FALSE)
  }else{
    SNP.loci[[1]] <- NULL
  }
  SNP.loci[[2]] <- SNP.loci.ori


  ## get data:
  DAT <- list(snps.sim = snps.sim.complete,
              snps.rec = snps.rec.complete,
              snps.sim.rec = snps.sim.rec.complete,
              phen.rec = phen.rec)

  ## get output:
  results <- list(dat=DAT,
                  vals=VALS,
                  thresh=thresholds,
                  res=RES,
                  treeWAS.combined=SNP.loci,
                  SCORE3=SCORE3)

  return(results)

} # end treeWAS




##############################################################################################
## legend (temporary?)
## Not necessarily needed--usually only a couple UNIQUE thresholds...
#     par(mfrow=c(1,2))
#     par(oma = c(5,4,0,0) + 0.1)
#     par(mar = c(0,0,1,1) + 0.1)
#     ## column 1:
#     midpoints1 <- barplot(rep(10, length(THRESH)/2),
#                         col = seasun(length(THRESH))[1:(length(THRESH)/2)],
#                         horiz=TRUE)
#     ##overlay names:
#     text(3, midpoints1, labels=names(THRESH)[1:(length(THRESH)/2)], cex=0.75, adj=0.3)
#
#     ## column 2:
#     midpoints2 <- barplot(rep(10, length(THRESH)/2),
#                           col = seasun(length(THRESH))[((length(THRESH)/2)+1):length(THRESH)],
#                           horiz=TRUE)
#     ##overlay names:
#     text(3, midpoints2, labels=names(THRESH)[((length(THRESH)/2)+1):length(THRESH)], cex=0.75, adj=0.3)
#
#     ## return to original par settings:
#     par(mfrow=c(1,1))
#
#     ## add title
#     title("Legend: Significance Thresholds")
#
#     par(oma=c(0,0,0,0))
#     par(mar=c(5,4,4,2)+0.1)
##############################################################################################



# ## only 2 unique sets of p.vals for each test
# ## (for n.snps & 10x n.snps):
# pv <- list()
# summ <- list()
# for(j in 1:length(results)){
#   pv[[j]] <- list()
#   summ[[j]] <- list()
# for(i in 1:length(results[[j]])){
#   pv[[j]][[i]] <- results[[j]][[i]]$p.vals
#   summ[[j]][[i]] <- summary(pv[[j]][[i]])
# }
# }
#
# length(unique(summ[[2]]))
# length(unique(summ[[2]][seq(2, length(summ[[3]]), 2)]))
# length(unique(summ[[2]][seq(2, length(summ[[3]]), 2)]))




###############
## CHECK!!!! ##
###############
## w assoc.prob == 100, still getting low/0(!) scores for "snps.assoc"
## get original phen for all terminal + internal nodes and edges...

## checking simultaneous score:

# str(foo)
# phen.str <- foo$phen.plot.col
# phen1 <- foo$phen
# phen2 <- phen.str$all.nodes
#
# head(phen1, 20)
# head(phen2, 20)
#
# phen2[which(phen2 == "blue")] <- "A"
# phen2[which(phen2 == "red")] <- "B"
# phen2 <- as.factor(phen2)
# names(phen2) <- paste("ind", 1:length(phen2), sep=".")
#
# phen.edges <- phen.str$edges
# phen.edges[which(phen.edges == "blue")] <- "A"
# phen.edges[which(phen.edges == "red")] <- "B"
# # phen.edges[which(phen.edges == "green")] <- "C"
# phen.edges <- as.factor(phen.edges)
# names(phen.edges) <- paste("ind", 1:length(phen.edges), sep=".")
# head(phen.edges)
#
# ## which edges should/do contain phen subs?
# which(phen.edges == "green")
# ## run relevant code in reconstruct to get phen.subs.edges list (to see which edges are identified):
# phen.subs.edges$total
#
# snps.diffs <- list()
# snps.assoc.index <- index[snps.assoc] ## GAK! -- all snps.assoc have the same index (915)! set.seed problem, for a start...
# for(i in snps.assoc.index){
#   snps.diffs[[i]] <- get.branch.diffs(var = snps.rec[,i],
#                                       edges = edges)
# }
# which(abs(snps.diffs[[i]]) == 1)
# ## ALSO--important to NOTE that a lot of the reconstructed edges are only off by one...
# ## (thus should the subsequent score not be doing much better????)
# length(which(which(abs(snps.diffs[[i]]) == 1) %in% phen.subs.edges$total)) ## SO Shouldn't the max corr.dat score2 be 11????

###############
## CHECK!!!! ##
###############
## The FPR for Score 2 (simultaneous) should NOT be that high-- something wrong with the threshold selection??
#
# snps.assoc.ori.ori <- snps.assoc
#
# sim.dat <- sim.dat.ori <- results$dat$simultaneous
# corr.dat <- corr.dat.ori <- sim.dat$corr.dat
# corr.sim <- corr.sim.ori <- sim.dat$corr.sim
# p.vals <- p.vals.ori <- sim.dat$p.vals
#
# sim.res <- sim.res.ori <- results$res$simultaneous
# str(sim.res[[1]])
# str(sim.res[[32]])
# thresh.ori <- thresh <- sim.res[[1]]$sig.thresh
# thresh.ori <- thresh <- sim.res[[32]]$sig.thresh
#
# table(corr.sim)
#
# hist(corr.sim)
# lines(density(corr.sim), col="red", lwd=2)
#
# str(density(corr.sim))
#
# ## with first 10000 only??
# corr.sim.ori <- corr.sim
# corr.sim <- corr.sim.ori[1:10000]
# table(corr.sim)



# ###
#
# t.corr.sim <- results$dat$terminal$corr.sim
#
# ########################
# ## get density curve: ##
# ########################
#
# ## SIMULTANEOUS SCORE: ##
# dat <- corr.sim[1:10000]
# d <- density(dat)
# # from=0 may be necessary for aligning polygon w hsit alon x-axis, BUT may cause problems for polygon drawing (try lines instead?)
# xmax <- max(hist(dat)$mids)+min(hist(dat)$mids)
# ymax <- ceiling(max(d$y))
# hist(dat, freq=F, xlim=c(0,xmax), ylim=c(0,ymax))
# # lines(d, col="red", lwd=2, xlim=c(0,1), ylim=c(0,ymax))
# polygon(d, col=transp("red", 0.25), border="red", lwd=2, xlim=c(0,xmax), ylim=c(0,ymax))
#
# ## TERMINAL SCORE ##
# dat <- t.corr.sim
# d <- density(dat, from=0)
# ymax <- ceiling(max(d$y))
# hist(dat, freq=F, xlim=c(0,1), ylim=c(0,ymax))
# # lines(d, col="red", lwd=2, xlim=c(0,1), ylim=c(0,ymax))
# polygon(d, col=transp("red", 0.25), border="red", lwd=2, xlim=c(0,1), ylim=c(0,ymax))


###




###
