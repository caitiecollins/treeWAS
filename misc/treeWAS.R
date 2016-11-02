
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
#                       assoc.prob = 95,
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
                    plot.null.dist = TRUE,
                    plot.dist = FALSE,
                    snps.reconstruction = "parsimony",
                    phen.reconstruction = "parsimony"){

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
  x <- snps
  n.snps <- ncol(snps)

  ## convert phenotype to factor
  phen <- as.factor(phen)
  y <- phen

  ## set n.ind:
  n.ind <- length(y)
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
    if(is.rooted(tree)==FALSE) tree <- midpoint(tree)


    if(plot.tree==TRUE){
      plot(tree)
      title("Phylogenetic tree (original)")
      axisPhylo()
    } # end plot.tree

  }# end tree...


  ##################################
  ## Check if COALESCENT or RTREE ## (is.ultrametric? any other tests required here???) ##
  ########################################################################################

  if(is.ultrametric(tree)){
    coaltree <- TRUE
  }else{
    coaltree <- FALSE
  }


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

    ## FOR NOW -- REPLACING tree$tip.label w c(1:N)
    ## WARNING -- ASSUMES that tree$edge cells containing 1:N CORRESPOND to snps rown 1:N!!!!!!!!

    tree$tip.label <- c(1:(tree$Nnode+1))

    #     if(all.is.numeric(tree$tip.label)){
    #       ## if we can convert to numeric, do so:
    #       tree$tip.label <- as.numeric(tree$tip.label)
    #     }else{
    #
    #       ## if we can remove "t" to get numeric, do so:
    #       prefix <- keepFirstN(tree$tip.label, 1)
    #       if(all(tolower(prefix) == "t")){
    #         temp <- removeFirstN(tree$tip.label, 1)
    #         if(all.is.numeric(temp)){
    #           tree$tip.label <- as.numeric(temp)
    #         }else{
    #           ## else, replace with numeric indices:
    #           # tree$tip.label <- c(1:length(tree$tip.label))
    #           warning("Site-wise parsimony scores (phangorn's
    #                   fitch parsimony function) may not be calculated correctly
    #                   when tip.labels are not numeric.
    #                   Please change tree$tip.label to numeric values.")
    #         }
    #       }
    #     }

    ## NODE labels ##
    if(all.is.numeric(tree$node.label)){
      tree$node.label <- as.numeric(tree$node.label)
    }else{
      ## if we can remove "NODE_" to get numeric, do so:
      prefix <- keepFirstN(tree$node.label, 4)
      if(all(tolower(prefix) == "node")){
        temp <- removeFirstN(tree$node.label, 5)
        if(all.is.numeric(temp)){
          tree$node.label <- as.numeric(temp)
        }else{
          ## else, replace with numeric indices:
          # tree$node.label <- c((n.ind+1):(n.ind+tree$Nnode))
          warning("Site-wise parsimony scores (phangorn's
                  fitch parsimony function) may not be calculated correctly
                  when node.labels are not numeric.
                  Please change tree$node.label to numeric values.")
        }
      }
    }
  }


  ###################
  ## HANDLE N.SUBS ##
  ###################

  ## if n.subs is a vector (ie. distribution) ##
  ## we use this distribution directly (but in proportion with the number of sites)
  ## to specify the n.subs per site. (Handled within snp.sim fn.)

  ## if n.subs is NULL ##
  ## we compute the distribution of the n.subs-per-site
  ## using the Fitch parsimony score calculation fns from phangorn.

  if(is.null(n.subs)){

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
                        snp.root = NULL,
                        n.snps.assoc = 0,
                        assoc.prob = 100,
                        tree = tree,
                        phen.loci = NULL,
                        heatmap = FALSE,
                        reconstruct = FALSE,
                        dist.dna.model = dist.dna.model,
                        row.names = row.names(snps),
                        coaltree = coaltree,
                        set = NULL,
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


  ## NOTE TO CHECK! ##   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

  ## Are we SLOWING things down significantly by identifying UNIQUE snps & snps.sim WITHIN the get.sig.snps fn??
  ## And could we identify unique snps/snps.sim HERE (AND add an extra INDEX argument to get.sig.snps)??
  ## (ie. get.sig.snps INPUT = UNIQUE snps, snps.sim + index & OUTPUT = results for ALL ORIGINAL/NON-unique sites...).

  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

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

    ## By PARSIMONY: ##
    if(snps.reconstruction == "parsimony"){
      ## Reconstruct REAL SNPs: ##
      snps.REC <- asr(var = snps, tree = tree, type = "parsimony")
      snps.rec <- snps.REC$var.rec

      ## Reconstruct SIMULATED SNPs: ##
      snps.sim.REC <- asr(var = snps.sim, tree = tree, type = "parsimony")
      snps.sim.rec <- snps.sim.REC$var.rec
    }

    ## By ACE: ##
    if(snps.reconstruction == "ace"){
      ## Reconstruct REAL SNPs: ##
      snps.REC <- asr(var = snps, tree = tree, type = "ace")
      snps.rec <- snps.REC$var.rec

      ## Reconstruct SIMULATED SNPs: ##
      snps.sim.REC <- asr(var = snps.sim, tree = tree, type = "ace")
      snps.sim.rec <- snps.sim.REC$var.rec
    }

    #######################
    ## Reconstruct phen: ##
    #######################

    ## By PARSIMONY: ##
    if(phen.reconstruction == "parsimony"){
      phen.REC <- asr(var = phen, tree = tree, type = "parsimony")
      phen.rec <- phen.REC$var.rec
    }

    ## By ACE: ##
    if(phen.reconstruction == "ace"){
      phen.REC <- asr(var = phen, tree = tree, type = "ace")
      phen.rec <- phen.REC$var.rec
    }

  } # end reconstruction for tests 2 & 3


  ###########################
  ## GET UNIQUE SNPS(.SIM) ##
  ###########################

  ## TO DO: ## !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ###
  ## NOTE--WOULD BE GOOD TO ADD SIMILAR SOLN FOR RECONSTRUCT ABOVE AS W GET.SIG.SNPS BELOW (ie.
  ## ALLOW FOR INPUT OF UNIQUE VAR AND INDEX AS ARGUMENTS).
  ## ONCE DONE--MOVE UNIQUE CODE BELOW TO ABOVE THE RECONSTRUCTION CODE SEGMENT...

  ## !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ###

  ## Get UNIQUE snps + index
  snps.ori <- snps
  temp <- get.unique.matrix(snps, MARGIN=2)
  snps.unique <- temp$unique.data
  snps.index <- temp$index

  ## Get UNIQUE snps.sim + index
  snps.sim.ori <- snps.sim
  temp <- get.unique.matrix(snps.sim, MARGIN=2)
  snps.sim.unique <- temp$unique.data
  snps.sim.index <- temp$index

  ## Get UNIQUE snps.reconstruction
  temp <- get.unique.matrix(snps.rec, MARGIN=2)
  snps.rec <- temp$unique.data
  snps.rec.index <- temp$index
  if(!identical(snps.rec.index, snps.index)){
    warning("Careful-- snps and snps.rec should have the same index when reduced
              to their unique forms!") ## SHOULD THIS BE A "STOP" INSTEAD? OR IS THIS ERROR NOT FATAL OR NOT POSSIBLE????
  }

  ## Get UNIQUE snps.sim.reconstruction
  temp <- get.unique.matrix(snps.sim.rec, MARGIN=2)
  snps.sim.rec <- temp$unique.data
  snps.sim.rec.index <- temp$index
  if(!identical(snps.sim.rec.index, snps.sim.index)){
    warning("Careful-- snps.sim and snps.sim.rec should have the same index when reduced
              to their unique forms!") ## SHOULD THIS BE A "STOP" INSTEAD? OR IS THIS ERROR NOT FATAL OR NOT POSSIBLE????
  }

  #######################
  ## identify sig.snps ##
  #######################

  ## Note: UNIQUE snps & snps.sim are identified WITHIN the get.sig.snps fn
  ## to reduce computational time, but results are identified on the basis of all
  ## ORIGINAL snps & snps.sim columns inputted.

  sig.list <- list()

  TEST <- as.list(test)

  ## Run get.sig.snps fn once for each association test:
  for(i in 1:length(TEST)){
    sig.list[[i]] <-  get.sig.snps(snps = snps,
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
                                   phen.reconstruction = phen.rec)
  }

  # sig <- list()
  # for(i in 1:length(sig.list[[1]])) sig[[i]] <- sig.list[[1]][[i]]$sig.snps


  ## DOUBLE CHECKING ##
  #   str(sig.list[[i]])
  #   sig.list[[i]]$sig.snps
  #   sig.list[[i]]$sig.corrs
  #   ## plot
  #   hist(sig.list[[i]]$corr.sim)
  #   hist(sig.list[[i]]$corr.dat)

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

  RES <- list()

  ## get results for each test run:
  for(i in 1:length(sig.list)){

  #############################################
  ## isolate elements of get.sig.snps output ##
  #############################################

  corr.dat <- sig.list[[i]]$corr.dat
  corr.sim <- sig.list[[i]]$corr.sim
  p.vals <- sig.list[[i]]$p.vals
  sig.snps.names <- sig.list[[i]]$sig.snps.names
  sig.snps <- sig.list[[i]]$sig.snps
  sig.corrs <- sig.list[[i]]$sig.corrs
  sig.p.vals <- sig.list[[i]]$sig.p.vals
  min.p <- sig.list[[i]]$min.p
  sig.thresh <- sig.list[[i]]$sig.thresh

  ########################################


  ##################################
  ## 4) (A) Plot the distribution ##
  ##################################

  plot.sig.snps(corr.dat, corr.sim, sig.corrs, sig.snps,
                sig.thresh=sig.thresh, test=TEST[[i]],
                plot.null.dist = plot.null.dist,
                plot.dist = plot.dist)


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
  results[[1]] <- corr.dat
  results[[2]] <- corr.sim
  results[[3]] <- p.vals
  results[[4]] <- sig.thresh
  results[[5]] <- df
  results[[6]] <- min.p

  names(results) <- c("corr.dat",
                      "corr.sim",
                      "p.vals",
                      "sig.thresh",
                      "sig.snps",
                      "min.p.value")

  RES[[i]] <- results
  } # end for loop

  names(RES) <- test
  results <- RES

  return(results)

} # end treeWAS


