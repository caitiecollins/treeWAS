
####################
## coalescent.sim ##
####################

## a function for simulating trees under a fully-linked coalescent model.
## optional simulation of a phenotype and phenotypically-associated SNPs is implemented.
## optional use of a distribution to guide the substitution rate of the non-associated SNPs is implemented.

## TO DO:
## 1) (Re-)implement associated SNP randomization procedure...
## want to implement procedures that combine the above options...
## 2) Allow phenotypically-associated SNPs simulation to be optionally guided
## by a user-inputted phenotype for the terminal nodes (--> would need to simulate
## phenotypic substitutions from the terminal nodes UP to the root, the reverse
## of the current procedure...)
## 3) Implement assoc.options (currently using deprecated "all" option without requiring argument,
## but would like to consider implementing alternative "model" option(s))


## ARGUMENTS ##
# n.ind <- 100 # n.genomes you want to end up with
# gen.size <- 1000000 # bases
# theta <- gen.size*2 # (if sim.by=="branch")# OR # 1*2 # (if sim.by=="locus")
# biallelic <- TRUE # if TRUE, select ONLY complementary nt; if FALSE,
#                select from 3 alternatives (ie. A/C/G/T-current nt)
# seed <- 1 # allow user to control randomization to get reproducible results.
# n.snps.assoc <- 5
# assoc.option <- c("all", "model") # deprecated (only "all" available)
# sim.by <- c("locus", "branch") # deprecated (only "locus" has all current protocols implemented)


## EXAMPLE ##
# out <- coalescent.sim(n.ind=100, gen.size=10000, sim.by="locus",
#                       theta=1*2, dist=NULL,
#                       theta_p=15, phen=NULL,
#                       n.snps.assoc=20, assoc.option="all", assoc.prob=90,
#                       haploid=TRUE, biallelic=TRUE, seed=1,
#                       plot=TRUE, heatmap=FALSE, plot2="UPGMA")

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param n.ind An integer specifying the number of individual genomes to simulate
#' (ie. the number of terminal nodes in the tree).
#' @param n.snps An integer specifying the number of genetic loci to simulate.
#' @param n.subs Either an integer or a vector (containing a distribution) that is
#' used to determine the number of substitutions
#' to occur on the phylogenetic tree for each genetic locus (see details).
#' @param n.snps.assoc An optional integer specifying the number of genetic loci
#' @param assoc.prob An optional integer (> 0, <= 100) specifying the strength of the
#' association between the n.snps.assoc loci and the phenotype (see details).
#' @param n.phen.subs An integer specifying the expected number of phenotypic
#' substitutions to occur on the phylogenetic tree (through the same process as
#' the n.subs parameter when n.subs is an integer (see details)).
#' @param phen An optional vector containing a phenotype for each of the
#' n.ind individuals if no phenotypic simulation is desired.
#' @param heatmap A logical indicating whether to produce a heatmap of the genetic distance
#' between the simulated genomes of the n.ind individuals.
#' @param reconstruct Either a logical indicating whether to attempt to reconstruct
#' a phylogenetic tree using the simulated genetic data, or one of c("UPGMA", "nj", "ml")
#' to specify that tree reconstruction is desired by one of these three methods
#' (Unweighted Pair Group Method with Arithmetic Mean, Neighbour-Joining, Maximum-Likelihood).
#' @param seed An optional integer controlling the pseudo-random process of simulation. Two
#' instances of coalescent.sim with the same seed and arguments will produce identical output.
#'
#' @details #### n.subs ####
#' If the value of the n.subs parameter is set to an integer, this integer is
#' used as the parameter of a Poisson distribution from which the number of substitutions to
#' occur on the phylogenetic tree is drawn for each of the n.snps simulated genetic loci.
#' If n.subs is a vector containing a distribution, this is used directly (in proportion to n.snps)
#' to define the number of substitutions per site. For example, if n.subs=c(3000, 900, 70, 20, 0, 10)
#' and n.snps=8000, then 6000 simulated sites will undergo exactly
#' one substitution somewhere on the phylogenetic tree, 1800 will undergo two,
#' 140 three, 40 four, 0 five, and 20 six.
#' #### assoc.prob ####
#' The assoc.prob parameter controls the strength of association through a process analagous to dilution.
#' All n.snps.assoc loci are initially simulated to undergo a substitution
#' every time the phenotype undergoes a substitution (ie. perfect association).
#' The assoc.prob parameter then acts like a dilution factor, removing (100 - assoc.prob)%
#' of the substitutions that occurred during simulation under perfect association.
#'
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

########################################################################

############
## NOTES: ##
############
## theta_p changed to n.phen.subs (and just n.subs in phen.sim.R)


## OLD ARGS: ##
# (n.ind=100, gen.size=10000, sim.by="locus",
#  theta=1*2, dist=NULL,
#  n.phen.subs=15, phen=NULL,
#  n.snps.assoc=5, assoc.option="all", assoc.prob=90,
#  haploid=TRUE, biallelic=TRUE, seed=NULL,
#  plot=TRUE, heatmap=FALSE, plot2="UPGMA")

## NEW ARGS: ##
# n.ind <- 100
# n.snps <- 10000
# n.subs <- dist_0.1
# n.snps.assoc <- 10
# assoc.prob <- 100
# n.phen.subs <- 15
# phen <- NULL
# plot <- TRUE
# heatmap <- FALSE
# reconstruct <- FALSE
# dist.dna.model <- "JC69"
# grp.min <- 0.25
# row.names <- NULL
# coaltree <- FALSE
# s <- 20
# af <- 10
# set <- 3 # NULL #
# filename <- NULL
# seed <- 77


## TOY EG for PRESENTATION:
# c.sim <- coalescent.sim(n.ind=15,
#                         n.snps=1000,
#                         n.subs=dist_0,
#                         n.snps.assoc=10,
#                         assoc.prob=90,
#                         n.phen.subs=5,
#                         phen=NULL,
#                         plot=TRUE,
#                         heatmap=FALSE,
#                         reconstruct=FALSE,
#                         dist.dna.model="JC69",
#                         grp.min = 0.25,
#                         row.names=NULL,
#                         coaltree = TRUE,
#                         s = 10,
#                         af = 5,
#                         filename = NULL,
#                         set=1,
#                         seed=5)

## EG:
# c.sim <- coalescent.sim(n.ind=100,
#                         n.snps=10000,
#                         n.subs=dist_0,
#                         n.snps.assoc=10,
#                         assoc.prob=100,
#                         n.phen.subs=15,
#                         phen=NULL,
#                         plot=TRUE,
#                         heatmap=FALSE,
#                         reconstruct=FALSE,
#                         dist.dna.model="JC69",
#                         grp.min = 0.25,
#                         row.names=NULL,
#                         coaltree = TRUE,
#                         s = 10,
#                         af = 5,
#                         filename = NULL,
#                         set=1,
#                         seed=78)


coalescent.sim <- function(n.ind=100,
                           n.snps=10000,
                           n.subs=1,
                           n.snps.assoc=0,
                           assoc.prob=100,
                           n.phen.subs=15,
                           phen=NULL,
                           plot=TRUE,
                           heatmap=FALSE,
                           reconstruct=FALSE,
                           dist.dna.model="JC69",
                           grp.min = NULL,
                           row.names=TRUE,
                           set=NULL,
                           coaltree = TRUE,
                           s = 20,
                           af = 10,
                           filename = NULL,
                           seed=NULL){
  ## load packages:
  # require(adegenet)
  # require(ape)
  # require(phangorn)

  sets <- NULL

  if(length(which(c(plot, heatmap, reconstruct)==TRUE))==1){
    par(ask=FALSE)
  }else{
    par(ask=TRUE)
  }

  ################################
  ## Simulate Phylogenetic Tree ##
  ################################
  if(coaltree == TRUE){
    tree <- coalescent.tree.sim(n.ind = n.ind, seed = seed)
  }else{
    set.seed(seed)
    tree <- rtree(n = n.ind)
  }
  tree <- reorder.phylo(tree, order="pruningwise")


  ########################
  ## Simulate Phenotype ##
  ########################
  if(set == 3){

    ############
    ## NEW Q: ##
    ############
    ## if Q contains RATES --> P contains probs
    ## QUESTION -- HOW TO INTERPRET/PREDICT THE RELATIVE EFFECTS OF ASSOC.FACTOR AND N.SUBS (+ BRANCH LENGTH) ON ASSOC STRENGTH, N.SUBS PER TREE ?

    if(is.null(s)) s <- 20 # n.subs
    if(is.null(af)) af <- 10 # association factor
    ## Modify s to account for sum(tree$edge.length):
    ## --> ~ s/10 (= 2) for coaltree
    ## OR --> ~ s/100 (= 0.2) for rtree
    s <- s/sum(tree$edge.length)


    Q.mat <- matrix(c(NA, 1*s, 1*s, 0,
                      1*af*s, NA, 0, 1*af*s,
                      1*af*s, 0, NA, 1*af*s,
                      0, 1*s, 1*s, NA),
                    nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))

    diag(Q.mat) <- sapply(c(1:nrow(Q.mat)), function(e) -sum(Q.mat[e, c(1:ncol(Q.mat))[-e]]))
    # Q <- Q.mat

    ## SAVE PANEL PLOT:
    if(!is.null(filename)){
      pdf(file=filename[[2]], width=7, height=11)
      # dev.copy(pdf, file=filename[[2]], width=7, height=11)
    }

    ## RUN SNP.SIM.Q: ##
    system.time(
    snps.list <- snp.sim.Q(n.snps = n.snps,
                           n.subs = n.subs,
                           snp.root = NULL,
                           n.snps.assoc = n.snps.assoc,
                           assoc.prob = assoc.prob,
                           ## dependent/corr' transition rate/prob mat:
                           Q = Q.mat,
                           tree = tree,
                           n.phen.subs = n.phen.subs,
                           phen.loci = NULL,
                           heatmap = FALSE,
                           reconstruct = FALSE,
                           dist.dna.model = "JC69",
                           grp.min = grp.min,
                           row.names = NULL,
                           set=set,
                           seed=seed)
    )
    snps <- snps.list$snps
    snps.assoc <- snps.list$snps.assoc
    sets <- NULL
    phen <- snps.list$phen
    phen.nodes <- snps.list$phen.nodes

    if(!is.null(filename)){
      ## end saving panel plot:
      dev.off()
    }

  }else{


    ##################
    ## SETS 1 and 2 ##
    ##################

  if(is.null(phen)){
    ## get list of phenotype simulation output
    phen.list <- phen.sim(tree, n.subs = n.phen.subs, grp.min = grp.min, seed = seed)

    ## get phenotype for terminal nodes only
    phen <- phen.list$phen

    ## get phenotype for all nodes,
    ## terminal and internal
    phen.nodes <- phen.list$phen.nodes

    ## get the indices of phen.subs (ie. branches)
    phen.loci <- phen.list$phen.loci
  }else{
    #############################
    ## User-provided Phenotype ##
    #############################
    phen.nodes <- phen
    phen.loci <- NULL
  }


  ###################
  ## Simulate SNPs ##
  ###################

  ## TO DO: #######################################
  ## CHECK SNP SIMULATION FOR COMPUTATIONAL SPEED! #########################################################################################
  #################################################
  ## 10 --> 53 --> 12.5
  ## Are the remaining extra 2.5 seconds still just a result of the while loop??
  ## Or have I slowed anything down in the post-processing steps as well??????????????????

  #   n.snps <- 10000 # 13.3
  #   n.snps <- 100000 # 153.8
  #   n.snps <- 1000000 # >> 1941.7 (stopped trying..)
  #
    # system.time(
  snps.list <- snp.sim(n.snps=n.snps,
                       n.subs=n.subs,
                       n.snps.assoc=n.snps.assoc,
                       assoc.prob=assoc.prob,
                       tree=tree,
                       phen.loci=phen.loci,
                       heatmap=heatmap,
                       reconstruct=reconstruct,
                       dist.dna.model=dist.dna.model,
                       row.names = NULL,
                       set=set,
                       seed=seed)
  # )

  snps <- snps.list$snps
  snps.assoc <- snps.list$snps.assoc
  sets <- snps.list$sets

  }

  #################################
  ## Plot Tree showing Phenotype ##
  #################################

  ## SAVE TREE PLOT:
  if(!is.null(filename)){
    pdf(file=filename[[1]], width=7, height=11)
    # dev.copy(pdf, file=filename[[1]], width=7, height=11)
  }

  if(plot==TRUE){
    if(class(try(plot.phen(tree = tree,
                           phen.nodes = phen.nodes,
                           plot = plot))) =="try-error"){
      plot(tree)
      warning("Oops-- something went wrong when trying to plot
              phenotypic changes on tree.")
    }else{
      phen.plot.col <- plot.phen(tree = tree,
                                 phen.nodes = phen.nodes,
                                 plot = plot, main.title=FALSE)
    }

    ##################
    ## SET 2 CLADES ##
    ##################
    ## plot set2 clade sets along tips:
    if(!is.null(sets)){
      ## Get CLADES:
      set1 <- names(sets)[which(sets == 1)]
      set2 <- names(sets)[which(sets == 2)]

      tip.labs <- tree$tip.label
      if(coaltree == FALSE) tip.labs <- removeFirstN(tip.labs, 1) ## assuming tip.labs are prefaced w/ "t" for all rtrees...
      cladeCol <- rep(NA, length(tip.labs))
      cladeCol <- replace(cladeCol, which(tip.labs %in% set1), "black")
      cladeCol <- replace(cladeCol, which(tip.labs %in% set2), "grey")
      ## PLOT CLADES along tips:
      ## coaltree:
      if(coaltree == TRUE){
        tiplabels(text=NULL, cex=0.6, adj=c(0.65, 0.5), col=cladeCol, pch=15) # adj=c(0.65, 0.75) ## NOT SURE WHY/WHEN ADJ WORKS/w WHAT VALUES?????
      }else{
        ## rtree:
        tiplabels(text=NULL, cex=0.75, adj=c(0.65, 0.75), col=cladeCol, pch=15) # adj=c(0.65, 0.75) ## NOT SURE WHY/WHEN ADJ WORKS/w WHAT VALUES?????
      }
    }

    if(!is.null(filename)){
      ## end saving tree plot:
      dev.off()
    }

  }

  ################
  ## Get Output ##
  ################
  out <- list(snps, snps.assoc, phen, phen.plot.col, tree, sets)
  names(out) <- c("snps", "snps.assoc", "phen", "phen.plot.col", "tree", "sets")
  return(out)

} # end coalescent.sim
