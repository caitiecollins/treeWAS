
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
# n.ind <- 10 # n.genomes you want to end up with
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
# n.subs <- 1
# n.snps.assoc <- 10
# assoc.prob <- 100
# n.phen.subs <- 15
# phen <- NULL
# plot <- TRUE
# heatmap <- FALSE
# reconstruct <- FALSE
# dist.dna.model <- "JC69"
# seed <- 1


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
                           row.names=NULL,
                           seed=1){
  ## load packages:
  require(adegenet)
  require(ape)
  require(phangorn)

  if(length(which(c(plot, heatmap, reconstruct)==TRUE))==1){
    par(ask=FALSE)
  }else{
    par(ask=TRUE)
  }

  ################################
  ## Simulate Phylogenetic Tree ##
  ################################
  tree <- coalescent.tree.sim(n.ind = n.ind, seed = seed)

  ########################
  ## Simulate Phenotype ##
  ########################
  if(is.null(phen)){
    ## get list of phenotype simulation output
    phen.list <- phen.sim(tree, n.subs = n.phen.subs)

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

  #################################
  ## Plot Tree showing Phenotype ##
  #################################
  if(plot==TRUE){
    phen.plot.col <- plot.phen(tree = tree,
                              phen.nodes = phen.nodes,
                              plot = plot)
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
                       seed=seed)
  # )

  snps <- snps.list$snps
  snps.assoc <- snps.list$snps.assoc

  ################
  ## Get Output ##
  ################
  out <- list(snps, snps.assoc, phen, phen.plot.col, tree)
  names(out) <- c("snps", "snps.assoc", "phen", "phen.plot.col", "tree")
  return(out)

} # end coalescent.sim
