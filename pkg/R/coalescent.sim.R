
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
# n.snp.assoc <- 5
# assoc.option <- c("all", "model") # deprecated (only "all" available)
# sim.by <- c("locus", "branch") # deprecated (only "locus" has all current protocols implemented)


## EXAMPLE ##
# out <- coalescent.sim(n.ind=100, gen.size=10000, sim.by="locus",
#                       theta=1*2, dist=NULL,
#                       theta_p=15, phen=NULL,
#                       n.snp.assoc=20, assoc.option="all", assoc.prob=90,
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
#' @param n.ind description.
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
#  n.snp.assoc=5, assoc.option="all", assoc.prob=90,
#  haploid=TRUE, biallelic=TRUE, seed=NULL,
#  plot=TRUE, heatmap=FALSE, plot2="UPGMA")


coalescent.sim <- function(n.ind=100,
                           n.snps=10000, n.subs=1,
                           n.snp.assoc=10, assoc.prob=100,
                           n.phen.subs=15, phen=NULL,
                           plot=TRUE,
                           heatmap=FALSE, reconstruct=FALSE,
                           seed=1){

  require(adegenet)
  require(ape)
  require(phangorn)

  #   ## load utils.R for selectBiallelicSNP switch.phen .is.integer0 .is.even .is.odd
  #   source("C:/Users/Caitlin/treeWAS/pkg/R/utils.R")

  if(plot==TRUE && heatmap==FALSE && reconstruct==FALSE){
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

    ## get phenotype for all edges in tree
    ## (ie. "A", "B", or "A"-and-"B")
    phen.edges <- phen.list$phen.edges

    ## get the indices of phen.subs (ie. branches)
    phen.loci <- phen.list$phen.loci
  }else{
    #############################
    ## User-provided Phenotype ##
    #############################
    phen.nodes <- phen
    phen.edges <- NULL
  }

  #################################
  ## Plot Tree showing Phenotype ##
  #################################

  if(plot==TRUE){
    phen.plot.col <- plot.phen(tree = tree,
                              phen.nodes = phen.nodes,
                              phen.edges = phen.edges,
                              plot = plot)
  }

  ###################
  ## Simulate SNPs ##
  ###################
  snps.list <- snp.sim(n.ind=n.ind,
                       n.snps=n.snps, n.subs=n.subs,
                       n.snp.assoc=n.snp.assoc, assoc.prob=assoc.prob,
                       tree=tree,
                       phen.loci=phen.loci,
                       heatmap=heatmap, reconstruct=reconstruct,
                       seed=seed)
  snps <- snps.list$snps
  snps.assoc <- snps.list$snps.assoc

  ################
  ## Get Output ##
  ################

  out <- list(snps, snps.assoc, phen, phen.plot.col, tree)
  names(out) <- c("snps", "snps.assoc", "phen", "phen.plot.col", "tree")
  return(out)

} # end coalescent.sim
