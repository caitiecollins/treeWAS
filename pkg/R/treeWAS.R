
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
#' @param x description.
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
#' @import adegenet ape phangorn ade4
#'
#' @export

########################################################################

## n.subs ## either an integer specifying Poisson parameter OR a distribution OR (if NULL?? or "fitch"??) compute w parsimony



treeWAS <- function(snps, phen, n.subs=NULL,
                    tree=c("UPGMA", "nj", "ml"), dist.dna.model="JC69", plot.tree=FALSE,
                    test="score", p.value=0.001, p.value.correct=c("bonf", "fdr"), p.value.by="count",
                    sim.n.snps=10000, n.reps=1,
                    plot.null.dist=TRUE, plot.dist=FALSE){

  ###################
  ## LOAD PACKAGES ##
  ###################
  require(adegenet)
  require(phangorn)
  require(ape)
  require(ade4) #?
  #require(adephylo)

  #####################################################################
  ## 0) HANDLE INPUT DATA #############################################
  #####################################################################

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
                             method=tree,
                             dist.dna.model=dist.dna.model,
                             plot=plot.tree)
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


  ## tree's tip.labels must be numeric (for Fitch parsimony step)
  tree.ori <- tree
  tree$tip.label <- c(1:length(tree$tip.label))
  tree$node.label <- c((n.ind+1):(n.ind+tree$Nnode))


  ###################
  ## HANDLE N.SUBS ##
  ###################

  ## if n.subs is an integer ##
  ## we use this as the parameter of a Poisson distribution
  ## from which the n.subs per site is drawn for each site.

  ## if n.subs is a vector (ie. distribution) ##
  ## we use this distribution directly (but in proportion with the number of sites)
  ## to specify the n.subs per site.

  ## if n.subs is NULL ##
  ## we compute the distribution of the n.subs-per-site
  ## using the Fitch parsimony score calculation fns from phangorn.

  if(is.null(n.subs)){
    n.subs <- get.fitch.n.mts(snps=snps, tree=tree)
    n.subs <- table(n.subs)
    noms <- as.numeric(names(n.subs))
    temp <- rep(0, max(noms))
    for(i in 1:max(noms)){
      if(i %in% noms) temp[i] <- n.subs[which(noms==i)]
    }
    n.subs <- temp
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
                        n.subs = n.subs, n.snps.assoc = 0, assoc.prob = 100,
                        tree = tree, phen.loci = NULL,
                        heatmap = FALSE,
                        reconstruct = FALSE, dist.dna.model = "JC69",
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
  #### Determine the PC-p-values for all SNP loci | ##############
  ##   null distributions of correlations from simulated data ####
  #### Synthesize results output: List of all significant SNPs, ##
  ##   their names/locations, their p-values for this phenotype ##
  ################################################################

  #######################
  ## identify sig.snps ##
  #######################

  sig.list <- get.sig.snps(snps=snps, snps.sim=snps.mat, phen=phen,
                           test=test,
                           p.value=p.value,
                           p.value.correct=p.value.correct,
                           p.value.by=p.value.by)

  #############################################
  ## isolate elements of get.sig.snps output ##
  #############################################

  corr.dat <- sig.list$corr.dat
  corr.sim <- sig.list$corr.sim
  sig.snps.names <- sig.list$sig.snps.names
  sig.snps <- sig.list$sig.snps
  sig.corrs <- sig.list$sig.corrs
  p.vals <- sig.list$p.vals
  min.p <- sig.list$min.p
  sig.thresh <- sig.list$sig.thresh

  ########################################


  ##################################
  ## 4) (A) Plot the distribution ##
  ##################################

  plot.sig.snps(corr.dat, corr.sim, sig.corrs, sig.snps,
                sig.thresh=sig.thresh, test,
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
    #                        which(dimnames(snps)[[2]] == sig.ssig.snps
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
                     p.vals,
                     sig.corrs,
                     S1P1, S0P0, S1P0, S0P1)
    names(df) <- c("SNP.locus",
                   "p.value",
                   "Test.statistic",
                   "S1P1", "S0P0", "S1P0", "S0P1")

    ## NOTE: Could return sig.snps.names somewhere here in addition to sig.snps loci ####    ####    ####    ####

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
  results[[3]] <- sig.thresh
  results[[4]] <- df
  results[[5]] <- min.p

  names(results) <- c("corr.dat",
                      "corr.sim",
                      "sig.thresh",
                      "sig.snps",
                      "min.p.value")

  return(results)

} # end treeWAS


