
######################
## tree.reconstruct ##
######################


########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param dna A matrix or DNAbin object containing genomes for (only)
#'                the terminal nodes of the tree to be reconstructed.
#'                Individuals should be in the rows and loci in the columns; rows and columns should be labelled.
#' @param method A character string specifying the method of phylogenetic reconstruction:
#'                one of \code{"NJ"}, \code{"BIONJ"} (the default), or \code{"parsimony"};
#'                or, if NAs are present in the distance matrix, one of: \code{"NJ*"} or \code{"BIONJ*"}.
#' @param dist.dna.model A character string specifying the type of model to use in
#'                          calculating the genetic distance between individual genomes (see ?dist.dna).
#' @param plot A logical specifying whether to plot the reconstructed phylogenetic tree.
#'
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' 
#'
#' @import ape
#' @importFrom phangorn as.phyDat
#' @importFrom phangorn midpoint
#' @importFrom phangorn pratchet
#' @importFrom phangorn acctran
#'
#'
#' @export

########################################################################
# @import phangorn
# @useDynLib phangorn as.phyDat midpoint optim.pml pml, .registration = TRUE
# @useDynLib phangorn, .registration = TRUE
# @importFrom phangorn optim.pml
# @importFrom phangorn pml
# @importFrom phangorn optim.parsimony

############
## TO DO: ##
############
## add all the options from hclust (stats) as available methods..?
## change all methods to either upper or lower case (or add to lower check).


tree.reconstruct <- function(dna,
                             method= c("BIONJ", "NJ", "parsimony", "BIONJ*", "NJ*"),
                             dist.dna.model="JC69",
                             plot=TRUE){

  ###################
  ## LOAD PACKAGES ##
  ###################
  # require(ape)
  # require(phangorn)

  ############
  ## CHECKS ##
  ############
  ## DNA ##
  if(class(dna) != "DNAbin"){
    # if(class(dna) == "genind"){
    #   dna <- dna@tab ## might be problems w ploidy...
    # }
    if(is.matrix(dna)){
      # dna <- as.DNAbin(dna)
      sp <- matrix(as.character(dna), nrow=nrow(dna), ncol=ncol(dna))
      rownames(sp) <- rownames(dna)
      colnames(sp) <- colnames(dna)

      ## Check/convert levels:
      levs <- unique(as.vector(unlist(sp)))
      nts <- c("a", "c", "g", "t")
      if(length(levs[!is.na(levs)]) > 4){
        stop("There must be no more than 4 unique values in dna, excluding NAs.")
      }
      if(!all(levs %in% nts)){
        for(i in 1:length(levs)){
          sp <- replace(sp, which(sp == levs[i]), nts[i])
        } # end for loop
      } # end levs conversion

      dna <- as.DNAbin(sp)
      rownames(dna) <- rownames(sp)
      colnames(dna) <- colnames(sp)
    }else{
      stop("dna should be of class DNAbin or matrix")
    }
  }
  ## TREE REC METHOD ##
  method <- tolower(method)
  if(method == "njs") method <- "nj*"
  if(method == "bionjs") method <- "bionj*"
  method <- match.arg(arg = method,
                      choices = c("bionj", "nj", "parsimony", "nj*", "bionj*"),
                      several.ok = FALSE)
  if(!any(c("nj", "bionj", "parsimony", "nj*", "bionj*") %in% method)){
    warning("method should be one of  'nj', 'bionj', 'parsimony', 'nj*', 'bionj*'. Choosing 'BIONJ'.")
    method <- "bionj"
  }else{
    ## use first arg if more than 1 present:
    if(length(method) > 1){
      method <- method[1]
    }
  }

  # if(method == "upgma"){
  #   warning("UPGMA enforces ultrametricity, which can bias treeWAS results.
  #           NJ or BIONJ may give more reliable results.")
  # }

  tree <- NULL

  ## TO DO: ADD PHENOTYPE-COLORING OPTIONS TO RECONSTRUCTED PHYLO PLOTS

  ## Get distance matrix:
  D <- dist.dna(dna, model = dist.dna.model)

  ## Handle MISSING data:
  ## NOTE: hclust not able to handle NAs/NaNs in D..
  ## NB: NAs are ok in dna, but NAs in D arise when dist.dna cannot find a dist btw. any 2 individuals,
  ## e.g., there is an NA at all loci in at least one of the 2 inds.
  ## --> Use phylo methods that can handle NAs in D (eg. NJ* and BIONJ*, from ape).
  if(any(is.na(D))){
    if(!method %in% c("nj*", "bionj*")){
      if(method == "nj"){
        method <- "nj*"
      }else{
        method <- "bionj*"
      }
      cat("NAs in distance matrix. Replacing method of phylo estimation with ", method, ".", sep="")
    }
  }

  ###################################
  ## Methods with NO missing data: ##     #####     #####     #####     #####     #####     #####     #####     #####     #####
  ###################################

  ########
  ## NJ ##
  ########
  if(method=="nj"){
    tree <- nj(D)
    #tree <- midpoint(ladderize(tree))
    ## Always work with tree in pruningwise order:
    tree <- reorder.phylo(tree, order="pruningwise")
    ## Trees must be rooted:
    if(!is.rooted(tree)) tree <- midpoint(tree)
    if(plot==TRUE){
      plot(tree, edge.width=2, cex=0.5)
      title("Neighbour-joining tree")
      axisPhylo()
    }
  }
  ###########
  ## BIONJ ##
  ###########
  if(method=="bionj"){
    tree <- bionj(D)
    #tree <- midpoint(ladderize(tree))
    ## Always work with tree in pruningwise order:
    tree <- reorder.phylo(tree, order="pruningwise")
    ## Trees must be rooted:
    if(!is.rooted(tree)) tree <- midpoint(tree)
    if(plot==TRUE){
      plot(tree, edge.width=2, cex=0.5)
      title("BIONJ tree")
      axisPhylo()
    }
  }

  ###########
  ## UPGMA ##
  ###########
  # if(method=="upgma"){
  #   tree <- hclust(D, method="average")
  #   tree <- as.phylo(tree)
  #   #tree <- midpoint(ladderize(tree))
  #   ## Always work with tree in pruningwise order:
  #   tree <- reorder.phylo(tree, order="pruningwise")
  #   ## Trees must be rooted:
  #   if(!is.rooted(tree)) tree <- midpoint(tree)
  #   if(plot==TRUE){
  #     plot(tree, edge.width=2, cex=0.5)
  #     title("UPGMA tree")
  #   }
  # }

  ###############
  ## parsimony ## ## a bit slow if many unique columns
  ###############
  if(method=="parsimony"){
    ## as.phyDat warns if NAs present (& doesn't include these...)
    dna4 <- suppressWarnings(as.phyDat(dna))
    ## get pars tree:
    # tre.ini <- nj(D)
    # tree <- optim.parsimony(tre.ini, dna4)
    tre.ini <- pratchet(dna4, trace=0) # better (can also return set of treeS)
                        # , maxit=400, k=8)
    ## add edge lengths w ACCTRAN:
    tree <- acctran(tre.ini, dna4) # edge lengths in n.subs (but relative lengths still fine).
    ## Always work with tree in pruningwise order:
    tree <- reorder.phylo(tree, order="pruningwise")
    ## Convert edge.lengths from parsimony cost to n.subs-per-site
    ## (s.t. parsimony lengths ~ lengths via NJ or UPGMA):
    tree$edge.length <- tree$edge.length/ncol(dna)
    ## Trees must be rooted:
    if(!is.rooted(tree)) tree <- midpoint(tree)
    if(plot==TRUE){
      plot(tree, edge.width=2, cex=0.5)
      title("Parsimony tree")
      axisPhylo()
    }
  }

  ########
  ## ML ##  ## discontinued: too slow
  ########
  # if(method=="ml"){
  #   dna4 <- suppressWarnings(as.phyDat(dna))
  #   tre.ini <- nj(D)
  #   fit.ini <- pml(tre.ini, dna4, k=nrow(dna))
  #   fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE,
  #                    optQ = TRUE, optGamma = TRUE)
  #
  #   ## NOTE--you may want to store these in a results.ml list
  #   ## and return it with your results instead of printing
  #   ## OR at least print a message
  #   ## (eg. "Printing maximum-likelihood calculations...")
  #   ## before printing these numbers...
  #
  #   #     anova(fit.ini, fit)
  #   #     AIC(fit.ini)
  #   #     AIC(fit)
  #
  #   tree <- fit$tree
  #   #tree <- midpoint(ladderize(tree))
  #   ## Always work with tree in pruningwise order:
  #   tree <- reorder.phylo(tree, order="pruningwise")
  #   ## Trees must be rooted:
  #   if(!is.rooted(tree)) tree <- midpoint(tree)
  #   if(plot==TRUE){
  #     plot(tree, show.tip=TRUE, edge.width=2)
  #     title("Maximum-likelihood tree")
  #     axisPhylo()
  #   }
  # }


  ######################################
  ## Methods with MISSING data (in D) ##     #####     #####     #####     #####     #####     #####     #####     #####     #####
  ######################################

  #########
  ## NJ* ##
  #########
  if(method=="nj*"){
    tree <- njs(D)
    #tree <- midpoint(ladderize(tree))
    ## Always work with tree in pruningwise order:
    tree <- reorder.phylo(tree, order="pruningwise")
    ## Trees must be rooted:
    if(!is.rooted(tree)) tree <- midpoint(tree)
    if(plot==TRUE){
      plot(tree, edge.width=2)
      title("Neighbour-joining* tree")
      axisPhylo()
    }
  }
  ############
  ## BIONJ* ##
  ############
  if(method=="bionj*"){
    tree <- bionjs(D)
    #tree <- midpoint(ladderize(tree))
    ## Always work with tree in pruningwise order:
    tree <- reorder.phylo(tree, order="pruningwise")
    ## Trees must be rooted:
    if(!is.rooted(tree)) tree <- midpoint(tree)
    if(plot==TRUE){
      plot(tree, edge.width=2)
      title("BIONJ* tree")
      axisPhylo()
    }
  }


  par(ask=FALSE)

  return(tree)
} # end tree.reconstruct

