
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
#' @param dna A binary matrix or DNAbin object containing genomes for (only)
#'                the terminal nodes of the tree to be reconstructed.
#'                Individuals should be in the rows and loci in the columns; rows and columns should be labelled.
#' @param dist.dna.model A character string specifying the type of model to use in
#'                          calculating the genetic distance between individual genomes (see ?dist.dna).
#' @param plot A logical specifying whether to plot the reconstructed phylogenetic tree.
#'
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
#' @import ape phangorn

########################################################################

############
## TO DO: ##
############
## add all the options from hclust (stats) as available methods..
## change all methods to either upper or lower case (or add to lower check).


tree.reconstruct <- function(dna,
                             method=c("UPGMA", "NJ", "ML"),
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
  if(class(dna) != "DNAbin"){
    # dna <- as.DNAbin(dna)
    snps <- dna
    sp <- as.character(snps)
    sp <- replace(sp, which(sp == "0"), "a")
    sp <- replace(sp, which(sp == "1"), "c")
    sp <- matrix(sp, ncol=ncol(snps), nrow=nrow(snps))
    dna <- as.DNAbin(sp)
    rownames(dna) <- rownames(snps)
    colnames(dna) <- colnames(snps)
  }
  method <- tolower(method)
  if(!any(c("upgma", "nj", "ml") %in% method)){
    warning("method should be one of 'UPGMA', 'NJ', 'ML'. Choosing 'UPGMA'.")
  }

  tree <- NULL
  ## TO DO: ADD PHENOTYPE-COLORING OPTIONS TO RECONSTRUCTED PHYLO PLOTS

  D <- dist.dna(dna, model = dist.dna.model)

  ###########
  ## UPGMA ##
  ###########
  if(method=="upgma"){
    tree <- hclust(D, method="average")
    tree <- as.phylo(tree)
    #tree <- midpoint(ladderize(tree))
    ## Always work with tree in pruningwise order:
    tree <- reorder.phylo(tree, order="pruningwise")
    ## Trees must be rooted:
    if(!is.rooted(tree)) tree <- midpoint(tree)
    if(plot==TRUE){
      plot(tree, main="")
      title("UPGMA tree")
    }
  }
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
      plot(tree, edge.width=2)
      title("Neighbour-joining tree")
      axisPhylo()
    }
  }
  ########
  ## ML ##
  ########
  if(method=="ml"){
    dna4 <- as.phyDat(dna)
    tre.ini <- nj(dist.dna(dna, model=dist.dna.model))
    fit.ini <- pml(tre.ini, dna4, k=n.ind)
    fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE,
                     optQ = TRUE, optGamma = TRUE)

    ## NOTE--you may want to store these in a results.ml list
    ## and return it with your results instead of printing
    ## OR at least print a message
    ## (eg. "Printing maximum-likelihood calculations...")
    ## before printing these numbers...

    #     anova(fit.ini, fit)
    #     AIC(fit.ini)
    #     AIC(fit)

    tree <- fit$tree
    #tree <- midpoint(ladderize(tree))
    ## Always work with tree in pruningwise order:
    tree <- reorder.phylo(tree, order="pruningwise")
    ## Trees must be rooted:
    if(!is.rooted(tree)) tree <- midpoint(tree)
    if(plot==TRUE){
      plot(tree, show.tip=TRUE, edge.width=2)
      title("Maximum-likelihood tree")
      axisPhylo()
    }
  }
  par(ask=FALSE)

  return(tree)
} # end tree.reconstruct
