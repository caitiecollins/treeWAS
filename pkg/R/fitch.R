
#####################
## get.fitch.n.mts ##
#####################
## phangorn-based fitch fn

########################################################################

###################
## DOCUMENTATION ##
###################

#' Caclulate parsimony scores.
#'
#' Determine parsimony scores for all genetic loci and a given tree.
#' An extension of the fitch function available in package phangorn.
#'
#' @param snps A numeric matrix containing two unique values with row.names matching tree tip.labels.
#' @param tree A phylo object.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @examples
#'
#' ## generate a tree
#' tree <- rtree(100)
#' ## geerate snps, a matrix of 0s and 1s
#' snps <- matrix(sample(c(0,1),100000,T), nrow=100)
#' row.names(snps) <- tree$tip.label
#'
#' ## run function
#' out <- get.fitch.n.mts(snps, tree)
#'
#' ## examine output
#' str(out)
#' table(out)
#' hist(out)
#'
#' @import phangorn
#'
#' @export

########################################################################

get.fitch.n.mts <- function(snps, tree){

  ## load packages
  # require(phangorn)

  ## checks
  if(!is.matrix(snps)) stop("snps must be of class matrix.")
  if(!is.numeric(snps) | length(unique(as.vector(snps)))!=2){
    stop("snps must be a numeric matrix with exactly two unique values.")
  }
  snps.levels <- unique(as.vector(snps))
  ## returns only unique patterns...
  snps.phyDat <- as.phyDat(as.matrix(snps),
                           type="USER", levels=snps.levels)
  ## get index of all original snps columns to map to unique pattern
  index <- attr(snps.phyDat, "index")

  fitch.phangorn <- phangorn::fitch
  ## get parsimony score for all unique patterns in snps
  ## NB: For fitch.phangorn, snps data must be of class phyDat
  fitch.unique <- fitch.phangorn(tree, snps.phyDat, site="site")
  # table(fitch.unique)

  ## get score for all original sites
  fitch.complete <- fitch.unique[index]
  return(fitch.complete)
} # end get.fitch.n.mts




