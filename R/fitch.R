
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
#' Determine parsimony scores for all genetic loci, or a phenotypic variable, along a given tree.
#' An extension of the fitch function available in package phangorn.
#'
#' @param x A numeric matrix or vector containing two unique values with row.names matching tree tip.labels.
#' @param tree A phylo object.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @examples
#'
#' ## generate a tree
#' tree <- ape::rtree(100)
#' ## generate snps, a matrix of 0s and 1s
#' snps <- matrix(sample(c(0,1),100000,TRUE), nrow=100)
#' row.names(snps) <- tree$tip.label
#'
#' ## run function
#' out <- get.fitch.n.mts(x=snps.mat, tree)
#'
#' ## examine output
#' str(out)
#' table(out)
#' hist(out)
#'
#' @importFrom phangorn fitch
#' @importFrom phangorn as.phyDat
#'
#' @export

########################################################################
# @useDynLib phangorn, .registration = TRUE


get.fitch.n.mts <- function(x, tree, snps=NULL){

  ## load packages
  # require(phangorn)

  ## Re-coding snps as x (to allow for phen/vectors).
  ## --> snps now deprecated:
  X <- NULL
  if(!missing(x)){
    X <- x
    if(!is.null(snps) & !is.null(x)){
      warning("As 'x' is specified, we ignore the 'snps' argument. \n
              (In get.fitch.n.mts the 'snps' argument has now been replaced by an argument named 'x'.)")
    }
  }else{
    if(!is.null(snps)){
      X <- snps
    }
  }
  ## If ONE of x or snps was specified, continue; else, stop:
  if(!is.null(X)){
    x <- X
  }else{
    stop("'x' must be specified.")
  }

  ## checks
  # if(!is.matrix(x)) stop("x must be a matrix.") ## no it mustn't... works for phen too.
  # levs <- unique(as.vector(unlist(x)))
  levs <- unique(as.vector(as.matrix(x)))
  # if(any(is.na(levs))) levs <- levs[-which(is.na(levs))]
  if((!is.numeric(x) & !is.logical(x)) | length(levs[!is.na(levs)])!=2){
    stop("x must be a numeric matrix or vector, with two unique values, excluding NAs
         (though we recommend that NAs be in the minority for each column).")
  }
  if(any(is.na(levs))){
    nnas <- sapply(c(1:ncol(snps)), function(e) length(which(is.na(snps[,e])))/nrow(snps))
    toRemove <- which(nnas > 0.5)
    if(length(toRemove) > 0){
      cat(length(toRemove), "snps columns are over 50% NAs.
          You may want to consider removing these columns as they are
          unlikely to be significant and can generate inappropriate inferences.")
    }
  }

  x.levels <- sort(levs, na.last = TRUE)
  ## returns only unique patterns...
  x.phyDat <- as.phyDat(as.matrix(x),
                           type="USER", levels=x.levels)
  ## get index of all original x columns to map to unique pattern
  index <- attr(x.phyDat, "index")

  fitch.phangorn <- phangorn::fitch
  ## get parsimony score for all unique patterns in x
  ## NB: For fitch.phangorn, x data must be of class phyDat
  fitch.unique <- fitch.phangorn(tree, x.phyDat, site="site")
  # table(fitch.unique)

  ## get score for all original sites
  fitch.complete <- fitch.unique[index]
  return(fitch.complete)
} # end get.fitch.n.mts


