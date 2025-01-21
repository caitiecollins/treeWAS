
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
#' \dontrun{
#'
#' ## generate a tree
#' tree <- ape::rtree(100)
#' ## generate snps, a matrix of 0s and 1s
#' snps <- matrix(sample(c(0,1),100000,TRUE), nrow=100)
#' row.names(snps) <- tree$tip.label
#'
#' ## run function
#' out <- get.fitch.n.mts(x=snps, tree)
#'
#' ## examine output
#' str(out)
#' table(out)
#' hist(out)
#' }
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
  ## do not include NA as a level:
  levs <- unique(as.vector(x[!is.na(x)]))
  if((!is.numeric(x) & !is.logical(x)) | length(levs[!is.na(levs)])!=2){
    stop("x must be a numeric matrix or vector, with two unique values, excluding NAs
         (though we recommend that NAs be in the minority for each column).\n")
  }
  # levs <- unique(as.vector(x))
  if(any(is.na(levs))){
    if(is.matrix(x)){
      nnas <- sapply(c(1:ncol(x)), function(e) length(which(is.na(x[,e])))/nrow(x))
      toRemove <- which(nnas > 0.5)
      if(length(toRemove) > 0){
        cat(length(toRemove), "snps columns are over 50% NAs.
            You may want to remove these columns as they are unlikely to be significant
            and can generate inappropriate inferences during ancestral state reconstruction.\n")
      }
    }else{
      nnas <- length(which(is.na(x)))/length(x)
      # toRemove <- which(nnas > 0.5)
      if(nnas > 0.5){
        cat("x is over 50% NAs.
            This may generate inappropriate inferences during ancestral state reconstruction.\n")
      }
    }
  }

  x.levels <- sort(levs, na.last = TRUE)
  ## returns only unique patterns...
  ## *use levels=states (eg. c(0,1)), but keep NAs in x and use ambiguity=NA
  ## to allow NAs without counting them twd parsimony score values.
  x.phyDat <- phangorn::as.phyDat(as.matrix(x),
                           type="USER", levels=x.levels, ambiguity=NA)
  ## get index of all original x columns to map to unique pattern
  index <- attr(x.phyDat, "index")

  ## get parsimony score for all unique patterns in x
  ## NB: For phangorn::fitch, x data must be of class phyDat
  fitch.unique <- phangorn::fitch(tree, x.phyDat, site="site")
  # table(fitch.unique)

  ## get score for all original sites
  fitch.complete <- fitch.unique[index]
  return(fitch.complete)
} # end get.fitch.n.mts


