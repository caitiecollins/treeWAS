##############
## ace.test ##
##############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param data A matrix or data.frame.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#'
#' ## example: ##
#'
#' # load data
#' data("snps.ace")
#' data("phen.ace")
#' data("tree.ace")
#' snps <- snps.ace
#' phen <- phen.ace
#' tree <- tree.ace
#'
#' ## run ace.test on subset of data (otherwise SLOW!)
#' out <- ace.test(snps[,1:10], phen, tree, method="discrete")
#'
#' ## examine output
#' round(out, 4)
#'
#' @import phangorn

########################################################################



ace.test <- function(snps, phen, tree, method="discrete"){

  edges <- tree$edge

  ###############################
  ## get unique SNPs patterns: ##
  ###############################
  tab.out <- matrix.table(t(snps))
  snps.unique <- t(as.matrix(tab.out$unique.data))
  index <- tab.out$index
  if(ncol(snps.unique) == ncol(snps)){
    all.unique <- TRUE
  }else{
    all.unique <- FALSE
  }

  ## work w only unique snps:
  snps.ori <- snps
  snps <- snps.unique

  ##################################
  ## get diffs for (unique) SNPS: ##
  ##################################
  diffs <- list()
  # system.time == 105.042 elapsed
  # w unique ncol(snps) == 2500
  for(i in 1:ncol(snps)){
    ## get variable
    var <- snps[,i]
    var.terminal <- var
    snps.ace.d <- ace(var, tree, type=method)
    var.internal <- snps.ace.d$lik.anc[,2]
    var <- c(var.terminal, var.internal)

    ## get differences for this variable's ace likelihoods
    diffs[[i]] <- get.ace.diffs(var, edges)
  }


  #########################
  ## get diffs for PHEN: ##
  #########################
  ## do we need to check phen is numeric??
  phen.terminal <- phen
  phen.ace.d <- ace(phen, tree, type=method)
  phen.internal <- phen.ace.d$lik.anc[,2]
  var <- c(phen.terminal, phen.internal)

  phen.diffs <-  get.ace.diffs(var, edges)

  ###########################################
  ## compare diffs for each SNPs vs. PHEN: ##
  ###########################################
  diff.corrs <- list()
  for(i in 1:ncol(snps)){
    snps.diffs <- diffs[[i]]
    sp.diffs <- snps.diffs * phen.diffs

    diff.corrs[[i]] <- sum(sp.diffs)
  }
  diff.corrs <- as.vector(unlist(diff.corrs))

  ############################################
  ## get values for duplicate snps columns: ##
  ############################################
  if(all.unique == TRUE){
    diff.corrs.complete <- diff.corrs
  }else{
    diff.corrs.complete <- rep(NA, ncol(snps.ori))
    for(i in 1:ncol(snps.unique)){
      diff.corrs.complete[which(index == i)] <- diff.corrs[i]
    }
  }

  ###################
  ## return output ##
  ###################
  ## Return vector containing measure of
  ## correlation btw SNPi and PHENi diffs across nodes in tree
  ## for all snps in original snps matrix.
  return(diff.corrs.complete)

} # end ace.test


########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param var A vector containing a variable whose change across edges we want to examine.
#' @param edges A 2-column matrix containing the upstream and downstream nodes
#'  in columns 1 and 2 of a tree's edge matrix, as found in a phylo object's tree$edge slot.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#'
#' ## basic use of fn
#' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
#'

########################################################################

###################
## get.ace.diffs ##
###################
## get diffs btw ace prob/likelihood at upstream vs. downstream node
## for a single variable for which you have a value for all terminal and internal nodes.
get.ace.diffs <- function(var, edges){
  ## CHECKS: ##
  ## var should be a vector or have 2 columns summing to 1 or 100:
  if(!is.null(dim(var))){
    if(ncol(var) > 2) warning("var contains more than one discrete variable;
                              selecting first variable.")
    var <- var[,2]
  }
  ## var and tree$edge should contain the same number of inds:
  if(length(var) != (nrow(edges)+1)) stop("var contains more individuals than tree$edge does.")

  ## FOR LOOP ##
  diffs <- var[edges[,1]] - var[edges[,2]]

  return(as.vector(diffs))
} # end get.ace.diffs



########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param data A matrix or data.frame.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#'
#' ## basic use of fn
#' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
#'

########################################################################

##################
## matrix.table ##
##################
## Get just the unique rows of a matrix,
## the pattern/index to map duplicate rows to,
## and a table counting repeats of unique rows.
matrix.table <- function(data){
  ## get df
  if(!is.data.frame(data)){
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  ## concatenate all rows into single elements
  dat <- do.call("paste", c(data, sep = "\r"))
  ## get unique rows
  unique.inds <- !duplicated(dat)
  ## keep only unique rows
  levels <- dat[unique.inds]
  ## make a factor in which each level
  ## is the smushed single-element vector
  ## of each unique row
  cat <- factor(dat, levels = levels)
  ## get n. unique levels
  n.levels <- length(levels(cat))
  ## get the unique index that
  ## each original ind/row should map to:
  map.to <- (as.integer(cat) - 1)
  map.to <- map.to[!is.na(map.to)]
  if(length(map.to)) map.to <- map.to + 1
  ## get the number of inds at each unique level:
  tab <- tabulate(map.to, n.levels)
  ## get output
  out <- list(index = map.to,
              table = tab,
              unique.data = data[unique.inds,])
  return(out)
} # end matrix.table

