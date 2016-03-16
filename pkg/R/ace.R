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
#' @param snps A matrix containing binary snps for all individuals.
#' @param phen A factor containing the phenotype for which to test for association.
#' @pram tree A phylo object containing the tree representing the ancestral relationships
#' between individuals for which snps and phen are known.
#' @param method A character string specifying the type of ACE method to implement.
#' @param snps.ace A logical indicating whether to run ACE on all snps (TRUE, slower)
#' or to run parsimony on snps instead (FALSE, faster; the default).
#'
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



ace.test <- function(snps, phen, tree, method="discrete", snps.ace=FALSE){

  edges <- tree$edge

  ###############################
  ## get unique SNPs patterns: ##
  ###############################
  tab.out <- table.matrix(t(snps))
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

  ##########################################
  ## RUN ACE on PHEN & PARSIMONY on SNPS: ##   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ##########################################
  if(snps.ace == FALSE){

    ## get n.subs per site:
    cost <- get.fitch.n.mts(snps, tree)

    #############################################
    ## get SNP states of all internal nodes ?? ##
    #############################################

    #########
    ## MPR ##
    #########

    ## pace == ancestral.pars
    #pa.MPR <- pace(tree, dna, type="MPR")

    #############
    ## ACCTRAN ##
    #############

    ## pace == ancestral.pars
    pa.ACCTRAN <- pace(tree, dna, type="ACCTRAN")

    ## pace  --> diff resuls w MPR vs. ACCTRAN
    diffs <- sapply(c(1:length(pa.ACCTRAN)), function(e) identical(pa.MPR[[e]], pa.ACCTRAN[[e]]))

    ## TO DO -- STILL NEED TO FIGURE OUT HOW TO HANDLE MPR/ACCTRAN RESULTS W DIFFERENT UPPER AND LOWER ESTIMATES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ###########################################
    ## convert reconstruction back to snps.. ##
    ###########################################
    ## each of the n.ind elements of pa is a matrix w n.snps rows and 4 columns, each for the 4 nts possible (acgt)

    # rec <- pa.MPR
    rec <- pa.ACCTRAN

    # str(rec)
    snps.rec <- list()
    for(i in 1:length(rec)){
      snps.rec[[i]] <- rec[[i]][,3]
    }
    snps.rec <- do.call("rbind", snps.rec)
    rownames(snps.rec) <- c(rownames(snps), c((nrow(snps)+1):((nrow(snps)*2)-1)))
    colnames(snps.rec) <- colnames(snps)
    # identical(snps.rec[1:nrow(snps),], snps)


    ##############################################
    ## get LOCATIONS (branches) of snps subs ?? ##
    ##############################################
    subs.edges <- rep(list(NULL), ncol(snps))
    for(i in 1:ncol(snps)){
      snp <- snps.rec[, i]
      subs.logical <- sapply(c(1:nrow(edges)),
                             function(e)
                               snp[edges[e,1]]
                             ==
                               snp[edges[e,2]])
      ## get indices of all edges containing a substitution
      subs.total <- which(subs.logical == FALSE)
      ## get df of states of ancestor and descendants nodes on these edges
      df <- data.frame(snp[edges[subs.total,1]], snp[edges[subs.total,2]])
      names(df) <- c("anc", "dec")
      ## get indices of all edges w a positive sub (0 --> 1)
      subs.pos <- subs.total[which(df$anc==0)]
      ## get indices of all edges w a negative sub (1 --> 0)
      subs.neg <- subs.total[which(df$anc==1)]

      ## get output list
      subs.edges[[i]] <- rep(list(NULL), 3)
      names(subs.edges[[i]]) <- c("total", "pos", "neg")
      if(length(subs.total) > 0) subs.edges[[i]][["total"]] <- subs.total
      if(length(subs.pos) > 0) subs.edges[[i]][["pos"]] <- subs.pos
      if(length(subs.neg) > 0) subs.edges[[i]][["neg"]] <- subs.neg
    }

    cost3 <- sapply(c(1:length(subs.edges)), function(e) length(subs.edges[[e]]))

    ## check:
    ## for SNP1, does it identify the correct/reasonable branches?
    ## see plot w edgeCol2--reasonable yes, but not the real simulated branches...
    # temp <- rep(0, nrow(edges))
    # temp <- replace(temp, subs.edges[[1]], 1)

    #######################################
    ## test for association w phen (ACE) ##
    #######################################
    ace.score <- list()
    ## NB: length(subs.edges) == ncol(snps.unique)
    for(i in 1:length(subs.edges)){
      ## get the absolute value of
      ## the sum of all the positive (0-->1) and negative (1-->0) subs
      ace.score[[i]] <- abs(sum(phen.diffs[subs.edges[[i]][["pos"]]])) +
        abs(sum(phen.diffs[subs.edges[[i]][["neg"]]]))
    }
    ace.score <- as.vector(unlist(ace.score))
  }

  ####################################
  ## RUN ACE on BOTH SNPS AND PHEN: ##   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ####################################

  if(snps.ace == TRUE){
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
    ace.score <- list()
    for(i in 1:ncol(snps)){
      snps.diffs <- diffs[[i]]
      sp.diffs <- snps.diffs * phen.diffs

      ace.score[[i]] <- sum(sp.diffs)
    }
    ace.score <- as.vector(unlist(ace.score))
  }


  ############################################
  ## get values for duplicate snps columns: ##
  ############################################
  if(all.unique == TRUE){
    ace.score.complete <- ace.score
  }else{
    ace.score.complete <- rep(NA, ncol(snps.ori))
    for(i in 1:ncol(snps.unique)){
      ace.score.complete[which(index == i)] <- ace.score[i]
    }
  }

  ###################
  ## return output ##
  ###################
  ## Return vector containing measure of
  ## correlation btw SNPi and PHENi diffs across nodes in tree
  ## for all snps in original snps matrix.
  return(ace.score.complete)

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
## table.matrix ##
##################
## Get just the unique rows of a matrix,
## the pattern/index to map duplicate rows to,
## and a table counting repeats of unique rows.
table.matrix <- function(data){
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
} # end table.matrix

