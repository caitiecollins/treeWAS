



#######################
## simultaneous.test ##
#######################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param snps.reconstruction A matrix containing the terminal and reconstructed
#' ancestral states of SNPs for all nodes in the tree.
#' @param phen.reconstruction A vector containing the terminal and reconstructed
#' ancestral states of the phenotype for all nodes in the tree.
#' @param tree A phylo object containing the tree representing the ancestral relationships
#' between the individuals for which snps and phen are known.
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#'
#' @importFrom scales rescale
#' @importFrom Hmisc all.is.numeric
#'
#' @export

########################################################################
# @useDynLib phangorn, .registration = TRUE
# @importFrom phangorn midpoint

simultaneous.test <- function(snps.reconstruction, # can be snps.REC OR snps.sim.REC matrix ## NOTE: subs.edges no longer required for any version of this test.
                              phen.reconstruction,
                              tree){

  snps.rec <- snps.reconstruction
  phen.rec <- phen.reconstruction
  rm(snps.reconstruction)
  rm(phen.reconstruction)

  ## Always work with tree in pruningwise order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  # if(!is.rooted(tree)) tree <- midpoint(tree) # require(phangorn)
  ## Get tree edges:
  edges <- tree$edge

  ####################################################################
  #####################
  ## Handle phen.rec ##
  #####################
  ## convert phenotype to numeric:
  phen.rec.ori <- phen.rec
  ## Convert to numeric (required for assoc tests):
  na.before <- length(which(is.na(phen.rec)))

  ## NB: can only be binary or continuous at this point...
  levs <- unique(as.vector(unlist(phen.rec)))
  n.levs <- length(levs[!is.na(levs)])
  if(n.levs == 2){
    if(!is.numeric(phen.rec)){
      if(all.is.numeric(phen.rec)){
        phen.rec <- as.numeric(as.character(phen.rec))
      }else{
        phen.rec <- as.numeric(as.factor(phen.rec))
      }
    }
  }else{
    if(!is.numeric(phen.rec)){
      if(all.is.numeric(phen.rec)){
        phen.rec <- as.numeric(as.character(phen.rec))
      }else{
        stop("phen.rec has more than 2 levels but is not numeric (and therefore neither binary nor continuous).")
      }
    }
  }
  ## ensure ind names not lost
  names(phen.rec) <- names(phen.rec.ori)

  ## Check that no errors occurred in conversion:
  na.after <- length(which(is.na(phen.rec)))
  if(na.after > na.before){
    stop("NAs created while converting phen.rec to numeric.")
  }
  ####################################################################

  ################################################
  ## RE-SCALE NON-BINARY VALUES (phen only ...) ##
  ################################################
  ## phen.rec (both Pa and Pd should be on same scale):
  # if(n.levs > 2)
  phen.rec <- rescale(phen.rec, to=c(0,1)) # require(scales)

  ###############################
  ## GET DIFFS ACROSS BRANCHES ##
  ###############################

  ## Get SNPs diffs: ##
  snps.diffs <- snps.rec[edges[,1], ] - snps.rec[edges[,2], ]

  ## Get phen diffs: ##
  phen.diffs <- phen.rec[edges[,1]] - phen.rec[edges[,2]]

  sp.diffs <- snps.diffs * phen.diffs

  ## Return with sign:
  score2 <- colSums(sp.diffs, na.rm=TRUE)
  # score2 <- abs(score2)
  names(score2) <- colnames(snps.rec)

  ## Return with and without sign?
  # score2.sign <- colSums(sp.diffs, na.rm=TRUE)
  # score2 <- abs(score2.sign)
  # names(score2) <- names(score2.sign) <- colnames(snps.rec)
  #
  # score2 <- list("score2" = score2,
  #                "score2.sign" = score2.sign)

  return(score2)

} # end simultaneous.test





# ######################
# ## get.branch.diffs ##
# ######################
# Not actually needed to get score 2!
# ########################################################################
#
# ###################
# ## DOCUMENTATION ##
# ###################
#
# #' Short one-phrase description.
# #'
# #' Longer proper discription of function...
# #'
# #' @param var A vector containing a variable whose change across edges we want to examine.
# #' @param edges A 2-column matrix containing the upstream and downstream nodes
# #'  in columns 1 and 2 of a tree's edge matrix, as found in a phylo object's tree$edge slot.
# #'
# #' @author Caitlin Collins \email{caitiecollins@@gmail.com}
# #' @export
# #'
#
# ########################################################################
#
# ## get diffs btw ace prob/likelihood at upstream vs. downstream node
# ## for a single variable for which you have a value for all terminal and internal nodes.
# get.branch.diffs <- function(var, edges){
#   ## CHECKS: ##
#   ## var should be a vector or have 2 columns summing to 1 or 100:
#   if(!is.null(dim(var))){
#     if(ncol(var) > 2) warning("var contains more than one discrete variable;
#                               selecting first variable.")
#     var <- var[,2]
#   }
#   ## var and tree$edge should contain the same number of inds:
#   if(length(var) != (nrow(edges)+1)) stop("var contains more
#                                           individuals than tree$edge does.")
#
#   ## ~ FOR LOOP ##
#   diffs <- var[edges[,1]] - var[edges[,2]]
#
#   return(as.vector(diffs))
# } # end get.branch.diffs



