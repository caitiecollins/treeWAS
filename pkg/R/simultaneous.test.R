



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
#' @export
#'

########################################################################

simultaneous.test <- function(snps.reconstruction, # can be snps.REC OR snps.sim.REC matrix ## NOTE: subs.edges no longer required for any version of this test.
                              phen.reconstruction,
                              tree){

  snps.rec <- snps.reconstruction
  phen.rec <- phen.reconstruction

  ## Get tree edges:
  edges <- tree$edge

  #########################
  ## Get UNIQUE snps.rec ##
  #########################
  temp <- get.unique.matrix(snps.rec, MARGIN=2)
  snps.rec.unique <- temp$unique.data
  index <- temp$index

  if(ncol(snps.rec.unique) == ncol(snps.rec)){
    all.unique <- TRUE
  }else{
    all.unique <- FALSE
  }

  ## work w only unique snps:
  snps.rec.ori <- snps.rec
  snps.rec <- snps.rec.unique

  ###############################
  ## GET DIFFS ACROSS BRANCHES ##
  ###############################

  ## Get SNPs diffs: ##
  snps.diffs <- snps.rec[edges[,1], ] - snps.rec[edges[,2], ]

  ## Get phen diffs: ##
  phen.diffs <- phen.rec[edges[,1]] - phen.rec[edges[,2]]

  sp.diffs <- snps.diffs * phen.diffs
  score <- abs(colSums(sp.diffs, na.rm=TRUE))

  ################################################
  ## get values for duplicate snps.rec columns: ##
  ################################################

  ## get reconstruction for all original sites
  if(all.unique == TRUE){
    score.complete <- score
  }else{
    score.complete <- score[index]
    names(score.complete) <- colnames(snps.rec.ori)
  }

  score <- score.complete

  return(score)

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



