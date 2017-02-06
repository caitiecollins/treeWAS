
#####################
## subsequent.test ## ## NEW ORIGINAL SCORE 3 (w integral score, no edge-length) ##
#####################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param tree A phylo object.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#'
#' ## basic use of fn
#' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
#'
#' @importFrom scales rescale
#'

########################################################################

subsequent.test <- function(snps.reconstruction,
                            phen.reconstruction,
                            tree){

  snps.rec <- snps.reconstruction
  phen.rec <- phen.reconstruction

  ## Always work with tree in pruningwise order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  if(!is.rooted(tree)) tree <- midpoint(tree)
  ## get tree edges:
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

  ################################################
  ## RE-SCALE NON-BINARY VALUES (phen only ...) ##
  ################################################
  ## phen.rec (both Pa and Pd should be on same scale):
  phen.rec <- rescale(phen.rec, to=c(0,1)) # require(scales)

  ###############################
  ## GET SCORE ACROSS BRANCHES ##
  ###############################

  Pa <- phen.rec[edges[,1]]
  Pd <- phen.rec[edges[,2]]
  Sa <- snps.rec[edges[,1], ]
  Sd <- snps.rec[edges[,2], ]
  bl <- tree$edge.length

  #################################################################     #####
  ###############
  ## SCORE 3.0 ##
  ###############
  ## ORIGINAL AND NEW INTEGRAL-BASED SCORE3 (without edge length):
  score3 <- get.score3(Pa = Pa, Pd = Pd, Sa = Sa, Sd = Sd, l = NULL)

  score3 <- abs(colSums(score3, na.rm=TRUE))
  names(score3) <- colnames(snps.rec)

  ################################################
  ## get values for duplicate snps.rec columns: ##
  ################################################

  ## get reconstruction for all original sites
  if(all.unique == TRUE){
    score3.complete <- score3
  }else{
    score3.complete <- score3[index]
    names(score3.complete) <- colnames(snps.rec.ori)
  }

  score3 <- score3.complete

  return(score3)

} # end subsequent.test








################
## get.score3 ##
################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param Pa A numeric value containing either the state,
#' or the probability of the state, of the phenotype at a given \emph{ancestral} node.
#' @param Pd A numeric value containing either the state,
#' or the probability of the state, of the phenotype at a given \emph{descendant} node.
#' @param Sa A numeric value containing either the state,
#' or the probability of the state, of SNPi at a given \emph{ancestral} node.
#' @param Sd A numeric value containing either the state,
#' or the probability of the state, of SNPi at a given \emph{descendant} node.
#' @param l A numeric value specifying the length of the branch in the phylogenetic tree
#' that joins the ancestral and descendant node.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#'
#' ## basic use of fn
#' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
#'

########################################################################

get.score3 <- function(Pa, Pd, Sa, Sd, l=NULL){

  score3 <- NULL

  if(!is.null(l)){
    ## NEW original integral-based score (WITHOUT edge-length!)...
    score3 <- l*(((4/3)*Pa*Sa) +
                   ((2/3)*Pa*Sd) +
                   ((2/3)*Pd*Sa) +
                   ((4/3)*Pd*Sd) -
                   Pa -
                   Pd -
                   Sa -
                   Sd +
                   1)
  }else{
    ## NEW original integral-based score (WITHOUT edge-length!)...
    score3 <- (((4/3)*Pa*Sa) +
                 ((2/3)*Pa*Sd) +
                 ((2/3)*Pd*Sa) +
                 ((4/3)*Pd*Sd) -
                 Pa -
                 Pd -
                 Sa -
                 Sd +
                 1)
  }

  return(score3)

} # end get.score3



#





#
