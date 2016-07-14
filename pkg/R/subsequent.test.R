
#####################
## subsequent.test ##
#####################

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


########################################################################

subsequent.test <- function(snps.reconstruction,
                            phen.reconstruction,
                            tree){

  snps.rec <- snps.reconstruction
  phen.rec <- phen.reconstruction

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

  ###############################
  ## GET SCORE ACROSS BRANCHES ##
  ###############################

  score <- list()

  for(i in 1:ncol(snps.rec)){

    score3 <- list()

    ##############
    ## FOR LOOP ##
    ##############
    for(e in 1:nrow(edges)){

      pa <- phen.rec[edges[e,1]]
      pd <- phen.rec[edges[e,2]]
      sa <- snps.rec[edges[e,1], i]
      sd <- snps.rec[edges[e,2], i]
      bl <- tree$edge.length[e]

      score3[[e]] <- get.score3(Pa = pa, Pd = pd, Sa = sa, Sd = sd, l = bl)

    } # end (e) for loop

    score[[i]] <- abs(sum(as.vector(unlist(score3))))

  } # end (i) for loop

  score <- as.vector(unlist(score))
  names(score) <- colnames(snps.rec)

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

get.score3 <- function(Pa, Pd, Sa, Sd, l){

  score3 <- NULL

  ## CHECKS:
  if(length(Pa) > 1) stop("Pa must be of length one
                                    (i.e., (the probability of)
                          the phenotypic state of ONE ancestral node.")
  if(length(Pd) > 1) stop("Pd must be of length one
                                    (i.e., (the probability of)
                          the phenotypic state of ONE descendant node.")
  if(length(Sa) > 1) stop("Sa must be of length one
                                    (i.e., (the probability of)
                          the SNPi state of ONE ancestral node.")
  if(length(Sd) > 1) stop("Sd must be of length one
                                    (i.e., (the probability of)
                          the SNPi state of ONE descendant node.")

  score3 <- l*(((4/3)*Pa*Sa) +
                 ((2/3)*Pa*Sd) +
                 ((2/3)*Pd*Sa) +
                 ((4/3)*Pd*Sd) -
                 Pa -
                 Pd -
                 Sa -
                 Sd +
                 1)

  return(score3)

} # end get.score3


