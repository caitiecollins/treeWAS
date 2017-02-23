
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

########################################################################

subsequent.test <- function(snps.reconstruction,
                            phen.reconstruction,
                            tree,
                            rec = "parsimony"){

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

  Pa <- phen.rec[edges[,1]]
  Pd <- phen.rec[edges[,2]]
  Sa <- snps.rec[edges[,1], ]
  Sd <- snps.rec[edges[,2], ]
  bl <- tree$edge.length

  #################################################################     #####

  ## NOTE: SCORE 3.1 will NOT work when reconstructions are done w ACE!
  if(rec == "ace"){
    score3.1 <- NULL
  }else{

    ###############
    ## SCORE 3.1 ##
    ###############
    ## SIMPLE UNWEIGHTED TALLY SCORE (1 point for each of "subsequent", "maintained", "simultaneous")

    score3.p <- c("00|00", "11|11", "00|11", "11|00", "01|00", "10|00", "01|11", "10|11")
    score3.n <- c("00|01", "00|10", "11|01", "11|10", "01|01", "10|10", "01|10", "10|01")

    score3.1 <- t(matrix(paste(Sa, Pa, "|", Sd, Pd, sep=""), nrow=ncol(snps.rec), byrow=T))
    score3.1 <- abs(sapply(c(1:ncol(score3.1)), function(e) length(which(score3.1[,e] %in% score3.p)) - length(which(score3.1[,e] %in% score3.n))))

    names(score3.1) <- colnames(snps.rec)
  }
  #################################################################     #####
  ###############
  ## SCORE 3.0 ##
  ###############
  ## ORIGINAL AND NEW INTEGRAL-BASED SCORE3 (with and without edge length):
  score3.L <- get.score3(Pa = Pa, Pd = Pd, Sa = Sa, Sd = Sd, l = bl)
  score3.NoL <- get.score3(Pa = Pa, Pd = Pd, Sa = Sa, Sd = Sd, l = NULL)

  score3.L <- abs(colSums(score3.L, na.rm=TRUE))
  names(score3.L) <- colnames(snps.rec)

  score3.NoL <- abs(colSums(score3.NoL, na.rm=TRUE))
  names(score3.NoL) <- colnames(snps.rec)

  ################################################
  ## get values for duplicate snps.rec columns: ##
  ################################################

  ## get reconstruction for all original sites
  if(all.unique == TRUE){
    score3.1.complete <- score3.1
    score3.L.complete <- score3.L
    score3.NoL.complete <- score3.NoL
  }else{
    ## only expand score3.1 if rec is by parsimony:
    if(!is.null(score3.1)){
      score3.1.complete <- score3.1[index]
      names(score3.1.complete) <- colnames(snps.rec.ori)
    }
    score3.L.complete <- score3.L[index]
    names(score3.L.complete) <- colnames(snps.rec.ori)
    score3.NoL.complete <- score3.NoL[index]
    names(score3.NoL.complete) <- colnames(snps.rec.ori)
  }

  score3.1 <- score3.1.complete
  score3.L <- score3.L.complete
  score3.NoL <- score3.NoL.complete

  ## only return score3.1 if present (ie. rec == parsimony)
  if(is.null(score3.1)){
    score <- list("score3.L" = score3.L,
                  "score3.NoL" = score3.NoL)
  }else{
    score <- list("score3.1" = score3.1,
                  "score3.L" = score3.L,
                  "score3.NoL" = score3.NoL)
  }
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
