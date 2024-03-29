
#####################
## subsequent.test ## ## SCORE 3 ##
#####################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Subsequent test
#'
#' Calculates treeWAS score 3, the subsequent test.
#'
#' @param tree A phylo object.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @examples
#'
#' ## basic use of fn
#' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
#'
#' @importFrom scales rescale
#' @importFrom Hmisc all.is.numeric
#' @export

########################################################################
# @useDynLib phangorn, .registration = TRUE
# @importFrom phangorn midpoint

subsequent.test <- function(snps.reconstruction,
                            phen.reconstruction,
                            tree,
                            correct.prop = FALSE,
                            categorical = FALSE){

  snps.rec <- snps.reconstruction
  phen.rec <- phen.reconstruction
  rm(snps.reconstruction)
  rm(phen.reconstruction)

  ## Always work with tree in pruningwise order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  # if(!is.rooted(tree)) tree <- midpoint(tree) # require(phangorn)
  ## get tree edges:
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
  if(!is.numeric(phen.rec)){
    if(all.is.numeric(phen.rec)){
      phen.rec <- as.numeric(as.character(phen.rec))
    }else{
      phen.rec <- as.numeric(as.factor(phen.rec))
      if(n.levs > 2){
        if(categorical != TRUE){
          warning("phen.rec has more than 2 levels but is not numeric.
                  Setting 'categorical' to TRUE.")
          categorical <- TRUE
        }
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
  if(categorical == FALSE){
    phen.rec <- rescale(phen.rec, to=c(0,1)) # require(scales)
  }

  ###############################
  ## GET SCORE ACROSS BRANCHES ##
  ###############################

  ## Get snps, phen values for all internal+terminal nodes:
  Sx <- snps.rec
  Px <- phen.rec

  #################################################################     #####
  ###############
  ## SCORE 3.0 ##
  ###############
  if(categorical == FALSE){
    if(correct.prop == FALSE){
      ## Get snps, phen values for ancestral & descendant nodes:
      Pa <- phen.rec[edges[,1]]
      Pd <- phen.rec[edges[,2]]
      Sa <- snps.rec[edges[,1], ]
      Sd <- snps.rec[edges[,2], ]
      bl <- tree$edge.length

      ## ORIGINAL INTEGRAL-BASED SCORE3 (without edge length):
      score3 <- get.score3(Pa = Pa, Pd = Pd, Sa = Sa, Sd = Sd, l = NULL)

      ## Return with sign:
      score3 <- colSums(score3, na.rm=TRUE)
    }else{
      ## MARGINAL-CORRECTED SCORE 1 (Phi):
      score3 <- ((colSums((1 - Px)*(1 - Sx), na.rm=TRUE)*colSums(Px*Sx, na.rm=TRUE)) -
                   (colSums((1 - Px)*Sx, na.rm=TRUE)*colSums(Px*(1 - Sx), na.rm=TRUE))) /
        (sqrt(colSums(1 - Sx, na.rm=TRUE)*colSums(Sx, na.rm=TRUE)*sum((1 - Px), na.rm=TRUE)*sum(Px, na.rm=TRUE)))
    }
  }else{
    ## CATEGORICAL SCORE 3 (Phi):
    score3 <- suppressWarnings(sqrt(sapply(c(1:ncol(Sx)), function(e)
      chisq.test(x=Px, y=Sx[,e], correct=F)$statistic)/length(Px)))
  }

  # score3 <- abs(score3)
  names(score3) <- colnames(snps.rec)

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
#' ## Example ##
#' \dontrun{
#' ## basic use of fn
#' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
#' }

########################################################################

get.score3 <- function(Pa, Pd, Sa, Sd, l=NULL){

  score3 <- NULL

  if(!is.null(l)){
    ## NEW integral-based score (WITH edge-length!)...
    score3 <- (l*(((4/3)*Pa*Sa) +
                   ((2/3)*Pa*Sd) +
                   ((2/3)*Pd*Sa) +
                   ((4/3)*Pd*Sd) -
                   Pa -
                   Pd -
                   Sa -
                   Sd +
                   1))/sum(l)
  }else{
    ## NEW integral-based score (WITHOUT edge-length!)...
    score3 <- (((4/3)*Pa*Sa) +
                 ((2/3)*Pa*Sd) +
                 ((2/3)*Pd*Sa) +
                 ((4/3)*Pd*Sd) -
                 Pa -
                 Pd -
                 Sa -
                 Sd +
                 1)/length(Pa)
  }

  return(score3)

} # end get.score3



###########################################################################################

###########################################################################################














#
