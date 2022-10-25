


###################
## terminal.test ## ## SCORE 1 ##
###################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Terminal test
#'
#' Calculates treeWAS score 1, the terminal test.
#'
#' @param tree A phylo object.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#' ## Example ##
#' \dontrun{
#' ## basic use of fn
#' out <- terminal.test(snps, phen)
#' }
#'
#' @importFrom scales rescale
#' @importFrom Hmisc all.is.numeric
#'

########################################################################


terminal.test <- function(snps,
                          phen,
                          correct.prop = FALSE,
                          categorical = FALSE){

  ####################################################################
  #################
  ## Handle phen ##
  #################
  ## convert phenotype to numeric:
  phen.ori <- phen
  ## Convert to numeric (required for assoc tests):
  na.before <- length(which(is.na(phen)))

  ## NB: can only be binary or continuous at this point...
  levs <- unique(as.vector(unlist(phen)))
  n.levs <- length(levs[!is.na(levs)])
  if(!is.numeric(phen)){
    if(all.is.numeric(phen)){
      phen <- as.numeric(as.character(phen))
    }else{
      phen <- as.numeric(as.factor(phen))
      if(n.levs > 2){
        if(categorical != TRUE){
          warning("phen has more than 2 levels but is not numeric.
                  Setting 'categorical' to TRUE.")
          categorical <- TRUE
        }
    }
    }
  }
  ## ensure ind names not lost
  names(phen) <- names(phen.ori)

  ## Check that no errors occurred in conversion:
  na.after <- length(which(is.na(phen)))
  if(na.after > na.before){
    stop("NAs created while converting phen to numeric.")
  }
  ####################################################################

  ##################################
  ## GET SCORE 1 @ TERMINAL NODES ##
  ##################################

  Pd <- phen # .rec[edges[,2]]
  Sd <- snps # .rec[edges[,2], ]

  ################################################
  ## RE-SCALE NON-BINARY VALUES (phen only (?)) ##
  ################################################
  Pd.ori <- Pd
  if(categorical == FALSE){
    Pd <- rescale(Pd, to=c(0,1))  ## require(scales)
  }

  #################################################################     #####
  #############
  ## SCORE 1 ##
  #############
  if(categorical == FALSE){
    if(correct.prop == FALSE){
      ## ORIGINAL TERMINAL SCORE 1:
      score1 <- (Pd*Sd - (1 - Pd)*Sd - Pd*(1 - Sd) + (1 - Pd)*(1 - Sd))  ## CALCULATE SCORE 1 EQUATION

      ## Return with sign:
      score1 <- colSums(score1, na.rm=TRUE)/length(Pd)
    }else{
      ## MARGINAL-CORRECTED SCORE 1 (Phi):
      score1 <- ((colSums((1 - Pd)*(1 - Sd), na.rm=TRUE)*colSums(Pd*Sd, na.rm=TRUE)) -
                   (colSums((1 - Pd)*Sd, na.rm=TRUE)*colSums(Pd*(1 - Sd), na.rm=TRUE))) /
        (sqrt(colSums(1 - Sd, na.rm=TRUE)*colSums(Sd, na.rm=TRUE)*sum((1 - Pd), na.rm=TRUE)*sum(Pd, na.rm=TRUE)))
    }
  }else{
    ## CATEGORICAL SCORE 1 (Phi):
    score1 <- suppressWarnings(sqrt(sapply(c(1:ncol(Sd)), function(e)
      chisq.test(x=Pd, y=Sd[,e], correct=F)$statistic)/length(Pd)))
  }

  # score1 <- abs(score1)
  names(score1) <- colnames(snps)

  return(score1)

} # end terminal.test
