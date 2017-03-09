


###################
## terminal.test ## ## ORIGINAL SCORE 1 ##
###################

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
#' out <- terminal.test(snps, phen)
#'
#' @importFrom scales rescale
#' @importFrom Hmisc all.is.numeric
#'

########################################################################


terminal.test <- function(snps,
                          phen){


  #########################
  ## Get UNIQUE snps.rec ##
  #########################
  temp <- get.unique.matrix(snps, MARGIN=2)
  snps.unique <- temp$unique.data
  index <- temp$index

  if(ncol(snps.unique) == ncol(snps)){
    all.unique <- TRUE
  }else{
    all.unique <- FALSE
  }

  ## work w only unique snps:
  snps.ori <- snps
  snps <- snps.unique

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
  if(n.levs == 2){
    if(!is.numeric(phen)){
      if(all.is.numeric(phen)){
        phen <- as.numeric(as.character(phen))
      }else{
        phen <- as.numeric(as.factor(phen))
      }
    }
  }else{
    if(!is.numeric(phen)){
      if(all.is.numeric(phen)){
        phen <- as.numeric(as.character(phen))
      }else{
        stop("phen has more than 2 levels but is not numeric (and therefore neither binary nor continuous).")
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

  ###############################
  ## GET SCORE ACROSS BRANCHES ##
  ###############################

  Pd <- phen # .rec[edges[,2]]
  Sd <- snps # .rec[edges[,2], ]

  ################################################
  ## RE-SCALE NON-BINARY VALUES (phen only (?)) ##
  ################################################
  Pd.ori <- Pd
  if(n.levs > 2) Pd <- rescale(Pd, to=c(0,1))  ## require(plotrix) OR require(scales)

  #################################################################     #####
  #############
  ## SCORE 1 ##
  #############
  ## ORIGINAL TERMINAL SCORE 1:

  score1 <- (Pd*Sd - (1 - Pd)*Sd - Pd*(1 - Sd) + (1 - Pd)*(1 - Sd))  ## CALCULATE SCORE 1 EQUATION

  score1 <- abs((colSums(score1, na.rm=TRUE)/length(Pd)))
  names(score1) <- colnames(snps)

  ############################################
  ## get values for duplicate snps columns: ##
  ############################################

  ## get reconstruction for all original sites
  if(all.unique == TRUE){
    score1.complete <- score1
  }else{
    score1.complete <- score1[index]
    names(score1.complete) <- colnames(snps.ori)
  }

  score1 <- score1.complete

  return(score1)

} # end terminal.test
