


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

  ###############################
  ## GET SCORE ACROSS BRANCHES ##
  ###############################

  Pd <- phen # .rec[edges[,2]]
  Sd <- snps # .rec[edges[,2], ]

  ################################################
  ## RE-SCALE NON-BINARY VALUES (phen only (?)) ##
  ################################################
  Pd.ori <- Pd
  Pd <- rescale(Pd, to=c(0,1))  ## require(plotrix) OR require(scales)

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
