


###################
## terminal.test ## 
###################


########################################################################

###################
## DOCUMENTATION ##
###################

#' Test for epistasis between genetic loci with Score 1.
#'
#' [*\emph{A work in progress; not curently integrated into treeWAS:}*]
#' Use the terminal.test (Score 1) to test for associations between genetic loci, 
#' which may indicate an epistatic interaction.
#' This function can be used either to test 
#' for pairwise association between all pairs of genetic loci
#' or for associations between a subset of snps and all other snps 
#' (recommended for large datasets; see details). 
#'
#' @param snps A matrix containing the states of SNPs (in columns) for all individuals (in rows).
#' @param snps.subset An optional vector (see details); else, NULL. 
#' The snps.subset vector can be a character vector, containing a subset of colnames(snps.rec), 
#' a logical vector, using TRUE or FALSE to indicate which columns are to be retained and excluded,
#' or an integer vector, specifying the column indices to be retained. 
#' 
#' 
#' @details The number of pairwise tests between all pairs of snps 
#' grows rapidly as the number of snps columns increases. 
#' As such, for datasets where ncol(snps.reconstruction) is large, we recommend that
#' the snps.subset argument is used to reduce the number of tests, by
#' indicating which snps to test for association with all other snps. 
#' The snps.subset index can be used to select any subset of snps of interest. 
#' For example, one may wish to test for interactions between all snps and a subset of snps that 
#' had been deemed significantly associated with a particular phenotype in a previous run of treeWAS.
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


terminal.test.epi <- function(snps,
                              snps.subset=NULL){

  ################################################
  ## GET SUBSET of SNPS (logical/names/indices) ##
  ################################################
  toKeep <- NULL
  if(!is.null(snps.subset)){
    # if(length(snps.subset) != ncol(snps)){
    #   stop("snps.subset must be of length ncol(snps).")
    # }else{
    if(is.logical(snps.subset)){
      toKeep <- snps.subset
    }else{
      if(all(snps.subset %in% colnames(snps))) toKeep <- which(colnames(snps) %in% snps.subset)
    }
    # snps.rec <- snps[,toKeep]
    toKeep <- which(colnames(snps) %in% snps.subset) # where (snps.subset = sig.snps.names)
    # }
  }

  ################################################
  ## RE-SCALE NON-BINARY VALUES (phen only (?)) ##
  ################################################
  # Pd <- phen
  # Pd.ori <- Pd
  # # if(n.levs > 2)
  # Pd <- rescale(Pd, to=c(0,1))  ## require(scales)

  #################################################################     #####
  #############
  ## SCORE 1 ##
  #############
  #######################
  ## ORIGINAL TERMINAL SCORE 1:
  # score1 <- (Pd*Sd - (1 - Pd)*Sd - Pd*(1 - Sd) + (1 - Pd)*(1 - Sd))  ## CALCULATE SCORE 1 EQUATION
  #######################

  ## Get snp1:snp2 diffs: ##
  s1s2 <- SCORE1 <- list()
  Sd <- snps
  ## If no snps.subset, run test over all columns...
  if(is.null(toKeep)) toKeep <- 1:ncol(snps)
  for(i in 1:length(toKeep)){
    Pd <- snps[,toKeep[i]]
    s1s2[[i]] <- (Pd*Sd - (1 - Pd)*Sd - Pd*(1 - Sd) + (1 - Pd)*(1 - Sd))  ## CALCULATE SCORE 1 EQUATION

    ## Return with sign:
    SCORE1[[i]] <- colSums(s1s2[[i]], na.rm=TRUE)/length(Pd)
    # SCORE1 <- abs(SCORE1)
    names(SCORE1[[i]]) <- paste(colnames(snps)[toKeep[i]], colnames(snps), sep="/")
  } # end for loop

  score1 <- unlist(SCORE1)
  #######################

  # names(score1) <- colnames(snps)
  #######################
  noms <-  strsplit(names(score1), "/")
  # str(noms)
  mat <- rep(NA, length(noms))
  mat <- cbind(mat, mat)
  for(i in 1:length(noms)){
    if(length(noms[[i]]) > 2){
      x <- noms[[i]][1]
      x[2] <- paste(noms[[i]][2:length(noms[[i]])], collapse="/")
      noms[[i]] <- x
    }
    mat[i,] <- noms[[i]]
  }  # end for loop
  noms <- do.call(rbind, noms)
  # str(noms)
  #######################
  # noms.ori <- names(score1)
  attr(score1, "snps1") <- noms[,1]
  attr(score1, "snps2") <- noms[,2]
  # names(score1) ## still there, just not visible w str(score1)
  #######################


  return(score1)

} # end terminal.test.epi
