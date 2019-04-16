



#######################
## simultaneous.test ##
#######################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Test for association between genetic loci with Score 2.
#'
#' [*\emph{A work in progress; not curently integrated into treeWAS:}*]
#' Use the simultaneous.test (Score 2) to test for associations between genetic loci, 
#' which may indicate an epistatic interaction.
#' This function can be used either to test 
#' for pairwise association between all pairs of genetic loci
#' or for associations between a subset of snps and all other snps 
#' (recommended for large datasets; see details). 
#'
#' @param snps.reconstruction A matrix containing the terminal and reconstructed
#' ancestral states of SNPs for all nodes in the tree.
#' @param tree A phylo object containing the tree representing the ancestral relationships
#' between the individuals for which snps and phen are known.
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
# @useDynLib phangorn, .registration = TRUE
# @importFrom phangorn midpoint

simultaneous.test.epi <- function(snps.reconstruction, # can be snps.REC OR snps.sim.REC matrix ## NOTE: subs.edges no longer required for any version of this test.
                                  tree,
                                  snps.subset=NULL){

  snps.rec <- snps.reconstruction
  rm(snps.reconstruction)
  

  ## Always work with tree in pruningwise order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  # if(!is.rooted(tree)) tree <- midpoint(tree) # require(phangorn)
  ## Get tree edges:
  edges <- tree$edge


  ################################################
  ## GET SUBSET of SNPS (logical/names/indices) ##
  ################################################
  toKeep <- NULL
  if(!is.null(snps.subset)){
    if(!is.vector(snps.subset)){
      stop("snps.subset must be a vector (either a logical or numerical index vector, 
           or a vector of snps.rec column names, indicating which columns are to be kept as a subset")
    }else{
      ## LOGICAL (snps.subset = T/F toKeep) ##
      if(is.logical(snps.subset)){
        toKeep <- which(snps.subset == TRUE)
      }else{
        ## NUMERIC (snps.subset = indices toKeep) ##
        if(is.numeric(snps.subset)){
          if(!all(snps.subset %in% c(1:ncol(snps.rec)))){
            stop("not all snps.subset correspond to indices in 1:ncol(snps.rec)")
          }else{
            toKeep <- snps.subset
          }
        }else{
          ## CHARACTER (snps.subset = colnames toKeep) ##
        if(!all(snps.subset %in% colnames(snps.rec))){
          stop("not all snps.subset are in colnames(snps.rec)")
        }else{
          toKeep <- which(colnames(snps.rec) %in% snps.subset)
        } 
        }
      }
      # snps.rec <- snps.rec[,toKeep]
      # toKeep <- which(colnames(snps.rec) %in% snps.subset) # where (snps.subset = sig.snps.names)
    }
  }

  ####################################################################

  ###############################
  ## GET DIFFS ACROSS BRANCHES ##
  ###############################

  ## Get SNPs diffs: ##
  snps.diffs <- snps.rec[edges[,1], ] - snps.rec[edges[,2], ]

  ## Get snp1:snp2 diffs: ##
  s1s2.diffs <- SCORE2 <- list()

  ## If no snps.subset, run test over all columns...
  if(is.null(toKeep)) toKeep <- 1:ncol(snps.diffs)
  for(i in 1:length(toKeep)){
    s1s2.diffs[[i]] <- snps.diffs[,toKeep[i]] * snps.diffs
    ## Return with sign:
    SCORE2[[i]] <- colSums(s1s2.diffs[[i]], na.rm=TRUE)
    # SCORE2 <- abs(SCORE2)
    names(SCORE2[[i]]) <- paste(colnames(snps.rec)[toKeep[i]], colnames(snps.rec), sep="/")
  } # end for loop

  #######################
  score2 <- unlist(SCORE2)
  #######################
  noms <-  strsplit(names(score2), "/")
  str(noms)
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
  # noms.ori <- names(score2)
  attr(score2, "snps1") <- noms[,1]
  attr(score2, "snps2") <- noms[,2]
  # names(score2) ## still there, just not visible w str(score2)
  #######################

  return(score2)

} # end simultaneous.test.epi
####################################################################################
####################################################################################






#






