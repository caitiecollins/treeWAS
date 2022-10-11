



#######################
## simultaneous.test ## ## SCORE 2 ##
#######################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Simultaneous test
#'
#' Calculates treeWAS score 2, the simultaneous test, as the number of 
#' substitutions or changes in genotype (\code{snps.reconstruction}) and phenotype 
#' (\code{phen.reconstruction}) that occur simultaneously on the same branches of the tree. 
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
#'
#'
#' @importFrom scales rescale
#' @importFrom Hmisc all.is.numeric
#' @importFrom utils combn
#'
#' @export

########################################################################
# @useDynLib phangorn, .registration = TRUE
# @importFrom phangorn midpoint

simultaneous.test <- function(snps.reconstruction,
                              phen.reconstruction,
                              tree,
                              categorical = FALSE){

  snps.rec <- snps.reconstruction
  phen.rec <- phen.reconstruction
  rm(snps.reconstruction)
  rm(phen.reconstruction)

  ## Always work with tree in pruningwise order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  # if(!is.rooted(tree)) tree <- midpoint(tree) # require(phangorn)
  ## Get tree edges:
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
  ## GET DIFFS ACROSS BRANCHES ##
  ###############################

  if(categorical == FALSE){
    ## ORIGINAL SCORE 2:
    ## Get SNPs diffs: ##
    snps.diffs <- snps.rec[edges[,1], ] - snps.rec[edges[,2], ]
    
    ## Get phen diffs: ##
    phen.diffs <- phen.rec[edges[,1]] - phen.rec[edges[,2]]
    
    sp.diffs <- snps.diffs * phen.diffs
    
    ## Return with sign:
    score2 <- colSums(sp.diffs, na.rm=TRUE)
    # score2 <- abs(score2)
    names(score2) <- colnames(snps.rec)
    
  }else{
    ## CATEGORICAL SCORE 2:
    
    ## Get SNPs diffs: ##
    snps.diffs <- snps.rec[edges[,1], ] - snps.rec[edges[,2], ]
    
    pairs <- t(combn(unique(phen.rec[!is.na(phen.rec)]), m=2))
    S2 <- list()
    for(p in 1:nrow(pairs)){
      
      ## Get phen diffs: ##
      pr <- phen.rec
      pr[which(!pr %in% pairs[p,])] <- NA
      pr <- as.numeric(as.factor(as.character(pr)))-1
      phen.diffs <- pr[edges[,1]] - pr[edges[,2]]
      
      sp.diffs <- snps.diffs * phen.diffs
      S2[[p]] <- colSums(sp.diffs, na.rm=TRUE)
    } # end for (p) loop
    
    s2 <- do.call(rbind, S2)
    score2 <- colSums(abs(s2), na.rm=TRUE)
    names(score2) <- colnames(snps.rec)
  }
  
  return(score2)

} # end simultaneous.test




