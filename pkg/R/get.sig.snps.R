
######################
## get.correlations ##
######################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Get significant SNPs, according to a given test of association.
#'
#' Identify which SNPs are deemed to be significantly associated with a phenotype,
#' according to a given test of association and p-value.
#' (Serves as the treeWAS association testing function;
#' runs the \code{assoc.test} function internally.)
#'
#' @param snps A matrix containing the real snps.
#' @param snps.sim A matrix or list of matrices containing simulated snps.
#' @param phen A factor or vector containing the phenotype (only allowed to contain two levels for now).
#' @param tree A phylo object containing a phylogenetic tree in which the number of tips is equal to the
#' length of \code{phen} and the number of rows of \code{snps} and \code{snps.sim}.
#' @param test A character string specifying the test to be used to measure correlations between the snps
#' (real and simulated) and the phenotype. Must be one of 'score', 'cor', 'fisher' (indicating, respectively,
#' the correlation score, standard correlation, fisher's exact test).
#' @param p.value A single number specifying the p.value below which correlations are deemed to be 'significant'.
#' @param p.value.correct Specify if/how to correct for multiple testing:
#' either FALSE, or one of 'bonf' or 'fdr' (indicating, respectively,
#' the Bonferroni and False Discovery Rate corrections).
#' @param p.value.by Specify how to determine the location of the p.value threshold:
#' either 'count' or 'density' (indicating, respectively, that the p.value threshold should
#' be determined by exact count or with the use of a density function).
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @examples
#'
#' ## load data
#' data(dist)
#' str(dist)
#'
#' ## basic use of fn
#'
#' fn(arg1, arg2)
#'
#' #' ## more elaborate use of fn
#' fn(arg1, arg2)
#'
#' @export

########################################################################

get.sig.snps <- function(snps, snps.sim,
                         phen, tree,
                         test=c("score", "cor", "fisher", "ace.cum", "ace.pagel"),
                         p.value=0.001,
                         p.value.correct=c("bonf", "fdr", FALSE),
                         p.value.by=c("count", "density")){

  #################
  ## HANDLE SNPS ##
  #################

  ##########
  ## snps ##
  ##########

  ## Check snps column names
  if(is.null(colnames(snps))) colnames(snps) <- c(1:ncol(snps))

  ## Get UNIQUE snps + index
  snps.ori <- snps
  temp <- get.unique.matrix(snps, MARGIN=2)
  snps.unique <- temp$unique.data
  snps.index <- temp$index
  snps <- snps.unique

  ## record whether all snps are unique or not for later:
  if(ncol(snps.unique) == ncol(snps.ori)){
    all.unique <- TRUE
  }else{
    all.unique <- FALSE
  }

  ##############
  ## snps.sim ##
  ##############

  ## Handle matrix/list input:
  if(class(snps.sim) == "list"){
    ## If list of length 1...
    if(length(snps.sim) == 1){
      ## keep matrix:
      snps.sim <- snps.sim[[1]]
    }else{
      ## If list of multiple matrices...
      ## merge all elements into one big matrix
      ## by pasting columns together:
      snps.sim <- do.call(cbind, snps.sim)
    }
  }

  ## Get UNIQUE snps.sim + index
  snps.sim.ori <- snps.sim
  temp <- get.unique.matrix(snps.sim, MARGIN=2)
  snps.sim.unique <- temp$unique.data
  snps.sim.index <- temp$index
  snps.sim <- snps.sim.unique

  ## record whether all snps are unique or not for later:
  if(ncol(snps.sim.unique) == ncol(snps.sim.ori)){
    all.unique.sim <- TRUE
  }else{
    all.unique.sim <- FALSE
  }


  #################
  ## HANDLE PHEN ##
  #################
  ## convert phenotype to numeric:
  phen.ori <- phen
  if(!is.numeric(phen)) phen <- as.numeric(phen)
  ## for ease of interpretation,
  ## if phen has 2 levels, 1 and 2,
  ## make these 0 and 1:
  if(length(unique(phen))!=2){
    stop("This function is only designed for phenotypes with two levels.")
  }else{
    if(length(phen[-c(which(phen==1), which(phen==2))])==0){
      phen <- replace(phen, which(phen==1), 0)
      phen <- replace(phen, which(phen==2), 1)
    }
  }

  #####################
  ## ACE-BASED TESTS ##
  #####################

  ## RUN ACE on all UNIQUE snps & snps.sim columns & on phen:
  if(test %in% c("ace.cum", "ace.pagel")){

    ## Get ACE for phen:
    ace.phen <- ace(x=phen, phy=tree, type="discrete",
                    method="ML", model="ER")

    ## Get ACE for SNPs (real):
    ace.snps <- list()
    for(i in 1:ncol(snps)){
      ace.snps[[i]] <- ace(x=snps[, i], phy=tree, type="discrete",
                           method="ML", model="ER")
    }

    ## Get ACE for SNPs (simulated):
    ace.snps.sim <- list()
    for(i in 1:ncol(snps.sim)){
      ace.snps.sim[[i]] <- ace(x=snps.sim[, i], phy=tree, type="discrete",
                              method="ML", model="ER")
    }

  } # end ACE for ace-based tests

  ######################
  ## ASSOCIATION TEST ##
  ######################

  ########################################################
  ## Calculate correlations btw REAL SNPs and phenotype ##
  ########################################################
  corr.dat <- assoc.test(snps=snps, phen=phen, test=test)


  #############################################################
  ## Calculate correlations btw SIMULATED SNPs and phenotype ##
  #############################################################
  corr.sim <- assoc.test(snps=snps.sim, phen=phen, test=test)


  ###################################
  ## HANDLE DUPLICATE SNPS COLUMNS ##
  ###################################

  ## Expand corr.dat (if not all snps columns unique):
  if(all.unique == FALSE){
    corr.dat.complete <- rep(NA, ncol(snps.ori))
    for(i in 1:ncol(snps.unique)){
      corr.dat.complete[which(snps.index == i)] <- corr.dat.unique[i]
    }
    corr.dat <- corr.dat.complete
  }

  ## Expand corr.sim (if not all snps.sim columns unique):
  if(all.unique.sim == FALSE){
    corr.sim.complete <- rep(NA, ncol(snps.sim.ori))
    for(i in 1:ncol(snps.sim.unique)){
      corr.sim.complete[which(snps.sim.index == i)] <- corr.sim.unique[i]
    }
    corr.sim <- corr.sim.complete
  }



  ############################
  ## HANDLE P.VALUE OPTIONS ##
  ############################

  #####################
  ## p.value.correct ##
  #####################
  if(p.value.correct == "bonf"){
    ##########
    ## bonf ##
    ##########
    p.value <- p.value/length(corr.dat)

    ################
    ## p.value.by ##
    ################
    if(p.value.by == "count") thresh <- quantile(corr.sim, probs=1-p.value)
    if(p.value.by == "density") thresh <- quantile(density(corr.sim)$x, probs=1-p.value)
  }


  if(p.value.correct == "fdr"){
    #########
    ## fdr ##
    #########

    p.vals <- sapply(c(1:length(corr.sim)),
                     function(e)
                       length(which(corr.sim > corr.sim[e]))
                     /length(corr.sim))
    p.vals <- sort(p.vals, decreasing=TRUE)
    p.fdr <- p.adjust(p.vals, method="fdr", n=length(p.vals))
    p.thresh <- quantile(p.fdr, probs=1-p.value)
    thresh <- quantile(corr.sim, probs=p.thresh)

    ################
    ## p.value.by ##
    ################
    if(p.value.by == "count") thresh <- quantile(corr.sim, probs=p.thresh)
    if(p.value.by == "density") thresh <- quantile(density(corr.sim)$x, probs=p.thresh)

    #     colnames(snps)[which(pval.fdr < p.thresh)]

    #     ## COMPARING THRESHOLDS ATTAINED BY DIFFERENT METHODS ##
    #
    #     ## get threshold by counts
    #     thresh.uncorr <- quantile(corr.sim, probs=1-p.value)
    #     thresh.bonf <- quantile(corr.sim, probs=1-(p.value/length(corr.sim)))
    #     thresh.fdr <- quantile(corr.sim, probs=p.thresh)
    #
    #     ## plot
    #     hist(corr.sim, breaks=50, xlim=c(0,1))
    #     abline(v=thresh.uncorr, col="black", lwd=2)
    #     abline(v=thresh.bonf, col="red", lwd=2)
    #     abline(v=thresh.fdr, col="blue", lwd=2)
    #
    #     ## get threshold by density
    #     thresh.uncorr.d <- quantile(density(corr.sim)$x, probs=1-p.value)
    #     thresh.bonf.d <- quantile(density(corr.sim)$x, probs=1-(p.value/length(corr.sim)))
    #     thresh.fdr.d <- quantile(density(corr.sim)$x, probs=p.thresh)
    #
    #     ## plot
    #     par(new=TRUE)
    #     plot(density(corr.sim), xlim=c(0,1), col="blue")
    #     abline(v=thresh.uncorr.d, col="black", lwd=2, lty=2)
    #     abline(v=thresh.bonf.d, col="red", lwd=2, lty=2)
    #     abline(v=thresh.fdr.d, col="blue", lwd=2, lty=2)

  } # end p.value.correct == "fdr"


  ##################################
  ## GET SIGNIFICANT CORRELATIONS ##
  ##################################
  if(test=="fisher"){
    ## Identify (real) SNPs w correlations > thresh:
    sig.snps <- which(corr.dat < thresh)
    p.vals <- sapply(c(1:length(sig.snps)),
                     function(e)
                       length(which(corr.sim < corr.dat[sig.snps[e]]))
                     /length(corr.sim))
  }else{
    ## Identify (real) SNPs w correlations > thresh:
    sig.snps <- which(corr.dat > thresh)
    p.vals <- sapply(c(1:length(sig.snps)),
                     function(e)
                       length(which(corr.sim > corr.dat[sig.snps[e]]))
                     /length(corr.sim))
  }


  ## 0 p.vals
  min.p <- paste("p-values listed as 0 are <",
                 1/length(corr.sim), sep=" ")

  ## get list of those correlation values
  sig.corrs <- corr.dat[sig.snps]
  ## get the list of those SNPs (ie. their locus names)
  # sig.snps <- dimnames(snps)[[2]][sig.snps]
  sig.snps.names <- dimnames(snps)[[2]][sig.snps]

  ## re-order list of sig.snps and sig.corrs by value of sig.corr
  if(test=="fisher"){
    NWO <- order(sig.corrs, decreasing=TRUE)
  }else{
    NWO <- order(sig.corrs, decreasing=FALSE)
  }
  sig.snps <- sig.snps[NWO]
  sig.corrs <- sig.corrs[NWO]
  p.vals <- p.vals[NWO]

  gc()

  #################
  ## GET RESULTS ##
  #################

  out <- list(corr.dat,
              corr.sim,
              thresh,
              sig.snps.names,
              sig.snps,
              sig.corrs,
              p.vals,
              min.p)

  names(out) <- c("corr.dat",
                  "corr.sim",
                  "sig.thresh",
                  "sig.snps.names",
                  "sig.snps",
                  "sig.corrs",
                  "p.vals",
                  "min.p")

  return(out)

} # end get.sig.snps







################
## assoc.test ##
################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Run a test of association between SNPs and a phenotype.
#'
#' Run one of five tests of association between each column of a SNPs matrix and a phenotype
#' (some tests only implemented for \emph{binary} SNPs and phenotype).
#'
#' @param snps A matrix containing the real snps.
#' @param phen A factor or vector containing the phenotype (only allowed to contain two levels for now).
#' @param test A character string specifying the test to be used to measure correlations between the snps
#' (real and simulated) and the phenotype. Must be one of 'score', 'cor', 'fisher', 'ace.cum', 'ace.pagel'
#' (indicating, respectively, the correlation score, standard correlation, fisher's exact test,
#' ancestral character estimation (ACE) based cumulative multiplicative test, and ACE-based Pagel-inspired test).
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @export

########################################################################

assoc.test <- function(snps, phen, test=c("score",
                                          "cor",
                                          "fisher",
                                          "ace.cum",
                                          "ace.pagel")){

  ###########
  ## SCORE ##
  ###########
  if(test=="score"){
    # ~ Correlation "SCORE" =
    # ((nS1P1 + nS0P0) - (nS1P0 + nS0P1) / (n.total))
    ## must be calculated for each SNP individually...
    corr.dat <- sapply(c(1:ncol(snps)), function(e)
      (((length(which(snps[which(phen==1),e]==1)) +
           length(which(snps[which(phen==0),e]==0)))
        - (length(which(snps[which(phen==1),e]==0)) +
             length(which(snps[which(phen==0),e]==1))))
       / nrow(snps)))
  } # end test score

  #################
  ## CORRELATION ##
  #################
  if(test=="cor"){
    corr.dat <- sapply(c(1:ncol(snps)), function(e)
      cor(snps[,e], phen)) # regular correlation...
  } # end test cor

  #########################
  ## FISHER'S EXACT TEST ##
  #########################
  if(test=="fisher"){
    corr.dat <- sapply(c(1:ncol(snps)),
                       function(e) fisher.test(snps[,e],
                                               y=phen, alternative="two.sided")$p.value)
    ## two.sided bc we want to know if inds w
    ## the phen have EITHER more 1s or 0s
  } # end test fisher



  #####################
  ## ACE-BASED TESTS ##
  #####################

  ## NOTE: For ace-based tests (cum & pagel), instead of inputting snps & phen as variables,
  ## input ace (or just ace$lik.anc (?)) for snps & phen...


  ## USE ABSOLUTE VALUE
  corr.dat <- abs(corr.dat)

  return(corr.dat)
} # end assoc.test
