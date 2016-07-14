
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
#' @param test A character string or vector containing one or more of the following available tests of association:
#' "terminal", "simultaneous", "subsequent", "cor", "fisher". By default, the terminal test is run
#' (note that within treeWAS, the first three tests are run in a loop by default).
#' See details for more information on what these tests do and when they may be appropriate.
#' @param n.tests An integer between 1 and 5 specifying the number of tests you are running on all loci,
#' to be used in appropriately correcting for multiple testing.
#' (i.e., the number of times you will be running the \code{get.sig.snps} function).
#' @param p.value A single number specifying the p.value below which correlations are deemed to be 'significant'.
#' @param p.value.correct Specify if/how to correct for multiple testing:
#' either FALSE, or one of 'bonf' or 'fdr' (indicating, respectively,
#' the Bonferroni and False Discovery Rate corrections). By default, 'bonf' is selected
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


##########
## FDR: ##
##########

## NOTES: ##
## FDR works by managing the number of FALSE discoveries, RELATIVE to the number of TOTAL discoveries.
## Maintains Q = FALSE discoveries / TOTAL discoveries.
## Eg. For Q = 0.05, both 5/100 and 50/1000 meet the criterion.
## Hence, FDR is described as being both adaptive and scalable.

## QUESTION! ##
## HOW SHOULD WE HANDLE MULTIPLE TESTING CORRECTION WITH FDR WHEN RUNNING MULITPLE TESTS OF ASSOC??????
## Eg. Bonf --> multiply divisor by 3
## BUT--if we just run FDR p-value correction by multiplying the "n.tests" by n.tests,
## would the result not be that we just increase the number of false positives accepted in each test?!
## Could we pool the test results somehow??
## Or is it OK to ignore the performance of multiple separate assoc tests when using FDR??

########################################################################

get.sig.snps <- function(snps,
                         snps.sim,
                         phen,
                         tree,
                         test = "terminal",
                         n.tests = 1,
                         p.value = 0.001,
                         p.value.correct = "bonf",
                         p.value.by = "count",
                         snps.reconstruction = snps.REC,
                         snps.sim.reconstruction = snps.sim.REC,
                         phen.reconstruction = phen.REC){

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
      snps.sim <- do.call("cbind", snps.sim)
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
  #################################
  ## NOTE--change above to below if moving UNIQUE check to OUTSIDE get.sig.snps fn (& add index arguments):
  #   ## record whether all snps.sim are unique or not for later:
  #   if(is.null(snps.sim.index)){
  #     all.unique.sim <- TRUE
  #   }else{
  #     all.unique.sim <- FALSE
  #   }


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
  ## ensure ind names not lost
  names(phen) <- names(phen.ori)


  ######################
  ## ASSOCIATION TEST ##
  ######################

  if(test != "simultaneous" & test != "subsequent"){

    ########################################################
    ## Calculate correlations btw REAL SNPs and phenotype ##
    ########################################################
    corr.dat <- assoc.test(snps=snps, phen=phen, tree=NULL, test=test)


    #############################################################
    ## Calculate correlations btw SIMULATED SNPs and phenotype ##
    #############################################################
    corr.sim <- assoc.test(snps=snps.sim, phen=phen, tree=NULL, test=test)

  }else{

    #####################################
    ## SIMULTANEOUS & SUBSEQUENT TESTS ##
    #####################################
    ## (run w/ RECONSTRUCTIONS) ##

    ########################################################
    ## Calculate correlations btw REAL SNPs and phenotype ##
    ########################################################
    corr.dat <- assoc.test(snps=snps.reconstruction, phen=phen.reconstruction, tree=tree, test=test)


    #############################################################
    ## Calculate correlations btw SIMULATED SNPs and phenotype ##
    #############################################################
    corr.sim <- assoc.test(snps=snps.sim.reconstruction, phen=phen.reconstruction, tree=tree, test=test)
  }

  ###################################
  ## HANDLE DUPLICATE SNPS COLUMNS ##
  ###################################

  ## Expand corr.dat (if not all snps columns unique):
  if(all.unique == FALSE){
    corr.dat.complete <- corr.dat[snps.index]
    names(corr.dat.complete) <- colnames(snps.ori)
    corr.dat <- corr.dat.complete
  }

  ## Expand corr.sim (if not all snps.sim columns unique):
  if(all.unique.sim == FALSE){
    corr.sim.complete <- corr.sim[snps.sim.index]
    names(corr.sim.complete) <- colnames(snps.sim.ori)
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

    p.value <- p.value/(length(corr.dat)*n.tests)

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


    # system.time( # 0.15 :)!
    p.vals <- .get.p.vals(corr.sim)
      # )

    #     # system.time( # 144 (for 100,000 sites)
    #     p.vals <- sapply(c(1:length(corr.sim)),
    #                      function(e)
    #                        length(which(corr.sim > corr.sim[e]))
    #                      /length(corr.sim))
    #     # )
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
                 1/(length(corr.sim)*n.tests), sep=" ") ## CHECK---IS THIS RIGHT? SHOULD WE BE MULTIPLYING THE DIVISOR BY N.TESTS ??????????

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
#' @param test A character string or vector containing one or more of the following available tests of association:
#' "terminal", "simultaneous", "subsequent", "cor", "fisher". By default, the first three tests are run.
#' See details for more information on what these tests do and when they may be appropriate.
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @export

########################################################################

assoc.test <- function(snps,
                       phen,
                       tree = NULL,
                       test = c("terminal",
                                "simultaneous",
                                "subsequent",
                                "cor",
                                "fisher")){

  ##########################################
  ## TERMINAL (test 1: correlation score) ##
  ##########################################
  if(test=="terminal"){
    # ~ Correlation "SCORE" =
    # ((nS1P1 + nS0P0) - (nS1P0 + nS0P1) / (n.total))
    ## must be calculated for each SNP individually...
    corr.dat <- sapply(c(1:ncol(snps)), function(e)
      (((length(which(snps[which(phen==1),e]==1)) +
           length(which(snps[which(phen==0),e]==0)))
        - (length(which(snps[which(phen==1),e]==0)) +
             length(which(snps[which(phen==0),e]==1))))
       / nrow(snps)))
  } # end test terminal

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


  #######################
  ## SIMULTANEOUS TEST ##
  #######################
  if(test == "simultaneous"){
    corr.dat <- simultaneous.test(snps.reconstruction = snps, phen.reconstruction = phen, tree = tree)
  } # end test simultaneous

  #####################
  ## SUBSEQUENT TEST ##
  #####################
  if(test == "subsequent"){
    corr.dat <- subsequent.test(snps.reconstruction = snps, phen.reconstruction = phen, tree = tree)
  } # end test subsequent



  ## USE ABSOLUTE VALUE
  corr.dat <- abs(corr.dat)

  return(corr.dat)
} # end assoc.test





#################
## .get.p.vals ##
#################
## NOTE: only used for FDR threshold calculation!
## get a p-value associated with every value of corr.sim

.get.p.vals <- function(corr.sim){
  ## faster with table:
  cs.tab <- table(corr.sim)
  cs.fac <- factor(corr.sim)

  p.vals.unique <- sapply(c(1:(length(cs.tab)-1)),
                          function(e)
                            sum(cs.tab[(e+1):length(cs.tab)])
                          /sum(cs.tab))
  ## need to add trailing 0 for max corr.sim:
  p.vals.unique <- c(p.vals.unique, 0)

  ## get the unique index that
  ## each original corr.sim should map to:
  map.to <- (as.integer(cs.fac) - 1)
  map.to <- map.to[!is.na(map.to)]
  if(length(map.to)) map.to <- map.to + 1

  ## get p.vals for all corr.sim:
  p.vals <- p.vals.unique[map.to]

  return(p.vals)
} # end .get.p.vals
