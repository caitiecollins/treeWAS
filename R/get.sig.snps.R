
######################
## get.assoc.scores ##
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
#'             length of \code{phen} and the number of rows of \code{snps} and \code{snps.sim}.
#' @param test A character string or vector containing one or more of the following available tests of association:
#'             "terminal", "simultaneous", "subsequent", "cor", "fisher". By default, the terminal test is run
#'             (note that within treeWAS, the first three tests are run in a loop by default).
#'             See details for more information on what these tests do and when they may be appropriate.
#' @param correct.prop A logical indicating whether the \code{"terminal"} and \code{"subsequent"} tests will be corrected for
#'                     phenotypic class imbalance. Recommended if the proportion of individuals varies significantly across
#'                     the levels of the phenotype (if binary) or if the phenotype is skewed (if continuous).
#'                     If \code{correct.prop} is \code{FALSE} (the default),
#'                     the original versions of each test will be run as described in our
#'                     \href{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005958}{PLOS Computational Biology paper}.
#'                     If \code{TRUE}, an alternate association metric (based on the phi correlation coefficient) is calculated
#'                     across the terminal and all (internal and terminal) nodes, respectively.
#' @param categorical A logical indicating whether \code{phen} should be treated as a nominal categorical variable
#'                    whose unique values should be treated as levels rather than as meaningful numbers.
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @export
#' @importFrom pryr mem_used
#'

########################################################################

get.assoc.scores <- function(snps,
                             snps.sim,
                             phen,
                             tree,
                             test = "terminal",
                             correct.prop = FALSE,
                             categorical = FALSE,
                             snps.reconstruction = NULL,
                             snps.sim.reconstruction = NULL,
                             phen.reconstruction = NULL,
                             unique.cols = FALSE){

  print(paste("Started running", test, "test @", Sys.time()))

  ########         ########         ########         ########         ########         ########         ########         ########

  #######################################
  ## HANDLE INPUT (UNIQUE or EXPANDED) ##
  #######################################

  ####################
  ## snps, snps.rec ##
  ####################
  if(unique.cols == TRUE){
    all.unique <- TRUE
  }else{
    #################################
    ## get unique column patterns: ##
    #################################
    if(!test %in% c("simultaneous", "subsequent")){
      ##########
      ## snps ##
      ##########
      temp <- get.unique.matrix(snps, MARGIN=2)
      snps.unique <- temp$unique.data
      index <- temp$index

      if(ncol(snps.unique) == ncol(snps)){
        all.unique <- TRUE
      }else{
        all.unique <- FALSE
        colnoms <- colnames(snps)
      }

      ## work w only unique snps:
      snps.ori <- snps
      snps <- snps.unique
    }else{
      #########################
      ## snps.reconstruction ##
      #########################
      temp <- get.unique.matrix(snps.reconstruction, MARGIN=2)
      snps.reconstruction.unique <- temp$unique.data
      index <- temp$index

      if(ncol(snps.reconstruction.unique) == ncol(snps.reconstruction)){
        all.unique <- TRUE
      }else{
        all.unique <- FALSE
        colnoms <- colnames(snps.reconstruction)
      }

      ## work w only unique snps:
      # snps.ori <- snps
      # snps <- snps.unique
      snps.reconstruction.ori <- snps.reconstruction
      snps.reconstruction <- snps.reconstruction.unique
    }
  }


  ############################
  ## snps.sim, snps.sim.rec ##
  ############################
  if(unique.cols == TRUE){
    all.unique.sim <- TRUE
  }else{
    #################################
    ## get unique column patterns: ##
    #################################
    if(!test %in% c("simultaneous", "subsequent")){
      ##############
      ## snps.sim ##
      ##############
      temp <- get.unique.matrix(snps.sim, MARGIN=2)
      snps.sim.unique <- temp$unique.data
      index.sim <- temp$index

      if(ncol(snps.sim.unique) == ncol(snps.sim)){
        all.unique.sim <- TRUE
      }else{
        all.unique.sim <- FALSE
        colnoms.sim <- colnames(snps.sim)
      }

      ## work w only unique snps:
      snps.sim.ori <- snps.sim
      snps.sim <- snps.sim.unique
    }else{
      #############################
      ## snps.sim.reconstruction ##
      #############################
      temp <- get.unique.matrix(snps.sim.reconstruction, MARGIN=2)
      snps.sim.reconstruction.unique <- temp$unique.data
      index.sim <- temp$index

      if(ncol(snps.sim.reconstruction.unique) == ncol(snps.sim.reconstruction)){
        all.unique.sim <- TRUE
      }else{
        all.unique.sim <- FALSE
        colnoms.sim <- colnames(snps.sim.reconstruction)
      }

      ## work w only unique snps:
      snps.sim.reconstruction.ori <- snps.sim.reconstruction
      snps.sim.reconstruction <- snps.sim.reconstruction.unique
    }
  }

  ########         ########         ########         ########         ########         ########         ########         ########

  #################
  ## HANDLE PHEN ##
  #################
  ## convert phenotype to numeric:
  phen.ori <- phen

  ## store ind names:
  noms <- names(phen)

  ## Check if phen is binary:
  levs <- unique(as.vector(unlist(phen)))
  levs <- levs[!is.na(levs)]
  n.levs <- length(levs)
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

  ## ensure ind names not lost:
  names(phen) <- noms

  ########         ########         ########         ########         ########         ########         ########         ########
  ###############################################################################################################################

  ######################
  ## ASSOCIATION TEST ##
  ######################

  if(test != "simultaneous" & test != "subsequent"){

    ########################################################
    ## Calculate correlations btw REAL SNPs and phenotype ##
    ########################################################
    corr.dat <- assoc.test(snps=snps, phen=phen, tree=NULL, test=test,
                           correct.prop=correct.prop, categorical=categorical)

    print(paste("Real data scores completed for", test, "test @", Sys.time()))


    #############################################################
    ## Calculate correlations btw SIMULATED SNPs and phenotype ##
    #############################################################
    corr.sim <- assoc.test(snps=snps.sim, phen=phen, tree=NULL, test=test,
                           correct.prop=correct.prop, categorical=categorical)

    print(paste("Simulated data scores completed for", test, "test @", Sys.time()))

  }else{

    #####################################
    ## SIMULTANEOUS & SUBSEQUENT TESTS ##
    #####################################
    ## (run w/ RECONSTRUCTIONS) ##

    ########################################################
    ## Calculate correlations btw REAL SNPs and phenotype ##
    ########################################################
    corr.dat <- assoc.test(snps=snps.reconstruction, phen=phen.reconstruction, tree=tree, test=test,
                           correct.prop=correct.prop, categorical=categorical)

    print(paste("Real data scores completed for", test, "test @", Sys.time()))


    #############################################################
    ## Calculate correlations btw SIMULATED SNPs and phenotype ##
    #############################################################
    corr.sim <- assoc.test(snps=snps.sim.reconstruction, phen=phen.reconstruction, tree=tree, test=test,
                           correct.prop=correct.prop, categorical=categorical)

    print(paste("Simulated data scores completed for", test, "test @", Sys.time()))
  }

  ###################################
  ## HANDLE DUPLICATE SNPS COLUMNS ##
  ###################################

  ## Expand corr.dat (if not all snps columns unique):
  if(all.unique == FALSE){
    corr.dat.complete <- corr.dat[index]
    names(corr.dat.complete) <- colnoms
    corr.dat <- corr.dat.complete
  }

  ## Expand corr.sim (if not all snps.sim columns unique):
  if(all.unique.sim == FALSE){
    corr.sim.complete <- corr.sim[index.sim]
    names(corr.sim.complete) <- colnoms.sim
    corr.sim <- corr.sim.complete
  }

  ## quick look at corr.sim & corr.dat
  # hist(corr.sim, xlim=c(0,1))
  # hist(corr.dat, xlim=c(0,1))


  out <- list("corr.dat" = corr.dat,
              "corr.sim" = corr.sim)

  return(out)

} # end get.assoc.scores


########################################################################



##################
## get.sig.snps ##
##################

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
#' @param corr.dat A vector containing the association score values, for a given association test, for the real data.
#' @param corr.sim A vector containing the association score values, for a given association test, for the simulated data.
#' @param snps.names The column names of the original \code{snps} matrix from which the association score values
#'                    in \code{corr.dat} were derived.
#' @param test A character string or vector containing one or more of the following available tests of association:
#'             "terminal", "simultaneous", "subsequent", "cor", "fisher". By default, the terminal test is run
#'             (note that within treeWAS, the first three tests are run in a loop by default).
#'             See details for more information on what these tests do and when they may be appropriate.
#' @param p.value A single number specifying the p.value below which correlations are deemed to be 'significant'.
#' @param p.value.correct Specify if/how to correct for multiple testing:
#'                        either FALSE, or one of 'bonf' or 'fdr' (indicating, respectively,
#'                        the Bonferroni and False Discovery Rate corrections). By default, 'bonf' is selected
#' @param p.value.by Specify how to determine the location of the p.value threshold:
#'                   either 'count' or 'density' (indicating, respectively, that the p.value threshold should
#'                   be determined by exact count or with the use of a density function).
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @export
#' @importFrom pryr mem_used
#'



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


get.sig.snps <- function(corr.dat,
                         corr.sim,
                         snps.names,
                         test = "terminal",
                         p.value = 0.01,
                         p.value.correct = "bonf",
                         p.value.by = "count"){

  ## Use ABS Vals:
  # corr.dat <- abs(corr.dat)
  # corr.sim <- abs(corr.sim)

  n.tests <- length(test)

  ############################
  ## HANDLE P.VALUE OPTIONS ##
  ############################
  out <- list()

  #####################
  ## p.value.correct ##
  #####################
  if(p.value.correct == "bonf" || p.value.correct == F) {
    ##########
    ## bonf ##
    ##########
    if (p.value.correct == "bonf") p.value <- p.value/(length(corr.dat)*n.tests)

    ################
    ## p.value.by ##
    ################
    if(test == "fisher"){
      if(p.value.by == "count") thresh <- quantile(abs(corr.sim), probs=p.value)
      if(p.value.by == "density") thresh <- quantile(density(abs(corr.sim))$x, probs=p.value)
      if(thresh < 0) thresh <- 0 # problem only for density approach
    }else{
      if(p.value.by == "count") thresh <- quantile(abs(corr.sim), probs=1-p.value)
      if(p.value.by == "density") thresh <- quantile(density(abs(corr.sim))$x, probs=1-p.value)
    }
  }

  # print(paste("Test flag #4; memory used:", as.character(round(as.numeric(as.character(mem_used()/1000000000)), 2)), "Gb @",  Sys.time()))
  # print(paste("Thresh =",thresh,"@", Sys.time()))

  if(p.value.correct == "fdr"){
    #########
    ## fdr ##
    #########
    p.vals <- .get.p.vals(abs(corr.sim),
                          corr.dat = NULL,
                          fisher.test = TRUE)
    p.vals <- sort(p.vals, decreasing=TRUE)
    p.fdr <- p.adjust(p.vals, method="fdr", n=length(p.vals)*n.tests) # CHECK--IS THIS THE CORRECT MT APPROACH FOR FDR???
    p.thresh <- quantile(p.fdr, probs=1-p.value)

    ################
    ## p.value.by ##
    ################
    if(test == "fisher"){
      if(p.value.by == "count") thresh <- quantile(abs(corr.sim), probs=1-p.thresh)
      if(p.value.by == "density") thresh <- quantile(density(abs(corr.sim))$x, probs=1-p.thresh)
    }else{
      if(p.value.by == "count") thresh <- quantile(abs(corr.sim), probs=p.thresh)
      if(p.value.by == "density") thresh <- quantile(density(abs(corr.sim))$x, probs=p.thresh)
    }

  } # end p.value.correct == "fdr"

  # print(paste("Test flag #5; memory used:", as.character(round(as.numeric(as.character(mem_used()/1000000000)), 2)), "Gb @",  Sys.time()))

  ##################################
  ## GET SIGNIFICANT CORRELATIONS ##
  ##################################
  if(test=="fisher"){
    ## Identify (real) SNPs w correlations > thresh:
    sig.snps <- which(abs(corr.dat) < thresh)
    p.vals <- .get.p.vals(corr.sim = abs(corr.sim),
                          corr.dat = abs(corr.dat),
                          fisher.test = TRUE)
    sig.p.vals <- p.vals[sig.snps]
  }else{
    ## Identify (real) SNPs w correlations > thresh:
    sig.snps <- which(abs(corr.dat) > thresh)
    p.vals <- .get.p.vals(corr.sim = abs(corr.sim),
                          corr.dat = abs(corr.dat),
                          fisher.test = FALSE)
    sig.p.vals <- p.vals[sig.snps]
  }

  # print(paste("Test flag #6; memory used:", as.character(round(as.numeric(as.character(mem_used()/1000000000)), 2)), "Gb @",  Sys.time()))

  ## 0 p.vals
  min.p <- paste("p-values listed as 0 are <",
                 1/length(corr.sim), sep=" ") ## CHECK---IS THIS RIGHT? SHOULD WE BE MULTIPLYING THE DIVISOR BY N.TESTS ?

  ## get list of those correlation values
  sig.corrs <- corr.dat[sig.snps]
  ## get the list of those SNPs (ie. their locus names)
  sig.snps.names <- snps.names[sig.snps]

  ## re-order list of sig.snps and sig.corrs by value of sig.corr
  if(test=="fisher"){
    NWO <- order(abs(sig.corrs), decreasing=FALSE)
  }else{
    NWO <- order(abs(sig.corrs), decreasing=TRUE)
  }
  sig.snps <- sig.snps[NWO]
  sig.corrs <- sig.corrs[NWO]
  sig.p.vals <- sig.p.vals[NWO]
  sig.snps.names <- sig.snps.names[NWO]
  # gc()

  # print(paste("Test flag #7; memory used:", as.character(round(as.numeric(as.character(mem_used()/1000000000)), 2)), "Gb @",  Sys.time()))

  #################
  ## GET RESULTS ##
  #################

  out <- list(corr.dat,
               corr.sim,
               p.vals,
               thresh,
               sig.snps.names,
               sig.snps,
               sig.corrs,
               sig.p.vals,
               min.p)

  names(out) <- c("corr.dat",
                 "corr.sim",
                 "p.vals",
                 "sig.thresh",
                 "sig.snps.names",
                 "sig.snps",
                 "sig.corrs",
                 "sig.p.vals",
                 "min.p")

  print(paste("Finished running", test, "test @", Sys.time()))


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
                                "subsequent"),
                       correct.prop=FALSE,
                       categorical=FALSE){

  #################
  ## CORRELATION ##
  #################
  if(test=="cor"){
    corr.dat <- as.vector(cor(snps[, 1:ncol(snps)], phen)) # regular correlation...
    # corr.dat <- as.vector(cor(snps, phen)) # ? identical ??

  } # end test cor

  #########################
  ## FISHER'S EXACT TEST ## (** SLOW! **)
  #########################
  if(test=="fisher"){
    ## (!) Consider using tryCatch to avoid stupid (unpredictable?)
    ## workspace (FEXACT, LDSTP) error? (soln: arg workspace=2e6 (slow)).
    # # if(class(tryCatch(
    #   corr.dat <- sapply(c(1:ncol(snps)),
    #                    function(e) fisher.test(as.factor(snps[,e]),
    #                                            y=as.factor(phen),
    #                                            alternative="two.sided")$p.value)
    #   # , silent=TRUE)) == "try-error") print(1)

    ## Must prevent errors w 1-level factors (eg. if only 0 or NA but no 1s in snps[,e]):
    levs <- sapply(c(1:ncol(snps)), function(e) length(which(!is.na(unique(snps[,e])))))
    toKeep <- which(levs >= 2)
    corr.dat <- rep(NA, ncol(snps))
    corr.dat[toKeep] <- sapply(c(1:ncol(snps[,toKeep])),
                               function(e) fisher.test(snps[,toKeep[e]],
                                                       y=phen,
                                                       alternative="two.sided")$p.value)
  } # end test fisher


  ################################################
  ## TERMINAL TEST (score 1: correlation score) ##
  ################################################
  if(test=="terminal"){
    corr.dat <- terminal.test(snps = snps, phen = phen,
                              correct.prop=correct.prop, categorical=categorical)
    # ~ Correlation "SCORE" =
    # ((nS1P1 + nS0P0) - (nS1P0 + nS0P1) / (n.total))
  } # end test terminal

  #######################
  ## SIMULTANEOUS TEST ##
  #######################
  if(test == "simultaneous"){
    corr.dat <- simultaneous.test(snps.reconstruction = snps, phen.reconstruction = phen, tree = tree, categorical=categorical)
  } # end test simultaneous

  #####################
  ## SUBSEQUENT TEST ##
  #####################
  if(test == "subsequent"){
    corr.dat <- subsequent.test(snps.reconstruction = snps, phen.reconstruction = phen, tree = tree,
                                correct.prop=correct.prop, categorical=categorical)
  } # end test subsequent


  ## USE ABSOLUTE VALUE?
  # corr.dat <- abs(corr.dat)

  return(corr.dat)
} # end assoc.test





#################
## .get.p.vals ##
#################
## get a p-value associated with every value of corr.sim/corr.dat given the distribution of corr.sim:

.get.p.vals <- function(corr.sim, corr.dat=NULL, fisher.test = FALSE){

  if(is.null(corr.sim)) stop("Cannot calculate p-values for null object.")
  if(is.null(corr.dat)) corr.dat <- corr.sim

  ## Use ABS Vals:
  corr.dat <- abs(corr.dat)
  corr.sim <- abs(corr.sim)

  # w cum dist fn(F) p-val (p) for a given value(T) is: p = 1 - F(T)
  cum.dist <- ecdf(corr.sim)

  if(fisher.test == FALSE){
    ## For all but Fisher test, we care about how may null/simulated/chance values are GREATER than corr.dat[e]
    p <- 1 - cum.dist(corr.dat)
  }else{
    ## Fisher test works w p-values itself.
    ## So we care about how many values of the dist are SMALLER than corr.dat[e]
    p <- cum.dist(corr.dat)
  }

  return(p)
}  # end .get.p.vals



