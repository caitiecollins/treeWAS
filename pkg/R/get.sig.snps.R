
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
                         snps.unique = NULL,
                         snps.index = NULL,
                         snps.sim,
                         snps.sim.unique = NULL,
                         snps.sim.index = NULL,
                         phen,
                         tree,
                         test = "terminal",
                         n.tests = 1,
                         p.value = 0.01,
                         p.value.correct = "bonf",
                         p.value.by = "count",
                         snps.reconstruction,
                         snps.sim.reconstruction,
                         phen.reconstruction,
                         rec = "parsimony"){

  #################
  ## HANDLE SNPS ##
  #################

  ##########
  ## snps ##
  ##########

  if(is.null(snps.unique)){
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
  }else{
    ## If snps.unique provided as well as snps:

    ## CHECK: is index provided as well??
    if(is.null(snps.index)){
      warning("if snps.unique is provided,
              snps.index must also be provided to indicate
              original mapping locations for all unique sites.
              Ignoring unique snps provided; working with snps only.")

      ## repeat above steps (as if no snps.unique was provided):
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

    }else{

      snps.ori <- snps
      snps <- snps.unique

      ## record whether all snps are unique or not for later:
      if(ncol(snps.unique) == ncol(snps.ori)){
        all.unique <- TRUE
      }else{
        all.unique <- FALSE
      }
    }
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
  if(is.null(snps.sim.unique)){
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
  }else{
    ## If snps.sim.unique provided as well as snps:

    ## CHECK: is index provided as well??
    if(is.null(snps.sim.index)){
      warning("if snps.sim.unique is provided,
              snps.sim.index must also be provided to indicate
              original mapping locations for all unique sites.
              Ignoring unique snps.sim provided; working with snps.sim only.")

      ## repeat above steps (as if no snps.unique was provided):
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
    }else{

      snps.sim.ori <- snps.sim
      snps.sim <- snps.sim.unique

      ## record whether all snps are unique or not for later:
      if(ncol(snps.sim.unique) == ncol(snps.sim.ori)){
        all.unique.sim <- TRUE
      }else{
        all.unique.sim <- FALSE
      }
    }
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

    ############################
    ## HANDLE RECONSTRUCTIONS ##
    ############################

    ## NOTE: If snps(.sim) are UNIQUE but snps(.sim).reconstruction are NOT, test will give INCORRECT OUTPUT!!!

    if(all.unique == FALSE){
      ## check if snsp.rec is already in UNIQUE form:
      if(ncol(snps.reconstruction) != ncol(snps.unique)){
        temp <- get.unique.matrix(snps.reconstruction, MARGIN=2)
        snps.reconstruction <- temp$unique.data
        snps.reconstruction.index <- temp$index
        if(!identical(snps.reconstruction.index, snps.index)){
          warning("Careful-- snps and snps.reconstruction should have the same index when reduced
                  to their unique forms!") ## SHOULD THIS BE A "STOP" INSTEAD? OR IS THIS ERROR NOT FATAL OR NOT POSSIBLE????
        }
        }else{
          snps.reconstruction.index <- snps.index
      }
    }

    if(all.unique.sim == FALSE){
      ## check if snsp.rec is already in UNIQUE form:
      if(ncol(snps.reconstruction) != ncol(snps.unique)){
        temp <- get.unique.matrix(snps.sim.reconstruction, MARGIN=2)
        snps.sim.reconstruction <- temp$unique.data
        snps.sim.reconstruction.index <- temp$index
        if(!identical(snps.sim.reconstruction.index, snps.sim.index)){
          warning("Careful-- snps.sim and snps.sim.reconstruction should have the same index when reduced
                  to their unique forms!") ## SHOULD THIS BE A "STOP" INSTEAD? OR IS THIS ERROR NOT FATAL OR NOT POSSIBLE????
        }
        }else{
          snps.sim.reconstruction.index <- snps.sim.index
      }
    }

    ########################################################
    ## Calculate correlations btw REAL SNPs and phenotype ##
    ########################################################
    corr.dat <- assoc.test(snps=snps.reconstruction, phen=phen.reconstruction, tree=tree, test=test, rec = rec)



    #############################################################
    ## Calculate correlations btw SIMULATED SNPs and phenotype ##
    #############################################################
    corr.sim <- assoc.test(snps=snps.sim.reconstruction, phen=phen.reconstruction, tree=tree, test=test, rec = rec)
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

  ## quick look at corr.sim & corr.dat
  # hist(corr.sim, xlim=c(0,1))
  # hist(corr.dat, xlim=c(0,1))


  ############################
  ## HANDLE P.VALUE OPTIONS ##
  ############################

  out <- nom <- list()
  corr.dat.ori <- corr.dat
  corr.sim.ori <- corr.sim


  ## n.snps.sim ##
  corr.sim <- corr.sim.ori


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
    p.vals <- .get.p.vals(corr.sim,
                          corr.dat = NULL,
                          fisher.test = TRUE)
    p.vals <- sort(p.vals, decreasing=TRUE)
    p.fdr <- p.adjust(p.vals, method="fdr", n=length(p.vals)*n.tests) # CHECK--IS THIS THE CORRECT MT APPROACH FOR FDR????????????????????????
    p.thresh <- quantile(p.fdr, probs=1-p.value)

    ################
    ## p.value.by ##
    ################
    if(p.value.by == "count") thresh <- quantile(corr.sim, probs=p.thresh)
    if(p.value.by == "density") thresh <- quantile(density(corr.sim)$x, probs=p.thresh)

  } # end p.value.correct == "fdr"


  ##################################
  ## GET SIGNIFICANT CORRELATIONS ##
  ##################################
  if(test=="fisher"){
    ## Identify (real) SNPs w correlations > thresh:
    sig.snps <- which(corr.dat < thresh)
    p.vals <- .get.p.vals(corr.sim = corr.sim,
                          corr.dat = corr.dat,
                          fisher.test = TRUE)
    sig.p.vals <- p.vals[sig.snps]
  }else{
    ## Identify (real) SNPs w correlations > thresh:
    sig.snps <- which(corr.dat > thresh)
    p.vals <- .get.p.vals(corr.sim = corr.sim,
                          corr.dat = corr.dat,
                          fisher.test = FALSE)
    sig.p.vals <- p.vals[sig.snps]
  }



  ## 0 p.vals
  min.p <- paste("p-values listed as 0 are <",
                 1/length(corr.sim), sep=" ") ## CHECK---IS THIS RIGHT? SHOULD WE BE MULTIPLYING THE DIVISOR BY N.TESTS ??????????

  ## get list of those correlation values
  sig.corrs <- corr.dat[sig.snps]
  ## get the list of those SNPs (ie. their locus names)
  # sig.snps <- dimnames(snps)[[2]][sig.snps]
  sig.snps.names <- dimnames(snps.ori)[[2]][sig.snps]

  ## re-order list of sig.snps and sig.corrs by value of sig.corr
  if(test=="fisher"){
    NWO <- order(sig.corrs, decreasing=FALSE)
  }else{
    NWO <- order(sig.corrs, decreasing=TRUE)
  }
  sig.snps <- sig.snps[NWO]
  sig.corrs <- sig.corrs[NWO]
  sig.p.vals <- sig.p.vals[NWO]
  sig.snps.names <- sig.snps.names[NWO]
  gc()

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
                                "fisher"),
                       rec = "parsimony"){

  #################
  ## CORRELATION ##
  #################
  if(test=="cor"){
    corr.dat <- as.vector(cor(snps[, 1:ncol(snps)], phen)) # regular correlation...
    # corr.dat <- as.vector(cor(snps, phen)) # ? identical ??

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


  ################################################
  ## TERMINAL TEST (score 1: correlation score) ##
  ################################################
  if(test=="terminal"){
    corr.dat <- terminal.test(snps = snps, phen = phen)
    # ~ Correlation "SCORE" =
    # ((nS1P1 + nS0P0) - (nS1P0 + nS0P1) / (n.total))
    ## must be calculated for each SNP individually...
    #     corr.dat <- sapply(c(1:ncol(snps)), function(e)
    #       (((length(which(snps[which(phen==1),e]==1)) +
    #            length(which(snps[which(phen==0),e]==0)))
    #         - (length(which(snps[which(phen==1),e]==0)) +
    #              length(which(snps[which(phen==0),e]==1))))
    #        / nrow(snps)))
  } # end test terminal

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
    corr.dat <- subsequent.test(snps.reconstruction = snps, phen.reconstruction = phen, tree = tree, rec = rec)
  } # end test subsequent



  ## USE ABSOLUTE VALUE
  # if(test != "subsequent"){
    corr.dat <- abs(corr.dat)
  # }

  return(corr.dat)
} # end assoc.test





#################
## .get.p.vals ##
#################
## NOTE: only used for FDR threshold calculation!
## get a p-value associated with every value of corr.sim

.get.p.vals <- function(corr.sim, corr.dat=NULL, fisher.test=FALSE){

  p.vals <- NULL

  ###################
  ## CORR.SIM ONLY ##
  ###################
  if(is.null(corr.dat)){
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
  }else{

    ###########################
    ## CORR.DAT vs. CORR.SIM ##
    ###########################
    ## get table for corr.sim:
    cs.tab <- table(corr.sim)

    ## get table for corr.dat:
    cd.tab <- table(corr.dat)
    cd.fac <- factor(corr.dat)

    ## get all possible p.vals | corr.sim:
    p.vals.unique <- sapply(c(1:(length(cs.tab)-1)),
                            function(e)
                              sum(cs.tab[(e+1):length(cs.tab)])
                            /sum(cs.tab))
    ## need to add trailing 0 for max corr.sim:
    p.vals.unique <- c(p.vals.unique, 0)

    ## get real p.vals using corr.dat for break points:
    p.vals.dat <- list()


    if(fisher.test == FALSE){
      for(i in 1:length(cd.tab)){
        tab.above <- which(as.numeric(names(cs.tab)) > as.numeric(names(cd.tab[i])))
        if(!.is.integer0(tab.above)){
          p.vals.dat[[i]] <- sum(cs.tab[tab.above])/sum(cs.tab)
        }else{
          p.vals.dat[[i]] <- 0
        }
      }
    }else{
      ## fisher.test == TRUE --> Reverse sign (> --> <) ##
      for(i in 1:length(cd.tab)){
        tab.below <- which(as.numeric(names(cs.tab)) < as.numeric(names(cd.tab[i])))
        if(!.is.integer0(tab.below)){
          p.vals.dat[[i]] <- sum(cs.tab[tab.below])/sum(cs.tab)
        }else{
          p.vals.dat[[i]] <- 0
        }
      }
    }

    p.vals.dat <- as.vector(unlist(p.vals.dat))

    ## get the unique index that
    ## each original corr.dat should map to:
    map.to <- (as.integer(cd.fac) - 1)
    map.to <- map.to[!is.na(map.to)]
    if(length(map.to)) map.to <- map.to + 1

    ## get p.vals for all corr.dat:
    p.vals <- p.vals.dat[map.to]
  }

  return(p.vals)
} # end .get.p.vals
