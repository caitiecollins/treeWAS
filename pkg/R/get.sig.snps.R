
######################
## get.correlations ##
######################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param snps A matrix containing the real snps.
#' @param snps.sim A matrix or list of matrices containing simulated snps.
#' @param phen A factor or vector containing the phenotype (only allowed to contain two levels for now).
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
#' @export
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
#' @import adegenet ape phangorn

########################################################################

get.sig.snps <- function(snps, snps.sim, phen,
                         test=c("score", "cor", "fisher", "Pagel"),
                         p.value=0.001,
                         p.value.correct=c("bonf", "fdr", FALSE),
                         p.value.by=c("count", "density"),
                         Pagel=NULL, Pagel.sim=NULL){
  if(test != "Pagel"){
  #################
  ## HANDLE PHEN ##
  #################
  ## convert phenotype to numeric:
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

  ##################
  ## get.corrs fn ##
  ##################
  get.corrs <- function(snps, phen){

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


    ## USE ABSOLUTE VALUE
    corr.dat <- abs(corr.dat)

    return(corr.dat)
  } # end get.corrs


  ########################################################
  ## Calculate correlations btw REAL SNPs and phenotype ##
  ########################################################

  corr.dat <- get.corrs(snps=snps, phen=phen)

  }else{
    ################
    ## Pagel test ##
    ################
    corr.dat <- Pagel
  }


  #############################################################
  ## Calculate correlations btw SIMULATED SNPs and phenotype ##
  #############################################################
  if(test != "Pagel"){
  if(class(snps.sim) == "matrix"){
    corr.sim <- get.corrs(snps=snps.sim, phen=phen)
  }else{
    corr.sim <- list()

    for(i in 1:length(snps.sim)){
      if(!is.null(snps.sim[[i]])){
        corr.sim[[i]] <- get.corrs(snps=snps.sim[[i]], phen=phen)
      }else{
        corr.sim[[i]] <- NULL
      }
    } # end for loop

    corr.sim <- as.vector(unlist(corr.sim))
  } # end corr.sim
  }else{

    ################
    ## Pagel test ##
    ################
    corr.sim <- Pagel.sim
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
  if(test=="fisher" | test=="Pagel"){
    if(test =="Pagel") thresh <- 1 - thresh
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

  out <- list(corr.dat, corr.sim, thresh, sig.snps.names, sig.snps, sig.corrs, p.vals, min.p)
  names(out) <- c("corr.dat", "corr.sim", "sig.thresh", "sig.snps.names", "sig.snps", "sig.corrs", "p.vals", "min.p")

  return(out)


} # end get.sig.snps
