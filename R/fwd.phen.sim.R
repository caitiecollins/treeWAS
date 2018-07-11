

##################
## fwd.phen.sim ##
##################

## TO DO ##
## CAREFUL--phen.sim seems not to be working with trees other than those
## produced with your coalescent.tree.sim fn (eg. rtree(100))!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param tree An phylo object.
#' @param n.subs An integer controlling the phenotypic substition rate (see details).
#'
#' @description The parameter n.subs controls the simulation of the phenotype by specifying
#' the expected value of the number of phenotypic substitions to occur on the tree provided.
#' The true number of phenotypic substitions is drawn from a Poisson distribution with parameter n.subs.
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @examples
#'
#' ## basic use of fn
#' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
#'
#' ## plot output
#' plot(tree)
#'
#' @export

########################################################################

## TO DO: ##
## Add arg continuous=FALSE --> continuous phen sim.
## Implement ASR for the snps.assoc only --> get phen for internal nodes too.


## OPTIONS: ##
## Cumulative probability (eg. if 7/10 SNPs, 70% chance of phen)--may allow for lots of noise...
## Threshold (eg. must have 7/10 SNPs to have phen)--will allow for lots of noise.
## Specific combinations (eg. Must have SNPs 1&2 OR 3&4)--will be very hard for treeWAS.
## Combination of above

fwd.phen.sim <- function(snps.assoc, p=1, tree=NULL){


  n.snps.assoc <- sapply(c(1:nrow(snps.assoc)),
                            function(e)
                              length(which(snps.assoc[e,] == 1)))
  ####################
  ## .get.phen.prob ##
  ####################
  .get.phen.prob <- function(n.snps.assoc, p){
    if(p == 1){
      ys <- n.snps.assoc/ncol(snps.assoc)
    }else{
      ys <- (1-p^n.snps.assoc)/(1-p^ncol(snps.assoc))
    }
    return(ys)
  } # end .get.phen.prob

  phen.prob <- sapply(c(1:length(n.snps.assoc)),
                      function(e)
                        .get.phen.prob(n.snps.assoc[e], p))
  phen <- as.factor(
                sapply(c(1:length(phen.prob)),
                 function(e)
                   sample(c("A", "B"),
                          size=1,
                          replace=TRUE,
                          prob=c(phen.prob[e], 1-phen.prob[e]))))

  #   ###############
  #   ## HISTOGRAM ##
  #   ###############
  #   hist(.get.phen.prob(n.snps.assoc=n.snps.assoc, p=p),
  #        breaks=10, col="blue", xlim=c(0,1),
  #        main=paste("Histogram of Pr(phen)
  #                   \n p = ", p, sep=""))
  #
  #   ################
  #   ## PROB CURVE ##
  #   ################
  #   plot_prob_phen(p=p, n.snps.assoc=ncol(snps.assoc))


  ###############
  ## w p = 0.6 ##
  ###############

  ###########
  ## TABLE ##
  ###########
  #table(phen)
  #   A  B
  #   95  5

  #########################
  ## CORRELATION (SCORE) ##
  #########################
  #abs(corr.dat[snps.assoc])
  #0.08 0.20 0.04 0.18 0.26 0.20 0.28 0.40 0.36 0.08


  ###############
  ## w p = 0.8 ##
  ###############

  ###########
  ## TABLE ##
  ###########
  #table(phen)
  #   A  B
  #   75 25

  #########################
  ## CORRELATION (SCORE) ##
  #########################
  #abs(corr.dat[snps.assoc])
  #0.16 0.08 0.08 0.14 0.10 0.12 0.36 0.24 0.32 0.08

  #############
  ## w p = 1 ##
  #############

  ###########
  ## TABLE ##
  ###########
  #table(phen)
  #   A  B
  #   49 51

  #########################
  ## CORRELATION (SCORE) ##
  #########################
  #abs(corr.dat[snps.assoc])
  #0.04 0.04 0.12 0.02 0.02 0.08 0.00 0.04 0.08 0.12

  ###############
  ## w p = 1.2 ##
  ###############

  ###########
  ## TABLE ##
  ###########
  #table(phen)
  #   A  B
  #   37 63

  #########################
  ## CORRELATION (SCORE) ##
  #########################
  #abs(corr.dat[snps.assoc])
  #0.08 0.00 0.12 0.10 0.06 0.04 0.04 0.28 0.16 0.04


  return(phen)

} # end fwd.phen.sim


