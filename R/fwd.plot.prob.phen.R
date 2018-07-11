
####################
## plot_prob_phen ##
####################

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
#' ## basic use of fn ##
#' ## compare probability of having phenotype with 10 SNPs at varying p:
#' plot_prob_phen(p=0.8, n.snps.assoc=10)
#' plot_prob_phen(p=0.5, n.snps.assoc=10)
#' plot_prob_phen(p=0.2, n.snps.assoc=10)
#'
#' @export

########################################################################

# plot_prob_phen(p=1.8, n.snps.assoc=10)
# plot_prob_phen(p=1.5, n.snps.assoc=10)
# plot_prob_phen(p=1.4, n.snps.assoc=10)
# plot_prob_phen(p=1, n.snps.assoc=10)
# plot_prob_phen(p=0.8, n.snps.assoc=10)
# plot_prob_phen(p=0.5, n.snps.assoc=10)
# plot_prob_phen(p=0.2, n.snps.assoc=10)


plot_prob_phen <- function(p=0.5, n.snps.assoc=10){

  xs <- 0:n.snps.assoc
  if(p == 1){
    ys <- xs/n.snps.assoc
  }else{
    ys <- (1-p^xs)/(1-p^10)
  }

  ## plot ##
  plot(xs,ys,xlim=c(0,10),ylim=c(0,1),
       main=paste("p = ", p, sep=""))

} # end plot_prob_phen


