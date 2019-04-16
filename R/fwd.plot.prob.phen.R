
####################
## plot.prob.phen ##
####################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Plot the probability of association, given \code{p} and \code{n.snps.assoc}.
#'
#' [*For use with the 'fwd.-.sim' functions:*] 
#' Plot the cumulative probability of association (Pr(phen=1)), with a given value of \code{p},
#' as the number of associated sites (SNPi=1) increases from i=0 to i=\code{n.snps.assoc}.
#'
#' @param p A numeric value indicating the probability of substitution, at each site, along the tree. 
#' @param n.snps.assoc An integer specifying the number of genetic loci that are associated with the phenotype.
#' 
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @examples
#' \dontrun{
#' ## basic use of fn ##
#' ## compare probability of having phenotype with 10 SNPs at varying p:
#' plot.prob.phen(p=0.8, n.snps.assoc=10)
#' plot.prob.phen(p=0.5, n.snps.assoc=10)
#' plot.prob.phen(p=0.2, n.snps.assoc=10)
#' }
#' @export

########################################################################


plot.prob.phen <- function(p=0.5, n.snps.assoc=10){

  xs <- 0:n.snps.assoc
  if(p == 1){
    ys <- xs/n.snps.assoc
  }else{
    ys <- (1-p^xs)/(1-p^10)
  }

  ## plot ##
  plot(xs,ys,xlim=c(0,10),ylim=c(0,1),
       main=paste("p = ", p, sep=""),
       xlab="Number of associated sites in state 1",
       ylab="Cumulative probability of association")

} # end plot.prob.phen




#################################
##  ENABLE ALTERNATE FN NAME:  ##
#################################
plot_prob_phen <- function(p, n.snps.assoc, ...){
  return(plot.prob.phen(p, n.snps.assoc,  ...))
} # end plot_prob_phen





