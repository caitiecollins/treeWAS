

####################
## heatmap.DNAbin ##
####################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param dna A DNAbin object.
#' @param dist.dna.model A character string specifying the type of model to use in
#' calculating the genetic distance between individual genomes (see ?dist.dna).
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#'
#' @rawNamespace import(ape, except = zoom)
#' @export

########################################################################

heatmap.DNAbin <- function(dna, dist.dna.model="JC69"){

  # require(ape)

  if(!"DNAbin" %in% class(dna)) dna <- as.DNAbin(dna)

  #############
  ## HEATMAP ##
  #############
  ## get a distance matrix between the genomes
  D <- dist.dna(dna, model = dist.dna.model)

  mat <- t(as.matrix(D))
  mat <- mat[,ncol(mat):1]
  par(mar=c(1,5,5,1))
  image(x=1:ncol(mat), y=1:ncol(mat), mat,
        col=rev(heat.colors(100)),
        xaxt="n", yaxt="n", xlab="", ylab="")
  axis(side=2, at=c(1:ncol(mat)),
       labels=rev(names(dna)), las=2, cex.axis=1)
  axis(side=3, at=c(1:ncol(mat)),
       labels=names(dna), las=1, cex.axis=1)
  ## return margin parameter to default:
  par(mar=c(5,4,4,2)+0.1)

} # end heatmap.DNAbin
