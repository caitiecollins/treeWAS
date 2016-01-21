
##############
## readCFML ##
##############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param prefix A character string containing the prefix of all file names to be read in.
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
#' fn(arg1, arg2)
#'
#' #' ## more elaborate use of fn
#' fn(arg1, arg2)
#'
#' @import ape

########################################################################




#This function reads the output of CFML and outputs the distribution of substitutions per site
readCFML <- function(prefix) {

  require(ape)

  tree<-read.tree(sprintf('%s.labelled_tree.newick',prefix))
  seqs<-read.dna(sprintf('%s.ML_sequence.fasta',prefix),format='fasta')
  mapping<-scan(sprintf('%s.position_cross_reference.txt',prefix),sep=',',quiet=T)
  l<-length(seqs[1,])

  #Modify the edge matrix so that it uses the same indices as the fasta file
  labs<-labels(seqs)
  edges<-tree$edge
  treelabs<-c(tree$tip.label,tree$node.label)
  for (i in 1:nrow(edges)) for (j in 1:ncol(edges)) edges[i,j]=which(labs==treelabs[edges[i,j]])

  #Count substitutions for each site
  subs=rep(0,l)#Number of substitutitions for patterns
  num=rep(0,l)#Number of times a given pattern is used
  for (i in 1:l) {
    if (length(unique(seqs[,i]))==2) num[i]=sum(mapping==i) #Only count biallelic sites
    for (b in 1:nrow(edges)) if (seqs[edges[b,1],i]!=seqs[edges[b,2],i]) subs[i]=subs[i]+1
  }

  #Build distribution
  dist=rep(0,max(subs))
  for (i in 1:max(subs)) dist[i]=sum(num[which(subs==i)])
  return(dist)
}
