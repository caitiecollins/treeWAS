


###############
## read.CFML ##
###############

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
#' @import ape adegenet

########################################################################

## Not to be confused with the original readCFML fn (which returns only dist),
## this is the NEW read.CFML...

## This function reads the output of CFML and returns a LIST containing:
## the distribution of substitutions per site (dist),
## the phylogenetic tree, converted to be binary and post-ordered (tree),
## the set of sequences containing no duplicate column patterns (seqs),
## the index of all the original sequences in the set of unique sequence columns (mapping).

read.CFML <- function(prefix, plot=TRUE) {

  # require(ape)
  # require(adegenet)

  tree<-read.tree(sprintf('%s.labelled_tree.newick',prefix))
  seqs<-read.dna(sprintf('%s.ML_sequence.fasta',prefix),format='fasta')
  mapping<-scan(sprintf('%s.position_cross_reference.txt',prefix),sep=',',quiet=T)
  l<-length(seqs[1,])

  ## check tree (violations cause problems, eg. in phen.sim)
  if(!is.binary.tree(tree)) tree <- multi2di(tree)
  if(!identical(tree, reorder.phylo(tree,"postorder"))) tree <- reorder.phylo(tree,"postorder")

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

  ## Plot distribution
  if(plot==TRUE){
    barplot(dist,
            main="Number of substitutions per site",
            names.arg=c(1:length(dist)),
            ylab="Frequency",
            col=transp("royalblue", alpha=0.5))
  }

  out <- list(tree = tree,
              seqs = seqs,
              index = mapping,
              dist = dist)

  return(out)

} # end read.CFML
















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
#' @import ape adegenet

########################################################################

######################################
# dist_0 <- readCFML(prefix="~/ClonalFrameML/src/CFML.R.0")
# str(dist_0) #
# barplot(dist_0, main="Number of substitutions per site \n R = 0",
#         names.arg=c(1:length(dist_0)), ylab="Frequency",  col=transp("royalblue", alpha=0.5))
# save(dist_0, file="~/ClonalFrameML/src/CFML_R_0_dist.Rdata")
######################################
# setwd("/media/caitiecollins/Seagate Backup Plus Drive/SimBac")
# dist_0.1 <- readCFML(prefix="./CFML_R_0.1")
# str(dist_0.1)
# save(dist_0.1, file="./CFML_R_0.1_dist.Rdata")
# barplot(dist_0.1, main="Number of substitutions per site \n R = 0.1",
#         names.arg=c(1:length(dist_0.1)), ylab="Frequency", col=transp("royalblue", alpha=0.5))
######################################
# data(dist)
# barplot(dist, main="Number of substitutions per site",
#         names.arg=c(1:length(dist)), ylab="Frequency", col=transp("royalblue", alpha=0.5))
# title(substitute(paste("(", italic("S. aureus"), " empirical dataset)", sep="")), line=0.5)
######################################
# set.seed(1)
# dist <- rpois(n=1000000, lambda=1)
# for(i in 1:length(dist)){
#   while(dist[i]==0){
#     dist[i] <- rpois(n=1, lambda=1)
#   }
# }
# dist <- table(dist) # NB: always check to make sure you have no missing n.subs in 1:max(dist)
# barplot(dist, main="Number of substitutions per site",
#         names.arg=c(1:length(dist)), ylab="Frequency", col=transp("royalblue", alpha=0.5))
# title(substitute(paste("(Poisson distribution, ", italic("lambda"), " = 1)", sep="")), line=0.5)
######################################


#This function reads the output of CFML and outputs the distribution of substitutions per site
readCFML <- function(prefix, plot=TRUE) {

  # require(ape)
  # require(adegenet)

  tree<-read.tree(sprintf('%s.labelled_tree.newick',prefix))
  seqs<-read.dna(sprintf('%s.ML_sequence.fasta',prefix),format='fasta')
  mapping<-scan(sprintf('%s.position_cross_reference.txt',prefix),sep=',',quiet=T)
  l<-length(seqs[1,])

  ## check tree (violations cause problems, eg. in phen.sim)
  if(!is.binary.tree(tree)) tree <- multi2di(tree)
  if(!identical(tree, reorder.phylo(tree,"postorder"))) tree <- reorder.phylo(tree,"postorder")

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

  ## Plot distribution
  if(plot==TRUE){
    barplot(dist,
            main="Number of substitutions per site",
            names.arg=c(1:length(dist)),
            ylab="Frequency",
            col=transp("royalblue", alpha=0.5))
  }


  return(dist)
}
