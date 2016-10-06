



############################################################################################################################################




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

read.CFML <- function(prefix, tree=NULL, plot=TRUE) {

  # require(ape)
  # require(adegenet)

  if(is.null(tree)){
    tree <- read.tree(sprintf('%s.labelled_tree.newick', prefix))
  }else{
    ## if tree is filename, read in:
    if(class(tree) == "character"){
      tree <- read.tree(tree)
    }
  }
  seqs <- read.dna(sprintf('%s.ML_sequence.fasta', prefix), format='fasta')
  mapping <- scan(sprintf('%s.position_cross_reference.txt', prefix), sep=',', quiet=T)
  l <- length(seqs[1,])

  ## check tree (violations cause problems, eg. in phen.sim)
  if(!is.binary.tree(tree)) tree <- multi2di(tree)
  # if(!identical(tree, reorder.phylo(tree,"postorder"))) tree <- reorder.phylo(tree,"postorder") ## not sure this makes any difference..

  ## Modify the edge matrix so that it uses the same indices as the fasta file
  labs <- labels(seqs)
  edges <- tree$edge
  ## node.label causing problems (?!) ###
  ## Not sure if order wrong, but think below is wrong bc., while tree$tip.lab is in order corresponding
  ## correctly to edge mat, internal node labs are not
  ## (Hyp-- edge mat expects internal nodes to be labelled in order, ie. 111 --> NODE_111 ??? ---> (BUT MAYBE NOT???!!!!!)).
  treelabs <- c(tree$tip.label, tree$node.label)
  # treelabs <- tree$tip.label
  # treelabs.ori <- treelabs
  # edges.ori <- edges

  ## If treelabs not of length labs (eg tree$node.label is NULL), fill in remainder w labs from seqs:
  # if(length(treelabs) < length(labs))  treelabs <- c(treelabs, labs[c((length(treelabs)+1):length(labs))])
  if(length(treelabs) > length(labs) | length(treelabs) < length(labs)){
      stop("The number of individuals (terminal and internal nodes) labelled in the tree
            exceeds the number of individuals (rows) labelled in the sequences.")
  }
  ## Double check that treelabs and labs match up (order does NOT matter yet):
  if(!all(labs %in% treelabs)){
    stop("Sequence labels do not match tree labels.")
  }
  ## Realign order of labs and treelabs within edge mat:
  for (i in 1:nrow(edges)) for (j in 1:ncol(edges)) edges[i,j] <- which(labs==treelabs[edges[i,j]])

  #Count substitutions for each site
  subs <- rep(0,l) # Number of substitutitions for patterns
  num <- rep(0,l) # Number of times a given pattern is used

  ## SLOW STEP!!!
  system.time(
  for (i in 1:l) {
    if (length(unique(seqs[,i]))==2){
      num[i] <- sum(mapping==i) # Only count biallelic sites
      for (b in 1:nrow(edges)) if (seqs[edges[b,1],i]!=seqs[edges[b,2],i]) subs[i] <- subs[i]+1
    } ## NOTE--if we're only counting biallelic sites for num, should do so fo subs as well...
  } # end for loop
  ) # end system.time

  ## Build distribution
  dist <- rep(0, max(subs)) # there should be no trailing zeros!
  ## trailing zeros still happening!
  ## bc.
  for (i in 1:max(subs)) dist[i] <- sum(num[which(subs==i)])
  ## assign names:
  names(dist) <- c(1:length(dist))

  ## Plot distribution
  if(plot==TRUE){
    barplot(dist,
            main="Number of substitutions per site",
            names.arg=names(dist),
            ylab="Frequency",
            col=transp("royalblue", alpha=0.5))
  }

  out <- list(tree = tree,
              seqs = seqs,
              index = mapping,
              dist = dist)

  return(out)

} # end read.CFML





#

##

####

########

#################################################################################################################################



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

##############
## readCFML ##
##############

## ORIGINAL readCFML FN-- RETURNS DIST ONLY w barplot.
## NOTE-- Use read.CFML (with a DOT) to get the LIST output!!

## This function reads the output of CFML and outputs the distribution of substitutions per site

readCFML <- function(prefix, plot=TRUE) {

  # require(ape)
  # require(adegenet)

  tree<-read.tree(sprintf('%s.labelled_tree.newick',prefix))
  seqs<-read.dna(sprintf('%s.ML_sequence.fasta',prefix),format='fasta')
  mapping<-scan(sprintf('%s.position_cross_reference.txt',prefix),sep=',',quiet=T)
  l<-length(seqs[1,])

  ## Check tree (violations cause problems, eg. in phen.sim)
  if(!is.binary.tree(tree)) tree <- multi2di(tree)
  if(!identical(tree, reorder.phylo(tree,"postorder"))) tree <- reorder.phylo(tree,"postorder")

  ## Modify the edge matrix so that it uses the same indices as the fasta file
  labs<-labels(seqs)
  edges<-tree$edge
  treelabs<-c(tree$tip.label,tree$node.label)
  for (i in 1:nrow(edges)) for (j in 1:ncol(edges)) edges[i,j]=which(labs==treelabs[edges[i,j]])

  ## Count substitutions for each site
  subs=rep(0,l) # Number of substitutitions for patterns
  num=rep(0,l) # Number of times a given pattern is used
  for (i in 1:l) {
    if (length(unique(seqs[,i]))==2) num[i]=sum(mapping==i) #Only count biallelic sites
    for (b in 1:nrow(edges)) if (seqs[edges[b,1],i]!=seqs[edges[b,2],i]) subs[i]=subs[i]+1
  }

  ## Build distribution
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
