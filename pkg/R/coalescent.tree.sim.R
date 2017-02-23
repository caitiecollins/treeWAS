

#########################
## coalescent.tree.sim ##
#########################


########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param n.ind An integer specifying the number of terminal nodes desired.
#' @param seed An optional integer controlling the pseudo-random process underlying the tree generation.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @examples
#'
#' ## basic use of fn
#' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
#'
#' ## plot output
#' plot(tree)
#'

########################################################################



coalescent.tree.sim <- function(n.ind=100, seed=NULL){

  if(!is.null(seed)) set.seed(seed)

  n.nodes <- n.ind + (n.ind-1) # total n.nodes (internal, external)
  inds <- c(1:n.ind) # terminal nodes
  nodes <- rev(c((n.ind+1):n.nodes)) # internal nodes
  tree.params <- list() # to store and update output


  for(i in 1:(length(inds)-1)){
    ## get inds.ori from last generation:
    if(i==1){
      inds.ori <- inds
    }else{
      inds.ori <- tree.params[[(i-1)]][["inds.remaining"]]
    }
    ## get N, the number of individuals (remaining at this generation)
    ## from which a random 2 are to be selected for coalescence
    N <- length(inds.ori)

    ###################
    ## BRANCH LENGTH ##
    ###################
    ## get lamda, the parameter of the exponential distribution,
    ## given the number of individuals at this generation
    lambda <- (N*(N-1)) / 2
    ## draw x, the length of time to coalescence at this generation
    x <- rexp(n=1, rate=lambda)

    #####################
    ## COALESCENT PAIR ##
    #####################
    ## get co.pair, the 2 inds to coalesce at this generation
    co.pair <- sample(inds.ori, 2)
    ## merge these 2 inds, replace with new internal node,
    ## update the list of inds to sample at the next generation
    inds.remaining <- c(inds.ori[-which(inds.ori %in% co.pair)], nodes[i])

    ############
    ## OUTPUT ##
    ############
    ## store the output in the ith element of our list tree.params:
    tree.params[[i]] <- list()
    tree.params[[i]][[1]] <- x
    tree.params[[i]][[2]] <- co.pair
    tree.params[[i]][[3]] <- inds.ori
    tree.params[[i]][[4]] <- inds.remaining
    names(tree.params[[i]]) <- c("Time", "co.pair", "inds.ori", "inds.remaining")
  } # end for loop




  ## get edge.list
  to <- as.vector(unlist(sapply(c(1:length(tree.params)),
                                function(e) tree.params[[e]][["co.pair"]])))
  from <- as.vector(unlist(sapply(c(1:length(nodes)),
                                  function(e) rep(nodes[e], 2))))
  edge.list <- data.frame(from,to)

  ## get edge lengths
  times <- as.vector(unlist(sapply(c(1:length(tree.params)),
                                   function(e) tree.params[[e]][["Time"]])))

  ## make empty edge.lengths vector to store output below:
  edge.lengths <- NA

  ## for all the edges in our edge.list data.frame:
  for(i in 1:nrow(edge.list)){
    if(edge.list$to[i] %in% inds){
      ## if downstream node = terminal, sum all time intervals til ancestor.
      edge.lengths[i] <- sum(times[1:which(nodes==edge.list$from[i])])
    }else{
      ## BUT, if the downstream node = internal, must subtract time btw.
      ## downstream node and final generation.
      length.total <- sum(times[1:which(nodes==edge.list$from[i])])
      length.toRemove <- sum(times[1:which(nodes==edge.list$to[i])])
      edge.lengths[i] <- length.total - length.toRemove
    }
  } # end for loop

  ## convert edge.list to matrix
  edge.list <- as.matrix(edge.list)
  colnames(edge.list) <- NULL
  dimnames(edge.list) <- NULL

  ## put output into tree list (phylo format):
  tree <- list()
  tree$edge <- edge.list
  tree$tip.label <- c(1:n.ind)
  tree$edge.length <- edge.lengths
  tree$Nnode <- as.integer(n.ind - 1)

  ## change class by force
  class(tree) <- "phylo"
  ## return tree in pruningwise order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## root tree:
  if(!is.rooted(tree)) tree <- midpoint(tree)

  return(tree)

} # end coalescent.tree.sim
