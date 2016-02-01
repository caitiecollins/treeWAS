

##############
## phen.sim ##
##############


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
#' The true number of phenotypic substitions is drawn from a Poisson distribution with parameter \texttt{n.subs}.
#'
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


phen.sim <- function(tree, n.subs=15){

  ####################################
  ## PHENOTYPE simulation procedure ## ~ sim.by.locus...
  ####################################

  ## simulate phenotype for root individual:
  if(!is.null(n.subs)){
    phen.root <- "A"
  }else{
    phen.root <- NULL
  }

  ## store the inputted desired number of phenotypic substitutions
  n.phen.subs <- n.subs

  ## make dummy variables in which to store the resulting n.mts variables:
  lambda_p <- n.subs <- NA

  ## ensure phen variables start as NULL
  phen.branch <- phen.nodes <- phen.leaves <- NULL



  #############################################################
  ## If the user has specified a "mt" rate for the phenotype ##
  #############################################################

  ## (indicating that they want to generate a NEW phenotype for the tree provided)
  if(!is.null(n.phen.subs)){

    ## draw the number of mutations to occur at each site:
    n.subs <- rpois(n=1, lambda=n.phen.subs)
    ## if n.subs==0 or ==1, re-sample
    while(n.subs <= 1){
      n.subs <- rpois(n=1, lambda=n.phen.subs)
    }

    ## draw the branches to which you will assign the
    ## n.subs to occur for the phenotype (~ branch length):
    phen.loci <- sample(c(1:length(tree$edge.length)),
                        n.subs, replace=FALSE, prob=tree$edge.length)
    ## rearrange phen.loci
    phen.loci <- sort(phen.loci, decreasing=TRUE)


    ###############################
    ## For Loop to get PHENOTYPE ##
    ###############################
    ## get phenotype for all branches/ nodes in tree
    ## (from root node (ie. tree$edge[nrow(tree$edge), 1]) down):
    phen.nodes <- phen.branch <- list()

    ## set phenotype for all branches and nodes to be phen.root:
    for(i in 1:length(tree$edge.length)) phen.branch[[i]] <- phen.root
    names(phen.branch) <- paste("e", c(1:length(phen.branch)), sep=".")
    for(i in 1:length(unique(as.vector(unlist(tree$edge))))) phen.nodes[[i]] <- phen.root
    names(phen.nodes) <- paste("n", c(1:length(phen.nodes)), sep=".")

    #############################################################################

    #############################################################################

    ## get the node names for all individuals (terminal and internal)
    all.inds <- sort(unique(as.vector(unlist(tree$edge))))

    x <- rev(c(1:nrow(tree$edge)))

    ## get phen of nodes
    for(i in 1:length(x)){
      if(x[i] %in% phen.loci){
        phen.nodes[[tree$edge[x[i],2]]] <- .switch.phen(phen.nodes[[tree$edge[x[i],1]]])
      }else{
        ## if no phen subs occur on branch i, set phen of
        ## downstram individual to be equal to ancestor's
        phen.nodes[[tree$edge[x[i],2]]] <- phen.nodes[[tree$edge[x[i], 1]]]
      }
    } # end for loop

    ## get phen of TERMINAL nodes (leaves)
    n.ind <- tree$Nnode+1
    phen.leaves <- as.factor(as.vector(unlist(phen.nodes[c(1:n.ind)])))
    names(phen.leaves) <- paste("ind", c(1:length(phen.leaves)), sep=".")

    ## get phen of branches
    for(i in 1:length(x)){
      ## Branches with ONE phenotype get labelled by that phenotype:
      if(length(unique(phen.nodes[tree$edge[x[i],]])) == 1){
        if("A" %in% phen.nodes[tree$edge[x[i],]]){
          phen.branch[[x[i]]] <- "A"
        }else{
          phen.branch[[x[i]]] <- "B"
        }
      }else{
        ## Branches with TWO phenotypes get labelled as such, in ORDER:
        temp <- as.vector(unlist(phen.nodes[tree$edge[x[i],]]))
        if(temp[1] == "A"){
          phen.branch[[x[i]]] <- c("A", "B")
        }else{
          phen.branch[[x[i]]] <- c("B", "A")
        }
      }
    } # end for loop

  } ## end PHEN sim procedure...

  ## convert phen.nodes to factor
  phen.nodes <- as.factor(as.vector(unlist(phen.nodes)))
  if(!is.null(names(phen.leaves))){
    names(phen.leaves) <- paste("ind",
                                c(1:length(phen.leaves)),
                                sep=".")
  }
  if(!is.null(names(phen.nodes))){
    names(phen.nodes) <- c(names(phen.leaves),
                           paste("node",
                                 c((length(phen.leaves)+1):length(phen.nodes)),
                                 sep="."))
  }
  ## make output list
  phen.list <- list(phen.leaves, phen.nodes, phen.branch, phen.loci)
  names(phen.list) <- c("phen", "phen.nodes", "phen.edges", "phen.loci")

  return(phen.list)
} # end phen.sim
