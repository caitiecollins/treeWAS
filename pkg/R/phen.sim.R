

##############
## phen.sim ##
##############

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
#' @param grp.min An optional numeric value < 0.5 specifying the minimum accepted proportion of terminal nodes
#' to be in the minor phenotypic group. It may be useful to specify a \code{grp.min} of,
#' for example, 0.2 (the default) to prevent excessive imbalance in the phenotypic group sizes. However,
#' it is important to note that (at least for the time being) \code{grp.min} values closer to
#' 0.5 are likely to cause the computational time of \code{phen.sim} to increase substantially,
#' as the function will run until acceptable group sizes are randomly generated.
#' @param seed An optional integer used to set the seed and control the pseudo-random process used in
#' \code{phen.sim}, enabling the repeatable regeneration of identical output.
#'
#' @description The parameter n.subs controls the simulation of the phenotype by specifying
#' the expected value of the number of phenotypic substitions to occur on the tree provided.
#' The true number of phenotypic substitions is drawn from a Poisson distribution with parameter n.subs.
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @examples
#'
#' ## basic use of fn
#' tree <- coalescent.tree.sim(n.ind = 100, seed = 1)
#'
#' ## plot output
#' plot(tree)
#'
#' @importFrom phangorn midpoint
#'
#' @export

########################################################################
# @useDynLib phangorn, .registration = TRUE

phen.sim <- function(tree,
                     n.subs = 15,
                     grp.min = 0.2,
                     # coaltree = TRUE,
                     seed = NULL){

  if(!is.null(seed)) set.seed(seed)

  ## HANDLE TREE: ##
  ## Always work with trees in "pruningwise" order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  if(!is.rooted(tree)) tree <- midpoint(tree)

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

    ## START WHILE LOOP HERE ###########

    toRepeat <- TRUE

    while(toRepeat == TRUE){

      ## draw the number of substitutions to occur:
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
      phen.branch[1:length(tree$edge.length)] <- phen.root
      names(phen.branch) <- paste("e", c(1:length(phen.branch)), sep=".")

      phen.nodes[1:length(unique(as.vector(unlist(tree$edge))))] <- phen.root
      names(phen.nodes) <- paste("n", c(1:length(phen.nodes)), sep=".")

      #############################################################################

      #############################################################################

      ## get the node INDICES for all individuals (terminal and internal)
      all.inds <- sort(unique(as.vector(unlist(tree$edge)))) # 1:(n.ind*2 - 1)

      ####################################################################
      ############################
      ## Get Anc-Des EDGE ORDER ##
      ############################
      ## Get sequence from lowest ("root", Nterm+1) to highest ancestral node:
      ix <- c(min(tree$edge[,1]):max(tree$edge[,1]))
      ## Get for loop index of rows in tree$edge[,1], in pairs, from lowest to highest:
      x <- as.vector(unlist(sapply(c(1:length(ix)), function(e) which(tree$edge[,1] == ix[e]))))
      ####################################################################



      ## get phen of nodes
      for(i in 1:length(x)){
        if(x[i] %in% phen.loci){
          phen.nodes[[tree$edge[x[i],2]]] <- .switch.phen(phen.nodes[[tree$edge[x[i],1]]])
        }else{
          ## if no phen subs occur on branch i, set phen of
          ## downstream individual to be equal to ancestor's
          phen.nodes[[tree$edge[x[i],2]]] <- phen.nodes[[tree$edge[x[i], 1]]]
        }
      } # end for loop

      ## get phen of TERMINAL nodes (leaves)
      n.ind <- tree$Nnode+1
      phen.leaves <- as.factor(as.vector(unlist(phen.nodes[c(1:n.ind)])))

      ## Assign names to phen as tree$tip.labs in original order ##
      ## If checks fail, individuals in phen.leaves will be named 1:N (not ideal)
      names(phen.leaves) <- c(1:length(phen.leaves))
      ## Check that tip.labs are not NULL:
      if(is.null(tree$tip.label)){
        warning("tree$tip.label was NULL.
                Assigning individuals names 1:N. Note that these may NOT match sequence labels!")
      }else{
        ## Check that tip.labs is of correct length:
        if(length(tree$tip.label) != length(phen.leaves)){
          warning("The length of tree$tip.label did not match
                  the number of terminal node phenotypes simulated.
                  Assigning individuals names 1:N. Note that these may NOT match sequence labels!")
        }else{
          ## If checks passed, assign tip.labs to be names of phen.leaves:
          names(phen.leaves) <- tree$tip.label
        }
      }

      ## CHECK THAT MIN GRP.SIZE >= THRESHOLD ##
      if(!is.null(grp.min)){
        tab <- table(phen.leaves)
        grp.thresh <- (tree$Nnode+1)*grp.min
        if(min(tab) < grp.thresh){
          toRepeat <- TRUE
        }else{
          toRepeat <- FALSE
        }
      }else{
        toRepeat <- FALSE
      }

    } # end WHILE LOOP #########


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


  ## Assign names to phen.nodes ##
  ## ... as c(tree$tip.labs, tree$node.labs) in original order:
  ## If checks fail, individuals in phen.leaves will be named 1:N, node.1:node.Ninternal

  ## Assign terminal nodes names same as phen.leaves (ideally = tree$tip.label)
  noms.term <- names(phen.leaves)

  ## Assign internal nodes either tree$node.label or node.1:node.N:
  int.inds <- c((length(phen.leaves)+1):length(phen.nodes))
  noms.int <- paste("node", int.inds, sep=".")
  ## Check that node.labs are not NULL:
  if(!is.null(tree$node.label)){
    ## Check that node.labs is of correct length:
    if(length(tree$node.label) == length(int.inds)){
      noms.int <- tree$node.label
    }
  }

  ## Assign these names to phen.nodes:
  names(phen.nodes) <- c(noms.term, noms.int)



  ## make output list
  phen.list <- list(phen.leaves, phen.nodes, phen.branch, phen.loci)
  names(phen.list) <- c("phen", "phen.nodes", "phen.edges", "phen.loci")

  return(phen.list)
} # end phen.sim
