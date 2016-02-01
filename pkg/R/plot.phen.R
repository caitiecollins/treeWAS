
###############
## plot.phen ##
###############

############
## TO DO: ##
############

## (1) Code to reconstruct phen of internal nodes (and branches) if
## only phen for terminal is provided with tree.


########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param phen.nodes A vector containing the phenotypic state of either
#' (i) only terminal nodes in \texttt{tree} or
#' (ii) all nodes, terminal and internal in \texttt{tree}.
#' @param phen.branch An optional vector containing the phenotypic
#' @param tree A phylo object.
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
#' @import adegenet ape phangorn

########################################################################


plot.phen <- function(tree, phen.nodes, phen.edges=NULL, plot=TRUE){

  require(phangorn)

  phen.branch <- phen.edges

  #############################################################################
  ######################## PLOT phylo with PHEN ###############################
  #############################################################################

  ## get number of terminal nodes
  n.ind <- tree$Nnode+1

  ## check if phen provided is for all nodes or only terminal nodes:
  if(length(phen.nodes == (n.ind + tree$Nnode))){

    ## if phen for edges is provided,
    ## get colours and plot:
    if(!is.null(phen.branch)){

      ## get COLOR for NODES
      nodeCol <- as.vector(phen.nodes)
      nodeCol <- replace(nodeCol, which(nodeCol == "B"), "red")
      nodeCol <- replace(nodeCol, which(nodeCol == "A"), "blue")
      nodeCol <- as.vector(unlist(nodeCol))
      ## get COLOR for LEAVES ONLY
      leafCol <- nodeCol[1:n.ind]
      ## get COLOR of INTERNAL nodes ONLY
      internalNodeCol <- nodeCol[(n.ind+1):length(nodeCol)]

      ## get COLOR for EDGES
      edgeCol <- phen.branch
      ## FOR NOW--all edges w > 1 phen (either c("A", "B") OR c("B", "A")) --> green...
      l_edgeCol <- which(sapply(c(1:length(edgeCol)),
                                function(e) length(edgeCol[[e]])) == 2)
      edgeCol <- replace(edgeCol, l_edgeCol, "green")
      edgeCol <- replace(edgeCol, which(edgeCol == "A"), "blue")
      edgeCol <- replace(edgeCol, which(edgeCol == "B"), "red")
      edgeCol <- as.vector(unlist(edgeCol))
      ## get COLOR for EDGE LABELS
      edgeLabCol <- rep("green", nrow(tree$edge))
      edgeLabCol <- replace(edgeLabCol, phen.loci, "green") # "mediumorchid1"

      ###############
      ## plot TREE ##
      ###############
      if(plot==TRUE){
        if(n.ind <= 20){
          plot(tree, show.tip=FALSE, edge.width=2, edge.color=edgeCol) # edgeCol
          title("Coalescent tree w/ phenotypic changes")
          axisPhylo()
          edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."),
                     cex=0.5, font=2, bg=transp(edgeLabCol, 0.3), adj=c(1,1))
          tiplabels(text=tree$tip.label, cex=0.6, adj=c(-0.5, 0), bg=transp(leafCol, 0.3))
          nodelabels(text=rev(unique(tree$edge[,1])), cex=0.5, bg=transp(internalNodeCol, 0.3))
        }else{
          plot(tree, show.tip=FALSE, edge.width=2, edge.color=edgeCol) # edgeCol
          title("Coalescent tree w/ phenotypic changes")
          axisPhylo()
          #edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."),
          #           cex=0.5, font=2, bg=transp(edgeLabCol, 0.3), adj=c(1,1))
          tiplabels(text=tree$tip.label, cex=0.6, adj=c(-0.5, 0), col=leafCol, frame="none")
          ## make sure this isn't backward...:
          #nodelabels(text=rev(unique(tree$edge[,1])), cex=0.5, bg=transp(internalNodeCol, 0.3))
          ## should be numbered s.t. the root node is n.term+1
          ## RECALL: terminal nodes are numbered 1:n.ind from bottom to top of plot of tree;
          ## edges are numbered 1:nrow(edges) by following the lowest trace on the plot
          ## (starting from the root down to the lowermost tips);
          ## thus, internal nodes are numbered (n.ind+1):(n.ind+(n.ind-1)),
          ## from root to the top-most internal node to be connected (ie. the highest in the plot)
        }
      } # end plot = TRUE
    } # end !is.null(phen.branch)
  }else{ # end phen for all nodes

    #######################################################################


    ####################################################################
    ## If the user has PROVIDED a phenotype (for terminal nodes only) ##
    ####################################################################
    if(!is.null(phen)){
      # phen <- as.factor(sample(c("A", "B", "C", "D"), 100, replace=TRUE))

      ## get COLOR for LEAVES
      n.levels <- length(levels(as.factor(phen)))
      ## if we have only 2 levels, use standard red and blue
      if(n.levels==2){
        leafCol <- c("red", "blue")
      }else{
        ## but if we have n levels (n != 2), use funky palette
        leafCol <- get("funky")(n.levels)
      }

      scheme <- as.numeric(phen)
      leafCol <- leafCol[scheme]

      ## get COLOR for EDGES
      edgeCol <- "black" ## for NOW...
      ########
      ## TO DO-- color ~ terminal edges red/blue same as phen of terminal node...
      #### ... UNTIL two edge colors meet at any internal node (then just black edges to root)

      ###############
      ## plot TREE ##
      ###############
      if(plot==TRUE){
        plot(tree, show.tip=FALSE, edge.width=2, edge.color=edgeCol)
        title("Coalescent tree w/ phenotypes of leaves")
        axisPhylo()
        tiplabels(text=tree$tip.label, cex=0.6, adj=c(-0.5, 0), bg=transp(leafCol, 0.3))
      } # end plot = TRUE
    } # end if(!is.null(phen)) ## ie. PROVIDED phenotype & plotting
  } # end if phen for terminal nodes only

  ## get all info relevant to plotting phenotype with colored phylo:
  phen.plot.colors <- list(edgeLabCol, edgeCol, nodeCol, internalNodeCol, leafCol)
  names(phen.plot.colors) <- c("edge.labels", "edges", "all.nodes", "internal.nodes", "tip.labels")

  return(phen.plot.colors)

} # end plot.phen


