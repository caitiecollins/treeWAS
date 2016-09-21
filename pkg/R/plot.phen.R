
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
#' @param tree A phylo object.
#' @param phen.nodes A vector containing the phenotypic state of either
#' (i) only terminal nodes in tree or
#' (ii) all nodes, terminal and internal in tree.
#' @param plot A logical specifying whether to display a plot
#' of the inputted phylogenetic tree with edges coloured to show the
#' simulated phenotypic substitution process.
#'
#'
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
#' @import adegenet phangorn

########################################################################

###############################
## SIMPLE PLOT FOR CHECKING: ##
###############################

# plot(tree, show.tip=FALSE, edge.width=2)
# tiplabels(text=tree$tip.label, cex=0.6, adj=c(-0.5, 0), col="black", frame="none")
# edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."),
#            cex=0.66, font=2, frame="none", col="blue", adj=c(1,1))
# ## coaltree??
# # nodelabels(text=rev(unique(tree$edge[,1])), cex=0.6, bg=transp("yellow", 0.7), col="black", frame="circle") # frame="none", col="red"
# ## rtree??
# nodelabels(text=unique(tree$edge[,1]), cex=0.6, bg=transp("yellow", 0.7), col="black", frame="circle") # frame="none", col="red"
# axisPhylo()



plot.phen <- function(tree, phen.nodes, snp.nodes=NULL, plot=TRUE, RTL=FALSE, main.title = TRUE, ...){

  # require(phangorn)
  # require(adegenet)

  #############################################################################
  ######################## PLOT phylo with PHEN ###############################
  #############################################################################

  ## SIDE-BY-SIDE PLOTS??
  if(!is.null(phen.nodes) & !is.null(snp.nodes)) par(mfrow=c(1,2))

  ## get number of terminal nodes
  n.ind <- tree$Nnode+1

  ## check if phen provided is for all nodes or only terminal nodes:
  if(length(phen.nodes) == (n.ind + tree$Nnode)){

    ## get COLOR for NODES
    nodeCol <- as.vector(phen.nodes)
    nodeCol <- as.character(nodeCol)
    # nodeCol <- replace(nodeCol, which(nodeCol == "B"), "red")
    # nodeCol <- replace(nodeCol, which(nodeCol == "A"), "blue")
    nodeCol <- replace(nodeCol, which(nodeCol %in% c("B", "1")), "red")
    nodeCol <- replace(nodeCol, which(nodeCol %in% c("A", "0")), "blue")
    nodeCol <- as.vector(unlist(nodeCol))
    ## get COLOR for LEAVES ONLY
    leafCol <- nodeCol[1:n.ind]
    ## get COLOR of INTERNAL nodes ONLY
    internalNodeCol <- nodeCol[(n.ind+1):length(nodeCol)]

    ## get COLOR for EDGES
    edgeCol <- rep("black", nrow(tree$edge))
    for(i in 1:nrow(tree$edge)){
      edgeCol[i] <- nodeCol[tree$edge[i,2]]
      if(nodeCol[tree$edge[i,1]] != nodeCol[tree$edge[i,2]]) edgeCol[i] <- "grey"
    }
    edgeLabCol <- edgeCol

    ###############
    ## plot TREE ##
    ###############
    if(plot==TRUE){
      if(n.ind <= 20){
        plot(tree, show.tip=FALSE, edge.width=2, edge.color=edgeCol, ...) # edgeCol
        ## Add title?
        if(main.title == TRUE){
          title("Coalescent tree w/ phenotypic changes")
        }else{
          if(!is.null(main.title)) title(main.title)
        }
        axisPhylo()
        edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."),
                   cex=0.5, font=2, bg=transp(edgeLabCol, 0.3), adj=c(1,1))
        tiplabels(text=tree$tip.label, cex=0.6, adj=c(-0.5, 0), bg=transp(leafCol, 0.3))
        nodelabels(text=rev(unique(tree$edge[,1])), cex=0.5, bg=transp(internalNodeCol, 0.3))

      }else{

        if(RTL == FALSE){
          plot(tree, show.tip=FALSE, edge.width=2, edge.color=edgeCol, ...) # edgeCol
        }else{
          plot(tree, show.tip=FALSE, edge.width=2, edge.color=edgeCol, direction = "leftwards", ...) # edgeCol
        }
        ## Add title?
        if(main.title == TRUE){
          title("Coalescent tree w/ phenotypic changes")
        }else{
          if(!is.null(main.title)) title(main.title)
        }
        axisPhylo()
        # edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."),
                  # cex=0.5, font=2, bg=transp(edgeLabCol, 0.3), adj=c(1,1))
        if(RTL == FALSE){
          tiplabels(text=tree$tip.label, cex=0.6, adj=c(-0.5, 0), col=leafCol, frame="none")
        }else{
          ## if direction="leftwards"
          tiplabels(text=tree$tip.label, cex=0.6, adj=c(1.5,0), col=leafCol, frame="none")
        }
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

    if(!is.null(snp.nodes)){

      phen.nodes <- snp.nodes # for convenience

      ## check if phen provided is for all nodes or only terminal nodes:
      if(length(phen.nodes) == (n.ind + tree$Nnode)){

        ## get COLOR for NODES
        nodeCol <- as.vector(phen.nodes)
        nodeCol <- as.character(nodeCol)
        # nodeCol <- replace(nodeCol, which(nodeCol == "B"), "red")
        # nodeCol <- replace(nodeCol, which(nodeCol == "A"), "blue")
        nodeCol <- replace(nodeCol, which(nodeCol %in% c("B", "1")), "red")
        nodeCol <- replace(nodeCol, which(nodeCol %in% c("A", "0")), "blue")
        nodeCol <- as.vector(unlist(nodeCol))
        ## get COLOR for LEAVES ONLY
        leafCol <- nodeCol[1:n.ind]
        ## get COLOR of INTERNAL nodes ONLY
        internalNodeCol <- nodeCol[(n.ind+1):length(nodeCol)]

        ## get COLOR for EDGES
        edgeCol <- rep("black", nrow(tree$edge))
        for(i in 1:nrow(tree$edge)){
          edgeCol[i] <- nodeCol[tree$edge[i,2]]
          if(nodeCol[tree$edge[i,1]] != nodeCol[tree$edge[i,2]]) edgeCol[i] <- "grey"
        }
        edgeLabCol <- edgeCol


        ###############
        ## plot TREE ##
        ###############
        if(plot==TRUE){
          if(n.ind <= 20){
            plot(tree, show.tip=FALSE, edge.width=2, edge.color=edgeCol, direction="leftwards") # edgeCol
            ## Add title?
            if(main.title == TRUE){
              title("Coalescent tree w/ phenotypic changes")
            }else{
              if(!is.null(main.title)) title(main.title)
            }
            axisPhylo()
            edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."),
                       cex=0.5, font=2, bg=transp(edgeLabCol, 0.3), adj=c(1,1))
            tiplabels(text=tree$tip.label, cex=0.6, adj=c(-0.5, 0), bg=transp(leafCol, 0.3))
            nodelabels(text=rev(unique(tree$edge[,1])), cex=0.5, bg=transp(internalNodeCol, 0.3))

          }else{

            plot(tree, show.tip=FALSE, edge.width=2, edge.color=edgeCol, direction="leftwards") # edgeCol
            ## Add title?
            if(main.title == TRUE){
              title("Coalescent tree w/ phenotypic changes")
            }else{
              if(!is.null(main.title)) title(main.title)
            }
            axisPhylo()
            # edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."),
            # cex=0.5, font=2, bg=transp(edgeLabCol, 0.3), adj=c(1,1))

            tiplabels(text=tree$tip.label, cex=0.6, adj=c(1.5, 0), col=leafCol, frame="none") # adj=c(1.5,0) # if direction="leftwards"

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
      } # end check for length snp.nodes
    } # end snp.nodes

    ## return plot par settings to 1 plot per window..
    if(!is.null(phen.nodes) & !is.null(snp.nodes)) par(mfrow=c(1,1))

  }else{ # end phen for all nodes

    #######################################################################

    ####################################################################
    ## If the user has PROVIDED a phenotype (for terminal nodes only) ##
    ####################################################################

    phen <- phen.nodes

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
      edgeCol <- edgeLabCol <- "black" ## for NOW...
      ########
      ## TO DO-- color ~ terminal edges red/blue same as phen of terminal node...
      #### ... UNTIL two edge colors meet at any internal node (then just black edges to root)

      ###############
      ## plot TREE ##
      ###############
      if(plot==TRUE){
        plot(tree, show.tip=FALSE, edge.width=2, edge.color=edgeCol, ...)
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



# i <- 7
# plot.phen(tree, phen.nodes=phen.nodes, snp.nodes=snps.assoc.nodes[,i])
# length(which(phen.nodes[1:100] == snps.assoc.nodes[,i][1:100]))

# i <- 6
# plot.phen(tree, phen.nodes=phen.rec, snp.nodes=snps.rec[,snps.assoc[i]])
# length(which(phen.rec[1:100] == snps.rec[,snps.assoc[i]][1:100]))

###############################
## plotting reconstructions: ##
###############################

# ## PLUS -- enable plotting of boxes of extra variables along tips..
#
# ## get rec (from simTest):
#
# tree <- out$tree[[1]]
#
# rec <- out$res[[1]]$dat$snps.rec
#
# snps.assoc <- out$performance[[1]]$snps.assoc
#
# var.rec <- rec[,snps.assoc[6]]
# var.rec <- round(var.rec)
#
# var.rec <- replace(var.rec, which(var.rec == 0), "A")
# var.rec <- replace(var.rec, which(var.rec == 1), "B")
#
# plot.phen(tree, phen.nodes = var.rec)
#
# ## compare to phen:
# phen.rec <- as.character(out$res[[1]]$dat$phen.rec)
# phen.rec <- replace(phen.rec, which(phen.rec == 0), "A")
# phen.rec <- replace(phen.rec, which(phen.rec == 1), "B")
#
# plot.phen(tree, phen.nodes=phen.rec)
#
# ## side-by-side?
# par(mfrow=c(1,2))
# ## phen:
# plot.phen(tree, phen.nodes=phen.rec)
# title("set1_31: phen.rec (left) vs.", line=-0.5)
# ## snp:
# plot.phen(tree, phen.nodes = var.rec, direction="leftwards") # tip adj=c(1.5,0)
# title("snps.rec 304 (right)", line=-0.5)








#

