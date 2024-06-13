
###############
## plot_phen ##
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

#' Plot the states of a phenotype or genotype along a phylogenetic tree.
#'
#' This function is designed to visualise the reconstructed ancestral states of a variable along a phylogenetic tree.
#' It uses colour to represent the states of terminal and internal nodes (if available),
#' indicating changes between states by grey branches (except in the case of truly continuous variables).
#'
#'
#' @param tree A phylo object.
#' @param phen.nodes A vector containing the phenotypic state of either
#' (i) only terminal nodes in tree or
#' (ii) all nodes, terminal and internal in tree.
#' @param snp.nodes An optional vector containing the states of
#' a second variable (e.g., a genotypic variable) for either
#' the terminal nodes or all nodes in the tree.
#' @param plot A logical specifying whether to display a plot
#' of the inputted phylogenetic tree with edges coloured to show the
#' simulated phenotypic substitution process.
#' @param RTL A logical variable indicating whether to plot the
#' first or only tree from right to left (TRUE),
#' or left to right (FALSE, the default).
#' @param LTR.snp A logical variable indicating whether to plot the
#' optional second tree from left to right (TRUE),
#' or right to left (FALSE, the default).
#' @param main.title Either NULL or a character vector specifying a main title for the plot.
#' @param align.tip.label A logical indicating whether to align tip labels with each other (TRUE) or
#' to place tip labels at terminal nodes (FALSE, the default).
#' @param show.axis A logical indicating whether to add an axis showing the scale of branch lengths
#' at the foot of the plot with \code{axisPhylo} (TRUE, the default) or not (FALSE).
#'
#'
#' @details Ancestral states must be inferred in advance, for example, using function \code{asr}.
#' States are then shown in the colour of terminal node labels and the colour of the edges of the tree.
#' If only terminal states are available, these can be plotted along the tips of the tree.
#' If desired, a second variable, for example, a particular SNP or genetic locus, can be shown along
#' a second phylogeny. In this case, the second variable will be shown on a toplogically identical tree,
#' which will be plotted from right to left, mirroring the first tree along the vertical axis of the plotting window.
#' The \code{RTL} and \code{LTR.snp} arguments can be used to change the
#' orientation/direction of the first and/or second tree.
#'
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @examples
#'
#' ## Example 1 ##
#' \dontrun{
#' ## load phylogenetic and phenotypic data:
#' data(tree)
#' data(phen)
#'
#' ## reconstruct phenotypic ancestral states:
#' phen.rec <- asr(var=phen, tree=tree, type="parsimony", method="discrete")
#'
#' ## plot phenotype along tree:
#' plot_phen(tree, phen.nodes=phen.rec)
#' }
#'
#'
#' ## Example 2 ##
#' \dontrun{
#' ## load phylogenetic and phenotypic data:
#' data(tree)
#' data(phen)
#'
#' ## load genotypic data:
#' data(snps)
#'
#' ## reconstruct phenotypic ancestral states:
#' phen.rec <- asr(var=phen, tree=tree, type="parsimony", method="discrete")
#'
#' ## reconstruct genotypic ancestral states:
#' snps.rec <- asr(var=snps, tree=tree, type="parsimony", method="discrete")
#'
#' ## plot both the phenotype and a genotype along tree:
#' plot_phen(tree, phen.nodes=phen.rec, snp.nodes=snps.rec[,1])
#' }
#'
#' @import adegenet
#' @rawNamespace import(ape, except = zoom)
#' @importFrom Hmisc all.is.numeric
#' @importFrom graphics mtext
#' @export

########################################################################

plot_phen <- function(tree, phen.nodes, snp.nodes=NULL, plot=TRUE, RTL=FALSE, LTR.snp=FALSE,
                      main.title = NULL, align.tip.label=FALSE, show.axis=TRUE, ...){


  ## HANDLE TREE: ##
  ## Note: plot_phen expects output from phen.sim as input, and phen.sim works w pruningwise trees...
  tree <- reorder.phylo(tree, order="pruningwise")

  ## get number of terminal nodes
  if(!is.null(tree$tip.label)){
    n.ind <- length(tree$tip.label)
  }else{
    n.ind <- tree$Nnode+1
  }

  ## REORDER SNPS/PHEN to match TREE LABELS:
  if(!is.null(names(phen.nodes))){
    if(length(phen.nodes) == max(tree$edge)){
      if(!is.null(tree$node.label) && all(names(phen.nodes) %in% c(tree$tip.label, tree$node.label))){
        if(length(unique(c(tree$tip.label, tree$node.label))) < length(c(tree$tip.label, tree$node.label))){
          warning("Unable to rearrange phen.nodes to match tree labels
                  because not all c(tree$tip.label, tree$node.label) are unique.
                  Order may be incorrect.")
        }else{
         if(length(unique(names(phen.nodes))) < length(names(phen.nodes))){
           warning("Unable to rearrange phen.nodes to match tree labels
                  because not all names(phen.nodes) are unique.
                  Order may be incorrect.")
         }else{
           ord <- match(c(tree$tip.label, tree$node.label), names(phen.nodes))
           phen.nodes <- phen.nodes[ord]
         }
        }
      }else{
        if(!is.null(tree$tip.label) && all(tree$tip.label %in% names(phen.nodes))){
          ord <- match(tree$tip.label, names(phen.nodes))
          ## If no node.labels, but order same for tip inds, assume order ok for node inds;
          ## else, warn:
          if(!identical(ord, c(1:length(tree$tip.label)))){
            warning("Unable to rearrange phen.nodes such that names(phen.nodes)
              match c(tree$tip.label, tree$node.label). Order may be incorrect.")
          }
        }else{
          warning("Unable to rearrange phen.nodes such that names(phen.nodes)
              match c(tree$tip.label, tree$node.label). Order may be incorrect.")
        }
      }
    }
    if(length(phen.nodes) == n.ind){
      if(!is.null(tree$tip.label) && all(tree$tip.label %in% names(phen.nodes))){
        ord <- match(tree$tip.label, names(phen.nodes))
        phen.nodes <- phen.nodes[ord]
      }else{
        warning("Unable to rearrange phen.nodes such that names(phen.nodes)
            match tree$tip.label. Order may be incorrect.")
      }
    }
  }


  ## PLOT MARGINS: ##
  mar.ori <- par()$mar
  # mar.ori <- c(5,4,4,1)+0.1
  if(RTL==FALSE){
    mar.new <- c(1.5,0.5,0.5,0.5)
    # mar.new <- c(4, 2, 4, 0.5) + 0.1
  }else{
    mar.new <- c(1.5,0.5,0.5,0.5)
    # mar.new <- c(5, 4, 4, 2) + 0.1
  }
  ## Set plot margins:
  par(mar=mar.new)

  if(is.null(main.title)) main.title <- FALSE

  edgeLabCol <- edgeCol <- nodeCol <- internalNodeCol <- leafCol <- NULL

  #############################################################################
  ######################## PLOT phylo with PHEN ###############################
  #############################################################################

  ## SIDE-BY-SIDE PLOTS??
  if(!is.null(phen.nodes) && !is.null(snp.nodes)){
    par(mfrow=c(1,2))
    main.title2 <- main.title
    main.title <- FALSE
  }



  nodeCol <- "grey"


  ## get COLOR for NODES
  # nodeCol <- "grey"
  if(Hmisc::all.is.numeric(phen.nodes[!is.na(phen.nodes)])){
    var <- as.numeric(as.character(phen.nodes))
    levs <- sort(unique(var[!is.na(var)]))
  }else{
    var <- as.character(phen.nodes)
    levs <- unique(var[!is.na(var)])
  }

  if(length(levs) == 2){
    ## binary:
    # myCol <- c("red", "blue")
    myCol <- c("blue", "red")
    nodeCol <- var
    ## for loop
    for(i in 1:length(levs)){
      nodeCol <- replace(nodeCol, which(nodeCol == levs[i]), myCol[i])
    } # end for loop
  }else{

    if(is.numeric(var)){
      ## numeric:
      myCol <- num2col(var, col.pal = seasun, na.col = NA)
      # myCol <- num2col(var, col.pal = grDevices::rgb(base,0,1-base), na.col = NA)
      nodeCol <- myCol
    }else{
      ## categorical...
      myCol <- funky(length(levs))
      nodeCol <- var
      ## for loop
      for(i in 1:length(levs)){
        nodeCol <- replace(nodeCol, which(nodeCol == levs[i]), myCol[i])
      } # end for loop
    }
  }


  nodeCol <- as.vector(unlist(nodeCol))
  ## get COLOR for LEAVES ONLY
  leafCol <- nodeCol[1:n.ind]


  ## get edgeCol if phen provided is for internal nodes:
  if(length(phen.nodes) == max(tree$edge)){

    ## get COLOR of INTERNAL nodes ONLY
    internalNodeCol <- nodeCol[(n.ind+1):length(nodeCol)]

    ## get COLOR for EDGES
    edgeCol <- rep("black", nrow(tree$edge))
    for(i in 1:nrow(tree$edge)){
      edgeCol[i] <- nodeCol[tree$edge[i,2]]
      if(is.na(nodeCol[tree$edge[i,1]]) | is.na(nodeCol[tree$edge[i,2]])){
        edgeCol[i] <- "grey"
      }else{
        ## No grey if truly continuous...
        if(length(levs) < length(tree$tip.label)/10){
          if(nodeCol[tree$edge[i,1]] != nodeCol[tree$edge[i,2]]) edgeCol[i] <- "grey"
        }
      }
    }
    edgeLabCol <- edgeCol

    ## replace NAs w grey:
    if(any(is.na(var))){
      nodeCol[which(is.na(var))] <- "grey"
    }

  }else{
    edgeCol <- "black"
    nodeCol <- internalNodeCol <- edgeCol
  }


    ###############
    ## plot TREE ##
    ###############
    if(plot==TRUE){

      if(RTL == FALSE){
        if(align.tip.label == FALSE){
          # suppressWarnings(plot(tree, show.tip=T, tip.col="transparent", cex=0.5, adj=c(-0.5, 0), edge.width=3, edge.color=edgeCol, ...))
          # Error in xx.tmp + lox :
          #   (converted from warning) longer object length is not a multiple of shorter object length ## from adj w 2 values...
          suppressWarnings(plot(tree, show.tip=T, tip.col="transparent", cex=0.5, adj=-0.5, edge.width=3, edge.color=edgeCol, ...))
          tiplabels(text=tree$tip.label, cex=0.5, adj=c(-0.5, 0), col=leafCol, frame="none")
        }else{
          plot(tree, show.tip=T, tip.col=leafCol, edge.width=3, edge.color=edgeCol, align.tip.label=T, cex=0.5, ...)
        }

      }else{
        plot(tree, show.tip=T, tip.col="white", cex=0.5, adj=c(-7.5), edge.width=3, edge.color=edgeCol, direction = "leftwards", ...)
        tiplabels(text=tree$tip.label, cex=0.5, adj=c(1.5, 0), col=leafCol, frame="none")
      }

      ## Add title?
      if(main.title == TRUE){
        # if(is.ultrametric(tree)) title("Coalescent tree w/ phenotypic changes")
        title("Tree with phenotypic changes", cex.main=1, line=-0.5)
      }else{
        if(!is.null(main.title)) if(main.title != FALSE) title(main.title, cex.main=1, line=-0.5)
      }

      ## Add axis?
      if(show.axis == TRUE){
        mar.new2 <- c(2,0.5,0,0.5)
        par(mar=mar.new2)
        axisPhylo(cex.axis=0.7)
      }

  } # end plot = TRUE

  #####################################################################################################

  if(!is.null(snp.nodes)){

    phen.nodes <- snp.nodes # for convenience

    ## Check if > 1 column of snps(.rec) has been provided:
    if(is.matrix(phen.nodes)){
      warning("More than one column has been provided to snp.nodes. Only plotting snp.nodes[,1].")
      phen.nodes <- phen.nodes[,1]
    }

    ## REORDER SNPS/PHEN to match TREE LABELS:
    if(!is.null(names(phen.nodes))){
      if(length(phen.nodes) == max(tree$edge)){
        if(!is.null(tree$node.label) && all(names(phen.nodes) %in% c(tree$tip.label, tree$node.label))){
          ord <- match(c(tree$tip.label, tree$node.label), names(phen.nodes))
          phen.nodes <- phen.nodes[ord]
        }else{
          if(!is.null(tree$tip.label) && all(tree$tip.label %in% names(phen.nodes))){
            ord <- match(tree$tip.label, names(phen.nodes))
            ## If no node.labels, but order same for tip inds, assume order ok for node inds;
            ## else, warn:
            if(!identical(ord, c(1:length(tree$tip.label)))){
              warning("Unable to rearrange snp.nodes such that names(snp.nodes)
              match c(tree$tip.label, tree$node.label). Order may be incorrect.")
            }
          }else{
              warning("Unable to rearrange snp.nodes such that names(snp.nodes)
                match c(tree$tip.label, tree$node.label). Order may be incorrect.")
          }
        }
      }
      if(length(phen.nodes) == n.ind){
        if(!is.null(tree$tip.label) && all(tree$tip.label %in% names(phen.nodes))){
          ord <- match(tree$tip.label, names(phen.nodes))
          phen.nodes <- phen.nodes[ord]
        }else{
          warning("Unable to rearrange snp.nodes such that names(snp.nodes)
            match tree$tip.label. Order may be incorrect.")
        }
      }
    }

    nodeCol <- "grey"


    ## get COLOR for NODES
    # nodeCol <- "grey"
    if(Hmisc::all.is.numeric(phen.nodes[!is.na(phen.nodes)])){
      var <- as.numeric(as.character(phen.nodes))
      levs <- sort(unique(var[!is.na(var)]))
    }else{
      var <- as.character(phen.nodes)
      levs <- unique(var[!is.na(var)])
    }

    if(length(levs) == 2){
      ## binary:
      # myCol <- c("red", "blue")
      myCol <- c("blue", "red")
      nodeCol <- var
      ## for loop
      for(i in 1:length(levs)){
        nodeCol <- replace(nodeCol, which(nodeCol == levs[i]), myCol[i])
      } # end for loop
    }else{

      if(is.numeric(var)){
        ## numeric:
        myCol <- num2col(var, col.pal = seasun, na.col = NA)
        nodeCol <- myCol
      }else{
        ## categorical...
        myCol <- funky(length(levs))
        nodeCol <- var
        ## for loop
        for(i in 1:length(levs)){
          nodeCol <- replace(nodeCol, which(nodeCol == levs[i]), myCol[i])
        } # end for loop
      }
    }


    nodeCol <- as.vector(unlist(nodeCol))
    ## get COLOR for LEAVES ONLY
    leafCol <- nodeCol[1:n.ind]


    ## get edgeCol if phen provided is for internal nodes:
    if(length(phen.nodes) == max(tree$edge)){

      ## get COLOR of INTERNAL nodes ONLY
      internalNodeCol <- nodeCol[(n.ind+1):length(nodeCol)]

      ## get COLOR for EDGES
      edgeCol <- rep("black", nrow(tree$edge))
      for(i in 1:nrow(tree$edge)){
        edgeCol[i] <- nodeCol[tree$edge[i,2]]
        if(is.na(nodeCol[tree$edge[i,1]]) | is.na(nodeCol[tree$edge[i,2]])){
          edgeCol[i] <- "grey"
        }else{
          ## No grey if truly continuous...
          if(length(levs) < length(tree$tip.label)/10){
            if(nodeCol[tree$edge[i,1]] != nodeCol[tree$edge[i,2]]) edgeCol[i] <- "grey"
          }
        }
      }
      edgeLabCol <- edgeCol

      ## replace NAs w grey:
      if(any(is.na(var))){
        nodeCol[which(is.na(var))] <- "grey"
      }

    }else{
      edgeCol <- "black"
      nodeCol <- internalNodeCol <- edgeCol
    }

      ###############
      ## plot TREE ##
      ###############
      if(plot==TRUE){

        ## Set plot margins:
        par(mar=mar.new)


        if(LTR.snp == FALSE){
          if(align.tip.label == FALSE){
            suppressWarnings(plot(tree, show.tip=T, tip.col="transparent", cex=0.5, adj=c(-0.5, 0), edge.width=3, edge.color=edgeCol, direction = "leftwards", ...))
            tiplabels(text=tree$tip.label, cex=0.5, adj=c(1.5, 0), col=leafCol, frame="none")
          }else{
            plot(tree, show.tip=T, tip.col=leafCol, edge.width=3, edge.color=edgeCol, align.tip.label=T, cex=0.5, direction = "leftwards", ...)
          }

        }else{
          plot(tree, show.tip=T, tip.col="white", cex=0.5, adj=c(-7.5), edge.width=3, edge.color=edgeCol, ...)
          tiplabels(text=tree$tip.label, cex=0.5, adj=c(-0.5, 0), col=leafCol, frame="none")
        }

        # # plot(tree, show.tip=T, tip.col=leafCol, edge.width=3, edge.color=edgeCol, align.tip.label=T, cex=0.5, direction="leftwards")
        # plot(tree, show.tip=FALSE, edge.width=3, edge.color=edgeCol), direction="leftwards", ...) # edgeCol
        # tiplabels(text=tree$tip.label, cex=0.6, adj=c(1.5, 0), col=leafCol, frame="none")

        ## Add title?
        if(main.title == TRUE){
          # if(is.ultrametric(tree)) title("Coalescent tree w/ phenotypic changes")
          title("Tree with phenotypic changes", cex.main=1, line=-0.5)
        }else{
          if(!is.null(main.title)) if(main.title != FALSE) title(main.title, cex.main=1, line=-0.5)
        }


        ## Add axis?
        if(show.axis == TRUE){
          mar.new2 <- c(2,0.5,0,0.5)
          par(mar=mar.new2)
          axisPhylo(cex.axis=0.7)
        }

      } # end plot = TRUE

  } # end snp.nodes

  #####################################################################################################

  ## One common main title:
  if(!is.null(phen.nodes) && !is.null(snp.nodes)){
    if(!is.null(main.title2)) if(main.title2 != FALSE){
      graphics::mtext(main.title2, cex=1, font=2, side=3, line=-1.1, outer=TRUE)
    }
    ## return plot par settings to 1 plot per window..
    par(mfrow=c(1,1))
  }


  ## Reset plot margins: ##
  par(mar=mar.ori)

  ## get all info relevant to plotting phenotype with colored phylo:
  phen.plot.colors <- list(edgeLabCol, edgeCol, nodeCol, internalNodeCol, leafCol)
  names(phen.plot.colors) <- c("edge.labels", "edges", "all.nodes", "internal.nodes", "tip.labels")

  invisible(phen.plot.colors)

} # end plot_phen



#################################
##  ENABLE ALTERNATE FN NAME:  ##
#################################
# plot.phen <- function(tree, phen.nodes, ...){
#   return(plot_phen(tree, phen.nodes,  ...))
# } # end plot.phen



#


