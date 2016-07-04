


########################
## get.ancestral.pars ##
########################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Ancestral sequence reconstruction via parsimony
#'
#' A wrapper for the \code{ancestral.pars} function from \emph{ape}. Can perform
#' parsimonious ASR for variables in matrix or vector form.
#'
#' @param var A matrix or vector containing a variable whose state at ancestral nodes we want to infer.
#' @param tree A phylo object containing a phylogenetic tree whose tips contain the same individuals as are
#' in the elements of \code{var}, if \code{var} is a vector,
#' or in the rows of \code{var}, if \code{var} is a matrix.
#'
#' @details Note that the (row)names of \code{var} should match the tip.labels of \code{tree}.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @import phangorn ape
#'
#' @export
#'

########################################################################

get.ancestral.pars <- function(var, tree){

  edges <- tree$edge

  #############################
  ## RUN PARSIMONY on MATRIX ##
  #############################

  if(is.matrix(var)){

    snps <- var

    ###############################
    ## get unique SNPs patterns: ##
    ###############################
    temp <- get.unique.matrix(snps, MARGIN=2)
    snps.unique <- temp$unique.data
    index <- temp$index

    if(ncol(snps.unique) == ncol(snps)){
      all.unique <- TRUE
    }else{
      all.unique <- FALSE
    }

    ## work w only unique snps:
    snps.ori <- snps
    snps <- snps.unique

    ##########################################
    ## get SNP states of all internal nodes ##
    ##########################################

    #############
    ## ACCTRAN ##
    #############

    ## pace == ancestral.pars ## (ape)
    ## parsimony version of ace ## based on the fitch algorithm.
    ## ACCTRAN = "ACCelerated TRANsformation"

    ## as.phyDat requires row/col names:
    if(is.null(row.names(snps))){
      if(!is.null(tree$tip.label)) row.names(snps) <- tree$tip.label
    }
    if(is.null(colnames(snps))) colnames(snps) <- c(1:ncol(snps))
    ## get levels (ie. 0, 1)
    snps.levels <- unique(as.vector(snps))
    ## returns only unique patterns...
    snps.phyDat <- as.phyDat(as.matrix(snps),
                             type="USER", levels=snps.levels)
    ## get index of all original snps columns to map to unique pattern
    index <- attr(snps.phyDat, "index")

    ## pace == ancestral.pars
    pa.ACCTRAN <- pace(tree, snps.phyDat, type="ACCTRAN")

    ## NOTE: pace  --> diff resuls w MPR vs. ACCTRAN
    # pa.MPR <- pace(tree, snps.phyDat, type="MPR")
    #diffs <- sapply(c(1:length(pa.ACCTRAN)), function(e) identical(pa.MPR[[e]], pa.ACCTRAN[[e]]))

    ###########################################
    ## convert reconstruction back to snps.. ##
    ###########################################
    ## each of the n.ind elements of pa is a matrix w n.snps rows and either:
    ## 2 columns, for the 2 binary SNP states, or
    ## 4 columns, each for the 4 nts possible (acgt)

    # rec <- pa.MPR
    rec <- pa.ACCTRAN

    snps.rec <- vector("list", length(rec))
    ## Handle terminal nodes
    for(i in 1:nrow(snps)){
      snps.rec[[which(row.names(snps) == tree$tip.label[i])]] <- rec[[i]][,2]
    }
    ## Handle internal nodes
    for(i in (nrow(snps)+1):length(rec)){
      snps.rec[[i]] <- rec[[i]][,2]
    }
    ## bind rows together
    snps.rec <- do.call("rbind", snps.rec)
    ## assign rownames for all terminal and internal nodes
    rownames(snps.rec) <- c(rownames(snps), c((nrow(snps)+1):((nrow(snps)*2)-1)))
    colnames(snps.rec) <- c(1:length(snps.phyDat[[1]]))


    ###########################################
    ## get LOCATIONS (branches) of snps subs ##
    ###########################################
    subs.edges <- rep(list(NULL), ncol(snps.rec))
    for(i in 1:ncol(snps.rec)){
      snp <- snps.rec[, i]
      subs.logical <- sapply(c(1:nrow(edges)),
                             function(e)
                               snp[edges[e,1]]
                             ==
                               snp[edges[e,2]])
      ## get indices of all edges containing a substitution
      subs.total <- which(subs.logical == FALSE)
      ## get df of states of ancestor and descendants nodes on these edges
      df <- data.frame(snp[edges[subs.total,1]], snp[edges[subs.total,2]])
      names(df) <- c("anc", "dec")
      ## get indices of all edges w a positive sub (0 --> 1)
      subs.pos <- subs.total[which(df$anc==0)]
      ## get indices of all edges w a negative sub (1 --> 0)
      subs.neg <- subs.total[which(df$anc==1)]

      ## get output list
      subs.edges[[i]] <- rep(list(NULL), 3)
      names(subs.edges[[i]]) <- c("total", "pos", "neg")
      if(length(subs.total) > 0) subs.edges[[i]][["total"]] <- subs.total
      if(length(subs.pos) > 0) subs.edges[[i]][["pos"]] <- subs.pos
      if(length(subs.neg) > 0) subs.edges[[i]][["neg"]] <- subs.neg
    }

    ####################
    ## PLOT to CHECK? ##
    ####################
    ## for SNP1, does it identify the correct/reasonable branches?
    #     edgeCol <- rep("black", nrow(edges))
    #     edgeCol <- replace(edgeCol, subs.edges[[1]][["total"]], "green")
    #
    #     ## plot the i'th character's reconstruction on the tree:
    #     #require(adegenet)
    #     plotAnc(tree, pa.ACCTRAN, i=1,
    #             col=transp(c("red", "royalblue"), 0.75),
    #             cex.pie=0.1, pos=NULL,
    #             edge.color=edgeCol, edge.width=2, use.edge.length=FALSE, type="c")

    ################
    ## Get output ##
    ################

    ## get reconstruction for all original sites
    if(ncol(snps.ori) == ncol(snps.rec)){
      snps.rec.complete <- snps.rec
    }else{
      snps.rec.complete <- matrix(NA, nrow=nrow(snps.rec), ncol=ncol(snps.ori))
      for(i in 1:ncol(snps.rec)){
        snps.rec.complete[, which(index == i)] <- snps.rec[, i]
      }
      rownames(snps.rec.complete) <- rownames(snps.rec)
      colnames(snps.rec.complete) <- colnames(snps.ori)
    }

    ## get sub locations on branches for all original sites
    snps.subs.edges <- subs.edges
    if(ncol(snps.ori) == ncol(snps.rec)){
      snps.subs.edges.complete <- snps.subs.edges
    }else{
      snps.subs.edges.complete <- vector("list", ncol(snps.ori))
      for(i in 1:length(snps.subs.edges)){
        loc <- which(index == i)
        if(length(loc) == 1){
          snps.subs.edges.complete[[loc]] <- snps.subs.edges[[i]]
        }else{
          for(j in 1:length(loc)){
            snps.subs.edges.complete[[loc[j]]] <- snps.subs.edges[[i]]
          }
        }
      }
    }


    ## CHECK-- compare cost from fitch and pace: ##
    ## get n.subs per site by fitch:
    cost <- get.fitch.n.mts(snps, tree)
    ## get n.subs per site by pace:
    cost2 <- sapply(c(1:length(snps.subs.edges.complete)),
                    function(e) length(snps.subs.edges.complete[[e]][["total"]]))
    ## NOTE: cost2 differs somewhat noticeably from original fitch cost
    ## (ie. parsimony shifts distribution toward 1/reduces the weight of the upper tail...)
    ## WHY? Which should we use to get n.subs????????????????????????????????????????????????????????????

    ## Get final output list:
    var.rec <- snps.rec.complete
    subs.edges <- snps.subs.edges.complete

    out <- list("var.rec" = var.rec,
                "subs.edges" = subs.edges)

  }else{ # end matrix (snps) parsimony

    ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###


    #############################
    ## RUN PARSIMONY on VECTOR ##
    #############################

    ## Eg. get PHEN states of all internal nodes

    phen <- var

    #############
    ## ACCTRAN ##
    #############

    phen.ori <- phen
    phen <- as.numeric(as.character(phen))
    ## as.phyDat requires names...
    if(is.null(names(phen))){
      if(!is.null(tree$tip.label)) names(phen) <- tree$tip.label
    }

    ## get levels (ie. 0, 1)
    phen.levels <- unique(phen)
    phen.phyDat <- as.phyDat(as.matrix(phen),
                             type="USER", levels=phen.levels)
    ## pace == ancestral.pars
    rec <- phen.pa.ACCTRAN <- pace(tree, phen.phyDat, type="ACCTRAN")

    ## get reconstruction:
    phen.rec <- list()
    for(i in 1:length(rec)){
      phen.rec[[i]] <- rec[[i]][,2]
    }
    phen.rec <- as.vector(unlist(phen.rec))
    names(phen.rec) <- c(names(phen), c((length(phen)+1):((length(phen)*2)-1)))

    ###########################################
    ## get LOCATIONS (branches) of phen subs ##
    ###########################################

    ## make empty output list
    phen.subs.edges <- rep(list(NULL), 3)
    names(phen.subs.edges) <- c("total", "pos", "neg")

    ## identify if subs occur on each branch:
    phen.subs.logical <- sapply(c(1:nrow(edges)),
                                function(e)
                                  phen.rec[edges[e,1]]
                                ==
                                  phen.rec[edges[e,2]])
    ## get indices of all edges containing a substitution
    phen.subs.total <- which(phen.subs.logical == FALSE)
    ## get df of states of ancestor and descendants nodes on these edges
    df <- data.frame(phen.rec[edges[phen.subs.total,1]], phen.rec[edges[phen.subs.total,2]])
    names(df) <- c("anc", "dec")
    ## get indices of all edges w a positive sub (0 --> 1)
    phen.subs.pos <- phen.subs.total[which(df$anc==0)]
    ## get indices of all edges w a negative sub (1 --> 0)
    phen.subs.neg <- phen.subs.total[which(df$anc==1)]

    ## get output list
    if(length(phen.subs.total) > 0) phen.subs.edges[["total"]] <- phen.subs.total
    if(length(phen.subs.pos) > 0) phen.subs.edges[["pos"]] <- phen.subs.pos
    if(length(phen.subs.neg) > 0) phen.subs.edges[["neg"]] <- phen.subs.neg

    ####################
    ## PLOT to CHECK? ##
    ####################
    ## for SNP1, does it identify the correct/reasonable branches?
    #     edgeCol <- rep("black", nrow(edges))
    #     edgeCol <- replace(edgeCol, phen.subs.edges[["total"]], "green")
    #
    #     ## plot the i'th character's reconstruction on the tree:
    #     #require(adegenet)
    #     plotAnc(tree, phen.pa.ACCTRAN, i=1,
    #             col=transp(c("red", "royalblue"), 0.75),
    #             cex.pie=0.1, pos=NULL,
    #             edge.color=edgeCol, edge.width=2, use.edge.length=FALSE, type="c")

    ################
    ## Get output ##
    ################
    var.rec <- phen.rec
    subs.edges <- phen.subs.edges
    out <- list("var.rec" = var.rec,
                "subs.edges" = subs.edges)

  } # end vector (phen) parsimony

  return(out)

} # end get.ancestral.pars








#################
## reconstruct ##
#################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param var Either a matrix or a vector containing the state of a variable (eg. SNPs or a phenotype)
#' for all individuals (ie. for all terminal nodes in the tree).
#' @param tree A phylo object containing the tree representing the ancestral relationships
#' between the individuals for which snps and phen are known.
#' @param type A character string specifying whether ancestral state reconstruction should be
#' performed by \code{parsimony} or \code{ace} (as performed in package \emph{ape}).
#' @param method A character string specifying the type of ACE method to implement (only used if
#' \code{type} is set to "ace").
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#'
#' @import ape

########################################################################

reconstruct <- function(var,
                        tree,
                        type = c("parsimony", "ace"),
                        method = "discrete"){

  require(ape)

  ## get tree edges
  edges <- tree$edge

  #######################
  ## MATRIX (eg. SNPs) ##
  #######################
  if(class(var) == "matrix"){

    snps <- var

    ###############################
    ## get unique SNPs patterns: ##
    ###############################
    temp <- get.unique.matrix(snps, MARGIN=2)
    snps.unique <- temp$unique.data
    index <- temp$index

    if(ncol(snps.unique) == ncol(snps)){
      all.unique <- TRUE
    }else{
      all.unique <- FALSE
    }

    ## work w only unique snps:
    snps.ori <- snps
    snps <- snps.unique

    ############################
    ## run PARSIMONY on SNPs: ##
    ############################
    if(type == "parsimony"){

      ## run get.ancestral.pars
      snps.pars <- get.ancestral.pars(var=snps, tree=tree)

      ## get elements of output
      snps.rec <- snps.pars$var.rec
      snps.subs.edges <- snps.pars$subs.edges

    } # end parsimony


    ######################
    ## run ACE on SNPs: ##
    ######################
    if(type == "ace"){

      snps.rec <- snps.ACE <- list()

      for(i in 1:ncol(snps)){

        ## get variable i
        var <- snps[,i]

        ## get terminal values
        var.terminal <- var

        ## get internal values (from ACE output for variable i)
        snps.ACE[[i]] <- ace(var, tree, type=method)
        var.internal <- snps.ACE[[i]]$lik.anc[,2]

        ## get reconstruction from terminal & internal values
        snps.rec[[i]] <- c(var.terminal, var.internal)
      }

      ## bind columns of snps.rec together
      snps.rec <- do.call("cbind", snps.rec)
      colnames(snps.rec) <- colnames(snps)

    } # end ace


    ############################################
    ## get values for duplicate snps columns: ##
    ############################################

    ## get reconstruction for all original sites
    if(all.unique == TRUE){
      var.rec <- snps.rec
    }else{
      var.rec <- matrix(NA, nrow=nrow(snps.rec), ncol=ncol(snps.ori))
      for(i in 1:ncol(snps.rec)){
        var.rec[, which(index == i)] <- snps.rec[, i]
      }
      rownames(var.rec) <- rownames(snps.rec)
      colnames(var.rec) <- colnames(snps.ori)
    }

    ## get sub locations on branches for all original sites
    if(all.unique == TRUE){
      subs.edges <- snps.subs.edges
    }else{
      subs.edges <- vector("list", ncol(snps.ori))
      for(i in 1:length(snps.subs.edges)){
        loc <- which(index == i)
        if(length(loc) == 1){
          subs.edges[[loc]] <- snps.subs.edges[[i]]
        }else{
          for(j in 1:length(loc)){
            subs.edges[[loc[j]]] <- snps.subs.edges[[i]]
          }
        }
      }
    }



  }else{ # end matrix (snps)   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

    #######################
    ## VECTOR (eg. phen) ##
    #######################
    phen <- var

    ############################
    ## run PARSIMONY on phen: ##
    ############################
    if(type == "parsimony"){
      ## run get.ancestral.pars
      phen.pars <- get.ancestral.pars(var=phen, tree=tree)

      ## get elements of output
      var.rec <- phen.pars$var.rec
      subs.edges <- phen.pars$subs.edges

    } # end parsimony


    ######################
    ## run ACE on phen: ##
    ######################
    if(type == "ace"){
      ## Do we need to check phen is numeric??

      ## get terminal values
      phen.terminal <- phen

      ## get internal values (from ACE output)
      phen.ACE <- ace(phen, tree, type=method)
      phen.internal <- phen.ACE$lik.anc[,2]

      ## get reconstruction from terminal & internal values
      var.rec <- c(phen.terminal, phen.internal)

    } # end ace

  } # end vector (phen)   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###


  ################
  ## GET OUTPUT ##
  ################
  if(type == "parsimony"){
    output <- list("var.rec" = var.rec,
                   "subs.edges" = subs.edges)
  }else{
    ## NOTE that ACE does NOT return subs.edges...
    output <- list("var.rec" = var.rec)
  }

  ## return output
  return(output)

} # end reconstruct


