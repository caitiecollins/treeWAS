


#########
## asr ##
#########

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
#' performed by \code{parsimony} or \code{ML} (as performed by the \code{ace} function in package \emph{ape}).
#' @param method A character string specifying the type of ASR method to implement (only used if
#' \code{type} is set to "ML").
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#'
#' @import ape
#' @importFrom Hmisc all.is.numeric

########################################################################

asr <- function(var,
                tree,
                type = c("parsimony", "ML", "ace"), ## keeping "ace", in case I missed any instances, though deprecated.
                method = c("discrete", "continuous")){


  ## HANDLE TREE: ##
  ## Always work with trees in "pruningwise" order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  if(!is.rooted(tree)) tree <- midpoint(tree)

  ## get tree edges
  edges <- tree$edge
  ord <- NULL

  ## Arg checks:
  if(length(type) > 1) type <- type[1]
  type <- tolower(type)
  if(type == "ace") type <- "ml"

  if(length(method) > 1) method <- method[1]
  method <- tolower(method)

  #####################################################
  ## CHECK: If non-binary, use "ml" and "continuous" ##
  #####################################################
  levs <- unique(as.vector(unlist(var)))
  levs <- levs[!is.na(levs)]
  if(length(levs) != 2){
    ## If NON-BINARY: ##
    if(type != "ml"){
      type <- "ml"
      cat("Variable is non-binary. Setting reconstruction type to 'ML'")
    }

    if(method != "continuous"){
      method <- "continuous"
      cat("Variable is non-binary. Setting reconstruction method to 'continuous'")
    }
  }else{
    ## If BINARY: ##
    if(!type %in% c("parsimony", "ml", "ace")){
      type <- "parsimony"
      cat("Reconstruction type must be one of 'parsimony' or 'ML'. Selecting 'parsimony' by default.")
    }
    if(!method %in% c("discrete", "continuous")){
      method <- "discrete"
      cat("Reconstruction method must be one of 'discrete' or 'continuous'. Selecting 'discrete' by default.")
    }
  }

  #######################
  ## MATRIX (eg. SNPs) ##
  #######################
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

    print("N snps.unique (reconstruction):"); print(ncol(snps.unique))

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
    ## run ML on SNPs: ##
    ######################
    if(type == "ml"){

      if(method == "continuous"){
        # Check snps is numeric:
        if(!is.numeric(snps)){
          if(!all.is.numeric(snps)){
            stop("For a continuous reconstruction, the matrix must be numeric.")
          }else{
            r.noms <- rownames(snps)
            c.noms <- colnames(snps)
            snps <- matrix(as.numeric(as.character(snps)), nrow=nrow(snps), ncol=ncol(snps))
            rownames(snps) <- r.noms
            colnames(snps) <- c.noms
          }
        }
      }

      ## Check for zero- or negative-length tree$edge.length:
      ## (For now??) just replace w 0.0000000001
      if(any(tree$edge.length <= 0)){
        toReplace <- which(tree$edge.length <= 0)
        tree$edge.length[toReplace] <- 1e-10
      }

      snps.rec <- snps.ML <- list()

      for(i in 1:ncol(snps)){

        ## get variable i
        var <- snps[,i]

        ## get terminal values
        var.terminal <- var

        ## get internal values (from ML output for variable i)
        snps.ML[[i]] <- ace(var, tree, type=method)
        if(method == "discrete"){
          var.internal <- snps.ML[[i]]$lik.anc[,2]
        }else{
          var.internal <- snps.ML[[i]]$ace
        }

        ## get reconstruction from terminal & internal values
        snps.rec[[i]] <- c(var.terminal, var.internal)
      }

      ## bind columns of snps.rec together
      snps.rec <- do.call("cbind", snps.rec)
      colnames(snps.rec) <- colnames(snps)#

    } # end ML


    ############################################
    ## get values for duplicate snps columns: ##
    ############################################

    ## get reconstruction for all original sites
    if(all.unique == TRUE){
      var.rec <- snps.rec
    }else{
      var.rec <- snps.rec[, index]
      rownames(var.rec) <- rownames(snps.rec)
      colnames(var.rec) <- colnames(snps.ori)
    }

    ## get sub locations on branches for all original sites
    if(type == "parsimony"){
      if(all.unique == TRUE){
        subs.edges <- snps.subs.edges
      }else{
        subs.edges <- snps.subs.edges[index]
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
    ## run ML on phen: ##
    ######################
    if(type == "ml"){

      if(method == "continuous"){
        ## Check phen is numeric:
        if(!is.numeric(phen)){
          if(!all.is.numeric(phen)){
            stop("For a continuous reconstruction, the vector must be numeric.")
          }else{
            noms <- names(phen)
            phen <- as.numeric(as.character(phen))
            names(phen) <- noms
          }
        }
      }

      ## Check for zero- or negative-length tree$edge.length:
      ## (For now??) just replace w 0.0000000001
      if(any(tree$edge.length <= 0)){
        toReplace <- which(tree$edge.length <= 0)
        tree$edge.length[toReplace] <- 1e-10
      }

      ## get terminal values
      phen.terminal <- phen

      ## get internal values (from ML output)
      phen.ML <- ace(phen, tree, type=method)
      if(method == "discrete"){
        phen.internal <- phen.ML$lik.anc[,2]
      }else{
        phen.internal <- phen.ML$ace
      }

      ## get reconstruction from terminal & internal values
      var.rec <- c(phen.terminal, phen.internal)

    } # end ML

  } # end vector (phen)   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###


  ################
  ## GET OUTPUT ##
  ################
  if(type == "parsimony"){
    output <- list("var.rec" = var.rec,
                   "subs.edges" = subs.edges)
  }else{
    ## NOTE that ML does NOT return subs.edges...
    output <- list("var.rec" = var.rec)
  }

  ## return output
  return(output)

} # end asr













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
#' @import ape
#' @importFrom phangorn as.phyDat
#' @importFrom phangorn phyDat
#' @importFrom phangorn pace
#'
#' @export
#'

########################################################################

get.ancestral.pars <- function(var, tree){

  ## HANDLE TREE: ##
  ## Always work with trees in "pruningwise" order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  if(!is.rooted(tree)) tree <- midpoint(tree)

  ord <- NULL
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
    ## RUN CHECKS TO ENSURE tree$tip.label and rownames(snps) CONTAIN SAME SET OF LABELS!
    if(is.null(tree$tip.label)) stop("Trees must have tip.labels corresponding to rownames(snps).")
    if(is.null(rownames(snps))) stop("SNPs must have rownames corresponding to tree$tip.label.")
    if(!all(tree$tip.label %in% rownames(snps))) stop("tree$tip.label and rownames(snps)
                                                      must contain the same set of labels
                                                      so that individuals can be correctly identified.")
    ## Assign colnames if NULL:
    if(is.null(colnames(snps))) colnames(snps) <- c(1:ncol(snps))

    ## get levels (ie. 0, 1)
    # snps.levels <- sort(unique(as.vector(snps)))
    snps.levels <- sort(unique(as.vector(snps)), na.last = TRUE)
    ## returns only unique patterns...
    snps.phyDat <- as.phyDat(as.matrix(snps),
                             type="USER", levels=snps.levels)
    ## get index of all original snps columns to map to unique pattern
    index.phyDat <- attr(snps.phyDat, "index")

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

    # REC <- rec

    ## If NAs are present, replace column 1 0s with NA values
    ## whenever column 3 (NA) has a 1 in it:
    if(any(is.na(snps.levels))){
      na.col <- which(is.na(snps.levels))
      for(i in 1:length(rec)){
        foo <- rec[[i]]
        toReplace <- which(foo[,na.col] == 1)
        if(length(toReplace) > 0){
          foo[toReplace,1] <- NA
          foo[toReplace,2] <- NA
        }
        rec[[i]] <- foo
      } # end for loop
    }

    ## NOTE: pace works with terminal SNPs in the order they appear in tree$tip.label
    ## First, check to ensure all row.names(snps) are matched in tree$tip.label
    if(all(row.names(snps) %in% tree$tip.label)){
      ord <- match(tree$tip.label, rownames(snps))
    }else{
      ord <- 1:length(rownames(snps))
      warning("rownames(snps) and tree$tip.label contain different labels.
              Careful-- we proceed by assuming snps rows and tree tips are labelled in the same order!")
    }
    ## Want rec (list) to be in order of tree$tip.label
    ## eg. if tree$tip.label[1] is "31", ord[1] should be 31 (assuming rownames(snps) are 1:nrow)
    ord <- c(ord, c(nrow(snps)+1):length(rec))
    snps.rec <- do.call(cbind, rec[ord])
    snps.rec <- t(snps.rec[, seq(2, ncol(snps.rec), length(snps.levels))])

    ## assign rownames for all terminal and internal nodes
    rownames(snps.rec) <- c(rownames(snps), c((nrow(snps)+1):((nrow(snps)*2)-1)))
    colnames(snps.rec) <- c(1:length(snps.phyDat[[1]]))


    ## Handle index.phyDat! ##
    ## get reconstruction for all pre-phyDat sites
    if(ncol(snps) != ncol(snps.rec)){
      snps.rec.complete <- snps.rec[, index.phyDat]
      rownames(snps.rec.complete) <- rownames(snps.rec)
      colnames(snps.rec.complete) <- 1:ncol(snps.rec.complete)
      snps.rec <- snps.rec.complete
    }






    ###########################################
    ## get LOCATIONS (branches) of snps subs ##
    ###########################################
    subs.edges <- rep(list(NULL), ncol(snps.rec))

    subs.logical <- matrix(snps.rec[edges[, 1], ] == snps.rec[edges[, 2], ], nrow=nrow(edges), byrow=F)

    ## get states of anc and dec:
    df.anc <- snps.rec[edges[, 1], ]
    # df.dec <- snps.rec[edges[, 2], ]


    ## get indices of all edges containing a substitution
    for(i in 1:ncol(snps.rec)){
      subs.total <- which(subs.logical[, i] == FALSE)

      ## get indices of all edges w a positive sub (0 --> 1)
      subs.pos <- subs.total[which(subs.total %in% which(df.anc[,i] == 0))]
      ## get indices of all edges w a negative sub (1 --> 0)
      subs.neg <- subs.total[which(subs.total %in% which(df.anc[,i] == 1))]

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


    ############################################
    ## get values for duplicate snps columns: ##
    ############################################

    ## get reconstruction for all original sites
    if(ncol(snps.ori) == ncol(snps.rec)){
      snps.rec.complete <- snps.rec
    }else{
      snps.rec.complete <- snps.rec[, index]
      rownames(snps.rec.complete) <- rownames(snps.rec)
      colnames(snps.rec.complete) <- colnames(snps.ori)
    }

    ## get sub locations on branches for all original sites
    snps.subs.edges <- subs.edges
    if(ncol(snps.ori) == ncol(snps.rec)){
      snps.subs.edges.complete <- snps.subs.edges
    }else{
      snps.subs.edges.complete <- snps.subs.edges[index]
    }


    ################
    ## Get output ##
    ################

    ## CHECK-- compare cost from fitch and pace: ##
    ###########
    ## get n.subs per site by fitch:
    # cost <- get.fitch.n.mts(snps, tree)

    ## get n.subs per site by pace:
    # cost2 <- sapply(c(1:length(snps.subs.edges.complete)),
    # function(e) length(snps.subs.edges.complete[[e]][["total"]]))

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
    phen <- as.numeric(as.factor(phen))
    names(phen) <- names(phen.ori)

    ## as.phyDat requires names...
    ## RUN CHECKS TO ENSURE tree$tip.label and rownames(snps) CONTAIN SAME SET OF LABELS!
    if(is.null(tree$tip.label)) stop("Trees must have tip.labels corresponding to names(phen).")
    if(is.null(names(phen))) stop("Phen must have names corresponding to tree$tip.label.")
    if(!all(tree$tip.label %in% names(phen))) stop("tree$tip.label and names(phen)
                                                      must contain the same set of labels
                                                      so that individuals can be correctly identified.")

    ## get levels (ie. 0, 1)
    phen.levels <- sort(unique(phen))
    phen.phyDat <- as.phyDat(as.matrix(phen),
                             type="USER", levels=phen.levels)
    ## pace == ancestral.pars
    rec <- phen.pa.ACCTRAN <- pace(tree, phen.phyDat, type="ACCTRAN")

    ## get reconstruction:
    if(is.null(ord)) ord <- 1:length(rec)
    phen.rec <- do.call(cbind, rec[ord])
    phen.rec <- as.vector(phen.rec[, seq(2, ncol(phen.rec), 2)])

    names(phen.rec) <- c(names(phen), c((length(phen)+1):((length(phen)*2)-1)))

    ###########################################
    ## get LOCATIONS (branches) of phen subs ##
    ###########################################

    ## make empty output list
    phen.subs.edges <- rep(list(NULL), 3)
    names(phen.subs.edges) <- c("total", "pos", "neg")

    ## identify if subs occur on each branch:
    phen.subs.logical <- phen.rec[edges[, 1]] == phen.rec[edges[, 2]]
    names(phen.subs.logical) <- 1:nrow(edges) ## WHY DID IT AUTOMATICALLY LABEL THE ROWS IN REVERSE ORDER (199:101) ?????
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



