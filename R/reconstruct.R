


#########
## asr ##
#########

########################################################################

###################
## DOCUMENTATION ##
###################

#' Ancestral state reconstruction
#'
#' Reconstruct the ancestral states of a vector or matrix object by using either
#' parsimony or maximum-likelihood methods to infer the states
#' at the internal nodes of a phylogenetic tree.
#'
#' @param var Either a matrix or a vector containing the state of a variable (eg. SNPs or a phenotype)
#' for all individuals (ie. for all terminal nodes in the tree).
#' @param tree A phylo object containing the tree representing the ancestral relationships
#' between the individuals for which snps and phen are known.
#' @param type A character string specifying whether ancestral state reconstruction should be
#' performed by \code{parsimony} or \code{ML} (as performed by the \code{ace} function in package \emph{ape}).
#' @param method A character string specifying the type of ASR method to implement,
#' either \code{'discrete'} or \code{'continuous'} (only used if \code{type} is set to "ML").
#' @param unique.cols A logical indicating whether only unique column patterns are present in \code{var},
#' if \code{var} is a matrix (if so (\code{TRUE}), a time-consuming step can be skipped);
#' by default, \code{FALSE}.
#'
#' @return Depending on the dimensions of the input \code{var} object,
#' either a matrix or a vector containing \emph{both} the known states
#' of the variable at the terminal nodes (in positions 1:Nterminal) and the
#' inferred states at internal nodes (in positions (Nterminal+1):Ntotal).
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#'
#' @rawNamespace import(ape, except = zoom)
#' @importFrom Hmisc all.is.numeric
#' @importFrom phytools anc.ML fastAnc
#' @importFrom phangorn pml ancestral.pml


########################################################################

asr <- function(var,
                tree,
                type = c("parsimony", "ML", "ace"), ## keeping "ace", deprecated.
                method = c("discrete", "continuous"),
                unique.cols = FALSE){


  ## HANDLE TREE: ##
  ## Always work with trees in "pruningwise" order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  # if(!is.rooted(tree)) tree <- midpoint(tree)

  ## get tree edges
  edges <- tree$edge
  ord <- NULL

  ## Arg checks:
  if(length(type) > 1) type <- type[1]
  type <- tolower(type)
  if(type == "ace") type <- "ml"

  if(length(method) > 1) method <- method[1]
  method <- tolower(method)

  #######################
  ## MATRIX (eg. SNPs) ##
  #######################
  if(is.matrix(var)){

    #################
    ## CHECK ARGS: ##
    #################
    ## BINARY (assumed for snps): ##
    if(!type %in% c("parsimony", "ml", "ace")){
      type <- "parsimony"
      cat("Reconstruction type must be one of 'parsimony' or 'ML'. Selecting 'parsimony' by default.\n")
    }
    if(!method %in% c("discrete", "continuous")){
      method <- "discrete"
      cat("Reconstruction method must be one of 'discrete' or 'continuous'. Selecting 'discrete' by default.\n")
    }

    ## Assign var to snps:
    snps <- var

    if(all(rownames(snps) %in% tree$tip.label)){
      ord <- match(tree$tip.label, rownames(snps))
      snps <- snps[ord, ]
    }else{
      warning("Unable to rearrange var such that rownames(var)
                   match tree$tip.label. Order may be incorrect")
    }

    if(unique.cols == TRUE){
      all.unique <- TRUE
    }else{
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
    }

    # print("N snps.unique (reconstruction):"); print(ncol(snps.unique))

    ############################
    ## run PARSIMONY on SNPs: ##
    ############################
    if(type == "parsimony"){

      ## run get.ancestral.pars
      snps.rec <- get.ancestral.pars(var=snps, tree=tree, unique.cols = TRUE)

    } # end parsimony


    ######################
    ## run ML on SNPs: ##
    ######################
    if(type == "ml"){

      if(method == "continuous"){
        # Check snps is numeric:
        if(!is.numeric(snps)){
          if(!all.is.numeric(snps[!is.na(snps)])){
            stop("For a continuous reconstruction, the matrix must be numeric.\n")
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

      ## Check if MISSING DATA in snps: ##
      na.var <- FALSE
      ## If ANY snps column contains NAs, we change rec fn
      ## for continuous recs for ALL columns (for consistency's sake).
      if(any(is.na(as.vector(unlist(snps))))) na.var <- TRUE


      ## With CONTINUOUS ML rec, use anc.ML (or fastAnc!): ##
      if(method == "continuous"){

        ######### (** SLOW! **) ##########
        ## MISSING & Continuous ML rec: ##
        ## (not yet allowed for snps..) ##
        ##################################

        snps.rec <- snps.ML <- list()

        lgc <- FALSE
        if(is.logical(snps)) lgc <- TRUE

        for(i in 1:ncol(snps)){

          ## get variable i
          var <- snps[,i]

          ## get terminal values
          var.terminal <- var

          ## With MISSING DATA in any columns: ##
          ## Continuous ML rec: ##
          ## get internal values (when (any) var contains NAs):
          var <- var[!is.na(var)]
          snps.ML[[i]] <- anc.ML(tree, var) ## require(phytools)
          var.internal <- snps.ML[[i]]$ace

          ## get reconstruction from terminal & internal values
          snps.rec[[i]] <- c(var.terminal, var.internal)
          if(lgc == TRUE){
            noms <- names(snps.rec[[i]])
            snps.rec[[i]] <- as.logical(snps.rec[[i]])
            names(snps.rec[[i]]) <- noms
          }
        } # end for loop

        ## bind columns of snps.rec together
        snps.rec <- do.call("cbind", snps.rec)
        colnames(snps.rec) <- colnames(snps)

      }else{

        ###############
        ## DISCRETE: ##
        ###############

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

        ## ML discrete reconstruction:
        ## Step 1:
        fit <- pml(tree, snps.phyDat)
        ## (Step 2 below -- if binary, get probs; if not, get states...)
        # rec <- rec.ml <- ancestral.pml(fit, type = "ml", return = "prob")
        # rec <- rec.ml <- ancestral.pml(fit, type = "ml", return = "phyDat")

        ###########################################
        ## convert reconstruction back to snps.. ##
        ###########################################
        ## each of the n.ind elements of rec is a matrix w n.snps rows and either:
        ## 2 columns, for the 2 binary SNP states, or
        ## 4 columns, each for the 4 nts possible (acgt)

        ## Want to KEEP rec list in order of tree$tip.label to match tree$edge!
        l <- max(tree$edge[,2])
        ord <- 1:l

        ## Binary:
        if(length(snps.levels[!is.na(snps.levels)]) == 2){

          ## ML discrete reconstruction:
          ## Step 2 (Binary --> get probs):
          # rec <- rec.ml <- ancestral.pml(fit, type = "ml", return = "prob")
          rec <- rec.ml <- ancestral.pml(fit, type = "ml", return = "phyDat")

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
          ## Bind list into matrix:
          snps.rec <- do.call(rbind, rec[ord])
          # snps.rec <- t(snps.rec[, seq(2, ncol(snps.rec), length(snps.levels))]) # only for return="prob"

        }else{

          ## Non-Binary (discrete):

          ## ML discrete reconstruction:
          ## Step 2 (Non-binary --> get states):
          # rec <- rec.ml <- ancestral.pml(fit, type = "ml", return = "prob")
          rec <- rec.ml <- ancestral.pml(fit, type = "ml", return = "phyDat")

          ## Print warning notice:
          cat("Reconstructing non-binary genetic data matrix.\n")

          ## Bind list & reorder elements:
          snps.rec <- do.call(rbind, rec[ord])

        } # end non-binary

        ## Convert to restore levels, variable class:
        r.noms <- rownames(snps.rec)
        c.noms <- colnames(snps.rec)
        if(is.logical(snps.ori)){
          snps.rec <- matrix(as.logical(snps.rec), nrow=nrow(snps.rec), ncol=ncol(snps.rec)) ## logical
        }else{
          ## replace level #s w levels:
          snps.rec <- matrix(attr(rec, "levels")[snps.rec], nrow=nrow(snps.rec), ncol=ncol(snps.rec))
        }
        rownames(snps.rec) <- r.noms ## names re-assigned below?
        colnames(snps.rec) <- c.noms

      } # end Discrete ML


      ## assign rownames for all terminal and internal nodes
      if(is.null(tree$node.label)){
        rownames(snps.rec) <- c(tree$tip.label, c((nrow(snps)+1):max(tree$edge[,2])))
      }else{
        rownames(snps.rec) <- c(tree$tip.label, tree$node.label)
      }
      colnames(snps.rec) <- c(1:length(snps.phyDat[[1]]))


      ## Handle index.phyDat! ##
      ## get reconstruction for all pre-phyDat sites
      if(ncol(snps) != ncol(snps.rec)){
        snps.rec.complete <- snps.rec[, index.phyDat]
        rownames(snps.rec.complete) <- rownames(snps.rec)
        colnames(snps.rec.complete) <- 1:ncol(snps.rec.complete)
        snps.rec <- snps.rec.complete
      }

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
    }
    if(ncol(var.rec) == ncol(snps.ori)){
      colnames(var.rec) <- colnames(snps.ori)
    }


  }else{ # end matrix (snps)   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

    #######################
    ## VECTOR (eg. phen) ##
    #######################

    #################
    ## CHECK ARGS: ##
    #################
    levs <- unique(as.vector(unlist(var)))
    levs <- levs[!is.na(levs)] # no NAs allowed in phen variables in GWAS.
    ## BINARY: ##
    if(length(levs) == 2){
      ## Choose parsimony if none:
      if(!type %in% c("parsimony", "ml")){
        type <- "parsimony"
        cat("Reconstruction type must be one of 'parsimony' or 'ML'. Selecting 'parsimony' by default.\n")
      }
      if(type == "ml"){
        if(!method %in% c("discrete", "continuous")){
          method <- "discrete"
        }
      }
    }else{
      ## DISCRETE or CONTINUOUS: ##
      if(!type %in% c("parsimony", "ml")){
        type <- "ml"
        cat("Reconstruction type must be one of 'parsimony' or 'ML'. Selecting 'ML' by default.\n")
      }
      if(!method %in% c("discrete", "continuous")){
        if(type == "ml"){
          method <- "continuous"
          cat("Reconstruction method must be one of 'discrete' or 'continuous'. Selecting 'continuous' by default.\n")
        }
      }
      if(method == "continuous"){
        if(type != "ml"){
          type <- "ml"
          cat("Reconstruction type must be 'ML' when variable is 'continuous'. Setting type to 'ML'.\n")
        }
      }

      ## Parsimony NOT available if >75% unique:
      if(length(levs)/length(var) > 0.75){
        if(type == "parsimony"){
          type <- "ml"
          cat("Parsimony not available for variables > 75% unique. Setting reconstruction type to 'ML'.\n")
        }
      }
    } # end checks


    ## Assign var to phen:
    phen <- var

    if(all(names(phen) %in% tree$tip.label)){
      ord <- match(tree$tip.label, names(phen))
      phen <- phen[ord]
    }else{
      warning("Unable to rearrange var such that rownames(var)
                   match tree$tip.label. Order may be incorrect")
    }

    ###############################
    ## *If phen is NOT variable: ##
    ###############################
    if(length(levs) == 1){
      phen.rec <- rep(levs, max(tree$edge[,2]))
    }else{

    ############################
    ## run PARSIMONY on phen: ##
    ############################
    if(type == "parsimony"){
      ## run get.ancestral.pars
      phen.rec <- get.ancestral.pars(var=phen, tree=tree)

    } # end parsimony

    ######################
    ## run ML on phen:  ##
    ######################
    if(type == "ml"){

      if(method == "continuous"){
        ## Check phen is numeric:
        if(!is.numeric(phen)){
          if(!all.is.numeric(phen)){
            stop("For a continuous reconstruction, phen must be numeric.\n")
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

      # ## get internal values (from ML output)
      # if(is.rooted(tree) & is.binary.tree(tree) & length(levs) == 2){
      #   phen.ML <- ace(phen, tree, type=method)
      #   if(method == "discrete"){
      #     phen.internal <- phen.ML$lik.anc[,2]
      #   }else{
      #     phen.internal <- phen.ML$ace
      #   }
      #
      #   ## get reconstruction from terminal & internal values
      #   phen.rec <- c(phen.terminal, phen.internal)
      #
      # }else{ # end ML for rooted & binary trees
      #

      ## ML for all (rooted/unrooted, binary/non-binary) trees: ##

      ## DISCRETE: ##
      if(method == "discrete"){

        ## get levels (ie. 0, 1)
        phen.levels <- sort(unique(phen), na.last = TRUE)
        ## returns only unique patterns...
        phen.phyDat <- as.phyDat(as.matrix(phen),
                                 type="USER", levels=phen.levels)
        ## get index of all original snps columns to map to unique pattern
        index.phyDat <- attr(phen.phyDat, "index")

        ## ML discrete reconstruction:
        ## Step 1:
        fit <- pml(tree, phen.phyDat)

        ## Want to KEEP rec list in order of tree$tip.label to match tree$edge!
        l <- max(tree$edge[,2])
        ord <- 1:l

        ## Step 2 (Non-binary --> get states):
        rec <- rec.ml <- ancestral.pml(fit, type = "ml", return = "phyDat")

        ## Bind list & reorder elements:
        phen.rec <- do.call(rbind, rec[ord])

        ## replace level #s w levels:
        phen.rec <- attr(rec, "levels")[phen.rec]

      }else{

        ## CONTINUOUS: ##

        ## get terminal values
        phen.terminal <- phen
        ox <- match(tree$tip.label, names(phen))
        phen.terminal <- phen.terminal[ox]

        ## With MISSING DATA in any columns: ##
        ## Continuous ML rec: ##
        ## get internal values:
        phen.internal <- fastAnc(tree, phen) ## require(phytools)

        ## get internal values (when (any) var contains NAs):
        # phen2 <- phen[!is.na(phen)]
        # phen.ML <- anc.ML(tree, phen2) ## require(phytools)
        # phen.internal <- phen.ML$ace

        ## get reconstruction from terminal & internal values
        phen.rec <- c(phen.terminal, phen.internal)

      }

    } # end ML

    } # end (rec | Nlevs > 1)

    ## assign rownames for all terminal and internal nodes
    if(is.null(tree$node.label)){
      names(phen.rec) <- c(tree$tip.label,  c((length(phen)+1):max(tree$edge[,2])))
    }else{
      names(phen.rec) <- c(tree$tip.label, tree$node.label)
    }

    var.rec <- phen.rec

  } # end vector (phen)   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###


  ################
  ## GET OUTPUT ##
  ################
  output <- var.rec

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
#' @rawNamespace import(ape, except = zoom)
#' @importFrom phangorn as.phyDat
#' @importFrom phangorn phyDat
#' @importFrom phangorn pace
#' @importFrom phangorn midpoint
#'
#' @export
#'

########################################################################
# @useDynLib phangorn, .registration = TRUE
# @useDynLib phangorn, as.phyDat,  phyDat, pace, midpoint, .registration = TRUE

get.ancestral.pars <- function(var, tree, unique.cols = FALSE){

  ## HANDLE TREE: ##
  ## Always work with trees in "pruningwise" order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  # if(!is.rooted(tree)) tree <- midpoint(tree)

  ord <- NULL
  edges <- tree$edge

  #############################
  ## RUN PARSIMONY on MATRIX ##
  #############################

  if(is.matrix(var)){

    snps <- var

    if(unique.cols == TRUE){
      all.unique <- TRUE
    }else{
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
    }

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
    rec <- pa.ACCTRAN <- pace(tree, snps.phyDat, type="ACCTRAN")

    ## NOTE: pace  --> diff resuls w MPR vs. ACCTRAN
    # rec <- pa.MPR <- pace(tree, snps.phyDat, type="MPR")
    #diffs <- sapply(c(1:length(pa.ACCTRAN)), function(e) identical(pa.MPR[[e]], pa.ACCTRAN[[e]]))

    ## ML discrete alternative:
    # fit <- pml(tree, snps.phyDat)
    # rec.ml <- ancestral.pml(fit, type = "ml")

    ###########################################
    ## convert reconstruction back to snps.. ##
    ###########################################
    ## each of the n.ind elements of pa is a matrix w n.snps rows and either:
    ## 2 columns, for the 2 binary SNP states, or
    ## 4 columns, each for the 4 nts possible (acgt)

    ## Want to KEEP rec list in order of tree$tip.label to match tree$edge!
    ord <- 1:length(rec)

    ## Binary:
    if(length(snps.levels[!is.na(snps.levels)]) == 2){

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

      snps.rec <- do.call(cbind, rec[ord])
      snps.rec <- t(snps.rec[, seq(2, ncol(snps.rec), length(snps.levels))])

      ## Convert to logical (but only for binary snps that were originally logical):
      if(is.logical(snps)){
        r.noms <- rownames(snps.rec)
        c.noms <- colnames(snps.rec)
        snps.rec <- matrix(as.logical(snps.rec), nrow=nrow(snps.rec), ncol=ncol(snps.rec)) ## logical
        rownames(snps.rec) <- r.noms
        colnames(snps.rec) <- c.noms
      }

    }else{

      ## Non-Binary (discrete):

      ## Print warning notice:
      cat("Reconstructing non-binary genetic data matrix.\n")

      ## Reorder elements:
      rec <- rec[ord]

      ## For each row, get values:
      snps.rec <- list()
      for(i in 1:length(rec)){
        ## Get reconstruction for this row for all sites:
        mat <- rec[i][[1]]
        ## Replace non-binary (uncertain) values w NA:
        mat <- replace(mat, which(!mat %in% c(0,1)), NA)
        ## Convert to logical:
        sr <- matrix(as.logical(mat), nrow=nrow(mat), ncol=ncol(mat))
        ## Replace any row containing NAs to ONE NA:
        toReplace <- which(is.na(rowSums(sr, na.rm = FALSE)))
        sr[toReplace,] <- FALSE
        sr[toReplace, 1] <- NA
        ## Get values:
        foo <- rep(snps.levels, nrow(sr))
        ## Append to list:
        snps.rec[[i]] <- foo[t(sr)]
      } # end for loop

      ## Bind list into matrix:
      snps.rec <- do.call(rbind, snps.rec)

    } # end non-binary


    ## assign rownames for all terminal and internal nodes
    if(is.null(tree$node.label)){
      rownames(snps.rec) <- c(tree$tip.label, c((nrow(snps)+1):max(tree$edge[,2])))
    }else{
      rownames(snps.rec) <- c(tree$tip.label, tree$node.label)
    }
    colnames(snps.rec) <- c(1:length(snps.phyDat[[1]]))


    ## Handle index.phyDat! ##
    ## get reconstruction for all pre-phyDat sites
    if(ncol(snps) != ncol(snps.rec)){
      snps.rec.complete <- snps.rec[, index.phyDat]
      rownames(snps.rec.complete) <- rownames(snps.rec)
      colnames(snps.rec.complete) <- 1:ncol(snps.rec.complete)
      snps.rec <- snps.rec.complete
    }

    ############################################
    ## get values for duplicate snps columns: ##
    ############################################

    ## get reconstruction for all original sites
    # if(ncol(snps.ori) == ncol(snps.rec)){
    if(all.unique == TRUE){
      snps.rec.complete <- snps.rec
    }else{
      snps.rec.complete <- snps.rec[, index]
      rownames(snps.rec.complete) <- rownames(snps.rec)
      colnames(snps.rec.complete) <- colnames(snps.ori)
    }

    ## Get Output:
    out <- snps.rec.complete


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
    ## Make phen numeric:
    if(!is.numeric(phen)){
      if(all.is.numeric(phen)){
        phen <- as.numeric(as.character(phen))
      }else{
        phen <- as.numeric(as.factor(phen))
        warning("Reconstructing non-numeric phen.")
      }
    }
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
    if(class(try(suppressWarnings(pace(tree, phen.phyDat, type="ACCTRAN")), silent=T)) == "try-error"){

      rec <- phen.pa.ACCTRAN <- pace(tree, phen.phyDat, type="MPR")

    }else{
      rec <- phen.pa.ACCTRAN <- pace(tree, phen.phyDat, type="ACCTRAN")

    }
    ## Want to KEEP rec list in order of tree$tip.label to match tree$edge!
    ord <- 1:length(rec)

    ## get reconstruction:

    ## Binary:
    if(length(phen.levels) == 2){
      ## Combine in order:
      phen.rec <- do.call(cbind, rec[ord])
      ## Remove redundant column:
      phen.rec <- as.vector(phen.rec[, seq(2, ncol(phen.rec), 2)])
    }else{

      ## Non-Binary (discrete):
      ## Get values:
      phen.rec <- do.call(rbind, rec)
      phen.rec <- phen.rec[ord,]

      # ## Replace non-binary (?!) values w NA:
      # phen.rec <- replace(phen.rec, which(!phen.rec %in% c(0,1)), NA)
      # ## Convert to logical matrix:
      # pr <- as.logical(phen.rec)
      # pr <- matrix(pr, nrow=nrow(phen.rec), ncol=ncol(phen.rec))
      # # Replace any row containing NAs to ONE NA:
      # toReplace <- which(is.na(rowSums(pr, na.rm = FALSE)))
      # pr[toReplace,] <- FALSE
      # pr[toReplace, 1] <- NA
      #####
      # ## Get values:
      # levs <- attr(rec, "levels")
      # foo <- rep(levs, nrow(pr))
      # phen.rec <- foo[t(pr)]

      ## Get values:
      pr <- phen.rec
      levs <- attr(rec, "levels")
      foo <- matrix(rep(levs, nrow(pr)), nrow=nrow(pr), byrow=TRUE)
      pr2 <-  rep(NA, nrow(pr))
      for(e in 1:nrow(foo)){
        val <- foo[e, as.logical(pr[e,])]
        ## if no val or ties, set to NA:
        if(length(val) != 1) val <- NA
        pr2[e] <- val
      }
      phen.rec <- pr2

    } # end non-binary

    if(is.null(tree$node.label)){
      names(phen.rec) <- c(tree$tip.label,  c((length(phen)+1):max(tree$edge[,2])))
    }else{
      names(phen.rec) <- c(tree$tip.label, tree$node.label)
    }

    ## Get Output
    out <- phen.rec

  } # end vector (phen) parsimony

  return(out)

} # end get.ancestral.pars












##############################################################################
#
# ###################
# ## OLD SNPS CODE ##
# ###################
#
# ###########################################
# ## get LOCATIONS (branches) of snps subs ##
# ###########################################
# subs.edges <- rep(list(NULL), ncol(snps.rec))
#
# subs.logical <- matrix(snps.rec[edges[, 1], ] == snps.rec[edges[, 2], ], nrow=nrow(edges), byrow=F)
#
# ## get states of anc and dec:
# df.anc <- snps.rec[edges[, 1], ]
# # df.dec <- snps.rec[edges[, 2], ]
#
#
# ## get indices of all edges containing a substitution
# for(i in 1:ncol(snps.rec)){
#   subs.total <- which(subs.logical[, i] == FALSE)
#
#   ## get indices of all edges w a positive sub (0 --> 1)
#   subs.pos <- subs.total[which(subs.total %in% which(df.anc[,i] == 0))]
#   ## get indices of all edges w a negative sub (1 --> 0)
#   subs.neg <- subs.total[which(subs.total %in% which(df.anc[,i] == 1))]
#
#   ## get output list
#   subs.edges[[i]] <- rep(list(NULL), 3)
#   names(subs.edges[[i]]) <- c("total", "pos", "neg")
#   if(length(subs.total) > 0) subs.edges[[i]][["total"]] <- subs.total
#   if(length(subs.pos) > 0) subs.edges[[i]][["pos"]] <- subs.pos
#   if(length(subs.neg) > 0) subs.edges[[i]][["neg"]] <- subs.neg
# }
#
# ####################
# ## PLOT to CHECK? ##
# ####################
# ## for SNP1, does it identify the correct/reasonable branches?
# #     edgeCol <- rep("black", nrow(edges))
# #     edgeCol <- replace(edgeCol, subs.edges[[1]][["total"]], "green")
# #
# #     ## plot the i'th character's reconstruction on the tree:
# #     #require(adegenet)
# #     plotAnc(tree, pa.ACCTRAN, i=1,
# #             col=transp(c("red", "royalblue"), 0.75),
# #             cex.pie=0.1, pos=NULL,
# #             edge.color=edgeCol, edge.width=2, use.edge.length=FALSE, type="c")
# ################
# ## Get output ##
# ################
#
# ## CHECK-- compare cost from fitch and pace: ##
# ###########
# ## get n.subs per site by fitch:
# # cost <- get.fitch.n.mts(snps, tree)
#
# ## get n.subs per site by pace:
# # cost2 <- sapply(c(1:length(snps.subs.edges.complete)),
# # function(e) length(snps.subs.edges.complete[[e]][["total"]]))
#
# ## NOTE: cost2 differs somewhat noticeably from original fitch cost
# ## (ie. parsimony shifts distribution toward 1/reduces the weight of the upper tail...)
# ## WHY? Which should we use to get n.subs????????????????????????????????????????????????????????????
#
# ## get sub locations on branches for all original sites
# snps.subs.edges <- subs.edges
# # if(ncol(snps.ori) == ncol(snps.rec)){
# if(all.unique == TRUE){
#   snps.subs.edges.complete <- snps.subs.edges
# }else{
#   snps.subs.edges.complete <- snps.subs.edges[index]
# }
# ## Get final output list:
# var.rec <- snps.rec.complete
# subs.edges <- snps.subs.edges.complete
#
# out <- list("var.rec" = var.rec,
#             "subs.edges" = subs.edges)
##############################################################################
#
# ###################
# ## OLD PHEN CODE ##
# ###################
#
# ###########################################
# ## get LOCATIONS (branches) of phen subs ##
# ###########################################
#
# ## make empty output list
# phen.subs.edges <- rep(list(NULL), 3)
# names(phen.subs.edges) <- c("total", "pos", "neg")
#
# ## identify if subs occur on each branch:
# phen.subs.logical <- phen.rec[edges[, 1]] == phen.rec[edges[, 2]]
# names(phen.subs.logical) <- 1:nrow(edges) ## WHY DID IT AUTOMATICALLY LABEL THE ROWS IN REVERSE ORDER (199:101) ?????
# ## get indices of all edges containing a substitution
# phen.subs.total <- which(phen.subs.logical == FALSE)
# ## get df of states of ancestor and descendants nodes on these edges
# df <- data.frame(phen.rec[edges[phen.subs.total,1]], phen.rec[edges[phen.subs.total,2]])
# names(df) <- c("anc", "dec")
# ## get indices of all edges w a positive sub (0 --> 1)
# phen.subs.pos <- phen.subs.total[which(df$anc==0)]
# ## get indices of all edges w a negative sub (1 --> 0)
# phen.subs.neg <- phen.subs.total[which(df$anc==1)]
#
# ## get output list
# if(length(phen.subs.total) > 0) phen.subs.edges[["total"]] <- phen.subs.total
# if(length(phen.subs.pos) > 0) phen.subs.edges[["pos"]] <- phen.subs.pos
# if(length(phen.subs.neg) > 0) phen.subs.edges[["neg"]] <- phen.subs.neg
#
# ####################
# ## PLOT to CHECK? ##
# ####################
# ## for SNP1, does it identify the correct/reasonable branches?
# #     edgeCol <- rep("black", nrow(edges))
# #     edgeCol <- replace(edgeCol, phen.subs.edges[["total"]], "green")
# #
# #     ## plot the i'th character's reconstruction on the tree:
# #     #require(adegenet)
# #     plotAnc(tree, phen.pa.ACCTRAN, i=1,
# #             col=transp(c("red", "royalblue"), 0.75),
# #             cex.pie=0.1, pos=NULL,
# #             edge.color=edgeCol, edge.width=2, use.edge.length=FALSE, type="c")
#
#
# ################
# ## Get output ##
# ################
# var.rec <- phen.rec
# subs.edges <- phen.subs.edges
# out <- list("var.rec" = var.rec,
#             "subs.edges" = subs.edges)
#
