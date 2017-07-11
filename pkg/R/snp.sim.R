

#############
## snp.sim ##
#############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param snps description.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#'
#' @examples
#' ## Example ##
#'
#' @import adegenet ape
#' @importFrom Hmisc all.is.numeric

########################################################################

## ARGUMENTS: ##

## n.subs <- either an integer or a vector containing a distribution of n.subs-per-site
## phen.loci <- a vector containing the indices of the edges on which phen subs occurred


snp.sim <- function(n.snps = 10000,
                    n.subs = 1,
                    snp.root = NULL,
                    n.snps.assoc = 0,
                    assoc.prob = 100,
                    tree = coalescent.tree.sim(100),
                    phen.loci = NULL,
                    heatmap = FALSE,
                    reconstruct = FALSE,
                    dist.dna.model = "JC69",
                    row.names = NULL,
                    set = NULL,
                    seed = 1){

  # require(adegenet)
  # require(ape)
  # require(phangorn)

  ##################
  ## HANDLE TREE: ##
  ##################
  ## Always work with trees in "pruningwise" order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  if(!is.rooted(tree)) tree <- midpoint(tree)

  ####################################################################
  ############################
  ## Get Anc-Des EDGE ORDER ##
  ############################
  ## Get sequence from lowest ("root", Nterm+1) to highest ancestral node:
  ix <- c(min(tree$edge[,1]):max(tree$edge[,1]))
  ## Get for loop index of rows in tree$edge[,1], in pairs, from lowest to highest:
  x <- as.vector(unlist(sapply(c(1:length(ix)), function(e) which(tree$edge[,1] == ix[e]))))
  ####################################################################


  ##################################
  ## GET MUTATIONS' branch & loci ##
  ##################################
  n.ind <- tree$Nnode+1
  gen.size <- n.snps
  rm(n.snps)
  edges <- tree$edge

  if(!is.null(seed)) set.seed(seed)

  ## Simulate genotype for root individual: ##

  ## For n.subs = n or = dist approaches:
  # gen.root <- sample(c("a", "c", "g", "t"), gen.size, replace=TRUE)
  gen.root <- sample(c(TRUE, FALSE), gen.size, replace=TRUE)

  ## get the sum of all branch lengths in the tree:
  time.total <- sum(tree$edge.length)

  ## make dummy variables in which to store the resulting n.mts variables:
  L <- lambda <- n.mts <- new.nt <- NA

  snps.assoc <- NULL

  ## if n.snps.assoc is neither NULL nor 0:
  if(is.null(n.snps.assoc)) n.snps.assoc <- 0
  if(n.snps.assoc != 0){

    ## get non.assoc gen.size
    gen.size.ori <- gen.size
    gen.size <- gen.size-n.snps.assoc

    ## assign snps.assoc to be the last n.snps.assoc snps columns
    snps.assoc <- c((gen.size+1):(gen.size+n.snps.assoc))
  }

  ###################
  ## Handle n.subs ##
  ###################

  ## Either an integer
  ## --> draw n.subs from a Poisson distribution w parameter n.subs
  ## OR a vector (containing a distribution)
  ## --> use this distribution to define n.subs-per-site

  if(length(n.subs)==1 & is.null(names(n.subs))){

    #####################
    ## NO DISTRIBUTION ##
    #####################
    ## if no distribution is inputted,
    ## use normal simulation procedure
    ## (ie. Poisson parameter 1):

    warning("Using n.subs as Poisson parameter because input n.subs was of length 1 and had no names.")

    ## draw the number of mutations to occur at each site:
    if(!is.null(seed)) set.seed(seed)
    n.mts <- rpois(n=gen.size, lambda=(n.subs))
    ## for any n.mts==0, re-sample
    if(any(n.mts == 0)){
      ## need to change seed or we'll get trapped in the while loop
      if(!is.null(seed)) seed.i <- seed
      for(i in 1:length(n.mts)){
        while(n.mts[i]==0){
          if(!is.null(seed)){
            seed.i <- seed.i+1
            set.seed(seed.i)
          }
          n.mts[i] <- rpois(n=1, lambda=(n.subs))
        }
      }
    }

    }else{

      ###############################################
      ## DISTRIBUTION or RATES (fitPagel Q matrix) ##
      ###############################################

      ##################
      ## DISTRIBUTION ##
      ##################

      ## if a distribution is provided by the user,
      ## we use this to determine the number of substitutions
      ## to occur at what proportion of sites (note that
      ## we may not be simulating the same number of sites)

      dist <- n.subs

      ## check for names first!
      if(!is.null(names(dist))){
        ## only modify if names are numeric
        if(all.is.numeric(names(dist))){
          noms <- as.numeric(names(dist))
          aligned <- sapply(c(1:length(dist)), function(e) noms[e] == e)
          ## if any names do not correspond to their index, add zeros where missing:
          if(any(aligned == FALSE)){
            dist.new <- rep(0, max(noms))
            dist.new[noms] <- dist
            names(dist.new) <- c(1:length(dist.new))
            dist <- dist.new
          }
        }
      } # end check for missing places

      ## get dist.prop, a distribution containing the counts
      ## of the number of SNPs to be simulated that will have
      ## i many substitutions
      dist.sum <- sum(dist)
      dist.prop <- round((dist/dist.sum)*gen.size)
      ## check that these counts sum to gen.size,
      ## else add the remainder to the largest n.subs count
      ## (OR should we just add these to the n.subs=1 set ??? ###
      ## likely to be the same thing, but could not be...)
      if(sum(dist.prop) != gen.size){
        m <- which.max(dist.prop)
        #m <- 1
        if(sum(dist.prop) < gen.size){
          dist.prop[m] <- dist.prop[m] + (gen.size - sum(dist.prop))
        }
        if(sum(dist.prop) > gen.size){
          dist.prop[m] <- dist.prop[m] - (sum(dist.prop) - gen.size)
        }
      }

      ## get rid of useless trailing 0s
      while(dist.prop[length(dist.prop)] == 0){
        dist.prop <- dist.prop[c(1:(length(dist.prop)-1))]
      }

      ## make n.mts, a vector of length ncol(snps)
      n.mts <- rep(1111, gen.size)
      loci.available <- c(1:gen.size)
      ## assign dist.prop[i] elements of n.mts
      ## to be the same as the n.subs
      ## indicated by i, the given element of dist.prop
      for(j in 1:length(dist.prop)){
        ## provided there are not 0 sites to have this number of substitutions...
        if(dist.prop[j] > 0){
          if(length(loci.available) > 1){
            if(!is.null(seed)) set.seed(seed)
            ## assign dist.prop[i] elements of n.mts to be i
            loci.selected <- sample(loci.available, dist.prop[j], replace = FALSE)
            loci.available <- loci.available[-which(loci.available %in% loci.selected)]
          }else{
            ## if there is only 1 (the last) loci available,
            ## we select this one:
            loci.selected <- loci.available
          }
          n.mts[loci.selected] <- j
        }
      }
      ## Remove unnecessary objects...
      rm(dist)
      rm(dist.prop)
      rm(dist.sum)
      # } # end dist
  } # end fitPagel or dist

  ############################
  ## Assign mts to branches ##
  ############################

  if(n.snps.assoc != 0){
    ## for snps.assoc (the last n.snps.assoc snps, for now),
    ## add n.mts == n.phen.loci s.t these sites mutate at each
    ## and every phen.loci (for now, to be diluted later
    ## according to assoc.prob if !=100)
    n.mts <- c(n.mts, rep(length(phen.loci), n.snps.assoc))
  }

  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####
  #############################################################################
  ## GENERATE ALL SNPs FIRST, THEN REPLACE ANY NON-POLYMORPHIC IN WHILE LOOP ##
  #############################################################################

  #############################
  ## GET NON-ASSOCIATED SNPS ##
  #############################

  ## for each site, draw the branches to which
  ## you will assign the mts for this site
  ## (~ branch length):

  l.edge <- length(tree$edge.length)
  ## Get vector of FALSEs of length tree$edge.length:
  null.vect <- rep(FALSE, l.edge)


  ## TO DO: Memory inefficient step... Improve if possible?
  if(max(n.mts) > l.edge){
    if(!is.null(seed)) set.seed(seed)
      for(e in 1:length(n.mts)){
        repTF <- FALSE
        if(n.mts[e] > l.edge) repTF <- TRUE
        snps.loci[[e]] <- replace(null.vect, sample(c(1:l.edge),
                                                    n.mts[e],
                                                    replace=repTF,
                                                    prob=tree$edge.length), TRUE)
      }
      snps.loci <- t(do.call(rbind, snps.loci))

  }else{
    if(!is.null(seed)) set.seed(seed)
    snps.loci <- sapply(c(1:length(n.mts)),
                        function(e)
                          replace(null.vect,
                          sample(c(1:l.edge),
                                 n.mts[e],
                                 replace=FALSE,
                                 prob=tree$edge.length), TRUE))
  }

  ## rearrange snps.loci s.t it becomes a
  ## list of length tree$edge.length,
  ## each element of which contains the
  ## locations of the mutations that will
  ## occur on that branch
  snps.loci <- sapply(c(1:nrow(snps.loci)),
                       function(e) which(snps.loci[e,] == TRUE))


  ## get the node names for all individuals (terminal and internal)
  all.inds <- sort(unique(as.vector(unlist(tree$edge))))
  # we will store the output in a list called snps:
  snps <- list()
  ## we start w all inds having same genotype as root:
  snps[all.inds] <- rep(list(gen.root), length(all.inds))

  ## store replacement nts in list new.nts:
  new.nts <- list()
  ## distinguish btw list of loci and unique list
  snps.loci.ori <- snps.loci
  ## will need to treat repeat loci differently...
  snps.loci.unique <- lapply(snps.loci, unique)


  #############################
  ## For Loop to get new nts ##
  #############################
  for(i in x){
    ## for all snps other than root, we mutate the
    ## genome of the node preceding it, according to snps.loci.
    ## Draw new nts for each locus selected for mutation:
    if(!.is.integer0(snps.loci.unique[[i]])){
      # new.nts[[i]] <- sapply(c(1:length(snps.loci.unique[[i]])), function(e)
      #   selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
      #                                                  %in% snps[[tree$edge[i,1]]]
      #                                                  [snps.loci.unique[[i]][e]])]))
      new.nts[[i]] <- !snps[[tree$edge[i,1]]][snps.loci.unique[[i]]]


      ## if any loci are selected for multiple mutations
      ## within their given branch length:
      if(length(snps.loci.ori[[i]]) != length(snps.loci.unique[[i]])){
        ## identify which loci are repeaters
        repeats <- table(snps.loci.ori[[i]])[which(table(snps.loci.ori[[i]])!=1)]
        ## how many times they repeat
        n.reps <- repeats - 1
        ## the positions of these loci in the vector of snps loci
        toRepeat <- which(snps.loci.unique[[i]] %in% names(repeats))
        ## run chain of re-sampling to end in our new nt for repeater loci:
        foo <- list()
        for(j in 1:length(toRepeat)){
          foo[[j]] <- new.nts[[i]][toRepeat[j]]
          for(k in 1:n.reps[j]){
            if(k==1){
              # foo[[j]][k] <- selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
              #                                                               %in% foo[[j]][1])])
              foo[[j]][k] <- !foo[[j]][1]

            }else{
              # foo[[j]][k] <- selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
              #                                                               %in% foo[[j]][k-1])])
              foo[[j]][k] <- !foo[[j]][k-1]
            }
          }
          ## retain only the last nt selected
          out <- sapply(c(1:length(foo)),
                        function(e) foo[[e]][length(foo[[e]])])
        }
        ## for the loci with repeated mts, replace these positions
        ## in new.nts with the corresponding elements of out, above.
        new.nts[[i]][toRepeat] <- out
      } # end of if statement for repeaters

      ## update ancestral genotype with new.nts:
      temp <- snps[[tree$edge[i,1]]]
      temp[snps.loci.unique[[i]]] <- new.nts[[i]]
      snps[[tree$edge[i,2]]] <- temp

    }else{
      ## if no mts occur on branch, set genotype of
      ## downstream individual to be equal to ancestor's
      snps[[tree$edge[i,2]]] <- snps[[tree$edge[i,1]]]
    }
  } # end of for loop selecting new nts at mutator loci

  ####################################################
  ## CHECK IF ALL LOCI ARE POLYMORPHIC (|polyThres) ##
  ####################################################

  ## temporarily assemble non-associated loci into matrix:
  # temp.ori <- do.call("rbind", snps)
  temp <- do.call("rbind", snps)

  ## keep only rows containing terminal individuals:
  # temp.ori <- temp.ori[1:n.ind, ]
  temp <- temp[1:n.ind, ]

  ######################################################################################################################################################################

  ## identify n.minor.allele required to meet polyThres:
  polyThres <- 0.01
  n.min <- n.ind*polyThres

  ## make a list of any NON-polymorphic loci:
  csum <- colSums(temp)
  toRepeat <- which(csum < n.min | csum > (nrow(temp) - n.min))







  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####
  #################################################################
  ## REPLACE ANY NON-POLYMORPHIC LOCI & GENERATE ASSOCIATED SNPS ##
  #################################################################
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####

  #################################################
  ## REPLACE NON-POLYMORPHIC NON-ASSOCIATED SNPS ##
  #################################################

  ## just working w rows containing inds (not internal nodes)

  ######################################################################################################################################################################
  ######################################################################################################################################################################

  ######################################################################################################################################################################
  ######################################################################################################################################################################

  ####################
  ## NEW while loop ##
  ####################

  counter <- 0
  if(!is.null(seed)) seed.i <- seed

  while(length(toRepeat) > 0){
    ## for each site, draw the branches to which
    ## you will assign the mts for this site
    ## (~ branch length):

    ## Get vector of FALSEs of length tree$edge.length:
    null.vect <- rep(FALSE, length(tree$edge.length))


    if(max(n.mts[toRepeat]) > length(tree$edge.length)){
      if(!is.null(seed)){
        seed.i <- seed.i+1
        set.seed(seed.i)
      }
      snps.loci <- sapply(c(1:length(n.mts[toRepeat])),
                          function(e)
                            replace(null.vect,
                                    sample(c(1:length(tree$edge.length)),
                                           n.mts[toRepeat][e],
                                           replace=TRUE,
                                           prob=tree$edge.length), TRUE))
    }else{
      if(!is.null(seed)){
        seed.i <- seed.i+1
        set.seed(seed.i)
      }
      snps.loci <- sapply(c(1:length(n.mts[toRepeat])),
                          function(e)
                            replace(null.vect,
                                    sample(c(1:length(tree$edge.length)),
                                           n.mts[toRepeat][e],
                                           replace=FALSE,
                                           prob=tree$edge.length), TRUE))
    }

    ## rearrange snps.loci s.t it becomes a
    ## list of length tree$edge.length,
    ## each element of which contains the
    ## locations of the mutations that will
    ## occur on that branch
    snps.loci <- sapply(c(1:nrow(snps.loci)),
                        function(e) which(snps.loci[e,] == TRUE))


    ## get the node names for all individuals (terminal and internal)
    all.inds <- sort(unique(as.vector(unlist(tree$edge))))
    # we will store the output in a list called snps:
    # snps <- list()
    ## we start w all inds having same genotype as root:
    # snps[all.inds][toRepeat] <- rep(list(gen.root[toRepeat]), length(all.inds))
    for(i in 1:length(all.inds)){
      snps[[all.inds[i]]][toRepeat] <- gen.root[toRepeat]
    }

    ## store replacement nts in list new.nts:
    new.nts <- list()
    ## distinguish btw list of loci and unique list
    snps.loci.ori <- snps.loci
    ## will need to treat repeat loci differently...
    snps.loci.unique <- lapply(snps.loci, unique) # (identical to snps.loci?)


    #############################
    ## For Loop to get new nts ##
    #############################
    for(i in x){
      ## for all snps other than root, we mutate the
      ## genome of the node preceding it, according to snps.loci.
      ## Draw new nts for each locus selected for mutation:
      if(!.is.integer0(snps.loci.unique[[i]])){
        # new.nts[[i]] <- sapply(c(1:length(snps.loci.unique[[i]])), function(e)
        #   selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
        #                                                  %in% snps[[tree$edge[i,1]]][toRepeat]
        #                                                  [snps.loci.unique[[i]][e]])]))
        new.nts[[i]] <- !snps[[tree$edge[i,1]]][toRepeat][snps.loci.unique[[i]]]

        ## if any loci are selected for multiple mutations
        ## within their given branch length:
        if(length(snps.loci.ori[[i]]) != length(snps.loci.unique[[i]])){
          ## identify which loci are repeaters
          repeats <- table(snps.loci.ori[[i]])[which(table(snps.loci.ori[[i]])!=1)]
          ## how many times they repeat
          n.reps <- repeats - 1
          ## the positions of these loci in the vector of snps loci
          toRep <- which(snps.loci.unique[[i]] %in% names(repeats))
          ## run chain of re-sampling to end in our new nt for repeater loci:
          foo <- list()
          for(j in 1:length(toRep)){
            foo[[j]] <- new.nts[[i]][toRep[j]]
            for(k in 1:n.reps[j]){
              if(k==1){
                # foo[[j]][k] <- selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
                #                                                               %in% foo[[j]][1])])
                foo[[j]][k] <- !foo[[j]][1]
              }else{
                # foo[[j]][k] <- selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
                #                                                               %in% foo[[j]][k-1])])
                foo[[j]][k] <- !foo[[j]][k-1]
              }
            }
            ## retain only the last nt selected
            out <- sapply(c(1:length(foo)),
                          function(e) foo[[e]][length(foo[[e]])])
          }
          ## for the loci with repeated mts, replace these positions
          ## in new.nts with the corresponding elements of out, above.
          new.nts[[i]][toRep] <- out
        } # end of if statement for repeaters

        ## update ancestral genotype with new.nts:
        toto <- snps[[tree$edge[i,1]]][toRepeat]
        toto[snps.loci.unique[[i]]] <- new.nts[[i]]
        snps[[tree$edge[i,2]]][toRepeat] <- toto

      }else{
        ## if no mts occur on branch, set genotype of
        ## downstream individual to be equal to ancestor's
        snps[[tree$edge[i,2]]][toRepeat] <- snps[[tree$edge[i,1]]][toRepeat]
      }
    } # end of for loop selecting new nts at mutator loci


    ## temporarily assemble non-associated loci into matrix:
    # temp.new <- do.call("rbind", snps)
    temp <- do.call("rbind", snps)

    ## keep only rows containing terminal individuals:
    # temp <- temp.new[1:n.ind, ]
    temp <- temp[1:n.ind, ]

    ##############################################################
    # temp <- temp.new

    ######################################
    ##### while loop CHECK here: #########
    ######################################
    ## CHECK IF ALL LOCI ARE POLYMORPHIC (|polyThres)

    ## identify n.minor.allele required to meet polyThres:
    polyThres <- 0.01
    n.min <- n.ind*polyThres

    ## make a list of any NON-polymorphic loci:
    toRepeat.ori <- toRepeat
    temp.toRepeat <- temp[, toRepeat.ori]

    ## make a vector of any NON-polymorphic loci:
    ## If ncol = 1:
    if(!is.matrix(temp.toRepeat)){
      csum <- sum(temp.toRepeat)
      if(csum < n.min | csum > (length(temp.toRepeat) - n.min)){
        toRepeat <- toRepeat.ori
      }else{
        toRepeat <- NULL
      }
    }else{
      ## if temp.toRepeat is a true matrix:
      if(ncol(temp.toRepeat) > 0){
        csum <- colSums(temp.toRepeat)
        toRepeat <- toRepeat.ori[which(csum < n.min | csum > (nrow(temp.toRepeat) - n.min))]
      }else{
        csum <- sum(temp.toRepeat)
        if(csum < n.min | csum > (length(temp.toRepeat) - n.min)){
          toRepeat <- toRepeat.ori
        }else{
          toRepeat <- NULL
        }
      }
    }

    counter <- counter+1
    print("COUNTER"); print(counter)
    print("toRepeat"); print(length(toRepeat))


  } # end NEW while loop...
  #######################################
  ### while loop ENDS here: #############
  #######################################

  colnames(temp) <- c(1:ncol(temp))

  ######################################################################################################################################################################
  ######################################################################################################################################################################

  ######################################################################################################################################################################
  ######################################################################################################################################################################

  # # toRepeat <- 1:length(n.mts)
  # # loci <- list()
  #
  # ######################################
  # ## while loop STARTS here: ###########
  # ######################################
  # ## AGAIN--NEED TO DOUBLE CHECK: No problems with seed? #############
  #
  # counter <- 0
  # if(!is.null(seed)) seed.i <- seed
  # system.time(
  # while(length(toRepeat) > 0){
  #
  # for(i in toRepeat){
  #   ## get the lth element of n.mts to work with:
  #   n.mt <- n.mts[i]
  #
  #   ## for each site, draw the branches to which
  #   ## you will assign the mts for this site
  #   ## (~ branch length):
  #
  #   if(n.mt > length(tree$edge.length)){
  #     if(!is.null(seed)){
  #       seed.i <- seed.i+1
  #       set.seed(seed.i)
  #     }
  #     subs.edges <- sample(c(1:length(tree$edge.length)),
  #                          n.mt,
  #                          replace=TRUE,
  #                          prob=tree$edge.length)
  #   }else{
  #     if(!is.null(seed)){
  #       seed.i <- seed.i+1
  #       set.seed(seed.i)
  #     }
  #     subs.edges <- sample(c(1:length(tree$edge.length)),
  #                          n.mt,
  #                          replace=FALSE,
  #                          prob=tree$edge.length)
  #   }
  #
  #   ## TO DO: COULD REPLACE (all instances!!!) LATER WITH:
  #   ## Get vector of FALSEs of length tree$edge.length:
  #   #     null.vect <- rep(FALSE, length(tree$edge.length))
  #   #     subs.edges <- replace(null.vect,
  #   #                           sample(c(1:length(tree$edge.length)),
  #   #                                  n.mt,
  #   #                                  replace=FALSE,
  #   #                                  prob=tree$edge.length), TRUE)
  #
  #
  #   ## get nt for root at this locus:
  #   root.nt <- gen.root[i]
  #
  #   ## get nt for each individual at this locus
  #   temp[,i] <- .get.locus(subs.edges = subs.edges,
  #                           root.nt = root.nt,
  #                           tree = tree)[1:n.ind]
  #
  # } # end FOR LOOP for NON-associated SNPs
  #
  # ######################################
  # ##### while loop CHECK here: #########
  # ######################################
  # ## CHECK IF ALL LOCI ARE POLYMORPHIC (|polyThres)
  #
  # ## identify n.minor.allele required to meet polyThres:
  # polyThres <- 0.01
  # n.min <- n.ind*polyThres
  #
  # ## make a list of any NON-polymorphic loci:
  # toRepeat.ori <- toRepeat
  # temp.toRepeat <- temp[, toRepeat.ori]
  #
  #
  # toRepeat <- list()
  # ## if temp.toRepeat is a true matrix:
  # if(!is.matrix(temp.toRepeat)){
  #   if(any(table(temp.toRepeat) < n.min) | length(table(temp.toRepeat)) == 1){
  #     toRepeat[[length(toRepeat)+1]] <- toRepeat.ori
  #   }
  # }else{
  #   if(ncol(temp.toRepeat) > 0){
  #     for(i in 1:ncol(temp.toRepeat)){
  #       if(any(table(temp.toRepeat[,i]) < n.min) | length(table(temp.toRepeat[,i])) == 1){
  #         toRepeat[[length(toRepeat)+1]] <- toRepeat.ori[i]
  #       }
  #     }
  #   }else{
  #     if(any(table(temp.toRepeat) < n.min) | length(table(temp.toRepeat)) == 1){
  #       toRepeat[[length(toRepeat)+1]] <- toRepeat.ori
  #     }
  #   }
  # }
  # if(length(toRepeat) > 0){
  #   toRepeat <- as.vector(unlist(toRepeat))
  # }
  #
  # counter <- counter+1
  # print("COUNTER"); print(counter)
  # print("toRepeat"); print(length(toRepeat))
  #
  # } # end of while loop
  # )
  # # ######################################
  # # ## while loop ENDS here: #############
  # # ######################################
  #
  #
  # ## GET ALL NON-UNIQUE SNPS COLUMNS: ##
  #
  # if(all.unique == TRUE){
  #   temp.complete <- temp
  # }else{
  #   temp.complete <- temp[, index]
  # }






  #########################
  ## GET ASSOCIATED SNPS ##
  #########################

  ## Need to treat ASSOCIATED SNPs differently:
  ## (non.assoc.snps do NOT need to pass the "while" check;
  ## they just need to match phen.loci at this point.)
  if(n.snps.assoc != 0){
    ## get snps.loci for the ASSOCIATED snps (ie. set to phen.loci) ##
    for(i in 1:n.snps.assoc){
      ## recall: phen.loci contains the tree EDGES on which phen subs occur
      subs.edges <- phen.loci

      ## get nt for root at this locus:
      root.nt <- gen.root[snps.assoc[i]]

      ## get nt for each individual at this locus
      ## assign to (and replace) the snps.assoc elements of loci
      temp[, snps.assoc[i]] <- .get.locus(subs.edges = subs.edges,
                                          root.nt = root.nt,
                                          tree = tree)[1:n.ind]
    }
  } # end of snps.assoc generation


  ###########################################
  ## GET COMPLETE SNPS MATRIX ("snps"): ##
  ###########################################

  ## Remove unnecessary objects...
  rm(snps)
  ## Create snps matrix:
  snps <- temp
  ## Remove unnecessary objects...
  rm(temp)


  ## keep only rows containing terminal individuals:
  snps <- snps[1:n.ind, ]

  ## clear memory
  gc()


  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####


  ###############################################
  ## MODIFY SNPS.ASSOC ACCORDING TO ASSOC.PROB ##
  ###############################################

  if(n.snps.assoc != 0){
    ## if we have any imperfect associations... ##
    if(any(assoc.prob != 100)){
      ## check length
      if(length(assoc.prob) != n.snps.assoc){
        ## if only 1 prob value given...
        if(length(assoc.prob) == 1){
          ## ... assume uniform assoc.prob;
          assoc.prob <- rep(assoc.prob, n.snps.assoc)
          ## no warning needed
        }else{
          ## BUT if assoc.prob of random length:
          ## repeat until of length n.snps.assoc
          assoc.prob <- rep(assoc.prob, length.out=n.snps.assoc)
          ## and print warning (only if not of length n.snps.assoc OR 1)
          warning("assoc.prob not of length n.snps.assoc;
                  sequence will be repeated until correct length is reached.")
        }
        } # end checks

      ## for each associated SNP,
      ## we undo some associations | assoc.prob for that snp.assoc
      for(i in 1:n.snps.assoc){

        ## re-pseudo-randomise seed:
        if(!is.null(seed)){
          seed.i <- seed*i*10
          set.seed(seed.i)
        }

        prob <- assoc.prob[i]
        ## only if the association is imperfect
        if(prob != 100){
          ## draw snps to change at snps.assoc[i]
          n.toChange <- round(nrow(snps)*(1 - (prob/100)))
          toChange <- sample(c(1:nrow(snps)), n.toChange)

          ## change those snps at rows toChange, loci snps.assoc[i]
          for(j in 1:length(toChange)){
            snps[toChange[j], snps.assoc[i]] <-
              selectBiallelicSNP(snps[toChange[j], snps.assoc[i]])
          } # end for loop
        }
      } # end for loop
    } # end any assoc.prob != 100
  } # end modification | assoc.prob


  ##############################
  ## PLOTS & TREECONSTRUCTION ##
  ##############################
  tree.reconstructed <- NULL
  if(heatmap == TRUE || reconstruct!=FALSE){

    ## CONVERT TO CHARACTER: ##
    dna <- snps
    dna <- replace(dna, which(dna == TRUE), "a")
    dna <- replace(dna, which(dna == "FALSE"), "t")

    ## Get DNAbin object:
    dna <- as.DNAbin(dna)
    # rownames(dna) <- c(1:nrow(snps))


    #############
    ## HEATMAP ##
    #############
    if(heatmap==TRUE){
      heatmap.DNAbin(dna=dna,
                     dist.dna.model=dist.dna.model)
    }

    ##########################################
    ## PLOT 2: RECONSTRUCTING THE PHYLOGENY ##
    ##########################################
    if(reconstruct!=FALSE){
      tree.reconstructed <- tree.reconstruct(dna, # [1:n.ind,]
                                           method=reconstruct,
                                           dist.dna.model=dist.dna.model,
                                           plot=TRUE)
    }

    ## Remove unnecessary object:
    rm(dna)

  } # end heatmap, treeconstruction

  ##################
  ## CONVERT SNPS ##
  ##################

  ## CONVERT TO NUMERIC: ##
  ## Convert from logical to binary SNPs (for terminal nodes only):
  snps <- replace(snps, which(snps == TRUE), 1)

  ## Keep only rows containing terminal individuals?:
  ## (NOTE -- Consider moving this to AFTER snps.assoc assoc.prob section!)
  # snps <- snps[1:n.ind, ]
  #
  # ## NO LONGER NEED THIS SECTION? ##
  # gen.size <- ncol(snps)
  #
  ## get snps as DNAbin ## DON'T THINK WE NEED THIS ANYMORE...
  # ploidy <- 1
  # snps <- as.DNAbin(snps, ploidy=ploidy)
  # ## get snps as genind
  # #source("C:/Users/caitiecollins/adegenet/R/sequences.R")
  # snps.gen <- DNAbin2genind(snps, polyThres=0.01)
  # ## get snps as binary matrix
  # snps <- snps.gen@tab
  #
  # ## correct genind for ploidy:
  # snps <- snps[,seq(1, ncol(snps), 2)]
  #
  # ## assign snps row and column names
  # ## (TEMPORARILY?? -- DO WE NEED To DO ThIS HERE??) ### ###### #############     #
  # rownames(snps) <- 1:nrow(snps)
  # colnames(snps) <- 1:ncol(snps)


  ## Reassort snps.assoc to new columns:
  if(!is.null(snps.assoc)){

    ## update snps.assoc to reflect true loci
    gen.size.final <- ncol(snps)
    snps.assoc.loci.ori <- c((gen.size.final-(n.snps.assoc-1)):gen.size.final)

    #########################################
    ## RANDOMIZE SNPS.ASSOC LOCI POSITIONS ##
    #########################################

    ## Re-enabled snps.assoc loci "randomization" by
    ## just drawing indices and shuffling the columns accordingly...
    ## draw which SNPs will be associated to the phenotype
    if(!is.null(seed)) set.seed(seed)
    snps.assoc.loci <- sort(sample(c(1:gen.size.final),
                                   n.snps.assoc,
                                   replace=FALSE))

    snps.indices <- c(1:gen.size.final)
    snps.ori <- snps

    snps.non.assoc <- snps[,c(1:(gen.size.final-n.snps.assoc))]
    snps.assoc <- snps[,snps.assoc.loci.ori]
    snps.new <- matrix(99, nrow=nrow(snps), ncol=gen.size.final)
    snps.new[,snps.indices[-snps.assoc.loci]] <- snps.non.assoc
    snps.new[,snps.assoc.loci] <- snps.assoc
    snps <- snps.new
    snps.assoc <- snps.assoc.loci

  } # end snps.assoc randomization

  ###############################
  ## Assign row & column names ##
  ###############################

  ## assign/generate row.names
  if(!is.null(row.names)){
    ## If row.names have been provided in args list, assign them:
    if(length(row.names) == nrow(snps)){
      ## match tree$tip.label?
      if(!is.null(tree$tip.label)){
        if(all(row.names %in% tree$tip.label) & all(tree$tip.label %in% row.names)){
          ## REORDER to match tree$tip.labs if possible:
          if(!identical(row.names, tree$tip.label)){
            ord <- match(tree$tip.label, row.names)
            row.names <- row.names[ord]
          }
        }
      }
      rownames(snps) <- row.names
    }else{
      if(is.null(rownames(snps))) rownames(snps) <- c(1:nrow(snps))
    }
  }else{
    ## Otherwise, try to assign rownames(snps) to match tree$tip.label:
    if(!is.null(tree$tip.label)){
      if(length(tree$tip.label) == nrow(snps)){
        rownames(snps) <- tree$tip.label
      }else{
        rownames(snps) <- 1:nrow(snps)
        warning("The length of tree$tip.label was not equal to nrow(snps) being simulated;
              rownames(snps) have been set to 1:N and will not match tree$tip.label.")
      }
    }
  }

  ## generate column names:
  colnames(snps) <- 1:ncol(snps)


  #################################################
  ## SIM SET 2 (complementary clade-wise assoc): ##
  #################################################
  sets <- NULL
  if(!is.null(snps.assoc)){
    if(!is.null(set)){
      if(set == 2){

        ## Want to divide tree into 2 sets of clades btw 1/3:2/3 and 1/2:1/2
        clades <- tab <- grp.options <- sets.complete <- list()

        min.size <- ceiling((tree$Nnode+1)*(1/3))
        max.size <- floor((tree$Nnode+1)*(2/3))
        grp1 <- tree$Nnode+1

        ## get 2 sets of clades:

        ##########################################
        ## Pick SETS for ANY TREE (pretty sure) ##
        ##########################################
        dec <- grp <- sets.temp <- sets.complete <- list()

        inds <- c(1:(tree$Nnode+1))
        # new.root <- tree$edge[1,1] # initial root
        new.root <- tree$Nnode+2 # initial root

        counter <- 0
        #######################################
        ## WHILE LOOP to get size of clades: ##
        #######################################
        while(grp1 < min.size | grp1 > max.size){

          ## get all descendants of root node:
          all.dec <- .getDescendants(tree, node=new.root)

          ## get all descendants in first 2 major clades:
          dec[[1]] <- .getDescendants(tree, node=all.dec[1])
          dec[[2]] <- .getDescendants(tree, node=all.dec[2])

          ## get terminal inds only:
          sets.temp[[1]] <- dec[[1]][which(dec[[1]] %in% inds)]
          sets.temp[[2]] <- dec[[2]][which(dec[[2]] %in% inds)]

          grp[[1]] <- length(sets.temp[[1]])
          grp[[2]] <- length(sets.temp[[2]])

          max.grp <- which.max(c(grp[[1]], grp[[2]]))
          new.root <- all.dec[max.grp]

          set1 <- sets.temp[[max.grp]]

          sets <- rep(2, length(inds))
          sets <- replace(sets, set1, 1)
          names(sets) <- tree$tip.label

          counter <- counter+1

          grp1 <- grp[[max.grp]]

        } # end while loop


        set1 <- names(sets)[which(sets == 1)]
        set2 <- names(sets)[which(sets == 2)]
        ###########

        ########################
        ## MODIFY SNPS.ASSOC: ##
        ########################
        snps.assoc.set1 <- 1:round(length(snps.assoc)/2)
        snps.assoc.set2 <- (round(length(snps.assoc)/2)+1):length(snps.assoc)

        ## replace set1 snps with 0 at all inds in clade.set1:
        for(e in 1:length(snps.assoc.set1)){
          snps[which(rownames(snps) %in% set1), snps.assoc[snps.assoc.set1[e]]] <- 0
        }
        ## replace set2 snps with 0 at all inds in clade.set2:
        for(e in 1:length(snps.assoc.set2)){
          snps[which(rownames(snps) %in% set2), snps.assoc[snps.assoc.set2[e]]] <- 0
        }
      }
    }
  } # end sim set 2


  ##################
  ## get RESULTS: ##
  ##################
  out <- list(snps, snps.assoc, tree.reconstructed, sets)
  names(out) <- c("snps", "snps.assoc", "tree.reconstructed", "sets")

  return(out)

} # end snp.sim














################
## .get.locus ##
################
.get.locus <- function(subs.edges, root.nt, tree){

  ##################
  ## HANDLE TREE: ##
  ##################
  ## Always work with trees in "pruningwise" order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  if(!is.rooted(tree)) tree <- midpoint(tree)

  ####################################################################
  ############################
  ## Get Anc-Des EDGE ORDER ##
  ############################
  ## Get sequence from lowest ("root", Nterm+1) to highest ancestral node:
  ix <- c(min(tree$edge[,1]):max(tree$edge[,1]))
  ## Get for loop index of rows in tree$edge[,1], in pairs, from lowest to highest:
  x <- as.vector(unlist(sapply(c(1:length(ix)), function(e) which(tree$edge[,1] == ix[e]))))
  ####################################################################


  ## convert subs.edges into appropriate format:
  snps.loci <- list()
  snps.loci[[1]] <- subs.edges

  ## rearrange snps.loci s.t it becomes a
  ## list of length tree$edge.length,
  ## each element of which contains the
  ## locations of the mutations that will
  ## occur on that branch
  snps.loci <- sapply(c(1:length(tree$edge.length)),
                      function(f)
                        seq_along(snps.loci)[sapply(snps.loci,
                                                    function(e) f %in% e)])


  # we will store the output in a list called locus:
  locus <- list()
  ## get the node names for all individuals (terminal and internal)
  all.inds <- sort(unique(as.vector(unlist(tree$edge))))
  ## we start w all inds having same genotype as root:
  for(j in all.inds){
    locus[[j]] <- root.nt
  }

  ## store replacement nts in list new.nts:
  new.nts <- list()
  ## distinguish btw list of loci and unique list
  snps.loci.ori <- snps.loci
  ## will need to treat repeat loci differently...
  snps.loci.unique <- lapply(snps.loci, unique)


  #############################
  ## For Loop to get new nts ##
  #############################
  for(i in x){
    ## for all locus other than root, we mutate the
    ## genome of the node preceding it, according to snps.loci.
    ## Draw new nts for each locus selected for mutation:
    if(!.is.integer0(snps.loci.unique[[i]])){
      new.nts[[i]] <- sapply(c(1:length(snps.loci.unique[[i]])), function(e)
        selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
                                                       %in% locus[[tree$edge[i,1]]]
                                                       [snps.loci.unique[[i]][e]])]))
      ## if any loci are selected for multiple mutations
      ## within their given branch length:
      if(length(snps.loci.ori[[i]]) != length(snps.loci.unique[[i]])){
        ## identify which loci are repeaters
        repeats <- table(snps.loci.ori[[i]])[which(table(snps.loci.ori[[i]])!=1)]
        ## how many times they repeat
        n.reps <- repeats - 1
        ## the positions of these loci in the vector of snps loci
        toRepeat <- which(snps.loci.unique[[i]] %in% names(repeats))
        ## run chain of re-sampling to end in our new nt for repeater loci:
        foo <- list()
        for(j in 1:length(toRepeat)){
          foo[[j]] <- new.nts[[i]][toRepeat[j]]
          for(k in 1:n.reps[j]){
            if(k==1){
              foo[[j]][k] <- selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
                                                                            %in% foo[[j]][1])])

            }else{
              foo[[j]][k] <- selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
                                                                            %in% foo[[j]][k-1])])
            }
          }
          ## retain only the last nt selected
          out <- sapply(c(1:length(foo)),
                        function(e) foo[[e]][length(foo[[e]])])
        }
        ## for the loci with repeated mts, replace these positions
        ## in new.nts with the corresponding elements of out, above.
        new.nts[[i]][toRepeat] <- out
      } # end of if statement for repeaters

      ## update ancestral genotype with new.nts:
      temp <- locus[[tree$edge[i,1]]]
      temp[snps.loci.unique[[i]]] <- new.nts[[i]]
      locus[[tree$edge[i,2]]] <- temp

    }else{
      ## if no mts occur on branch, set genotype of
      ## downstream individual to be equal to ancestor's
      locus[[tree$edge[i,2]]] <- locus[[tree$edge[i,1]]]
    }
  } # end of for loop selecting new nts at mutator loci

  ## turn locus into a vector for easier post-handling
  locus <- as.vector(unlist(locus))

  return(locus)

} # end .get.locus
