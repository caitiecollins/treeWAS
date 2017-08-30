
## WARNING:
## SOME LINES OF SNP.SIM.Q.R HAVE FALLEN OUT OF DATE w SNP.SIM.R.
## BEFORE USING SNP.SIM.Q, UPDATE (at least!):
## sapply, add n.mts > l.edge check,
## convert snps to logical, for loop (selectBiallelicSNP --> !l[[i]])

#############
## snp.sim ##
#############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Aternative SNPs simulation fn.
#'
#' NOT currently in use. Please use the regular snp.sim function to simulate genetic data.
#'
#' @param n.snps An integer specifying the number of snps columns to be simulated.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#'
#' @examples
#' ## Example ##
#'
#' @import adegenet ape
#' @importFrom Hmisc all.is.numeric
#' @importFrom phangorn midpoint
#'
#' @export

########################################################################
#  @useDynLib phangorn, .registration = TRUE

## ARGUMENTS: ##

## n.subs <- either an integer or a vector containing a distribution of n.subs-per-site
## phen.loci <- a vector containing the indices of the edges on which phen subs occurred

#########
# tree <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_31_tree.Rdata"))

# ## ARGS (Eg): ##
# n.snps = 1000
# n.subs = 1
# snp.root = NULL
# n.snps.assoc = 10
# assoc.prob = 100
# # ## dependent/corr' transition rate/prob mat:
# # # Q = matrix(c(2, 0.75, 0.75, 1, 3, 0.5, 0.25,  3, 3, 0.25, 0.5, 3, 1, 0.75, 0.75, 2),
# # #            nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
# ### tree = coalescent.tree.sim(100)
# phen.loci = NULL
# n.phen.subs <- 15
# heatmap = FALSE
# reconstruct = FALSE
# dist.dna.model = "JC69"
# grp.min <- 0.25
# row.names = NULL
# set=3
# seed=1

snp.sim.Q <- function(n.snps = 10000,
                      n.subs = 1,
                      snp.root = NULL,
                      n.snps.assoc = 10,
                      assoc.prob = 100,
                      ## dependent/corr' transition rate/prob mat:
                      Q = matrix(c(2, 0.75, 0.75, 1, 3, 0.5, 0.25,  3, 3, 0.25, 0.5, 3, 1, 0.75, 0.75, 2),
                                 nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2)),
                      tree = coalescent.tree.sim(100),
                      n.phen.subs = 15,
                      phen.loci = NULL,
                      heatmap = FALSE,
                      reconstruct = FALSE,
                      dist.dna.model = "JC69",
                      grp.min = 0.25,
                      row.names = NULL,
                      set=3,
                      seed=1){

  # require(adegenet)
  # require(ape)


  ## HANDLE TREE: ##
  ## Always work with trees in "pruningwise" order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  # if(!is.rooted(tree)) tree <- midpoint(tree)

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
  n.ind <- min(tree$edge[,1])-1 # tree$Nnode+1
  gen.size <- n.snps
  edges <- tree$edge

  if(!is.null(seed)) set.seed(seed)

  ## Simulate genotype (& phen) for root individual: ##

  ## if snp.root given:
  if(!is.null(snp.root)){
    if(length(snp.root) == 1){
      ## select only root state --> different SNP sim method (???)
      if(snp.root %in% c(0, FALSE)) gen.root <- rep(FALSE, gen.size)
      if(snp.root %in% c(1, TRUE)) gen.root <- rep(TRUE, gen.size)
    }else{
      ## if snp.root provided for all loci:
      if(length(snp.root) == gen.size){
        if(length(unique(snp.root[!is.na(snp.root)])) == 2){
          gen.root <- as.logical(as.numeric(as.factor(snp.root))-1)
        }else{
          warning("snp.root must be binary; ignoring.")
        }
      }else{
        warning("snp.root should either be of length 1 or length n.snps; ignoring.")
      }
    }
  }

  ## For n.subs = n or = dist approaches:
  gen.root <- sample(c(TRUE, FALSE), gen.size, replace=TRUE)
  phen.root <- sample(c(TRUE, FALSE), 1)

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
    n.mts <- c(n.mts, rep(1, n.snps.assoc))
    # n.mts <- c(n.mts, rep(length(phen.loci), n.snps.assoc))
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
              foo[[j]][k] <- !foo[[j]][1]

            }else{
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

    # print("LENGTH toREPEAT"); print(length(toRepeat))
    # print("toRepeat NAs?"); print(length(which(is.na(toRepeat))))
    # print("n.mts NAs?"); print(length(which(is.na(n.mts))))
    # print("str(n.mts)"); print(str(n.mts))
    # print("n.mts[toRepeat]"); print(n.mts[toRepeat])
    # print("max(n.mts[toRepeat])"); print(max(n.mts[toRepeat]))
    # print("length(tree$edge.length)"); print(length(tree$edge.length))

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
                foo[[j]][k] <- !foo[[j]][1]
              }else{
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
    # print("COUNTER"); print(counter)
    # print("toRepeat"); print(length(toRepeat))


  } # end NEW while loop...
  #######################################
  ### while loop ENDS here: #############
  #######################################

  colnames(temp) <- c(1:ncol(temp))

  # print("WHILE LOOP DONE; toRepeat?"); print(length(toRepeat))
  ######################################################################################################################################################################
  ######################################################################################################################################################################



  #########################
  ## GET ASSOCIATED SNPS ##
  #########################

  ## Need to treat ASSOCIATED SNPs differently:
  ## (non.assoc.snps do NOT need to pass the "while" check;
  ## they just need to match phen.loci at this point.)

  snps.assoc.nodes <- phen.nodes <- NULL

  if(n.snps.assoc != 0){

    snps.assoc.nodes <- list()
    N.OVERLAP <- list()

    Qt <- list()

    ## get snps.loci for FIRST ASSSOCIATED snp WITH PHEN.loci:
    i <- 1

    edges <- tree$edge

    root.nt <- gen.root[snps.assoc[i]]


    ## Need to update tree$edge.length probs to mirror sampling without replacement...
    ## BUT only remove branches from the calculation once they have had a sub on them...?
    N.SUBS.TOTAL <- 15
    # N.SUBS.TOTAL <- rpois(1, n.phen.subs)  ############

    ################
    ## WHILE LOOP ##
    ################
    toRepeat <- TRUE
    ## WHILE LOOP STARTS HERE:
    while(toRepeat){ ############################

      snp.node <- phen.node <- as.list(rep(NA, length(unique(as.vector(edges)))))

      snp.node[[edges[x[1], 1]]] <- root.nt
      phen.node[[edges[x[1], 1]]] <- phen.root

      # probs.e <- round((tree$edge.length/sum(tree$edge.length)),3)
      # probs.e <- (tree$edge.length/sum(tree$edge.length))
      probs.e <- NULL

      N.SUBS.COUNTER <- 0

      ##################
      ## FOR (e) LOOP ##
      ##################
      ## go from last to first edge in edges:
      for(e in x){

        # print("E"); print(e)

        probs <- NULL
        toKeep <- toSub <- NULL

        ## WANT -- 15 subs/tree (thus 183 no-sub branches): ##

        ####################################################
        ## get conditional probs for each edge w matexpo! ##
        ####################################################
        ## (run within code...)
        Qt[[e]] <- matexpo(Q*tree$edge.length[e])

        P <- Qt[[e]]
        rownames(P) <- rownames(Qt[[e]]) <- rownames(Q)
        colnames(P) <- colnames(Qt[[e]]) <- colnames(Q)

        ##############################
        ## SNP.anc = 0 Phen.anc = 0 ##
        ##############################
        if(snp.node[[edges[e, 1]]] == FALSE & phen.node[[edges[e, 1]]] == FALSE){
          probs <- P[1,]
        }
        ##############################
        ## SNP.anc = 0 Phen.anc = 1 ##
        ##############################
        if(snp.node[[edges[e, 1]]] == FALSE & phen.node[[edges[e, 1]]] == TRUE){
          probs <- P[2,]
        }
        ##############################
        ## SNP.anc = 1 Phen.anc = 0 ##
        ##############################
        if(snp.node[[edges[e, 1]]] == TRUE & phen.node[[edges[e, 1]]] == FALSE){
          probs <- P[3,]
        }
        ##############################
        ## SNP.anc = 1 Phen.anc = 1 ##
        ##############################
        if(snp.node[[edges[e, 1]]] == TRUE & phen.node[[edges[e, 1]]] == TRUE){
          probs <- P[4,]
        }

        SP.dec <- sample(colnames(Q), 1, prob = probs)

        S.dec <- as.logical(as.numeric(keepFirstN(SP.dec, 1)))
        P.dec <- as.logical(as.numeric(keepLastN(SP.dec, 1)))
        names(S.dec) <- names(P.dec) <- NULL

        snp.node[[edges[e, 2]]] <- S.dec
        phen.node[[edges[e, 2]]] <- P.dec

        ## Did a PHEN sub occur?
        phen.sub <- FALSE
        if(phen.node[[edges[e, 1]]] != phen.node[[edges[e, 2]]]) phen.sub <- TRUE
        # SP.dec %in% colnames(Q)[toSub] ## Did EITHER a phen sub AND/OR a snps sub occur?

        # print("EDGE"); print(e)
        ## If a sub has occurred on branch e, add it to the n.subs counter
        if(phen.sub){
          N.SUBS.COUNTER <- N.SUBS.COUNTER+1
          # print("TRUE")
        }
        # print("N.SUBS.COUNTER"); print(N.SUBS.COUNTER)

      } # end for (e) loop


      ## STORE FIRST SNPS.ASSOC:
      snps.assoc.nodes[[i]] <- as.vector(unlist(snp.node))

      ## STORE PHEN (FOR ALL NODES):
      phen.nodes <- as.vector(unlist(phen.node))
      names(phen.nodes) <- c(1:length(phen.nodes))

      phen <- phen.nodes[1:n.ind]

      # ## TEMP -- CHECK w PLOT:
      # plot.phen(tree, phen.nodes=phen.nodes, snp.nodes=snps.assoc.nodes[[i]])
      # cor(as.numeric(phen.nodes[1:n.ind]), as.numeric(snps.assoc.nodes[[i]][1:n.ind]))
      N.overlap <- length(which(phen.nodes[1:n.ind] == snps.assoc.nodes[[i]][1:n.ind]))
      N.overlap <- max(N.overlap, (n.ind-N.overlap))
      N.OVERLAP[[i]] <- N.overlap
      #######################################################################


      ## CHECK THAT MIN GRP.SIZE >= THRESHOLD ##
      if(!is.null(grp.min)){
        tab <- table(phen)
        grp.thresh <- (n.ind)*grp.min
        if(min(tab) < grp.thresh){
          toRepeat <- TRUE
        }else{
          toRepeat <- FALSE
        }
      }else{
        toRepeat <- FALSE
      }

    } # end WHILE LOOP #########

    print("N.SUBS.COUNTER"); print(N.SUBS.COUNTER)
    print("N.overlap"); print(N.overlap[[i]])

    ### TEMP -- CHECK:
    # phen.node.ori <- phen.node
    # phen.rec <- as.vector(unlist(phen.node))
    # snp.rec <- as.vector(unlist(snps.assoc.nodes[[3]]))
    # plot.phen(tree, phen.nodes=phen.rec, snp.nodes=snp.rec)
    # title("set3_1 phen vs. snps.assoc3", line=0)
    #
    # cor(as.numeric(snp.rec[1:n.ind]), as.numeric(phen.rec[1:n.ind])) # 0.46 0.63 0.48
    # length(which(as.numeric(snp.rec[1:n.ind]) == as.numeric(phen.rec[1:n.ind])))/n.ind # 0.74 0.81 0.75
    ######

    ## get snps.loci for the REMAINING ASSOCIATED snps (ie. 2:n.assoc, conditional on phen) ##
    for(i in 2:n.snps.assoc){
      ## get nt for root at this locus:
      root.nt <- gen.root[snps.assoc[i]]

      snp.node <- as.list(rep(NA, length(unique(as.vector(edges)))))

      snp.node[[edges[x[1], 1]]] <- root.nt

      ## go from last to first edge in edges:
      for(e in x){

        probs <- NULL

        P <- Qt[[e]]

        if(snp.node[[edges[e, 1]]] == FALSE & phen.node[[edges[e, 1]]] == FALSE) probs <- P[1,]
        if(snp.node[[edges[e, 1]]] == FALSE & phen.node[[edges[e, 1]]] == TRUE) probs <- P[2,]
        if(snp.node[[edges[e, 1]]] == TRUE & phen.node[[edges[e, 1]]] == FALSE) probs <- P[3,]
        if(snp.node[[edges[e, 1]]] == TRUE & phen.node[[edges[e, 1]]] == TRUE) probs <- P[4,]

        ## Now we KNOW, A PRIORI, the phen.node for the DESCENDANT!
        probs.mod <- replace(probs, which(!as.logical(as.numeric(keepLastN(colnames(Q), 1))) == phen.node[[edges[e, 2]]]), 0)
        SP.dec <- sample(colnames(Q), 1, prob = probs.mod)

        S.dec <- as.logical(as.numeric(keepFirstN(SP.dec, 1)))
        names(S.dec) <- NULL

        snp.node[[edges[e, 2]]] <- S.dec

      } # end for (e) loop

      ## STORE SNPS.ASSOC (FOR ALL NODES):
      snps.assoc.nodes[[i]] <- as.vector(unlist(snp.node))

      ## Get proportion overlap btw phen and snps.assoc.i:
      N.overlap <- length(which(phen.nodes[1:n.ind] == snps.assoc.nodes[[i]][1:n.ind]))
      N.overlap <- max(N.overlap, (n.ind-N.overlap))
      N.OVERLAP[[i]] <- N.overlap

    } # end for loop

    ## Bind SNPs.ASSOC into matrix:
    snps.assoc.nodes <- do.call("cbind", snps.assoc.nodes)


    ###################################################
    ## TEMP -- COMPARE PHEN & ALL SNPS.ASSOC w PLOT: ##
    ###################################################
    par(mfrow=c(2,6))
    plot.phen(tree, phen.nodes=phen.nodes, main.title="phen")
    for(i in 1:5){
      plot.phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE,
                main.title=paste("snp.assoc", i, sep=" "))
      title(N.OVERLAP[[i]], line=0, font.main=1)
    }
    plot.phen(tree, phen.nodes=phen.nodes, main.title="phen")
    for(i in 6:10){
      plot.phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE,
                main.title=paste("snp.assoc", i, sep=" "))
      title(N.OVERLAP[[i]], line=0, font.main=1)
    }


    par(mfrow=c(1,1)) # end temp panel plot

    gc()

  } # end of snps.assoc generation



  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###


  ###########################################
  ## GET COMPLETE SNPS MATRIX ("genomes"): ##
  ###########################################

  ## Create snps matrix:
  snps <- temp
  rm(temp)

  ## Attach snps.assoc loci to last column:
  if(!is.null(snps.assoc.nodes)){
    snps <- cbind(snps[,c(1:(ncol(snps)-(n.snps.assoc)))], snps.assoc.nodes[1:n.ind,])
  }


  ## keep only rows containing terminal individuals:
  snps <- snps[1:n.ind, ]

  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####


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

        min.size <- ceiling((n.ind)*(1/3))
        max.size <- floor((n.ind)*(2/3))
        grp1 <- n.ind

        ## get 2 sets of clades:

        ##########################################
        ## Pick SETS for ANY TREE (pretty sure) ##
        ##########################################
        dec <- grp <- sets.temp <- sets.complete <- list()

        inds <- c(1:(n.ind))
        # new.root <- tree$edge[1,1] # initial root
        new.root <- n.ind+1 # initial root

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
  out <- list(snps, snps.assoc, tree.reconstructed, sets, phen, phen.nodes)
  names(out) <- c("snps", "snps.assoc", "tree.reconstructed", "sets", "phen", "phen.nodes")

  return(out)

} # end snp.sim.Q







#####################
## .getDescendants ##
#####################
## get clades --> numbers of descendants in successive major clades:
## phytools (?) fn from Revell blog:
## http://blog.phytools.org/2012/01/function-to-get-descendant-node-numbers.html
.getDescendants<-function(tree,node,curr=NULL){

  if(is.null(curr)) curr <- vector()

  daughters <- tree$edge[which(tree$edge[,1]==node),2]
  curr <- c(curr,daughters)

  w <- which(daughters>=length(tree$tip))
  if(length(w) > 0) for(i in 1:length(w))
    curr <- .getDescendants(tree, daughters[w[i]], curr)

  return(curr)
} # end .getDescendants




#
#
# ##################
# ## .get.locus01 ##
# ##################
#
# # ########################################################################
# #
# # ###################
# # ## DOCUMENTATION ##
# # ###################
# #
# # #' Short one-phrase description.
# # #'
# # #' Longer proper discription of function...
# # #'
# # #' @param snps description.
# # #'
# # #' @author Caitlin Collins \email{caitiecollins@@gmail.com}
# # #' @export
# # #'
# # #' @examples
# # #' ## Example ##
# #
# # ########################################################################
#
# .get.locus01 <- function(subs.edges, root.nt, tree){
#
#
#   ##################
#   ## HANDLE TREE: ##
#   ##################
#   ## Always work with trees in "pruningwise" order:
#   tree <- reorder.phylo(tree, order="pruningwise")
#   ## Trees must be rooted:
#   if(!is.rooted(tree)) tree <- midpoint(tree)
#
#   ####################################################################
#   ############################
#   ## Get Anc-Des EDGE ORDER ##
#   ############################
#   ## Get sequence from lowest ("root", Nterm+1) to highest ancestral node:
#   ix <- c(min(tree$edge[,1]):max(tree$edge[,1]))
#   ## Get for loop index of rows in tree$edge[,1], in pairs, from lowest to highest:
#   x <- as.vector(unlist(sapply(c(1:length(ix)), function(e) which(tree$edge[,1] == ix[e]))))
#   ####################################################################
#
#
#   ## convert subs.edges into appropriate format:
#   snps.loci <- list()
#   snps.loci[[1]] <- subs.edges
#
#   ## rearrange snps.loci s.t it becomes a
#   ## list of length tree$edge.length,
#   ## each element of which contains the
#   ## locations of the mutations that will
#   ## occur on that branch
#   snps.loci <- sapply(c(1:length(tree$edge.length)),
#                       function(f)
#                         seq_along(snps.loci)[sapply(snps.loci,
#                                                     function(e) f %in% e)])
#
#
#   # we will store the output in a list called locus:
#   locus <- list()
#   ## get the node names for all individuals (terminal and internal)
#   all.inds <- sort(unique(as.vector(unlist(tree$edge))))
#   ## we start w all inds having same genotype as root:
#   for(j in all.inds){
#     locus[[j]] <- root.nt
#   }
#
#   ## store replacement nts in list new.nts:
#   new.nts <- list()
#   ## distinguish btw list of loci and unique list
#   snps.loci.ori <- snps.loci
#   ## will need to treat repeat loci differently...
#   snps.loci.unique <- lapply(snps.loci, unique)
#
#
#   #############################
#   ## For Loop to get new nts ##
#   #############################
#   for(i in x){
#     ## for all locus other than root, we mutate the
#     ## genome of the node preceding it, according to snps.loci.
#     ## Draw new nts for each locus selected for mutation:
#     if(!.is.integer0(snps.loci.unique[[i]])){
#       new.nts[[i]] <- sapply(c(1:length(snps.loci.unique[[i]])), function(e)
#         selectBiallelicSNP(c("0", "1")[which(c("0", "1")
#                                              %in% locus[[tree$edge[i,1]]]
#                                              [snps.loci.unique[[i]][e]])]))
#       ## if any loci are selected for multiple mutations
#       ## within their given branch length:
#       if(length(snps.loci.ori[[i]]) != length(snps.loci.unique[[i]])){
#         ## identify which loci are repeaters
#         repeats <-table(snps.loci.ori[[i]])[which(table(snps.loci.ori[[i]])!=1)]
#         ## how many times they repeat
#         n.reps <- repeats - 1
#         ## the positions of these loci in the vector of snps loci
#         toRepeat <- which(snps.loci.unique[[i]] %in% names(repeats))
#         ## run chain of re-sampling to end in our new nt for repeater loci:
#         foo <- list()
#         for(j in 1:length(toRepeat)){
#           foo[[j]] <- new.nts[[i]][toRepeat[j]]
#           for(k in 1:n.reps[j]){
#             if(k==1){
#               foo[[j]][k] <- selectBiallelicSNP(c("0", "1")[which(c("0", "1")
#                                                                   %in% foo[[j]][1])])
#
#             }else{
#               foo[[j]][k] <- selectBiallelicSNP(c("0", "1")[which(c("0", "1")
#                                                                   %in% foo[[j]][k-1])])
#             }
#           }
#           ## retain only the last nt selected
#           out <- sapply(c(1:length(foo)),
#                         function(e) foo[[e]][length(foo[[e]])])
#         }
#         ## for the loci with repeated mts, replace these positions
#         ## in new.nts with the corresponding elements of out, above.
#         new.nts[[i]][toRepeat] <- out
#       } # end of if statement for repeaters
#
#       ## update ancestral genotype with new.nts:
#       temp <- locus[[tree$edge[i,1]]]
#       temp[snps.loci.unique[[i]]] <- new.nts[[i]]
#       locus[[tree$edge[i,2]]] <- temp
#
#     }else{
#       ## if no mts occur on branch, set genotype of
#       ## downstream individual to be equal to ancestor's
#       locus[[tree$edge[i,2]]] <- locus[[tree$edge[i,1]]]
#     }
#   } # end of for loop selecting new nts at mutator loci
#
#   ## turn locus into a vector for easier post-handling
#   locus <- as.vector(unlist(locus))
#
#   return(locus)
#
# } # end .get.locus01
















# snps <- get(load("C:/Cait 2016/Work/Xavier/Sims/set3/set3_9_snps.Rdata"))
# phen <- get(load("C:/Cait 2016/Work/Xavier/Sims/set3/set3_9_phen.Rdata"))
# perf <- get(load("C:/Cait 2016/Work/Xavier/Sims/set3/set3_9_performance.Rdata"))
# snps.assoc <- perf$snps.assoc
#
# res <- get(load("C:/Cait 2016/Work/Xavier/Sims/set3/set3_9_res.Rdata"))
# tree <- get(load("C:/Cait 2016/Work/Xavier/Sims/set3/set3_9_tree.Rdata"))
# tree <- get(load("C:/Cait 2016/Work/Xavier/Sims/set3/set3_Q.corr_unweighted/set3_9_tree.Rdata"))
# tree <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set3/08.2016/set3_9_tree.Rdata"))

# phen.rec <- res$dat$phen.rec
# snps.rec <- res$dat$snps.rec
#
# plot.phen(tree, phen.nodes=phen.rec, snp.nodes=snps.rec[,i])
#
# ## NOTE-- snps.rec has inverted snps numbering (ie. 0,1 --> 2,1)
# ## AND kept plink version w added columns:
# snps <- snps[, seq.int(1, ncol(snps), 2)]
# rownames(snps) <- c(1:nrow(snps))
# colnames(snps) <- c(1:ncol(snps))
# str(snps)

###########################################
## trying to get RATES out of Q probs... ##
###########################################
# SCORE3 <- get(load("C:/Cait 2016/Work/Xavier/Sims/set3/set3_9_score3.Rdata"))
#
# Q.ew <- SCORE3$corr.dat[[13]]
# Q.uw <- SCORE3$corr.dat[[14]]
#
# score3 <- SCORE3$corr.dat[[6]]
#
# i <- snps.assoc[3]
#
# #######################
# ## UN-edge-weighted: ##
# #######################
#
# Qi <- Q.uw[[i]]
#
# ## get NUMBERS of subs:
# Qi.subs <- Qi*length(score3[[i]][!is.na(score3[[i]])])
# # NOTICE much heavier diagonal :)
#
# ## re-weight by ROW?:
# temp <- t(sapply(c(1:nrow(Qi.subs)),
#                function(e)
#                  Qi.subs[e,]/sum(Qi.subs[e,])))
# colnames(temp) <- rownames(temp) <- colnames(Qi.subs)
# temp.uw <- temp
# temp.uw
#
# ## all rows should sum to 1:
# sum(temp.uw[1,])
#
#
# ####################
# ## EDGE-weighted: ##
# ####################
# # i <- snps.assoc[1]
#
# Qi <- Q.ew[[i]]
#
# ## get sum of edge lengths with these subs:
# Qi.ew <- Qi*sum(tree$edge.length[!is.na(score3[[i]])])
# # NOTICE much heavier diagonal :)
#
# sum(tree$edge.length)
# sum(tree$edge.length[!is.na(score3[[i]])])
# sum(Qi.ew)
#
# ## re-weight by ROW?:
# temp <- t(sapply(c(1:nrow(Qi.ew)),
#                  function(e)
#                    Qi.ew[e,]/sum(Qi.ew[e,])))
# colnames(temp) <- rownames(temp) <- colnames(Qi.ew)
# temp.ew <- temp
# temp.ew
#
# ## all rows should sum to 1:
# sum(temp.ew[1,])
#
# #############
# ## Q.corr: ##
# #############
#
# ## Q.corr (original): ##
# Q.corr <- matrix(c(2, 0.75, 0.75, 1, 3, 0.5, 0.25,  3, 3, 0.25, 0.5, 3, 1, 0.75, 0.75, 2),
#                  nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
#
# ## re-weight by ROW?:
# temp <- t(sapply(c(1:nrow(Q.corr)),
#                  function(e)
#                    Q.corr[e,]/sum(Q.corr[e,])))
# colnames(temp) <- rownames(temp) <- colnames(Q.corr)
# temp.corr <- temp
# temp.corr
#
# ## NEW Q.corr ( --> reduced n.subs??)
# Q <- matrix(c(0.6, 0.1, 0.1, 0.2, 0.4, 0.15, 0.05,  0.4, 0.4, 0.05, 0.15, 0.4, 0.2, 0.1, 0.1, 0.6),
#             nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
# Q
# Q.new1 <- Q.new
#
# ## 2nd NEW Q.corr ( --> reduced n.subs??)
# Q <- matrix(c(0.7, 0.05, 0.05, 0.2, 0.45, 0.1, 0.0,  0.45, 0.45, 0.0, 0.1, 0.45, 0.2, 0.05, 0.05, 0.7),
#             nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
# Q
# Q.new2 <- Q.new <- Q
#
############################

############
## Q.inst ##
############
## IGNORE !!
# Q.inst <- matrix(c(0, 0.5, 0.5, 0, 2, 0, 0,  2, 2, 0, 0, 2, 0, 0.5, 0.5, 0),
#                  nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
# Q.inst
# Q <- Q.inst
#################################

##########################
## Q.new3 (w/ matexpo!) ##
##########################
## if Q contains RATES --> P contains probs
# Q.new3 <- matrix(c(0, 0.5, 0.5, 1,
#                    2, 0, 0.25,  2,
#                    2, 0.25, 0, 2,
#                    1, 0.5, 0.5, 0),
#                  nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
#
# Q.mat <- matrix(c(-1*s, 0.5*s, 0.5*s, 0,
#                   2*s, -4*s, 0, 2*s,
#                   2*s, 0, -4*s, 2*s,
#                   0, 0.5*s, 0.5*s, -1*s),
#                 nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
#
# af <- 4
# s <- 0.5
# Q.mat <- matrix(c(NA, 1*s, 1*s, 0,
#                   1*af*s, NA, 0, 1*af*s,
#                   1*af*s, 0, NA, 1*af*s,
#                   0, 1*s, 1*s, NA),
#                 nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
#
# ## GET PROBABILITY MATRIX:
# diag(Q.mat) <- sapply(c(1:nrow(Q.mat)), function(e) -sum(Q.mat[e, c(1:ncol(Q.mat))[-e]]))
# Q <- Q.mat
# Q.new3 <- Q.mat
# Q.mat
#
# ## get conditional probs for each edge w matexpo! ##
# ## (run within code...)
# x <- rev(1:nrow(tree$edge))
# Qt <- list()
# for(e in x){
#   Qt[[e]] <- matexpo(Q*tree$edge.length[e])
# }
# e <- x[1]
# P <- Qt[[e]]
# rownames(P) <- rownames(Q)
# colnames(P) <- colnames(Q)
# P

## CHECK -- Does expm (pkg "Matrix") give same results as matexpo?? ## (YES)

##########################
## Q.new4 (w/ matexpo!) ##
##########################
## if Q contains RATES --> P contains probs

## ?? -- (WHY) DOES THE BL:TR DIAGONAL NEED TO BE ALL ZEROS???? ################################    ####        ####        ####    ??????????       ####
# s <- 1
# Q.new4 <- matrix(c(-2*s, 0.5*s, 0.5*s, 1*s,
#                    2*s, -4.25*s, 0.25*s,  2*s,
#                    2*s, 0.25*s, -4.25*s, 2*s,
#                    1*s, 0.5*s, 0.5*s, -2*s),
#                  nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))

# s <- 1
# Q.new4 <- matrix(c(-1*s, 0.5*s, 0.5*s, 0,
#                    2*s, -4*s, 0, 2*s,
#                    2*s, 0, -4*s, 2*s,
#                    0, 0.5*s, 0.5*s, -1*s),
#                  nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
#
# Q.new4
# Q <- Q.new4
#
# ## GET PROBABILITY MATRIX:
# ## get conditional probs for each edge w matexpo! ##
# ## (run within code...)
# x <- rev(1:nrow(tree$edge))
# Qt <- list()
# for(e in x){
#   Qt[[e]] <- matexpo(Q*tree$edge.length[e])
# }
# e <- x[2]
# P <- Qt[[e]]
# rownames(P) <- rownames(Q)
# colnames(P) <- colnames(Q)
# P




####################

#########################
## Q.xav (w/ matexpo!) ##
#########################
## if Q contains RATES --> P contains probs

# ## s = substitution rate
# s <- 1
# ## NOT SURE IF ROW/COLNAMES ARE RIGHT FOR THIS ONE:
# Q.xav <- matrix(c(-2*s, 1*s, 1*s, 0,
#                    2*s, -4*s, 0,  2*s,
#                    2*s, 0, -4*s, 2*s,
#                    0, 1*s, 1*s, -2*s),
#                  nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
# Q <- Q.xav
# Q.xav
#
# #########
#
# ##  NOT SURE HOW TO INTERPRET/PREDICT THE RELATIVE EFFECTS OF ASSOC.FACTOR AND N.SUBS
# ## (ie AF, S) ON ASSOC STRENGTH AND N.SUBS PER TREE (AND HOW SUB PROBS VARY W BRANCH LENGTH)
# ## ALSO -- (WHY) DO THE DIAGONALS NEED TO BE NEG AND 0?
# ## (CONSIDERING WE WANT TO SIM SOME SIMULTANEOUS SUBS..)
# ## (SEEMS LIKE THE 0-DIAGONAL HAS NO IMPACT WHETHER IT'S 0S OR OTHER VALUES?)
#
# af <- 6 # association factor
# s <- 20 # n.subs
#
# # Q.mat <- matrix(c(NA, 1*s, 1*s, 0,
# #                   1*af*s, NA, 0, 1*af*s,
# #                   1*af*s, 0, NA, 1*af*s,
# #                   0, 1*s, 1*s, NA),
# #                 nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
#
# Q.mat <- matrix(c(NA, 1*s, 1*s, 0,
#                   1*af*s, NA, 0, 1*af*s,
#                   1*af*s, 0, NA, 1*af*s,
#                   0, 1*s, 1*s, NA),
#                 nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
#
# diag(Q.mat) <- sapply(c(1:nrow(Q.mat)), function(e) -sum(Q.mat[e, c(1:ncol(Q.mat))[-e]]))
#
# Q <- Q.mat
# Q
#
#
# ## GET PROBABILITY MATRIX:
#
# ## get conditional probs for each edge w matexpo! ##
# ## (run within code...)
# x <- rev(1:nrow(tree$edge))
# Qt <- list()
# for(e in x){
#   Qt[[e]] <- matexpo(Q*tree$edge.length[e])
# }
# e <- x[112]
# P <- Qt[[e]]
# rownames(P) <- rownames(Q)
# colnames(P) <- colnames(Q)
# P


#####################

# ######################################################
# ## WANT -- 15 subs/tree (thus 183 no-sub branches): ##
# ######################################################
# # probs.e <- round((tree$edge.length/sum(tree$edge.length)),3)
# probs.e <- (tree$edge.length/sum(tree$edge.length))
#
# ## IDEA -- multiply all 3 non-stay-the-same cells of Q by above (or make it s.t. these cells are x% of the row sum)
# ## EG -- if you would expect 1.4/15 subs to occur on edge 1
#
# ## --> Need to find Pr(at least one sub occurs on branch e)
# ## bc we are/were working with sampling WITHOUT replacement w probs proportional to edge.length
#
# ## Pr(no sub on branch 198) =
# (1-probs.e[198])^15 # 0.022
# ## Therefore: Pr(at least one sub on branch 198) =
# 1 - (1-probs.e[198])^15 # 0.98

## Then -- need to re-work probs from Q s.t Pr(stay the same) = 0.02 and Pr(sub) sums to 0.98
