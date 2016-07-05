

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

## ARGUMENTS: ##

## n.subs <- either an integer or a vector containing a distribution of n.subs-per-site
## phen.loci <- a vector containing the indices of the edges on which phen subs occurred


snp.sim <- function(n.snps = 10000,
                    n.subs = 1,
                    snp.root = NULL,
                    n.snps.assoc = 10,
                    assoc.prob = 100,
                    tree = coalescent.tree.sim(100),
                    phen.loci = NULL,
                    heatmap = FALSE,
                    reconstruct = FALSE,
                    dist.dna.model = "JC69",
                    row.names = NULL,
                    seed=1){

  require(adegenet)
  require(ape)
  require(phangorn)

  ##################################
  ## GET MUTATIONS' branch & loci ##
  ##################################
  n.ind <- tree$Nnode+1
  gen.size <- n.snps

  if(!is.null(seed)) set.seed(seed)

  ## Simulate genotype for root individual:
  if(class(n.subs) != "matrix"){

    ## For n.subs = n or = dist approaches:
    gen.root <- sample(c("a", "c", "g", "t"), gen.size, replace=TRUE)

  }else{
    ## For ACE/Pagel-test approach:
    ## Input = either:
    ## one state (chosen by directly selecting the more likely state), or
    ## two likelihoods (taken from fit.iQ$lik.anc[1,] from w/in fitPagel).

    ## One state:
    if(!is.null(snp.root)){
      if(length(snp.root) == 1){
        ## select only root state --> different SNP sim method (???)
        if(snp.root == 0) gen.root <- "a"
        if(snp.root == 1) gen.root <- "t"

        ## select root state & assign this state to all nodes,
        ## to be changed later by a modifiction of the existing SNP sim method...
        #if(snp.root == 0) gen.root <- rep("a", gen.size)
        #if(snp.root == 1) gen.root <- rep("t", gen.size)
      }
    }

  }
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

  if(length(n.subs)==1){

    #####################
    ## NO DISTRIBUTION ##
    #####################

    ## if no distribution is inputted,
    ## use normal simulation procedure
    ## (ie. Poisson parameter 1):

    ## draw the number of mutations to occur at each site:
    n.mts <- rpois(n=gen.size, lambda=(n.subs))
    ## for any n.mts==0, re-sample
    for(i in 1:length(n.mts)){
      while(n.mts[i]==0){
        n.mts[i] <- rpois(n=1, lambda=(n.subs))
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
      # } # end dist
  } # end fitPagel or dist

  ############################
  ## Assign mts to branches ##
  ############################

  ## whether n.mts is chosen by Poisson or according to a Distribution...

  if(n.snps.assoc != 0){
    ## for snps.assoc (the last n.snps.assoc snps, for now),
    ## add n.mts == n.phen.loci s.t these sites mutate at each
    ## and every phen.loci (for now, to be diluted later
    ## according to assoc.prob if !=100)
    n.mts <- c(n.mts, rep(length(phen.loci), n.snps.assoc))
  }

  ## for each site, draw the branches to which
  ## you will assign the mts for this site
  ## (~ branch length):
  snps.loci <- sapply(c(1:length(n.mts)),
                      function(e)
                        sample(c(1:length(tree$edge.length)),
                               n.mts[e],
                               replace=FALSE,
                               prob=tree$edge.length))
  if(n.snps.assoc != 0){
    ## get snps.loci for the ASSOCIATED snps (ie. set to phen.loci) ##
    for(i in 1:n.snps.assoc){
      snps.loci[[snps.assoc[i]]] <- phen.loci
      #print(snps.loci[[snps.assoc[i]]])
    }
  }
  ## rearrange snps.loci s.t it becomes a
  ## list of length tree$edge.length,
  ## each element of which contains the
  ## locations of the mutations that will
  ## occur on that branch
  snps.loci <- sapply(c(1:length(tree$edge.length)),
                      function(f)
                        seq_along(snps.loci)[sapply(snps.loci,
                                                    function(e) f %in% e)])


  # we will store the output in a list called genomes:
  genomes <- list()
  ## get the node names for all individuals (terminal and internal)
  all.inds <- sort(unique(as.vector(unlist(tree$edge))))
  ## we start w all inds having same genotype as root:
  for(i in all.inds){
    genomes[[i]] <- gen.root
  }

  ## store replacement nts in list new.nts:
  new.nts <- list()
  ## distinguish btw list of loci and unique list
  snps.loci.ori <- snps.loci
  ## will need to treat repeat loci differently...
  snps.loci.unique <- lapply(snps.loci, unique)
  ## the last individual in the first column of tree$edge
  ## (ie. ind.length(tree$tip.label)+1 ) is our root individual:
  x <- rev(c(1:nrow(tree$edge)))


  #############################
  ## For Loop to get new nts ##
  #############################
  for(i in x){
    ## for all genomes other than root, we mutate the
    ## genome of the node preceding it, according to snps.loci.
    ## Draw new nts for each locus selected for mutation:
    if(!.is.integer0(snps.loci.unique[[i]])){
      new.nts[[i]] <- sapply(c(1:length(snps.loci.unique[[i]])), function(e)
        selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
                                                       %in% genomes[[tree$edge[i,1]]]
                                                       [snps.loci.unique[[i]][e]])]))
      ## if any loci are selected for multiple mutations
      ## within their given branch length:
      if(length(snps.loci.ori[[i]]) != length(snps.loci.unique[[i]])){
        ## identify which loci are repeaters
        repeats <-table(snps.loci.ori[[i]])[which(table(snps.loci.ori[[i]])!=1)]
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
      temp <- genomes[[tree$edge[i,1]]]
      temp[snps.loci.unique[[i]]] <- new.nts[[i]]
      genomes[[tree$edge[i,2]]] <- temp

    }else{
      ## if no mts occur on branch, set genotype of
      ## downstream individual to be equal to ancestor's
      genomes[[tree$edge[i,2]]] <- genomes[[tree$edge[i,1]]]
    }
  } # end of for loop selecting new nts at mutator loci

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
        prob <- assoc.prob[i]
        ## only if the association is imperfect
        if(prob != 100){
          ## draw genomes to change at snps.assoc[i]
          n.toChange <- round(length(genomes)*(1 - (prob/100)))
          toChange <- sample(c(1:length(genomes)), n.toChange)

          ## change those genomes at snps.assoc[i]
          for(j in 1:length(toChange)){
            genomes[[toChange[j]]][snps.assoc[i]] <-
              selectBiallelicSNP(genomes[[toChange[j]]][snps.assoc[i]])
          } # end for loop
        }
      } # end for loop
    } # end any assoc.prob != 100
  } # end modification | assoc.prob

  ##############################
  ## PLOTS & TREECONSTRUCTION ##
  ##############################

  if(heatmap == TRUE || reconstruct!=FALSE){
    dna <- as.DNAbin(genomes)
    names(dna) <- c(1:length(genomes))
  }

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
  tree.reconstructed <- NULL
  if(reconstruct!=FALSE){
    if(reconstruct==TRUE){
      warning("reconstruct should be one of 'UPGMA', 'nj', 'ml'. Choosing 'UPGMA'.")
    }

    tree.reconstructed <- tree.reconstruct(dna[1:n.ind],
                                         method=reconstruct,
                                         dist.dna.model=dist.dna.model,
                                         plot=TRUE)
  }

  ##################
  ## CONVERT SNPS ##
  ##################

  ## keep and return ONLY genomes for TERMINAL individuals
  genomes <- genomes[1:n.ind]

  x <- genomes
  gen.size <- length(x[[1]])
  ## working with snps in matrix form
  snps <- do.call("rbind", x)
  ## get snps as DNAbin
  ploidy <- 1
  snps.bin <- as.DNAbin(snps, ploidy=ploidy)
  ## get snps as genind
  #source("C:/Users/caitiecollins/adegenet/R/sequences.R")
  snps.gen <- DNAbin2genind(snps.bin, polyThres=0.01)
  ## get snps as binary matrix
  snps <- snps.gen@tab

  ## correct genind for ploidy:
  snps <- snps[,seq(1, ncol(snps), 2)]
  colnames(snps) <- 1:ncol(snps)
  # colnames(snps) <- gsub("[:.:]|.$*", "", colnames(snps))

  if(!is.null(snps.assoc)){
    #############################################
    ## Ensure snps.assoc loci are ###############
    ## correctly labelled | n.columns retained ##
    ## in snps.gen object #######################
    #############################################
    ## identify any columns of snps.bin that
    ## do NOT meet DNAbin2genind polyThres
    ## or are NOT SNPs
    x <- snps.bin
    if(is.list(x)) x <- as.matrix(x)
    if(is.null(colnames(x))) colnames(x) <- 1:ncol(x)

    temp <- lapply(1:ncol(x), function(i)
      .getFixed(x[,i], i)) # process all loci, return a list
    fixed.loci <- which(temp==TRUE) ## identify loci that are NOT SNPs

    ## update snps.assoc to reflect true loci
    gen.size.final <- ncol(snps)
    snps.assoc.loci.ori <- c((gen.size.final-(n.snps.assoc-1)):gen.size.final)

    ## Re-enabled snps.assoc loci "randomization" by
    ## just drawing indices and shuffling the columns accordingly...
    ## draw which SNPs will be associated to the phenotype
    snps.assoc.loci <- sort(sample(c(1:gen.size.final),
                                   n.snps.assoc,
                                   replace=FALSE))

    snps.indices <- c(1:gen.size.final)
    snps.ori <- snps
    #     snps.names.ori <- ind.names.ori <- NULL
    #     if(!is.null(dimnames(snps))){
    #       snps.names.ori <- dimnames(snps)[[2]]
    #       ind.names.ori <- dimnames(snps)[[1]]
    #     }
    snps.non.assoc <- snps[,c(1:(gen.size.final-n.snps.assoc))]
    snps.assoc <- snps[,snps.assoc.loci.ori]
    snps.new <- matrix(99, nrow=nrow(snps), ncol=gen.size.final)
    # snps.new <- matrix(99, nrow=100, ncol=gen.size.final)
    snps.new[,snps.indices[-snps.assoc.loci]] <- snps.non.assoc
    snps.new[,snps.assoc.loci] <- snps.assoc
    snps <- snps.new
    if(!is.null(snps.names.ori)) colnames(snps) <- snps.names.ori
    if(!is.null(ind.names.ori)) rownames(snps) <- ind.names.ori
    snps.assoc <- snps.assoc.loci

    ## update names of snps.assoc (??? REMOVED; to change L00NN --> NN 05/07/2016)
    #     gen.size <- ncol(snps)
    #     names(snps.assoc) <- as.vector(unlist(sapply(c(1:length(snps.assoc)),
    #                                                  function(e)
    #                                                    paste("L",
    #                                                          paste(rep(0, (nchar(gen.size)-
    #                                                                          nchar(snps.assoc[e]))),
    #                                                                sep="",
    #                                                                collapse=""),
    #                                                          snps.assoc[e], sep=""))))
  }

  ## Assign row & col names...

  ## assign/generate row.names
  if(!is.null(row.names)){
    if(length(row.names) == nrow(snps)){
      rownames(snps) <- row.names
    }else{
      if(is.null(rownames(snps))) rownames(snps) <- c(1:nrow(snps))
    }
  }else{
    if(is.null(rownames(snps))) rownames(snps) <- c(1:nrow(snps))
  }

  ## ensure colnames are complete & sequential
  ## CAREFUL----------STILL NEED TO CHECK THAT DROPPING SNPS IS NOT AFFECTING THE SNPS.ASSOC LOCI NAMES BEING REPORTED WITH OUTPUT!!!!!!
  #################################################################################################################################################################
  colnames(snps) <- 1:ncol(snps)

  ##################
  ## get RESULTS: ##
  ##################
  out <- list(snps, snps.assoc, tree.reconstructed)
  names(out) <- c("snps", "snps.assoc", "tree.reconstructed")

  return(out)

} # end snp.sim
