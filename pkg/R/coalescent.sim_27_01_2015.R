
####################
## coalescent.sim ##
####################

## a function for simulating trees under a fully-linked coalescent model.
## optional simulation of a phenotype and phenotypically-associated SNPs is implemented.
## optional use of a distribution to guide the substitution rate of the non-associated SNPs is implemented.

## TO DO:
## 1) (Re-)implement associated SNP randomization procedure...
## want to implement procedures that combine the above options...
## 2) Allow phenotypically-associated SNPs simulation to be optionally guided
## by a user-inputted phenotype for the terminal nodes (--> would need to simulate
## phenotypic substitutions from the terminal nodes UP to the root, the reverse
## of the current procedure...)
## 3) Implement assoc.options (currently using deprecated "all" option without requiring argument,
## but would like to consider implementing alternative "model" option(s))


## ARGUMENTS ##
# n.ind <- 10 # n.genomes you want to end up with
# gen.size <- 1000000 # bases
# theta <- gen.size*2 # (if sim.by=="branch")# OR # 1*2 # (if sim.by=="locus")
# biallelic <- TRUE # if TRUE, select ONLY complementary nt; if FALSE, select from 3 alternatives (ie. A/C/G/T-current nt)
# seed <- 1 # allow user to control randomization to get reproducible results.
# n.snp.assoc <- 5
# assoc.option <- c("all", "model") # deprecated (only "all" available)
# sim.by <- c("locus", "branch") # deprecated (only "locus" has all current protocols implemented)


## EXAMPLE ##
# out <- coalescent.sim(n.ind=100, gen.size=10000, sim.by="locus",
#                       theta=1*2, dist=NULL,
#                       theta_p=15, phen=NULL,
#                       n.snp.assoc=20, assoc.option="all", assoc.prob=90,
#                       haploid=TRUE, biallelic=TRUE, seed=1,
#                       plot=TRUE, heatmap=FALSE, plot2="UPGMA")

########################################################################

###################
## DOCUMENTATION ##
###################

#' Short one-phrase description.
#'
#' Longer proper discription of function...
#'
#' @param n.ind description.
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




coalescent.sim <- function(n.ind=100, gen.size=10000, sim.by="locus",
                           theta=1*2, dist=NULL,
                           theta_p=15, phen=NULL,
                           n.snp.assoc=5, assoc.option="all", assoc.prob=90,
                           haploid=TRUE, biallelic=TRUE, seed=NULL,
                           plot=TRUE, heatmap=FALSE, plot2="UPGMA"){

  require(adegenet)
  require(ape)
  require(phangorn)

  #   ## load utils.R for selectBiallelicSNP switch.phen .is.integer0 .is.even .is.odd
  #   source("C:/Users/Caitlin/treeWAS/pkg/R/utils.R")

  if(plot==TRUE && heatmap==FALSE && plot2==FALSE){
    par(ask=FALSE)
  }else{
    par(ask=TRUE)
  }

  if(missing(theta)) theta <- NULL
  if(is.null(theta)){
    if(sim.by=="branch"){
      theta <- gen.size*2
    }
    if(sim.by=="locus"){
      theta <- 1*2
    }
  }


  if(!is.null(seed)) set.seed(seed)

  n.nodes <- n.ind + (n.ind-1) # total n.nodes (internal, external)
  inds <- c(1:n.ind) # terminal nodes
  nodes <- rev(c((n.ind+1):n.nodes)) # internal nodes
  tree.params <- list() # to store and update output


  for(i in 1:(length(inds)-1)){
    ## get inds.ori from last generation:
    if(i==1){
      inds.ori <- inds
    }else{
      inds.ori <- tree.params[[(i-1)]][["inds.remaining"]]
    }
    ## get N, the number of individuals (remaining at this generation)
    ## from which a random 2 are to be selected for coalescence
    N <- length(inds.ori)

    ###################
    ## BRANCH LENGTH ##
    ###################
    ## get lamda, the parameter of the exponential distribution,
    ## given the number of individuals at this generation
    lambda <- (N*(N-1)) / 2
    ## draw x, the length of time to coalescence at this generation
    x <- rexp(n=1, rate=lambda)

    #####################
    ## COALESCENT PAIR ##
    #####################
    ## get co.pair, the 2 inds to coalesce at this generation
    co.pair <- sample(inds.ori, 2)
    ## merge these 2 inds, replace with new internal node,
    ## update the list of inds to sample at the next generation
    inds.remaining <- c(inds.ori[-which(inds.ori %in% co.pair)], nodes[i])

    ############
    ## OUTPUT ##
    ############
    ## store the output in the ith element of our list tree.params:
    tree.params[[i]] <- list()
    tree.params[[i]][[1]] <- x
    tree.params[[i]][[2]] <- co.pair
    tree.params[[i]][[3]] <- inds.ori
    tree.params[[i]][[4]] <- inds.remaining
    names(tree.params[[i]]) <- c("Time", "co.pair", "inds.ori", "inds.remaining")
  } # end for loop




  ## get edge.list
  to <- as.vector(unlist(sapply(c(1:length(tree.params)),
                                function(e) tree.params[[e]][["co.pair"]])))
  from <- as.vector(unlist(sapply(c(1:length(nodes)),
                                  function(e) rep(nodes[e], 2))))
  edge.list <- data.frame(from,to)

  ## get edge lengths
  times <- as.vector(unlist(sapply(c(1:length(tree.params)),
                                   function(e) tree.params[[e]][["Time"]])))

  ## make empty edge.lengths vector to store output below:
  edge.lengths <- NA

  ## for all the edges in our edge.list data.frame:
  for(i in 1:nrow(edge.list)){
    if(edge.list$to[i] %in% inds){
      ## if downstream node = terminal, sum all time intervals til ancestor.
      edge.lengths[i] <- sum(times[1:which(nodes==edge.list$from[i])])
    }else{
      ## BUT, if the downstream node = internal, must subtract time btw.
      ## downstream node and final generation.
      length.total <- sum(times[1:which(nodes==edge.list$from[i])])
      length.toRemove <- sum(times[1:which(nodes==edge.list$to[i])])
      edge.lengths[i] <- length.total - length.toRemove
    }
  } # end for loop

  ## convert edge.list to matrix
  edge.list <- as.matrix(edge.list)
  colnames(edge.list) <- NULL
  dimnames(edge.list) <- NULL

  ## put output into tree list (phylo format):
  tree <- list()
  tree$edge <- edge.list
  tree$tip.label <- c(1:n.ind)
  tree$edge.length <- edge.lengths
  tree$Nnode <- as.integer(n.ind - 1)

  ## change class by force
  class(tree) <- "phylo"

  #######################################################################
  #################### ^ PHYLOGENETIC TREE CREATED ^ #################### # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #######################################################################


  ####################################
  ## PHENOTYPE simulation procedure ## ~ sim.by.locus...
  ####################################

  ## simulate phenotype for root individual:
  if(!is.null(theta_p)){
    phen.root <- "A"
  }else{
    phen.root <- NULL
  }

  ## make dummy variables in which to store the resulting n.mts variables:
  lambda_p <- n.subs <- NA

  ## ensure phen variables start as NULL
  phen.branch <- phen.nodes <- phen.leaves <- NULL



  #############################################################
  ## If the user has specified a "mt" rate for the phenotype ##
  #############################################################

  ## (indicating that they want to generate a NEW phenotype for the tree provided)
  if(!is.null(theta_p)){

    ## draw the number of mutations to occur at each site:
    n.subs <- rpois(n=1, lambda=theta_p)
    ## if n.subs==0 or ==1, re-sample
    while(n.subs <= 1){
      n.subs <- rpois(n=1, lambda=theta_p)
    }

    ## draw the branches to which you will assign the n.subs to occur for the phenotype (~ branch length):
    phen.loci <- sample(c(1:length(tree$edge.length)), n.subs, replace=FALSE, prob=tree$edge.length)
    ## rearrange phen.loci
    phen.loci <- sort(phen.loci, decreasing=TRUE)


    ###############################
    ## For Loop to get PHENOTYPE ##
    ###############################
    ## get phenotype for all branches/ nodes in tree (from root node (ie. tree$edge[nrow(tree$edge), 1]) down):
    phen.nodes <- phen.branch <- list()

    ## set phenotype for all branches and nodes to be phen.root:
    for(i in 1:length(tree$edge.length)) phen.branch[[i]] <- phen.root
    names(phen.branch) <- paste("e", c(1:length(phen.branch)), sep=".")
    for(i in 1:length(unique(as.vector(unlist(tree$edge))))) phen.nodes[[i]] <- phen.root
    names(phen.nodes) <- paste("n", c(1:length(phen.nodes)), sep=".")

    #############################################################################

    #############################################################################

    ## get the node names for all individuals (terminal and internal)
    all.inds <- sort(unique(as.vector(unlist(tree$edge))))

    x <- rev(c(1:nrow(tree$edge)))

    ## get phen of nodes
    for(i in 1:length(x)){
      if(x[i] %in% phen.loci){
        phen.nodes[[tree$edge[x[i],2]]] <- .switch.phen(phen.nodes[[tree$edge[x[i],1]]])
      }else{
        ## if no phen subs occur on branch i, set phen of
        ## downstram individual to be equal to ancestor's
        phen.nodes[[tree$edge[x[i],2]]] <- phen.nodes[[tree$edge[x[i], 1]]]
      }
    } # end for loop

    ## get phen of TERMINAL nodes (leaves)
    phen.leaves <- as.factor(as.vector(unlist(phen.nodes[c(1:n.ind)])))
    names(phen.leaves) <- paste("ind", c(1:length(phen.leaves)), sep=".")

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

    ## end PHEN sim procedure...

    #############################################################################
    ######################## PLOT phylo with PHEN ###############################
    #############################################################################

    ## get COLOR for NODES
    nodeCol <- phen.nodes
    nodeCol <- replace(nodeCol, which(nodeCol == "B"), "red")
    nodeCol <- replace(nodeCol, which(nodeCol == "A"), "blue")
    nodeCol <- as.vector(unlist(nodeCol))
    ## get COLOR for LEAVES ONLY
    leafCol <- nodeCol[1:n.ind]
    ## get COLOR of INTERNAL nodes ONLY
    internalNodeCol <- nodeCol[(n.ind+1):length(nodeCol)]

    ## get COLOR for EDGES
    edgeCol <- phen.branch
    ## FOR NOW--all edges w > 1 phen (either c("A", "B") OR c("B", "A")) --> purple...
    l_edgeCol <- which(sapply(c(1:length(edgeCol)), function(e) length(edgeCol[[e]])) == 2)
    edgeCol <- replace(edgeCol, l_edgeCol, "purple")
    edgeCol <- replace(edgeCol, which(edgeCol == "A"), "blue")
    edgeCol <- replace(edgeCol, which(edgeCol == "B"), "red")
    edgeCol <- as.vector(unlist(edgeCol))
    ## get COLOR for EDGE LABELS
    edgeLabCol <- rep("green", nrow(tree$edge))
    edgeLabCol <- replace(edgeLabCol, phen.loci, "mediumorchid1")

    ###############
    ## plot TREE ##
    ###############
    if(plot==TRUE){
      if(n.ind <= 20){
        plot(tree, show.tip=FALSE, edge.width=2, edge.color=edgeCol) # edgeCol
        title("Coalescent tree w/ phenotypic changes")
        axisPhylo()
        edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."), cex=0.5, font=2, bg=transp(edgeLabCol, 0.3), adj=c(1,1))
        tiplabels(text=tree$tip.label, cex=0.6, adj=c(-0.5, 0), bg=transp(leafCol, 0.3))
        nodelabels(text=rev(unique(tree$edge[,1])), cex=0.5, bg=transp(internalNodeCol, 0.3))
      }else{
        plot(tree, show.tip=FALSE, edge.width=2, edge.color=edgeCol) # edgeCol
        title("Coalescent tree w/ phenotypic changes")
        axisPhylo()
        #edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."), cex=0.5, font=2, bg=transp(edgeLabCol, 0.3), adj=c(1,1))
        tiplabels(text=tree$tip.label, cex=0.6, adj=c(-0.5, 0), col=leafCol, frame="none")
        #nodelabels(text=rev(unique(tree$edge[,1])), cex=0.5, bg=transp(internalNodeCol, 0.3)) ## make sure this isn't backward...
        ## should be numbered s.t. the root node is n.term+1
        ## RECALL: terminal nodes are numbered 1:n.ind from bottom to top of plot of tree;
        ## edges are numbered 1:nrow(edges) by following the lowest trace on the plot
        ## (starting from the root down to the lowermost tips);
        ## thus, internal nodes are numbered (n.ind+1):(n.ind+(n.ind-1)),
        ## from root to the top-most internal node to be connected (ie. the highest in the plot)
      }
    } # end plot = TRUE
  } # end if(!is.null(theta_p)) ## ie. SIMULATED phenotype & plotting


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


  #######################################################################

  #######################################################################


  if(!is.null(seed)) set.seed(seed)
  ## simulate genotype for root individual:
  gen.root <- sample(c("a", "c", "g", "t"), gen.size, replace=TRUE)
  ## get the sum of all branch lengths in the tree:
  time.total <- sum(tree$edge.length)

  ## make dummy variables in which to store the resulting n.mts variables:
  L <- lambda <- n.mts <- new.nt <- NA

  ##################################
  ## GET MUTATIONS' branch & loci ##
  ##################################

  ###############
  ## BY BRANCH ## ## NOTE: ## for sim.by BRANCH, MT.RATE ~= the EXPECTED n.mts TOTAL (in the WHOLE GENOME) (ie. gen.size(*2), eg. 1000(*2))!
  ###############             #######################################################################################################
  if(sim.by=="branch"){

    ## generate mts for each of n.ind genomes:
    ## NOTE: We're working from the root down,  so from the bottom row of
    ## tree$edge up to the top row --> Mts accumulate, w/ more Mts on longer branches...
    ## NOTE2: The placement of the root is arbitrary/ irrelevant to the outcome w.r.t the genetic
    ## distance/ relationship between the nodes...

    ##############################
    ## OLD simulation procedure ##
    ##############################

    snps.loci <- list()

    ## draw n.mts per branch and their loci
    for(i in rev(1:length(tree$edge.length))){
      ## get branch length i
      L[i] <- rev(tree$edge.length)[i]

      ## get lambda, the mean of the poisson dist
      ## from which we draw the n.mts to occur on branch i
      ## (proportional to branch length i out of the total sum of all branch lengths in the tree)
      lambda[i] <- ((theta/2) * (L[i] / time.total))

      ## get n.mts from poisson dist:
      n.mts[i] <- rpois(n=1, lambda[i])

      ## draw n.mts locations in the genome
      ## (NO LONGER allowing for replacement)
      snps.loci[[i]] <- sample(c(1:gen.size), n.mts[i], replace=FALSE)
    } # end for loop

  } # end sim.by branch


  ##############
  ## BY LOCUS ## ## NOTE: ## for sim.by LOCUS, MT.RATE ~= the EXPECTED n.mts AT EACH SITE (eg. 1(*2)) !
  ##############             #######################################################################
  if(sim.by=="locus"){

    #########################################
    ## NEW SNP-BY-SNP simulation procedure ##
    #########################################


    #############################
    ##  ## ASSOCIATED SNPS ##  ##
    #############################

    ##########################################
    ## OPTION 1: ALL ASSOCIATED SNPS CHANGE WITH PHENOTYPE*
    ## * (optionally modified later in code with assoc.prob to get probabilistic relationship rather than absolute)
    #### REPLACE n.mts[snps.assoc] w length(phen.loci)
    #### THEN, INSTEAD of drawing branches for those phen-associated SNPs, deliberately choose the phen.loci branches

    ## OPTION 2: SOME ASSOCIATED SNPS CHANGE WITH PHENOTYPE (ALL CHANGE, BUT NOT ALL TOGETHER)
    #### --ie. some phen caused by snps.assoc[1:3] and remaining phen branches caused by snps.assoc[4:n.snp.assoc], eg.

    ## OPTION 3: ASSOCIATED SNPS CHANGE ACCORDING TO A MODEL
    #### --ie. INDIVIDUAL snps.assoc have some LOW probability of existing on NON-phen branches,
    ####  BUT on PHEN branches, ALL (or almost all, with some HIGH prob) snps.assoc change to (or remain) set as 1s
    ##########################################

    snps.assoc <- NULL

    ## if n.snp.assoc is neither NULL nor 0:
    if(is.null(n.snp.assoc)) n.snp.assoc <- 0
    if(n.snp.assoc != 0){

      ## get non.assoc gen.size
      gen.size.ori <- gen.size
      gen.size <- gen.size-n.snp.assoc

      ## assign snps.assoc to be the last n.snp.assoc snps columns
      snps.assoc <- c((gen.size+1):(gen.size+n.snp.assoc))

    }


    if(is.null(dist)){

      #####################
      ## NO DISTRIBUTION ##
      #####################

      ## if no distribution is inputted,
      ## use normal simulation procedure
      ## (ie. Poisson parameter 1):

      ## draw the number of mutations to occur at each site:
      n.mts <- rpois(n=gen.size, lambda=(theta/2))
      ## for any n.mts==0, re-sample
      for(i in 1:length(n.mts)){
        while(n.mts[i]==0){
          n.mts[i] <- rpois(n=1, lambda=(theta/2))
        }
      }

    }else{
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

      #set.seed(1)
      ## make n.mts, a vector of length ncol(snps)
      n.mts <- rep(1111, gen.size)
      loci.accounted <- c(1:gen.size)
      ## assign dist.prop[i] elements of n.mts
      ## to be the same as the n.subs
      ## indicated by i, the given element of dist.prop
      for(j in 1:length(dist.prop)){
        #print("J"); print(j)
        ## provided there are
        if(dist.prop[j] > 0){

          #################################################################################################################
          #################################################################################################################
          ## For the love of GOD I have NO idea why, but if I do not set the seed here
          ## (which results in the same sampling every time thus is problematic, obviously)
          ## then it chooses ONE element of n.mts AT RANDOM and appears to switch it BACK to its original state.
          ## What the heck is going on????
          #################################################################################################################
          #################################################################################################################
          set.seed(1)
          ## assign dist.prop[i] elements of n.mts to be i
          loci.selected <- sample(loci.accounted, dist.prop[j], replace = FALSE)
          #print("LOCI SELECTED, dist.prop[j]"); print(length(loci.selected)); print(dist.prop[j])
          loci.accounted <- loci.accounted[-which(loci.accounted %in% loci.selected)]
          #print("LOCI ACCOUNTED"); print(length(loci.accounted))
          n.mts[loci.selected] <- j
          #print("NMTS"); print(length(n.mts[loci.selected]))
        }
      }
      #tab <- table(n.mts)
      #tab
      #which(tab[1:length(tab)-1] != dist.prop)

    } # end dist

    ############################
    ## Assign mts to branches ##
    ############################
    ## whether n.mts is chosen by Poisson or according to a Distribution...

    if(n.snp.assoc != 0){
      ## for snps.assoc (the last n.snp.assoc snps, for now),
      ## add n.mts == n.phen.loci s.t these sites mutate at each
      ## and every phen.loci (for now, to be diluted later
      ## according to assoc.prob if !=100)
      n.mts <- c(n.mts, rep(length(phen.loci), n.snp.assoc))
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
    if(n.snp.assoc != 0){
      ## get snps.loci for the ASSOCIATED snps (ie. set to phen.loci) ##
      for(i in 1:n.snp.assoc){
        snps.loci[[snps.assoc[i]]] <- phen.loci
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


  } # end sim.by locus


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
      if(biallelic==FALSE){
        new.nts[[i]] <- sapply(c(1:length(snps.loci.unique[[i]])), function(e)
          sample(c("a", "c", "g", "t")[-which(c("a", "c", "g", "t")
                                              %in% genomes[[tree$edge[i,1]]]
                                              [snps.loci.unique[[i]][e]])], 1))
      }else{
        new.nts[[i]] <- sapply(c(1:length(snps.loci.unique[[i]])), function(e)
          selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
                                                         %in% genomes[[tree$edge[i,1]]]
                                                         [snps.loci.unique[[i]][e]])]))
      }
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
              if(biallelic==FALSE){
                foo[[j]][k] <- sample(c("a", "c", "g", "t")
                                      [-which(c("a", "c", "g", "t")
                                              %in% foo[[j]][1])], 1)
              }else{
                foo[[j]][k] <- selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
                                                                              %in% foo[[j]][1])])
              }
            }else{
              if(biallelic==FALSE){
                foo[[j]][k] <- sample(c("a", "c", "g", "t")
                                      [-which(c("a", "c", "g", "t")
                                              %in% foo[[j]][k-1])], 1)
              }else{
                foo[[j]][k] <- selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t")
                                                                              %in% foo[[j]][k-1])])
              }
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
  if(n.snp.assoc != 0){
    if(assoc.option == "all"){

      #       ## to compare to at end
      #       temp <- t(sapply(c(1:length(genomes)), function(i) genomes[[i]][snps.assoc]))
      #       if(nrow(temp) > 1){
      #       rownames(temp) <- paste("edge", c(1:length(genomes)), sep=".")
      #       temp.ori <- temp[c((n.ind+1):length(genomes), 1:n.ind),]
      #       }else{
      #         names(temp) <- paste("edge", c(1:length(genomes)), sep=".")
      #         temp.ori <- temp[c((n.ind+1):length(genomes), 1:n.ind)]
      #       }

      ## CHECKS ##
      ## if assoc.prob not set,
      ## set it to be 100% associated for all snps.assoc
      if(missing(assoc.prob)) assoc.prob <- NULL
      if(is.null(assoc.prob)) assoc.prob <- rep(100, n.snp.assoc)

      ## if we have any imperfect associations... ##
      if(any(assoc.prob != 100)){
        ## check length
        if(length(assoc.prob) != n.snp.assoc){
          ## if only 1 prob value given...
          if(length(assoc.prob) == 1){
            ## ... assume uniform assoc.prob;
            assoc.prob <- rep(assoc.prob, n.snp.assoc)
            ## no warning needed
          }else{
            ## BUT if assoc.prob of random length:
            ## repeat until of length n.snp.assoc
            assoc.prob <- rep(assoc.prob, length.out=n.snp.assoc)
            ## and print warning (only if not of length n.snp.assoc OR 1)
            warning("assoc.prob not of length n.snp.assoc;
                    sequence will be repeated until correct length is reached.")
          }
          } # end checks


        ## for each associated SNP,
        ## we undo some associations | assoc.prob for that snp.assoc
        for(i in 1:n.snp.assoc){
          prob <- assoc.prob[i]
          ## only if the association is imperfect
          if(prob != 100){
            ## draw genomes to change at snps.assoc[i]
            n.toChange <- round(length(genomes)*(1 - (prob/100)))
            toChange <- sample(c(1:length(genomes)), n.toChange)

            ## change those genomes at snps.assoc[i]
            for(j in 1:length(toChange)){
              genomes[[toChange[j]]][snps.assoc[i]] <- selectBiallelicSNP(genomes[[toChange[j]]][snps.assoc[i]])
            } # end for loop
          }
        } # end for loop
        } # end any assoc.prob != 100

      ## check (compare to temp from before for loop)
      ## to compare to at end
      #       temp <- t(sapply(c(1:length(genomes)), function(i) genomes[[i]][snps.assoc]))
      #       if(nrow(temp) > 1){
      #         rownames(temp) <- paste("edge", c(1:length(genomes)), sep=".")
      #         temp <- temp[c((n.ind+1):length(genomes), 1:n.ind),]
      #       }else{
      #         names(temp) <- paste("edge", c(1:length(genomes)), sep=".")
      #         temp <- temp[c((n.ind+1):length(genomes), 1:n.ind)]
      #       }
      #       length(which(temp == temp.ori))

    }
  } # end modification | assoc.prob



  if(heatmap == TRUE || plot2!=FALSE){
    dna <- as.DNAbin(genomes)
    names(dna) <- c(1:length(genomes))
  }


  #############
  ## HEATMAP ##
  #############
  if(heatmap==TRUE){
    ## get a distance matrix between the genomes
    D <- dist.dna(dna, model="JC69")

    mat <- t(as.matrix(D))
    mat <- mat[,ncol(mat):1]
    par(mar=c(1,5,5,1))
    image(x=1:ncol(mat), y=1:ncol(mat), mat,
          col=rev(heat.colors(100)),
          xaxt="n", yaxt="n", xlab="", ylab="")
    axis(side=2, at=c(1:ncol(mat)),
         lab=rev(names(dna)), las=2, cex.axis=1)
    axis(side=3, at=c(1:ncol(mat)),
         lab=names(dna), las=1, cex.axis=1)
    ## return margin parameter to default:
    par(mar=c(5,4,4,2)+0.1)
  }


  ####################
  ## PLOT (regular) ##
  ####################
  ## NOTE-- we only print the plot here IF we have NOT generated a phenotype OR read one in
  ########## (thereby having printed it earlier in the phen section)
  if(plot == TRUE && is.null(theta_p) && is.null(phen)){
    if(sim.by=="branch"){
      plot(tree, show.tip=FALSE, edge.width=2)
      title("Coalescent tree")
      axisPhylo()
      tiplabels(text=tree$tip.label, cex=1, adj=-.5)
      nodelabels(text=rev(unique(tree$edge[,1])), cex=0.75)
      edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."), cex=0.66)
      edgelabels(text=rev(n.mts), col="red", frame="none", cex=1.1, adj=c(1,-0.5))
    } # end sim.by branch

    if(sim.by=="locus"){
      plot(tree, show.tip=FALSE, edge.width=2)
      title("Coalescent tree")
      axisPhylo()
      tiplabels(text=tree$tip.label, cex=1, adj=-.5)
      nodelabels(text=rev(unique(tree$edge[,1])), cex=0.75)
      edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."), cex=0.66)
      edgelabels(text=sapply(c(1:length(snps.loci)), function(e)
        length(snps.loci[[e]])), col="red", frame="none", cex=1.1, adj=c(1,-0.5))
    } # end sim.by locus

  }

  ###################
  ## PLOT (simple) ##
  ###################
  if(plot=="simple"){
    plot(tree, show.tip=TRUE, edge.width=2, cex=0.6, adj=.25)
    title("Coalescent tree")
    axisPhylo()
  }



  ##########################################
  ## PLOT 2: RECONSTRUCTING THE PHYLOGENY ##
  ##########################################

  ## TO DO: ADD PHENOTYPE-COLORING OPTIONS TO RECONSTRUCTED PHYLO PLOTS!
  if(plot2!=FALSE){
    D <- dist.dna(dna[1:n.ind], model="JC69")
  }
  if(plot2=="nj"){
    tree1 <- nj(D)
    #tree1 <- midpoint(ladderize(tree1))
    tree1 <- midpoint(tree1)
    plot(tree1, edge.width=2)
    title("Neighbour-joining tree")
    axisPhylo()
  }

  if(plot2=="UPGMA"){
    tree2 <- hclust(D, method="average")
    tree2 <- as.phylo(tree2)
    #tree2 <- midpoint(ladderize(tree2))
    tree2 <- midpoint(tree2)
    plot(tree2, main="")
    title("UPGMA tree")
  }

  if(plot2=="ml"){
    dna4 <- as.phyDat(dna[1:n.ind])
    tre.ini <- nj(dist.dna(dna[1:n.ind], model="JC69"))
    fit.ini <- pml(tre.ini, dna4, k=n.ind)
    fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE,
                     optQ = TRUE, optGamma = TRUE)

    ## NOTE--you may want to store these in a results.ml list and return it with your results instead of printing
    ## OR at least print a message (eg. "Printing maximum-likelihood calculations...") before printing these numbers...
    #     anova(fit.ini, fit)
    #     AIC(fit.ini)
    #     AIC(fit)

    tree3 <- fit$tree
    #tree3 <- midpoint(ladderize(tree3))
    tree3 <- midpoint(tree3)
    plot(tree3, show.tip=TRUE, edge.width=2)
    title("Maximum-likelihood tree")
    axisPhylo()
  }

  par(ask=FALSE)


  ##################
  ## CONVERT SNPS ##
  ##################

  ## keep and return ONLY genomes for TERMINAL individuals
  genomes <- genomes[1:n.ind]

  x <- genomes
  if(is.null(names(x))) names(x) <- paste("ind", c(1:length(x)), sep=".")
  gen.size <- length(x[[1]])
  ## working with snps in matrix form
  snps <- do.call("rbind", x)
  ## get snps as DNAbin
  if(haploid==TRUE){
    ploidy <- 1
  }else{
    ploidy <- 2
  }
  snps.bin <- as.DNAbin(snps, ploidy=ploidy)
  ## get snps as genind
  #source("C:/Users/caitiecollins/adegenet/R/sequences.R")
  snps.gen <- DNAbin2genind(snps.bin, polyThres=0.01)
  ## get snps as binary matrix
  snps <- snps.gen@tab

  ## correct genind for ploidy if haploid:
  if(haploid==TRUE){
    snps <- snps[,seq(1, ncol(snps), 2)]
    colnames(snps) <- gsub("[:.:]|.$*", "", colnames(snps))
  }

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

    getFixed <- function(locus, posi,
                         exp.char=c("a","t","g","c"),
                         polyThres=0.01){
      vec <- as.character(locus)
      vec[!vec %in% exp.char] <- NA
      N <- sum(!is.na(vec)) # N: number of sequences
      if(N==0 || sum(table(vec)/N >= polyThres )<2){
        return(TRUE) # escape if untyped locus or no SNPs
      }else{
        return(FALSE)
      }
    } # end getFixed

    temp <- lapply(1:ncol(x), function(i) getFixed(x[,i], i)) # process all loci, return a list
    fixed.loci <- which(temp==TRUE) ## identify loci that are NOT SNPs

    ## update snps.assoc to reflect true loci
    gen.size.final <- ncol(snps)
    snps.assoc.loci.ori <- c((gen.size.final-(n.snp.assoc-1)):gen.size.final)

    ## Re-enabled snps.assoc loci "randomization" by
    ## just drawing indices and shuffling the columns accordingly...

    ## draw which SNPs will be associated to the phenotype
    snps.assoc.loci <- sort(sample(c(1:gen.size.final), n.snp.assoc, replace=FALSE))

    snps.indices <- c(1:gen.size.final)
    snps.ori <- snps
    snps.names.ori <- ind.names.ori <- NULL
    if(!is.null(dimnames(snps))){
      snps.names.ori <- dimnames(snps)[[2]]
      ind.names.ori <- dimnames(snps)[[1]]
    }
    snps.non.assoc <- snps[,c(1:(gen.size.final-n.snp.assoc))]
    snps.assoc <- snps[,snps.assoc.loci.ori]
    snps.new <- matrix(99, nrow=100, ncol=gen.size.final)
    snps.new[,snps.indices[-snps.assoc.loci]] <- snps.non.assoc
    snps.new[,snps.assoc.loci] <- snps.assoc
    snps <- snps.new
    if(!is.null(snps.names.ori)) colnames(snps) <- snps.names.ori
    if(!is.null(ind.names.ori)) rownames(snps) <- ind.names.ori
    snps.assoc <- snps.assoc.loci

  }

  ## update gen.size
  if(haploid==FALSE){
    ## update "gen.size" to reflect ACTUAL dim of snps output
    gen.size <- as.numeric(gsub("L|.$|[:.:]", "", colnames(snps)[ncol(snps)]))
  }else{
    gen.size <- ncol(snps)
  }



  ##################
  ## get RESULTS: ##
  ##################

  if(is.null(phen.branch)){

    ## If NO phen sim: ##
    #     out <- list(genomes, tree)
    #     names(out) <- c("genomes", "tree")
    out <- list(snps, tree)
    names(out) <- c("snps", "tree")

  }else{

    ## If PHEN HAS been SIMULATED: ##

    ## get all info relevant to plotting phenotype with colored phylo:
    phen.plot.colors <- list(edgeLabCol, edgeCol, nodeCol, internalNodeCol, leafCol)
    names(phen.plot.colors) <- c("edge.labels", "edges", "all.nodes", "internal.nodes", "tip.labels")

    ## get all info on simulated phenotype:
    node.noms <- names(phen.nodes)
    phen.nodes <- as.vector(unlist(phen.nodes))
    names(phen.nodes) <- node.noms

    ## update names of snps.assoc
    names(snps.assoc) <- as.vector(unlist(sapply(c(1:length(snps.assoc)),
                                                 function(e)
                                                   paste("L",
                                                         paste(rep(0, (nchar(gen.size)-nchar(snps.assoc[e]))), sep="", collapse=""),
                                                         snps.assoc[e], sep=""))))

    phen.loci.ori <- phen.loci
    phen.vars <- list(phen.loci, phen.loci.ori, phen.branch, phen.nodes, phen.leaves, snps.assoc)
    names(phen.vars) <- c("phen.sub.branches_effective", "phen.sub.branches", "edges", "nodes", "leaves", "snps.assoc")

    ## TO DO: change output to get rid of old phen.loci.ori (now always == phen.loci)
    ## .. want to check this wont mess up any indexing elsewhere...
    #     phen.vars <- list(phen.loci,  phen.branch, phen.nodes, phen.leaves, snps.assoc)
    #     names(phen.vars) <- c(phen.sub.branches", "edges", "nodes", "leaves", "snps.assoc")

    phenotype <- list(phen.plot.colors, phen.vars)
    names(phenotype) <- c("plot.colors", "variables")

    ## add phen info to results list
    #     out <- list(snps, snps.gen, tree, phenotype)
    #     names(out) <- c("snps", "snps.gen", "tree", "phenotype")
    out <- list(snps, tree, phenotype)
    names(out) <- c("snps", "tree", "phenotype")
  }


  return(out)

} # end coalescent.sim


##########################################################################
