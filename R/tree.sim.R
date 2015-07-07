 

##############
## tree.sim ##
##############


# n.ind <- 10 # n.genomes you want to end up with
# gen.size <- 1000000 # bases
# theta <- gen.size*2 # (if sim.by=="branch")# OR # 1*2 # (if sim.by=="locus") 
# biallelic <- TRUE # if TRUE, select ONLY complementary nt; 
###########   if FALSE, select from 3 alternatives (ie. A/C/G/T-current nt)
# seed <- 1 # allow user to control randomization to get reproducible results.
# sim.by <- c("branch", "locus")

#tree <- tree.indian
#out <- tree.sim(n.ind=length(tree$tip.label), 
# tree=tree, gen.size=10000, theta=100*2, 
# biallelic=FALSE, seed=2, 
# plot=TRUE, heatmap=FALSE, plot2="ml")

tree.sim <- function(n.ind=length(tree$tip.label), 
                     tree=tree, gen.size=10000, 
                     sim.by="locus", theta=NULL, 
                     theta_p=NULL, phen=NULL, 
                     biallelic=TRUE, seed=NULL, plot=TRUE,
                           heatmap=FALSE, plot2="UPGMA"){
  
  require(adegenet)
  require(ape)
  require(phangorn)
  
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
  ## simulate genotype for root individual:
  gen.root <- sample(c("a", "c", "g", "t"), 
                     gen.size, replace=TRUE)
  
  ## simulate phenotype for root individual:
  if(!is.null(theta_p)){
    phen.root <- "A"
  }else{
    phen.root <- NULL
  }
  
  ## get the sum of all branch lengths in the tree:
  time.total <- sum(tree$edge.length)
  
  ## make dummy variables in which to store 
  ## the resulting n.mts variables:
  L <- lambda <- lambda_p <- 
    n.subs <- n.mts <- new.nt <- NA
  
  
  ####################################
  ## PHENOTYPE simulation procedure ## ~ sim.by.locus...
  ####################################
  
  ## ensure phen variables start as NULL
  phen.branch <- phen.nodes <- phen.leaves <- NULL
    
  ## If the user has specified a 
  ## "mt" rate for the phenotype
  ## (indicating that they want to 
  ## generate a NEW phenotype for the tree provided)
  if(!is.null(theta_p)){  
    
    ## draw the number of mutations to occur at each site:    
    n.subs <- rpois(n=1, lambda=theta_p)
    ## if n.subs==0 or ==1, re-sample    
    while(n.subs <= 1){
      n.subs <- rpois(n=1, lambda=theta_p)
    }  
    
    ## draw the branches to which you will 
    ## assign the n.subs to occur for the phenotype 
    ## (~ branch length):
    phen.loci <- sample(c(1:length(tree$edge.length)), 
                        n.subs, 
                        replace=TRUE, 
                        prob=tree$edge.length)
    ## rearrange phen.loci
    phen.loci <- sort(phen.loci, 
                      decreasing=TRUE)
    
    #####################
    ## .switch.phen fn ##
    #####################
    .switch.phen <- function(x){
      out <- NULL
      if(x == "A") out <- "B"
      if(x == "B") out <- "A"
      return(out)
    } # end .switch.phen
    
    #################################################
    ## deal with multiple phen subs on same branch ##
    
    ## NOTE--If you want to keep ALL switches 
    ## (even if some have no effect on nodes), 
    ## you could do so by:
    ## 1) using this next part to get a 
    ## COUNT of changes, WITHOUT keeping 
    ## only changes that affect nodes
    ## 2) change phen.branch "A" --> c("B", "A") 
    ## part to reflect the count for that edge
    ## 3) Use part II of the next bit to 
    ## KEEP ONLY changes that affect nodes...
    #### .../ OR / ONLY look at 
    ## phen.branch[[phen.loci[i]]][1] and "[last]...
    #### ... to assign phen to nodes in 
    ## FROM and TO columns of relevant tree$edge row
    
    ########################
    ## Parity-testing fns ##
    ########################
    .is.even <- function(x) x %% 2 == 0 
    .is.odd <- function(x) x %% 2 != 0 
    
    ## if phen sub branches occur in EVEN 
    ## numbered multiples, REMOVE; if ODD, keep ONE:
    phen.loci.ori <- phen.loci
    for(i in 1:length(unique(phen.loci))){
      n <- length(which(phen.loci == unique(phen.loci)[i]))
      #if(.is.odd(n)) phen.loci[which(phen.loci == 
      #     unique(phen.loci)[i])] 
      #     <- unique(phen.loci)[i]
      if(.is.odd(n)){
        ## replace odd-numbered reps of a 
        ## phen.loci with ONE instance of that phen.loci
        phen.loci <- replace(phen.loci, 
                             which(phen.loci == 
                                     unique(phen.loci)[i]), 
                             c(unique(phen.loci)[i], 
                               rep(NA, 
                                   (length(which(phen.loci == 
                                     unique(phen.loci)[i]))-1)
                               )))
      } 
      if(.is.even(n)) phen.loci[which(phen.loci == 
                                      unique(phen.loci)[i])] <- NA
      phen.loci <- phen.loci[!is.na(phen.loci)]
    } # end for loop #################################
    
    
    ###############################
    ## For Loop to get PHENOTYPE ##
    ###############################
    ## get phenotype for all branches/ nodes in tree 
    ## (from root node (ie. tree$edge[nrow(tree$edge), 1]) down):
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
        phen.nodes[[tree$edge[x[i],2]]] <- 
          .switch.phen(phen.nodes[[tree$edge[x[i],1]]])
      }else{
        ## if no phen subs occur on branch i, set phen of 
        ## downstram individual to be equal to ancestor's
        phen.nodes[[tree$edge[x[i],2]]] <- 
          phen.nodes[[tree$edge[x[i], 1]]]
      }
    } # end for loop
    
    ## get phen of TERMINAL nodes (leaves)
    phen.leaves <- 
      as.factor(as.vector(unlist(phen.nodes[c(1:n.ind)])))
    names(phen.leaves) <- 
      paste("ind", c(1:length(phen.leaves)), sep=".")
    
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
    ## FOR NOW--all edges w > 1 phen 
    ## (either c("A", "B") OR c("B", "A")) --> purple... 
    l_edgeCol <- which(sapply(c(1:length(edgeCol)), 
                              function(e) 
                                length(edgeCol[[e]])) == 2)
    edgeCol <- replace(edgeCol, 
                       l_edgeCol, 
                       "purple")
    edgeCol <- replace(edgeCol, 
                       which(edgeCol == "A"), 
                       "blue")
    edgeCol <- replace(edgeCol, 
                       which(edgeCol == "B"), 
                       "red")
    edgeCol <- as.vector(unlist(edgeCol))
    ## get COLOR for EDGE LABELS
    edgeLabCol <- rep("green", 
                      nrow(tree$edge))
    edgeLabCol <- replace(edgeLabCol, 
                          phen.loci, 
                          "mediumorchid1")
    
    ###############
    ## plot TREE ##
    ###############
    if(plot==TRUE){
      if(n.ind <= 20){
        plot(tree, show.tip=FALSE, 
             edge.width=2, 
             edge.color=edgeCol) # edgeCol
        title("Coalescent tree w/ phenotypic changes")
        axisPhylo()
        edgelabels(text=paste("e", 
                              c(1:nrow(tree$edge)), 
                              sep="."), 
                   cex=0.5, font=2, 
                   bg=transp(edgeLabCol, 0.3), adj=c(1,1))
        tiplabels(text=tree$tip.label, 
                  cex=0.6, adj=c(-0.5, 0), 
                  bg=transp(leafCol, 0.3)) 
        nodelabels(text=rev(unique(tree$edge[,1])), 
                   cex=0.5, 
                   bg=transp(internalNodeCol, 0.3))  
      }else{
        plot(tree, show.tip=FALSE, 
             edge.width=2, 
             edge.color=edgeCol) # edgeCol
        title("Coalescent tree w/ phenotypic changes")
        axisPhylo()
        #edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."), 
        #   cex=0.5, font=2, bg=transp(edgeLabCol, 0.3), adj=c(1,1))
        tiplabels(text=tree$tip.label, 
                  cex=0.6, adj=c(-0.5, 0), 
                  col=leafCol, frame="none")
        #nodelabels(text=rev(unique(tree$edge[,1])), 
        #  cex=0.5, bg=transp(internalNodeCol, 0.3)) 
      }      
    } # end plot = TRUE
    
  } # end if(!is.null(theta_p)) 
  ## ie. SIMULATED phenotype & plotting
  
  #######################################################################
  
  ## If the user has PROVIDED a phenotype
  ## (for terminal nodes only)
  if(!is.null(phen)){
    # phen <- as.factor(sample(c("A", "B", "C", "D"), 
    #            100, replace=TRUE))
    
    ## get COLOR for LEAVES
    n.levels <- length(levels(as.factor(phen)))
    leafCol <- get("funky")(n.levels)
    scheme <- as.numeric(phen)
    leafCol <- leafCol[scheme]
    
    ## get COLOR for EDGES
    edgeCol <- "black" ## for NOW...
    ########
    ## TO DO-- color ~ terminal edges red/blue 
    ## same as phen of terminal node...
    #### ... UNTIL two edge colors meet at any 
    ## internal node (then just black edges to root)
    
    ###############
    ## plot TREE ##
    ###############
    if(plot==TRUE){
      plot(tree, show.tip=FALSE, 
           edge.width=2, 
           edge.color=edgeCol) 
      title("Coalescent tree w/ phenotypes of leaves")
      axisPhylo()
      tiplabels(text=tree$tip.label, 
                cex=0.6, adj=c(-0.5, 0), 
                bg=transp(leafCol, 0.3)) 
    } # end plot = TRUE
  } # end if(!is.null(phen)) 
  ## ie. PROVIDED phenotype & plotting
  
  
  #######################################################################
  
  #######################################################################
  
  
  
  ##################################
  ## GET MUTATIONS' branch & loci ##
  ##################################
  
  ###############
  ## BY BRANCH ## 
  ###############
  ## NOTE: ## for sim.by BRANCH, MT.RATE ~= 
  ## the EXPECTED n.mts TOTAL 
  ## (in the WHOLE GENOME) 
  ## (ie. gen.size(*2), eg. 1000(*2))!
  #####################################
  if(sim.by=="branch"){
  
  ## generate mts for each of n.ind genomes: 
  ## NOTE: We're working from the root down, 
    ## so from the bottom row of 
  ## tree$edge up to the top row --> Mts accumulate,
    ## w/ more Mts on longer branches... 
  ## NOTE2: The placement of the root is 
    ## arbitrary/ irrelevant to the outcome w.r.t the genetic
  ## distance/ relationship between the nodes...

  ##################################
  ## REGULAR simulation procedure ##
  ##################################
  
  snps.loci <- list()
  
  ## draw n.mts per branch and their loci 
  for(i in rev(1:length(tree$edge.length))){  
    ## get branch length i 
    L[i] <- rev(tree$edge.length)[i]
    
    ## get lambda, the mean of the poisson dist 
    ## from which we draw the n.mts to occur on branch i 
    ## (proportional to branch length i out of 
    ## the total sum of all branch lengths in the tree)
    lambda[i] <- ((theta/2) * (L[i] / time.total))
    
    ## get n.mts from poisson dist: 
    n.mts[i] <- rpois(n=1, lambda[i])
    
    ## draw n.mts locations in the genome 
    ## (allowing for replacement)
    snps.loci[[i]] <- sample(c(1:gen.size), n.mts[i], replace=TRUE)  
  } # end for loop

  } # end sim.by branch
  
  
  
  
  
  ##############
  ## BY LOCUS ## 
  ##############
  ## NOTE: ## for sim.by LOCUS, MT.RATE ~= 
  ## the EXPECTED n.mts AT EACH SITE (eg. 1(*2)) !
  ###################################################
  if(sim.by=="locus"){
    
    #####################################
    ## SNP-BY-SNP simulation procedure ##
    #####################################
    
    ## draw the number of mutations to occur at each site:    
    n.mts <- rpois(n=gen.size, lambda=(theta/2))
    ## for any n.mts==0, re-sample
    for(i in 1:length(n.mts)){      
      while(n.mts[i]==0){
        n.mts[i] <- rpois(n=1, lambda=(theta/2))
      }
    }      
    
    ## for each site, draw the branches to which 
    ## you will assign the mts for this site 
    ## (~ branch length):
    snps.loci <- sapply(c(1:length(n.mts)), 
                        function(e) 
                          sample(c(1:length(tree$edge.length)), 
                                 n.mts[e], 
                                 replace=TRUE, 
                                 prob=tree$edge.length))
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
  
  
  ## mini fn testing for output==integer(0) ##
  .is.integer0 <- function(x){
    is.integer(x) && length(x) == 0L
  } # end .is.integer0()
  
  
  # we will store the output in a list called genomes:
  genomes <- list()
  ## get the node names for all individuals 
  ## (terminal and internal)
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
        new.nts[[i]] <- sapply(c(1:length(snps.loci.unique[[i]])), 
                               function(e) 
          sample(c("a", "c", "g", "t")[-which(c("a", "c", "g", "t") 
                                    %in% genomes[[tree$edge[i,1]]]
                                    [snps.loci.unique[[i]][e]])], 1))
      }else{
        new.nts[[i]] <- sapply(c(1:length(snps.loci.unique[[i]])), 
                               function(e) 
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
                foo[[j]][k] <- selectBiallelicSNP(c("a", "c", "g", "t")
                                                  [which(c("a", "c", "g", "t") 
                                                            %in% foo[[j]][1])])
              }
            }else{
              if(biallelic==FALSE){
                foo[[j]][k] <- sample(c("a", "c", "g", "t")
                                      [-which(c("a", "c", "g", "t") 
                                              %in% foo[[j]][k-1])], 1)
              }else{
                foo[[j]][k] <- selectBiallelicSNP(c("a", "c", "g", "t")
                                                  [which(c("a", "c", "g", "t") 
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
  
  
  
  
    
  if(heatmap == TRUE || plot2!=FALSE){
    dna <- as.DNAbin(genomes)
    names(dna) <- c(1:length(genomes))
  }
  
  
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
  ## NOTE-- we only print the plot here IF we have 
  ## NOT generated a phenotype OR read one in
  ########## (thereby having printed it earlier 
  ## in the phen section)
  if(plot == TRUE && is.null(theta_p) && is.null(phen)){
    if(sim.by=="branch"){
      plot(tree, show.tip=FALSE, edge.width=2)
      title("Coalescent tree")
      axisPhylo()
      tiplabels(text=tree$tip.label, 
                cex=1, adj=-.5) 
      nodelabels(text=rev(unique(tree$edge[,1])), 
                 cex=0.75)
      edgelabels(text=paste("e", 
                            c(1:nrow(tree$edge)), 
                            sep="."), 
                            cex=0.66)
      edgelabels(text=rev(n.mts), 
                 col="red", 
                 frame="none", 
                 cex=1.1, 
                 adj=c(1,-0.5))
    } # end sim.by branch
    
    if(sim.by=="locus"){
      plot(tree, show.tip=FALSE, edge.width=2)
      title("Coalescent tree")
      axisPhylo()
      tiplabels(text=tree$tip.label, cex=1, adj=-.5) 
      nodelabels(text=rev(unique(tree$edge[,1])), cex=0.75)
      edgelabels(text=paste("e", 
                            c(1:nrow(tree$edge)), 
                            sep="."), cex=0.66)
      edgelabels(text=sapply(c(1:length(snps.loci)), 
                             function(e) 
                               length(snps.loci[[e]])), 
                               col="red", frame="none", 
                               cex=1.1, adj=c(1,-0.5))
    } # end sim.by locus
    
  } # end plot (regular) = TRUE (for non-phen-sim runs)
  
  ###################
  ## PLOT (simple) ##
  ###################
  if(plot=="simple"){
    plot(tree, show.tip=TRUE, edge.width=2, 
         cex=0.6, adj=.25)
    title("Coalescent tree") 
  } # end plot = simple
  
  
  
  
  ##########################################
  ## PLOT 2: RECONSTRUCTING THE PHYLOGENY ##
  ##########################################
  
  ## TO DO: ADD PHENOTYPE-COLORING OPTIONS TO 
  ## RECONSTRUCTED PHYLO PLOTS!
  if(plot2!=FALSE){
    D <- dist.dna(dna[1:n.ind], model="JC69")
  }

  if(plot2=="nj"){
    tree1 <- nj(D)    
    tree1 <- midpoint(tree1)
    plot(tree1, edge.width=2)
    title("Neighbour-joining tree")
    axisPhylo()
  }
  
  if(plot2=="UPGMA"){
    tree2 <- hclust(D, method="average") 
    tree2 <- as.phylo(tree2)
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
    
    ## NOTE--you may want to store these in a 
    ## results.ml list and return it with 
    ## your results instead of printing
    ## OR at least print a message 
    ## (eg. "Printing maximum-likelihood calculations...") 
    ## before printing these numbers... 
    #     anova(fit.ini, fit)
    #     AIC(fit.ini)
    #     AIC(fit) 
    
    tree3 <- fit$tree
    tree3 <- midpoint(tree3)
    plot(tree3, show.tip=TRUE, edge.width=2)
    title("Maximum-likelihood tree")
    axisPhylo()
  }
  
  par(ask=FALSE)
  
  ## keep and return ONLY genomes for 
  ## TERMINAL individuals
  genomes <- genomes[1:n.ind]  
  
  ##################
  ## CONVERT SNPS ##
  ##################
  x <- genomes
  if(is.null(names(x))) names(x) <- 
    paste("ind", c(1:length(x)), sep=".")
  gen.size <- length(x[[1]])
  ## working with snps in matrix form
  snps <- do.call("rbind", x)
  ## get snps as DNAbin 
  snps.bin <- as.DNAbin(snps)
  ## get snps as genind
  #source("C:/Users/caitiecollins/adegenet/R/sequences.R")
  snps.gen <- DNAbin2genind(snps.bin, 
                            polyThres=0.01)
  ## get snps as binary matrix
  snps <- snps.gen@tab  
  
  ## update "gen.size" to reflect 
  ## ACTUAL dim of snps output
  gen.size <- as.numeric(gsub("L|.$|[:.:]", "", 
                          colnames(snps)[ncol(snps)]))


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
    phen.plot.colors <- list(edgeLabCol, 
                             edgeCol, 
                             nodeCol, 
                             internalNodeCol, 
                             leafCol)
    names(phen.plot.colors) <- c("edge.labels",
                                 "edges",
                                 "all.nodes",
                                 "internal.nodes",
                                 "tip.labels")
    
    ## get all info on simulated phenotype:
    node.noms <- names(phen.nodes)
    phen.nodes <- as.vector(unlist(phen.nodes))
    names(phen.nodes) <- node.noms
    phen.vars <- list(phen.loci,
                      phen.loci.ori, 
                      phen.branch, 
                      phen.nodes, 
                      phen.leaves)
    names(phen.vars) <- c("phen.sub.branches_effective",
                          "phen.sub.branches",
                          "edges",
                          "nodes",
                          "leaves")
    
    phenotype <- list(phen.plot.colors, 
                      phen.vars)
    names(phenotype) <- c("plot.colors",
                          "variables")
    
    ## add phen info to results list
    #     out <- list(genomes, tree, phenotype)
    #     names(out) <- c("genomes", "tree", "phenotype")
    out <- list(snps, tree, phenotype)
    names(out) <- c("snps", "tree", "phenotype")
    
  }
  
  
  return(out)
  
} # end coalescent.sim





#########################
## selectBiallelicSNP: ##
#########################

## fn that returns the alternative nt| the nt inut

## NOTE-- while this is not inherently the 
## definition of a biallelic SNP, 
## it is currently suiting my purposes
#### by fulfilling the function of 
## ensuring that sites always revert 
## back and for between one state and ONE other,
#### hence never creating any 
## triallelic sites or tetralellic sites 
## --> binary encoding guaranteed to work fine. 

selectBiallelicSNP <- function(x, DNA=TRUE){
  
  out <- NULL
  
  if(!is.null(x)){
    
    x <- as.character(x)
    
    ## for DNA encoding:
    if(DNA==TRUE){
      if(x=="a") out <- "t"
      if(x=="t") out <- "a"
      if(x=="c") out <- "g"
      if(x=="g") out <- "c"
      
      if(x=="A") out <- "T"
      if(x=="T") out <- "A"
      if(x=="C") out <- "G"
      if(x=="G") out <- "C"
    } # end DNA
    
    ## for RNA encoding: 
    if(DNA==FALSE){  
      if(x=="a") out <- "u"
      if(x=="u") out <- "a"
      if(x=="c") out <- "g"
      if(x=="g") out <- "c"
      
      if(x=="A") out <- "U"
      if(x=="U") out <- "A"
      if(x=="C") out <- "G"
      if(x=="G") out <- "C"
    } # end RNA
  } # end !is.null(x)
  
  return(out)
  
} # end selectBiallelicSNP
#######################################################