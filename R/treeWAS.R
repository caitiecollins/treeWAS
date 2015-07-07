#############
## treeWAS ##
#############


## Re-implementation of GWAS method developed in Sheppard et al 2014,
## based on the generation of a phylogenetically-correct p-value
## derived by comparing correlation btw SNPs
## and a phenotype of interest to a null distribution 
## (Poisson dist w parameter 1) OR 
## (modification based on similar method developed in Farhat et al 2013)
## based on permutation of empirically-derived
## n.mts per site (except we determine this with the Fitch algorithm). 

########################################################################


treeWAS <- function(x, y, tree=NULL, 
                    mt.rate=NULL, 
                    p.value=0.001, mt.correct=FALSE, p.value.option="count",
                    plot.null.dist=TRUE, plot.dist=FALSE, plot.Fitch=TRUE, 
                    plot.tree=FALSE, tree.method="ml", 
                    n.reps=1, sim.gen.size=10000, sim.by="locus", 
                    method="Didelot", test="score", corr.dat.fisher=NULL){  
  
  
  ## load packages:
  require(adegenet)
  require(phangorn)
  require(adephylo)
  require(ape) 
  require(ade4) #?
  
  ## load my simulation fn:
  source("C:/Cait 2012/Work/Xavier/treeWAS Method/tree.sim.R")
  ## load fitch simulation/randomization fn:
  source("C:/Cait 2012/Work/Xavier/treeWAS Method/fitch.sim.R")  
  
  #####################################################################
  ## 0) Handle input data #############################################
  #####################################################################
  ## get & check snps
  if(!is.matrix(x)) x <- as.matrix(x)
  snps <- x
  
  ## convert phenotype to factor
  y <- as.factor(y)
  
  ## set n.ind:
  n.ind <- length(y)
  inds <- c(1:n.ind)

  
  #   if(is.null(tree)){
  #     ## reconstruct phylogenetic tree if necessary:
  #     
  #     ## if no tree.method selected, use maximum-likelihood as default:
  #     if(is.null(tree.method)) tree.method <- "ml"
  #     
  #     ## working with data in DNAbin format...
  #     ## get distance matrix for terminal node individuals
  #     D <- dist.dna(snps.bin[1:n.ind,], model="JC69")
  #     
  #     ## MAXIMUM-LIKELIHOOD TREE
  #     if(tree.method=="ml"){
  #       dna4 <- as.phyDat(snps.bin[1:n.ind,])
  #       tre.ini <- nj(D)
  #       fit.ini <- pml(tre.ini, dna4, k=n.ind)
  #       fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, 
  #                        optQ = TRUE, optGamma = TRUE)
  #       ## STORE RESULTS for maximum-likelihood tree reconstruction:
  #       ml.results <- list(anova(fit.ini, fit), 
  #                          AIC(fit.ini), 
  #                          AIC(fit))
  #       names(ml.results) <- 
  #         c("ANOVA for initial vs. optimal fit", 
  #           "AIC for initial fit", 
  #           "AIC for optimal fit")
  #       
  #       tree <- fit$tree
  #       tree <- midpoint(ladderize(tree))
  #       
  #       ## plot tree
  #       if(plot.tree==TRUE){
  #         plot(tree, show.tip=TRUE, edge.width=2)        
  #         title("Maximum-likelihood tree")
  #         axisPhylo()
  #       } # end plot.tree ml
  #     } # end max-likelihood tree
  #     
  #     
  #     ## NEIGHBOUR-JOINING TREE
  #     if(tree.method=="nj"){
  #       tree <- nj(D)
  #       tree <- midpoint(ladderize(tree))
  #       ## plot tree
  #       if(plot.tree==TRUE){
  #         plot(tree, edge.width=2)
  #         title("Neighbour-joining tree")
  #         axisPhylo()
  #       } # end plot.tree nj
  #     } # end neighbour-joining tree
  #     
  #     
  #     ## UPGMA TREE
  #     if(tree.method=="UPGMA"){
  #       tree <- hclust(D, method="average") 
  #       tree <- as.phylo(tree)
  #       tree <- midpoint(ladderize(tree))
  #       ## plot tree
  #       if(plot.tree==TRUE){
  #         plot(tree, main="")
  #         title("UPGMA tree")
  #       } # end plot.tree UPGMA
  #     } # end UPGMA tree    
  #     
  #     # end (optional) tree reconstruction 
  #     # (and (optional) plotting of reconstructed tree)
  #   }else{
  #     
    ## If user has already submitted a tree as input:    
    ## Work with a centered phylo tree for 
    ## consistency and visualisation's sake:
    if(class(tree) != "phylo") tree <- as.phylo(tree)
    ## if the tree is not already rooted, root it:
    if(is.rooted(tree)==FALSE) tree <- midpoint(tree)
    
    
    if(plot.tree==TRUE){
      plot(tree)
      title("Phylogenetic tree (original)")
      axisPhylo()
    } # end plot.tree
  #   } # end tree...
  
  ## tree's tip.labels must be numeric (for Fitch parsimony step)
  tree.ori <- tree
  tree$tip.label <- c(1:length(tree$tip.label))
  tree$node.label <- c((n.ind+1):(n.ind+tree$Nnode))


  ## MT.RATE ##
  
  ## set default mt rate
  if(missing(mt.rate)) mt.rate <- NULL
  
  if(sim.by=="branch"){
    ## USE FITCH.N.MTS AS THETA (MT RATE) IN treeWAS...
    if(is.null(mt.rate)){ # mt.rate <- 1  
      snps.phyDat <- as.phyDat(as.matrix(snps), 
                               type="USER", levels=c(0,1))
      mt.rate <- parsimony(tree, snps.phyDat, method="fitch") 
      ## --> single-number parsimony score 
      ## (for indian crypto = 1658889)
    }
  } # end sim.by branch
  
  if(sim.by=="locus"){
    if(is.null(mt.rate)) mt.rate <- 1
    } # end sim.by locus
  
  
  
  
  #####################################################
  ## 1) Simulate multiple trees to compare your real/ #
  ## reconstructed original tree to ###################
  #####################################################
  
  #####################################################
  #### A) Sheppard Method: ############################
  ##- Generate N comparator datasets using Mt.rate = ##
  ## theta/2 = 1/(sum of branch lengths) --> n.mts = ##
  ## Poisson dist'd w parameter 1######################
  #####################################################
  
  if(method=="Didelot"){
    
    out <- list()
    genomes <- list()
    trees <- list()
    snps.mat <- list()
    
    for(i in 1:n.reps){
      if(is.null(sim.gen.size)) sim.gen.size <- 
        gen.size
      if(sim.by=="branch"){
        if(sim.gen.size != gen.size) mt.rate <- 
          mt.rate*(sim.gen.size/gen.size)
      }
      
      ## SIMULATE A DATASET | your tree ##
      out[[i]] <- tree.sim(n.ind=n.ind, tree=tree, 
                           gen.size=sim.gen.size, 
                           sim.by=sim.by, theta=mt.rate*2, 
                           theta_p=NULL, phen=NULL,
                           biallelic=TRUE, seed=NULL, 
                           plot=FALSE, heatmap=FALSE, plot2=FALSE)
      
      
      genomes[[i]] <- out[[i]][[1]]
      trees[[i]] <- out[[i]][[2]]
      
      ## Modify genomes/snps matrices
      if(!is.null(genomes[[i]])){
        snps.mat[[i]] <- genomes[[i]]
        ## keeping only every other column of simulate SNPs (bc. all binary)
        ## NOTE: handling this in coalescent.sim now w arg haploid:
        #snps.mat[[i]] <- snps.mat[[i]][,seq(1, ncol(snps.mat[[i]]), 2)] 
        
      }else{
        snps.mat[[i]] <- NULL
      }
      
      gc()
      
    } # end for loop
    
  } # end Didelot method
  
  ## check avg n.mts
  #mean(sapply(c(1:length(snps.mat)), 
  #function(e) ncol(snps.mat[[e]]))) # /2 # 
  
  ####################################################################
  #### B) Farhat Method: NOT UP TO DATE!!!!!##########################
  ###### - Modify the method implemented in Farhat et al 2013 ########
  ###### - Determine n.mts per site* with the Fitch algorithm ########
  ############ *(ie. The site-specific n.mts) ########################
  ###### - Assign this n.mts to new branches at random ###############
  ###### (with Prob corr. branch length) (repeat Nx to get N trees)###
  ####################################################################
  
  if(method=="Fitch"){
    
    #########################################################
    ## get n.mts per site with fitch's parsimony algorith: ##
    #########################################################
    snps.phyDat <- list()
    fitch.n.mts <- list()
    
    toKeep <- seq(1, ncol(snps), 2)
    
    if(length(toKeep)!=0){
      for(i in 1:length(toKeep)){      
        snps.phyDat[[i]] <- as.phyDat(as.matrix(snps[,toKeep[i]]), 
                                      type="USER", levels=c(0,1))
        fitch.n.mts[[i]] <- parsimony(tree, snps.phyDat[[i]], 
                                      method="fitch")
        ###########################################################
        ## length(toKeep) 113382 --> system.time = 15 hours! :(  ##
        ###########################################################
        gc()
        
      } # end for loop
    }
    
    fitch.n.mts <- as.vector(unlist(fitch.n.mts))
    
    ## inspect distribution of n.mts:
    #table(fitch.n.mts)
    
    
    ####################################################
    ## Randomize this n.mts per site onto ##############
    ## different branches with prob ~ branch lenghts: ##
    ####################################################
    
    out <- list()
    genomes <- list()
    trees <- list()
    snps.mat <- list()
    
    for(i in 1:n.reps){
      if(is.null(sim.gen.size)) sim.gen.size <- 
        gen.size
      out[[i]] <- fitch.sim(snps, fitch.n.mts, 
                            tree, gen.size=sim.gen.size,
                            biallelic=TRUE, seed=NULL)
      
      genomes[[i]] <- out[[i]][[1]]
      
      ## Convert genomes into snps matrices
      dat <- genomes[[i]]
      dat <- do.call("rbind", dat)
      dat.bin <- as.DNAbin(dat)
      dat.gen <- DNAbin2genind(dat, polyThres=0.01)
      if(!is.null(dat.gen)){
        snps.mat[[i]] <- dat.gen@tab
        noms <- sapply(c(1:length(dat.gen@loc.names)), function(e) 
          paste(dat.gen@loc.names[e], dat.gen@all.names[[e]], sep="."))    
        colnames(snps.mat[[i]]) <- noms
      }else{
        snps.mat[[i]] <- NULL
      }
      
      gc()
      
    } # end for loop
    
  } # end Fitch/ modified-Farhat method
  
  
  ################################################################
  ## 3) Get results:##############################################
  #### Determine the PC-p-values for all SNP loci | ##############
  ##   null distributions of correlations from simulated data ####
  #### Synthesize results output: List of all significant SNPs, ##
  ##   their names/locations, their p-values for this phenotype ##
  ################################################################
  
  if(method=="Didelot"){
    ## Calculate correlation between simulated SNPs and phenotype, 
    ## then between real SNPs and phenotype:
    
    ## convert phenotype factor(y) to numeric:
    phen <- as.numeric(y)
    ## for ease of interpretation, 
    ## if phen has 2 levels, 1 and 2, 
    ## make these 0 and 1:
    if(length(unique(phen))==2){
      if(length(phen[-c(which(phen==1), which(phen==2))])==0){
        phen <- replace(phen, which(phen==1), 0)
        phen <- replace(phen, which(phen==2), 1)
      }
    }  
    
    ################################################
    ## calculate a correlation btw the (simulated) #
    ## SNPs in matrix[[i]] and the phenotype  ######
    ################################################
    
    ## create empty list to store correlation results 
    ## for each of the list of simulated SNPs matrices
    corr.sim <- list()
    
    ## for loop to get correlations:
    for(i in 1:length(snps.mat)){
      ## isolate this element of the list of 
      ## the simulated SNPs matrices
      foo <- snps.mat[[i]]
      
      ## get correlation ##      
      if(!is.null(foo)){
        
        ## CORRELATION ##
        if(test=="cor"){
          # regular correlation...
          all.corr <- cor(foo, phen)   
        } # end test cor
        
        ## SCORE ##
        if(test=="score"){
          # ~ Correlation "SCORE" = 
          # ((nS1P1 + nS0P0) - (nS1P0 + nS0P1) / (n.total)) 
          ## must be calculated for each SNP individually...   
          all.corr <- sapply(c(1:ncol(foo)), function(e) 
            (((length(which(foo[which(phen==1),e]==1)) + 
                 length(which(foo[which(phen==0),e]==0)))
              - (length(which(foo[which(phen==1),e]==0)) + 
                   length(which(foo[which(phen==0),e]==1))))
             / nrow(foo))) 
          ## hist(all.corr) 
          ## should be symmetric about 0, and in [0,1]...
        } # end test score
        
        ## FISHER'S EXACT TEST ##
        if(test=="fisher"){
          # fisher's exact test? ##
          all.corr <- sapply(c(1:ncol(foo)), function(e) 
            fisher.test(foo[,e], y=phen, 
                        alternative="two.sided")$p.value) 
          ## two.sided bc we want to know if inds w 
          ## the phen have EITHER more 1s or 0s 
        } # end test fisher
        
        corr.sim[[i]] <- abs(all.corr)
        
      }else{
        corr.sim[[i]] <- NULL
      }
    } # end for loop
    
    
    
    ############################################################
    ## Calculate correlations between real SNPs and phenotype ##
    ############################################################
    
    ## CORRELATION ##
    if(test=="cor"){
      corr.dat <- sapply(c(1:ncol(snps)), function(e) 
        cor(snps[,e], phen)) # regular correlation...
    } # end test cor
    
    ## SCORE ##
    if(test=="score"){
      # ~ Correlation "SCORE" = 
      # ((nS1P1 + nS0P0) - (nS1P0 + nS0P1) / (n.total)) 
      ## must be calculated for each SNP individually...   
      corr.dat <- sapply(c(1:ncol(snps)), function(e) 
        (((length(which(snps[which(phen==1),e]==1)) + 
             length(which(snps[which(phen==0),e]==0)))
          - (length(which(snps[which(phen==1),e]==0)) + 
               length(which(snps[which(phen==0),e]==1)))) 
         / nrow(snps)))
    } # end test score
    
    ## FISHER'S EXACT TEST ##
    if(test=="fisher"){
      corr.dat <- sapply(c(1:ncol(snps)), 
                         function(e) fisher.test(snps[,e], 
                         y=phen, alternative="two.sided")$p.value) 
      ## two.sided bc we want to know if inds w 
      ## the phen have EITHER more 1s or 0s      
    } # end test fisher
    
    
    ## USE ABSOLUTE VALUE
    corr.dat <- abs(corr.dat)
    
    
    ## get empirical p-value based on corr.sim ##
    ## INTERPRETATION: you want the threshold above which 
    ## only 0.001 proportion (ie. 0.1%) of the data lies:
    
    corr.sim <- as.vector(unlist(corr.sim))
    
    ## handle p-value input
    if(mt.correct == TRUE) if(p.value != "max") p.value <- p.value/length(corr.dat)
    
    
    ## CORRELATION ##
    if(test=="cor"){
      ## Identify threshold of significance
      if(p.value=="max"){
        if(p.value.option == "count") thresh <- max(corr.sim)
        if(p.value.option == "density") thresh <- max(density(corr.sim)$x)
      }else{
        if(p.value.option == "count") thresh <- quantile(corr.sim, probs=1-p.value)
        if(p.value.option == "density") thresh <- quantile(density(corr.sim)$x, probs=1-p.value)
      }
      ## Identify (real) SNPs w correlations > thresh:
      sig.snps <- which(corr.dat > thresh)
      p.vals <- sapply(c(1:length(sig.snps)), 
                       function(e) length(which(corr.sim > 
                       corr.dat[sig.snps[e]]))/length(corr.sim))
      #p.vals <- sapply(c(1:length(sig.snps)), 
      #           function(e) 
      #           ecdf(corr.sim)(corr.dat[sig.snps[e]]))
    } # end test cor
    
    ## SCORE ##
    if(test=="score"){
      ## Identify threshold of significance
      if(p.value=="max"){
        if(p.value.option == "count") thresh <- max(corr.sim)
        if(p.value.option == "density") thresh <- max(density(corr.sim)$x)
      }else{
        if(p.value.option == "count") thresh <- quantile(corr.sim, probs=1-p.value)
        if(p.value.option == "density") thresh <- quantile(density(corr.sim)$x, probs=1-p.value)
      }
      ## Identify (real) SNPs w correlations > thresh:
      sig.snps <- which(corr.dat > thresh)
      p.vals <- sapply(c(1:length(sig.snps)), 
                       function(e) length(which(corr.sim > 
                       corr.dat[sig.snps[e]]))/length(corr.sim))
      #p.vals <- sapply(c(1:length(sig.snps)), 
      #           function(e) 
      #           ecdf(corr.sim)(corr.dat[sig.snps[e]]))
    } # end test score
    
    ## FISHER'S EXACT TEST ##
    if(test=="fisher"){  
      ## Identify threshold of significance
      if(p.value=="max"){
        
        ## NOTE: FOR FISHER:
        ## SHOULD "MAX" --> min(corr.sim)
        ## & SHOULD PROBS = P.VALUE (instead of 1 - p.value???)
        
        if(p.value.option == "count") thresh <- max(corr.sim)
        if(p.value.option == "density") thresh <- max(density(corr.sim)$x)
      }else{
        if(p.value.option == "count") thresh <- quantile(corr.sim, probs=1-p.value)
        if(p.value.option == "density") thresh <- quantile(density(corr.sim)$x, probs=1-p.value)
      }
      ## Identify (real) SNPs w correlations > thresh:
      sig.snps <- which(corr.dat < thresh)
      p.vals <- sapply(c(1:length(sig.snps)), 
                       function(e) length(which(corr.sim < 
                       corr.dat[sig.snps[e]]))/length(corr.sim))
      #p.vals <- sapply(c(1:length(sig.snps)), 
      #           function(e) 
      #           ecdf(corr.sim)(corr.dat[sig.snps[e]]))
    } # end test fisher
    
    
    
    ## get list of those correlation values
    sig.corrs <- corr.dat[sig.snps]
    ## get the list of those SNPs (ie. their locus names)
    sig.snps <- dimnames(snps)[[2]][sig.snps]
    ## re-order list of sig.snps and sig.corrs by value of sig.corr
    if(test=="fisher"){
      NWO <- order(sig.corrs, decreasing=TRUE)
    }else{
      NWO <- order(sig.corrs, decreasing=FALSE)
    }          
    sig.snps <- sig.snps[NWO]
    sig.corrs <- sig.corrs[NWO]
    p.vals <- p.vals[NWO]
    
    
    gc()
    
    ##
    
  } # end  method Didelot
  
  
  ########################################
  
  
  
  
  ## for Fitch, we want to calculate the site-specific 
  ## correlations between the SNP at each locus 
  ## and the trait, getting
  ## the distributions, across all reps, 
  ## of the correlations for each site. 
  
  if(method=="Fitch"){
    ## Calculate correlation between simulated SNPs and phenotype
    ## convert phenotype factor(y) to numeric:
    phen <- as.numeric(y)
    ## create empty list to store correlation results 
    ## for each of the list of simulated SNPs matrices
    corr.sim <- list()
    ## and an empty list, each element of which will contain 
    ## all the SNPs in column i from 
    ## each snp set made (ie from each rep)
    foo <- list()
    
    ## for loop to get correlations:
    for(i in 1:ncol(snps.mat[[1]])){ 
      ## all SNPs matrices wil have same ncol as 
      ## they are generated by permutation...
      ## for each site in the SNPs matrix, get the distribution 
      ## of correlations across all reps:
      ## make foo a set of all the SNP column i's from 
      ## each element of the list, snps.mat
      foo[[i]] <- sapply(c(1:length(snps.mat)), 
                         function(e) snps.mat[[e]][,i])
      ## get the distribution of the correlations between 
      ## those snp columns and phen:
      if(!is.null(foo[[i]])){
        # regular correlation..
        all.corr <- cor(foo[[i]], phen)  
        corr.sim[[i]] <- abs(all.corr)
      }else{
        corr.sim[[i]] <- NULL
      }
    } # end for loop
    
    
    ## Calculate correlations between real SNPs and phenotype
    # regular correlation...
    corr.dat <- sapply(c(1:ncol(snps)), 
                       function(e) cor(snps[,e], phen)) 
    ## USE ABSOLUTE VALUE??
    corr.dat <- abs(corr.dat)
    
    ## get empirical p-value based on corr.sim FOR EACH SNP LOCUS
    ## make a list to store the threshold sig cut-offs for 
    ## each SNP's correlation w the phenotype
    thresh <- list()
    corr.sim.old <- corr.sim
    corr.sim <- list()
    are.sig.snps <- list()
    sig.snps <- NA
    sig.corrs <- NA
    ## we only need to calculate for every other site:
    toKeep <- seq(1, ncol(snps), 2)
    for(i in toKeep){
      ## INTERPRETATION: you want the threshold above which 
      ## only 0.001 proportion (ie. 0.1%) of the data lies:
      round <- which(toKeep==i)
      corr.sim[[round]] <- as.vector(unlist(corr.sim.old[[i]]))
      thresh[[round]] <- quantile(corr.sim[[round]], probs=1-p.value)
    } # end for loop
    
    for(i in 1:length(thresh)){
      ## Identify (real) SNPs w correlations > thresh:
      are.sig.snps[[i]] <- corr.dat[i] > thresh[[i]]
      ## get list of those correlation values
      if(are.sig.snps[[i]]==TRUE){
        if(is.na(sig.snps[1])){
          sig.snps[1] <- 
            dimnames(snps)[[2]][toKeep[i]]
        }else{
          sig.snps[(length(sig.snps)+1)] <- 
            dimnames(snps)[[2]][toKeep[i]]
        }
        if(is.na(sig.corrs[1])){
          sig.corrs[1] <- corr.dat[i]
        }else{
          sig.corrs[(length(sig.corrs)+1)] <- 
            corr.dat[i]
        }  
      }
    } # end for loop
    
    ## re-order list of sig.snps and 
    ## sig.corrs by value of sig.corr  
    keepMe <- which(dimnames(snps)[[2]][toKeep]
                    %in% sig.snps)
    thresh <- as.vector(unlist(thresh[keepMe]))
    NWO <- order(sapply(c(1:length(sig.corrs)),
                        function(e) thresh[e]/sig.corrs[e]), 
                 decreasing=TRUE)
    sig.snps <- sig.snps[NWO]
    sig.corrs <- sig.corrs[NWO]
    thresh <- thresh[NWO]
    
    
    gc()
    
    ##
    
  } # end method Fitch
  
  
  
  
  
  
  ##################################
  ## 4) (A) Plot the distribution ##
  ##################################
  
  if(method=="Didelot"){
    
    
    ## plot.dist ###################
    
    
    #################################################################
    ## plot histogram of correlations btw real SNPs and phenotype: ##
    #################################################################
    
    #################
    ## CORRELATION ##
    #################
    if(test=="cor"){
      ## plot.dist ##
      
      #################################################################
      ## plot histogram of correlations btw real SNPs and phenotype: ##
      #################################################################
      ## get histogram of correlations btw real SNPs and phenotype:
      h <- hist(corr.dat, plot=FALSE)
      ## get X and Y coords for labelling positions of all sig SNPs
      X <- sig.corrs
      ## get the number of sig SNPs in each bin of the histogram
      ## counting correlations==upper limit of each bin as falling 
      ## within that bin...
      sig.counts <- sapply(c(1:(length(h$breaks)-1)), 
                           function(e) length(X[which(X[which(X <= 
                           h$breaks[(e+1)])] > h$breaks[e])]))
      ## keep only counts > 0
      sig.counts <- sig.counts[which(sig.counts > 0)]
      ## get label heights:
      Y <- list()
      ## get average height of labels
      Y.avg <- max(h$counts)/4
      ## for bins with > 1 sig SNP, adjust height  
      if(length(sig.counts)!=0){
        for(i in 1:length(sig.counts)){
          if(!.is.integer0(sig.counts[i])){
            if(sig.counts[i]==1){
              Y[[i]] <- Y.avg
            }else{
              ## divide up the space on the y-axis into 
              ## increments (adding 1 to the denomenator 
              ## s.t y-max not exceeded)
              increment <- max(h$counts)/(sig.counts[i]+1)
              Y[[i]] <- increment*c(1:sig.counts[i])
            }  
          }
        } # end for loop
      }
      Y <- as.vector(unlist(Y))
      
      if(plot.dist==TRUE){
        
        ## plot histogram of correlations btw real 
        ## SNPs and phenotype: ##
        plot(h, main="Distribution of SNP-trait 
             correlations \n (with significant SNPs indicated if present)")
        ## ADD threshold line in red on x-axis where thresh hits... 
        abline(v=thresh, col="grey", lwd=2, lty=2)
        ## label threshold line(?)
        text(x=thresh, y=(max(h$counts)), 
             labels="significance threshold", 
             col="grey", pos=2, font=4)
        
        ## Only ask to draw arrows if at least 1 significant SNP:
        if(length(sig.snps) > 0){
          ## ADD arrows pointing from each 
          ## label to position on X-axis:
          arrows(x0=X , y0=(Y-(max(h$counts)/50)) , 
                 x1=X , y1=0 , col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs > 
          ## threshold at their location on the x-axis:
          text(x=X, y=Y, labels=sig.snps, col="red", font=2)
        }else{
          text(x=thresh, y=(max(h$counts)/4)*3, 
               labels="no significant SNPs found", 
               col="red", font=2, pos=2)
        }
      } # end plot.dist
    } # end test cor
    
    
    ###########
    ## SCORE ##
    ###########
    if(test=="score"){
      ## plot.dist ##
      
      #################################################################
      ## plot histogram of correlations btw real SNPs and phenotype: ##
      #################################################################
      ## get histogram of correlations btw real SNPs and phenotype:
      h <- hist(corr.dat, plot=FALSE)
      ## get X and Y coords for labelling positions of all sig SNPs
      X <- sig.corrs
      ## get the number of sig SNPs in each bin of the histogram
      ## counting correlations==upper limit of each bin 
      ## as falling within that bin...
      sig.counts <- sapply(c(1:(length(h$breaks)-1)), 
                           function(e) length(X[which(X[which(X <= 
                           h$breaks[(e+1)])] > h$breaks[e])]))
      ## keep only counts > 0
      sig.counts <- sig.counts[which(sig.counts > 0)]
      ## get label heights:
      Y <- list()
      ## get average height of labels
      Y.avg <- max(h$counts)/4
      ## for bins with > 1 sig SNP, adjust height  
      if(length(sig.counts)!=0){
        for(i in 1:length(sig.counts)){
          if(!.is.integer0(sig.counts[i])){
            if(sig.counts[i]==1){
              Y[[i]] <- Y.avg
            }else{
              ## divide up the space on the y-axis into 
              ## increments (adding 1 to the denomenator 
              ## s.t y-max not exceeded)
              increment <- max(h$counts)/(sig.counts[i]+1)
              Y[[i]] <- increment*c(1:sig.counts[i])
            }  
          }
        } # end for loop
      }
      Y <- as.vector(unlist(Y))
      
      if(plot.dist==TRUE){
        
        ## plot histogram of correlations btw real SNPs and phenotype: ##
        plot(h, main="Distribution of SNP-trait correlations 
             \n (with significant SNPs indicated if present)")
        ## ADD threshold line in red on x-axis where thresh hits... 
        abline(v=thresh, col="grey", lwd=2, lty=2)
        ## label threshold line(?)
        text(x=thresh, y=(max(h$counts)), 
             labels="significance threshold", 
             col="grey", pos=2, font=4)
        
        ## Only ask to draw arrows if at least 1 significant SNP:
        if(length(sig.snps) > 0){
          ## ADD arrows pointing from each label to position on X-axis:
          arrows(x0=X , y0=(Y-(max(h$counts)/50)) , 
                 x1=X , y1=0 , 
                 col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs > 
          ## threshold at their location on the x-axis:
          text(x=X, y=Y, labels=sig.snps, 
               col="red", font=2, pos=2)
        }else{
          text(x=thresh, y=(max(h$counts)/4)*3, 
               labels="no significant SNPs found", 
               col="red", font=2, pos=2)
        }
      } # end plot.dist
    } # end test score
    
    
    #########################
    ## FISHER'S EXACT TEST ##
    #########################
    if(test=="fisher"){  
      ## plot.dist ##
      
      ############################################################
      ## plot histogram of correlations btw real SNPs and phen: ##
      ############################################################
      
      ## get histogram of correlations btw real SNPs and phenotype:
      h <- hist(corr.dat, plot=FALSE)
      ## get X and Y coords for labelling positions of all sig SNPs
      X <- thresh
      if(length(sig.corrs) > 0) X <- sig.corrs
      ## get the number of sig SNPs in each bin of the histogram
      ## counting correlations==upper limit of each bin as 
      ## falling within that bin...
      sig.counts <- sapply(c(1:(length(h$breaks)-1)), 
                           function(e) length(X[which(X[which(X <= 
                           h$breaks[(e+1)])] > h$breaks[e])]))
      ## keep only counts > 0
      sig.counts <- sig.counts[which(sig.counts > 0)]
      ## get label heights:
      Y <- list()
      ## get average height of labels
      Y.avg <- max(h$counts)/4
      ## for bins with > 1 sig SNP, adjust height  
      if(length(sig.counts)!=0){
        for(i in 1:length(sig.counts)){
          if(!.is.integer0(sig.counts[i])){
            if(sig.counts[i]==1){
              Y[[i]] <- Y.avg
            }else{
              ## divide up the space on the y-axis into 
              ## increments (adding 1 to the denomenator s.t y-max not exceeded)
              increment <- max(h$counts)/(sig.counts[i]+1)
              Y[[i]] <- increment*c(1:sig.counts[i])
            }  
          }
        } # end for loop
      }
      Y <- as.vector(unlist(Y))
      
      if(plot.dist==TRUE){
        
        ## plot histogram of correlations btw real SNPs and phenotype: ##
        plot(h, main="Distribution of SNP-trait correlations 
             \n (with significant SNPs indicated if present)")
        ## ADD threshold line in red on x-axis where thresh hits... 
        abline(v=thresh, col="grey", lwd=2, lty=2)
        ## label threshold line(?)
        text(x=thresh, y=(max(h$counts)), 
             labels="significance threshold", 
             col="grey", pos=4, font=4)
        
        ## Only ask to draw arrows if at least 1 significant SNP:
        if(length(sig.snps) > 0){
          ## ADD arrows pointing from each label to 
          ## position on X-axis:
          arrows(x0=X , y0=(Y-(max(h$counts)/50)) , 
                 x1=X , y1=0 , col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs > threshold at their location on the x-axis:
          text(x=X, y=Y, labels=sig.snps, 
               col="red", font=2, pos=4)
        }else{
          text(x=thresh, y=(max(h$counts)/4)*3, 
               labels="no significant SNPs found", 
               col="red", font=2, pos=4)
        }
      } # end plot.dist
    } # end test fisher
    
    
    
    
    
    ## plot.null.dist #####################################
    
    #######################################################
    ## OR--add SNP annotations to histogram of null dist ##
    #######################################################
    
    #################
    ## CORRELATION ##
    #################
    if(test=="cor"){
      ## plot.null.dist ##
      
      ###########################################################
      ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
      ###########################################################
      ## plot correlations btw simulated SNPs and phenotype:
      h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)
      
      ## get alternate (null dist) label heights:
      X <- thresh
      if(length(sig.corrs) > 0) X <- sig.corrs
      Y.null <- runif(n=length(sig.snps), min=0.00001, 
                      max=max(h.null$counts))
      
      
      if(plot.null.dist==TRUE){
        
        ## if the true correlation value for SNP i is > 
        ## max bin, then extend the x-axis of the plot 
        ## to accommodate annotation:
        if(max(X) > max(h.null$breaks)){
          ## plot histogram of correlations btw real 
          ## SNPs and phenotype: ##
          ## EXTENDING THE X-AXIS
          plot(h.null, 
               main="Null distribution of correlations 
               \n (with significant SNPs indicated)",
               xlab="Correlation",
               xlim=c(min(h.null$breaks), max(X)+.05))
        }else{
          ## plot histogram of correlations btw real 
          ## SNPs and phenotype: ##
          ## WITHOUT EXTENDING THE X-AXIS
          plot(h.null, 
               main="Null distribution of correlations 
               \n (with significant SNPs indicated)",
               xlab="Correlation"
          )           
        }
        
        ## ADD threshold line in red on x-axis where thresh hits... 
        abline(v=thresh, col="grey", lwd=2, lty=2)
        ## label threshold line(?)
        text(x=thresh, y=(max(h.null$counts)), 
             labels="significance threshold", pos=2, 
             col="grey", font=4)
        if(length(sig.snps) > 0){
          ## ADD arrows pointing from each label to 
          ## position on X-axis:
          arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) , 
                 x1=X , y1=0 , col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs > 
          ## threshold at their location on the x-axis:
          text(x=X, y=Y.null, labels=sig.snps, 
               col="red", font=2)
        }else{
          text(x=thresh, y=(max(h.null$counts)/4)*3, 
               labels="no significant SNPs found", 
               col="red", font=2, pos=2)
        }
      } # end plot.null.dist
      } # end test cor
    
    
    ###########
    ## SCORE ##
    ###########
    if(test=="score"){
      ## plot.null.dist ##
      
      ###########################################################
      ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
      ###########################################################
      ## plot correlations btw simulated SNPs and phenotype:
      h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)
      
      ## get alternate (null dist) label heights:
      X <- thresh
      if(length(sig.corrs) > 0) X <- sig.corrs
      Y.null <- runif(n=length(sig.snps), min=0.00001, 
                      max=max(h.null$counts))
      
      
      if(plot.null.dist==TRUE){    
        ## if the true correlation value for SNP i is > 
        ## max bin, then extend the x-axis of the plot 
        ## to accommodate annotation:
        if(max(X) > max(h.null$breaks)){
          ## plot histogram of correlations btw real 
          ## SNPs and phenotype: ##
          ## EXTENDING THE X-AXIS
          plot(h.null, 
               main="Null distribution of correlation scores 
               \n (with significant SNPs indicated)",
               xlab="(Correlation) Score",
               xlim=c(min(h.null$breaks), max(X)+.05))
        }else{
          ## plot histogram of correlations btw real 
          ## SNPs and phenotype: ##
          ## WITHOUT EXTENDING THE X-AXIS
          plot(h.null, 
               main="Null distribution of correlation scores 
               \n (with significant SNPs indicated)",
               xlab="(Correlation) Score"
          )           
        }
        
        ## ADD threshold line in red on x-axis where thresh hits... 
        abline(v=thresh, col="grey", lwd=2, lty=2)
        ## label threshold line(?)
        text(x=thresh, y=(max(h.null$counts)), 
             labels="significance threshold", pos=2, 
             col="grey", font=4)
        if(length(sig.snps) > 0){
          ## ADD arrows pointing from each label to 
          ## position on X-axis:
          arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) , 
                 x1=X , y1=0 , col="blue", 
                 length=0.1, lwd=1)
          ## add annotation text labelling SNPs > 
          ## threshold at their location on the x-axis:
          text(x=X, y=Y.null, labels=sig.snps, 
               col="red", font=2, pos=2)
        }else{
          text(x=thresh, y=(max(h.null$counts)/4)*3, 
               labels="no significant SNPs found", 
               col="red", font=2, pos=2)
        }
      } # end plot.null.dist
      } # end test score
    
    
    #########################
    ## FISHER'S EXACT TEST ##
    #########################
    if(test=="fisher"){ 
      ## plot.null.dist ##
      
      ###########################################################
      ## plot (null) hist of corr's btw SIMULATED SNPs n phen: ##
      ###########################################################
      ## plot correlations btw simulated SNPs and phenotype:
      h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)
      
      ## get alternate (null dist) label heights:
      X <- thresh
      if(length(sig.corrs) > 0) X <- sig.corrs
      Y.null <- runif(n=length(sig.snps), min=0.00001, 
                      max=max(h.null$counts))
      
      if(plot.null.dist==TRUE){    
        ## if the true correlation value for SNP i is < 
        ## min bin, then extend the x-axis of the plot to 
        ## accommodate annotation:
        if(min(X) < min(h.null$breaks)){
          ## plot histogram of correlations btw real 
          ## SNPs and phenotype: ##
          ## EXTENDING THE X-AXIS
          plot(h.null, 
               main="Null distribution of Fisher's exact test p-values 
               \n (with significant SNPs indicated)",
               xlab="p-value",
               xlim=c((min(X)-.05), max(h.null$breaks)))
        }else{
          ## plot histogram of correlations btw real SNPs and phenotype: ##
          ## WITHOUT EXTENDING THE X-AXIS
          plot(h.null, 
               main="Null distribution of Fisher's exact test p-values
               \n (with significant SNPs indicated)",
               xlab="p-value"
          )           
        }
        
        ## ADD threshold line in red on x-axis where thresh hits... 
        abline(v=thresh, col="grey", lwd=2, lty=2)
        ## label threshold line(?)
        text(x=thresh, y=(max(h.null$counts)), 
             labels="significance threshold", pos=4, 
             col="grey", font=4)
        ## only ask to draw arrows if sig snps exist
        if(length(sig.snps) > 0){
          ## ADD arrows pointing from each label to 
          ## position on X-axis:
          arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) , 
                 x1=X , y1=0 , col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs >
          ## threshold at their location on the x-axis:
          text(x=X, y=Y.null, labels=sig.snps, 
               col="red", font=2, pos=4)
        }else{
          text(x=thresh, y=(max(h.null$counts)/4)*3, 
               labels="no significant SNPs found", 
               col="red", font=2, pos=4)
        }
      } # end plot.null.dist
      } # end test fisher
    
    
    
    
    
  } # end Didelot method
  
  
  ########################################################
  #### B) Plot the distributions* ########################
  ## (ONLY for the SIGNIFICANT SNPs) of the  #############
  ## correlations btw the phen and the randomized SNPs ###
  ########## *as these will be specific to each site...? #
  ########################################################
  
  if(method=="Fitch"){
    
    ## plot Fitch ##
    
    ## Only call to plot if there are some significant SNPs:
    if(length(sig.snps) > 0){
      
      ## we only need to calculate for every other site:
      toKeep <- seq(1, ncol(snps), 2)
      ## get indices of sig snps
      snpsToKeep <- which(dimnames(snps)[[2]][toKeep] 
                          %in% sig.snps)
      
      ## for each of the snps in snpsToKeep, 
      ## plot the null distribution of the 
      ## correlations btw that SNP and the phen
      if(length(snpsToKeep) > 1) par(ask=TRUE)
      for(i in snpsToKeep[NWO]){
        round <- which(snpsToKeep[NWO]==i)
        ####################################################
        ## plot hist of corr's btw simulated SNPs n phen: ##
        ####################################################
        ## get histogram of correlations btw simulated 
        ## SNPs and phenotype:
        h <- hist(corr.sim[[i]], plot=FALSE)
        ## get X and Y coords for labelling position of sig SNP
        ## get position on the x-axis of the labels
        X <- sig.corrs[round]
        ## get height of labels
        Y <- max(h$counts)/4
        
        if(plot.Fitch==TRUE){
          
          ## if the true correlation value for SNP i is > 
          ## max bin, then extend the x-axis of the 
          ## plot to accommodate annotation:
          if(X > max(h$breaks)){
            ## plot histogram of correlations btw real 
            ## SNPs and phenotype: ##
            ## EXTENDING THE X-AXIS
            plot(h, 
                 main=paste("Distribution of SNP-trait 
                            ## correlations for SNP", 
                            dimnames(snps)[[2]][toKeep][i], 
                            "\n (with true correlation indicated)"),
                 xlab="Correlation",
                 xlim=c(min(h$breaks), X+.05))
          }else{
            ## plot histogram of correlations btw real 
            ## SNPs and phenotype: ##
            ## WITHOUT EXTENDING THE X-AXIS
            plot(h, 
                 main=paste("Distribution of SNP-trait 
                            correlations for SNP", 
                            dimnames(snps)[[2]][toKeep][i], 
                            "\n (with true correlation indicated)"))           
          }
          
          ## ADD threshold line in red on x-axis where thresh hits... 
          abline(v=thresh[[round]], col="grey", lwd=2, lty=2)
          ## label threshold line(?)
          text(x=thresh[[round]], y=(max(h$counts)), 
               labels="significance threshold", 
               col="grey", pos=2, font=4)
          ## ADD arrows pointing from each label to position on X-axis:
          arrows(x0=X , y0=Y , x1=X , y1=0 , 
                 col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs > threshold at their location on the x-axis:
          text(x=X, y=Y, labels=sig.snps[round], 
               col="red", font=2)
        } # end plot.Fitch
      } # end for loop
      }
    
  } # end Fitch/ modified-Farhat method
  
  par(ask=FALSE)
  
  ########################################
  ## 5) Return results list ##############
  ########################################
  
  if(length(sig.snps)==0) sig.snps <- sig.corrs <- NULL
  
  ###########
  ## make a data.frame containing all relevant output for sig.snps
  if(length(sig.snps) > 0){
    
    ## Get counts for n.sig.snps in each cell of the contingency table: 
    toKeep <- sapply(c(1:length(sig.snps)), 
                     function(e) 
                       which(dimnames(snps)[[2]] == sig.snps[e]))
    snps.toKeep <- snps[,toKeep]
    
    ## 
    S1P1 <- sapply(c(1:ncol(snps.toKeep)), 
                   function(e) 
                     length(which(snps.toKeep[which(phen==1),e]==1)))
    S0P0 <- sapply(c(1:ncol(snps.toKeep)), 
                   function(e) 
                     length(which(snps.toKeep[which(phen==0),e]==0)))
    S1P0 <- sapply(c(1:ncol(snps.toKeep)), 
                   function(e) 
                     length(which(snps.toKeep[which(phen==0),e]==1)))
    S0P1 <- sapply(c(1:ncol(snps.toKeep)), 
                   function(e) 
                     length(which(snps.toKeep[which(phen==1),e]==0)))
    
    df <- data.frame(sig.snps, 
                     p.vals, 
                     sig.corrs, 
                     S1P1, S0P0, S1P0, S0P1)
    names(df) <- c("SNP.locus", 
                   "p.value", 
                   "Test.statistic", 
                   "S1P1", "S0P0", "S1P0", "S0P1")
  }else{
    df <- "No significant SNPs found."
  }
  
  ## 0 p.vals 
  min.p <- paste("p-values listed as 0 are <",
                 1/length(corr.sim), sep=" ")
  
  ## TO DO: 
  ## ADD MANHATTAN PLOT
  
  if(method=="Didelot"){
    results <- list()
    results[[1]] <- corr.dat
    results[[2]] <- corr.sim
    results[[3]] <- thresh
    results[[4]] <- df
    results[[5]] <- min.p
    
    
    names(results) <- c("Correlation values for 
                        empirical SNP-trait associations",
                        "Correlation values for 
                        simulated SNP-trait associations",
                        "Significance threshold",
                        "Significant SNPs",
                        "Note on p-values"
    )
    
    
  } # end method Didelot
  ##
  
  
  
  if(method=="Fitch"){
    results <- list()
    results[[1]] <- sig.snps
    results[[2]] <- sig.corrs
    results[[3]] <- thresh
    names(results) <- c("SNPs significantly associated with 
                        trait (lowest-highest excess of
                        corr' > thresh)",
                        "Correlations between (significant) 
                        SNPs and phenotype/trait",
                        "Significance threshold")
    
  } # end method Fitch
  
  
  return(results)
  
  
} # end treeWAS

##

###################################################################
###################################################################
###################################################################

##################
## .is.integer0 ##
##################

## fn testing for output==integer(0)
.is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
} # end .is.integer0




##