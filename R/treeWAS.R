
#############
## treeWAS ##
#############


#'treeWAS
#': tree-based GWAS.
#'
#'This function runs genome-wide association studies (GWAS) for bacteria by comparing SNP-trait correlations in the data 
#'with those derived from simulated genomic datasets based on phylogenetic trees taken from the data. 
#' 
#' @param x A list, matrix, DNAbin, or genind object containing whole-genome or SNP datasets with individuals in rows and genetic variables in columns.
#' @param y A vector or factor containing the phenotype or trait status of all individuals in the dataset \code{x}.
#' @param snps.gen An optional genind object containing SNPs (to save time required for conversion).
#' @param tree An optional phylo or hclust object containing a phylogenetic tree to be used to 
#' represent the real data and as a basis for simulating the comparator datasets.
#' @param mt.rate An optional integer specifying the mutation rate to be used when simulating comparator datasets.
#'  If \code{mt.rate = NULL}, when \code{sim.by = "locus"}, \code{mt.rate} will be set to 1*2; 
#'  and when \code{sim.by = "branch"}, \code{mt.rate} will be estimated based on the data using the function parsimony from package phangorn.
#' @param p.value A number specifying the p-value to be used when setting the threshold for determining which SNP-trait correlations
#'  are deemed to be statistically significant.
#' @param plot.null.dist A logical indicating whether to plot the significant SNPs identified (if any) on top of the null distribution
#'  of correlations between the simulated SNPs and the phenotype, \code{y}. (TRUE by default.)
#' @param plot.dist A logical indicating whether to plot the significant SNPs identified (if any) on top of the distribution of correlations
#'  between the SNPs provided in \code{x} and the phenotype, \code{y}. (FALSE by default.)
#' @param plot.Fitch A logical indicating whether to produce a plot, if \code{method = "Fitch"}. (TRUE by default.)
#' @param plot.tree A logical indicating whether to produce a plot of the phylogenetic tree. (FALSE by default.)
#' @param tree.method The (optional) method used to reconstruct a phylogenetic tree if \code{tree = NULL}.
#'  Must be one of c("ml", "nj", "UPGMA"), specifying maximum-likelihood, neighbour-joining, and UPGMA (hclust, method = "average"), respectively.
#' @param n.reps An integer specifying the number of repititions desired when generating the comparator simulated datasets (ie. the number
#'  of simulated datasets desired). (\code{n.reps = 1} by default.)
#' @param sim.gen.size An optional integer specifying the desired genome size for the simulated SNP datasets, if different from the number of
#'  SNPs in the real data provided in \code{x}.(\code{sim.gen.size = 10,000} by default.)
#' @param sim.by The method by which the simulated data is generated (when \code{method = "Didelot"}).
#' \itemize{
#'  \item If \code{sim.by = "locus"}, then \code{mt.rate} denotes the expected number of mutations to occur at each locus.
#'  \item If \code{sim.by = "branch"}, then \code{mt.rate} denotes the expected number of mutations to occur in the whole genome.
#'  }
#' @param method The method by which SNP simulation and association testing proceeds.
#' \itemize{
#'  \item If \code{method = "Didelot"}, SNPs are simulated either by locus or by branch (see parameter \code{sim.by}), and association testing
#'  occurs by comparing the distribution of correlations between real SNPs, \code{x}, and the phenotype, \code{y}, across all loci.
#'  \item If \code{method = "Fitch"}, SNPs are simulated according to a site-specific mutation rate, determined by Fitch's parsimony, 
#'  and association testing occurs site-by-site, comparing, at each locus, the correlation between the real SNP and the phenotype and the
#'  correlation between the simulated SNP and the phenotype.
#'  }
#' @param test The association test statistic to be used.
#' \itemize{
#'  \item If \code{test = "cor"}, the standard correlation between SNP and phenotype will be used. ("cor" is not recommended, as it is intended
#'  for use with continuous variables.)
#'  \item If \code{test = "score"}, a correlation score (currently for binary SNPs with binary phenotypes only)
#'  is calculated based on the SNP-phenotype contingency table at each locus, using the folowing
#'  formula to calculate the score for SNPi: \cr
#'  Scorei = ((S1P1 + S0P0) - (S1P0 + S0P1)) / N , \cr 
#'  where S1 denotes SNPi=1, S0 denotes SNPi=0, P1 denotes Phenotype=1, P0 denotes Phenotype=0, and N denotes the number of individuals 
#'  (ie. the total number of rows in the dataset \code{x}). 
#'  \item If \code{test = "fisher"}, the p-values returned by a Fisher's exact test between each SNP and the phenotype are used. 
#'  }
#' 
#' @author Caitlin Collins \email{caitlin.collins12@@imperial.ac.uk}
#' @export
#' @examples
#' \dontrun{
#' #############
#' ## EXAMPLE ##
#' #############
#' 
#' ## 1) Simulate data with your coalescent fn (see coalescent.R for commented fn BUT load coalescent.sim.R for actual fn)
#' ## 2) Identify SNPs* to test--no sense in testing each site in the genome if only some vary enough... 
#' #### *determine threshold to define a "SNP"... (eg. >= 1% of sites have alternate nt)
#' 
#' # simulate data
#' out <- coalescent.sim(n.ind=10, gen.size=10000, theta=100*2, biallelic=TRUE, seed=2, plot=TRUE, heatmap=FALSE, plot2="ml")
#' ## isolate elements of out-->
#' ## get list of genomes
#' x <- snps <- out[[1]]
#' ## get tree
#' tree <- out[[2]]
#' ## take a look at the phylogeny:
#' plot(tree)
#' 
#' ## Map phenotype perfectly to SNP 1
#' ## convert genomes to get SNPs as: (1) data.frame, (2) genind object, (3) SNPs matrix:
#' if(is.null(names(x))) names(x) <- paste("ind", c(1:length(x)), sep=".")
#' gen.size <- length(x[[1]])
#' ## working with snps in matrix form
#' snps <- do.call("rbind", x)
#' ## get snps as DNAbin 
#' snps.bin <- as.DNAbin(snps)
#' snps.gen <- DNAbin2genind(snps.bin, polyThres=0.01)
#' snps <- snps.gen@@tab
#' 
#' R <- which(snps[,1]==1)
#' y <- c(rep("S", nrow(snps)))
#' y <- as.factor(replace(y, R, "R"))
#' ## snps (list)
#' x <- snps <- out[[1]]
#' 
#' result <- treeWAS(x, y, snps.gen=snps.gen, tree=tree, mt.rate=NULL, p.value=0.001, 
#'                   plot.null.dist=TRUE, plot.dist=FALSE, plot.Fitch=TRUE, plot.tree=FALSE, tree.method=NULL, 
#'                   n.reps=1, sim.gen.size=1000,
#'                   sim.by="locus", method="Didelot", test="score")
#' 
#' ## inspect results
#' str(result)
#' 
#' ## inspect the data.frame containing information about the SNPs identified as "significant"
#' result[[5]]
#' 
#' ## NOTE: the name of SNP i is not L.i, so this is the best way to use the index of a SNP
#' ## to test whether that SNP has been identified as significant in your results
#' dimnames(x)[[2]][1] %in% result[[5]]$SNP.locus
#' }
#' @import adegenet phangorn ape ade4
#'



######################################################################################################################################
######################################################################################################################################
######################################################################################################################################



#############
## treeWAS ##
#############


## Re-implementation of GWAS method developed in Sheppard et al 2014 (Xavier's Campylobacter GWAS paper),
## based on the generation of a phylogenetically-correct p-value derived by comparing correlation btw SNPs
## and a phenotype of interest to a null distribution (Poisson dist w parameter 1) OR 
## (modification based on similar method developed in Farhat et al 2013) based on permutation of empirically-derived
## n.mts per site (except we will determine this with the Fitch algorithm...). 

########################################################################

## tree.method = c("nj", "UPGMA", "ml") +/ "parsimony"??, "Bayesian"??
## method = c("Didelot", "Fitch")

## mt.rate <- 1(*2) # IF sim.by=="locus" # ie. poisson mean of 1 for each site (simulating SNP by SNP) ## OR # gen.size*2 # IF sim.by=="branch" #
## sim.gen.size <- 10000 by default (if NULL --> use actual genome size of data supplied (ie. ncol(x)))
## test <- c("cor", "score", "fisher")

########################################################################


treeWAS <- function(x, y, snps.gen=NULL, tree=NULL, mt.rate=NULL, p.value=0.001, 
                    plot.null.dist=TRUE, plot.dist=FALSE, plot.Fitch=TRUE, plot.tree=FALSE, tree.method="ml", 
                    n.reps=1, sim.gen.size=10000,
                    sim.by="locus", method="Didelot", test="score"){  
  
  ## load packages:
  require(adegenet)
  require(phangorn)
  require(adephylo)
  require(ade4) #?
  
  ## load my simulation fn:
  source("tree.sim.R")
  ## load fitch simulation/randomization fn:
  source("fitch.sim.R")
  
  #   ## load my simulation fn:
  #   source("C:/Cait 2012/Work/Xavier/treeWAS Method/tree.sim.R")
  #   ## load fitch simulation/randomization fn:
  #   source("C:/Cait 2012/Work/Xavier/treeWAS Method/fitch.sim.R")  
  
  ##################################################################################################################################
  ## 0) Handle input data ##########################################################################################################
  ##################################################################################################################################
  
  ## convert genomes to get SNPs as: (1) data.frame, (2) genind object, (3) SNPs matrix:
  if(class(x)=="list"){
    if(is.null(names(x))) names(x) <- paste("ind", c(1:length(x)), sep=".")
    gen.size <- length(x[[1]])
    ## working with snps in matrix form
    snps <- do.call("rbind", x)
    ## get snps as DNAbin 
    snps.bin <- as.DNAbin(snps)
    if(is.null(snps.gen)){
      snps.gen <- DNAbin2genind(snps.bin, polyThres=0.01)
      snps <- snps.gen@tab
    }else{
      snps <- snps.gen@tab
    }
  }
  
  if(class(x)=="matrix"){
    if(is.null(row.names(x))) row.names(x) <- paste("ind", c(1:length(x)), sep=".")
    gen.size <- ncol(x)
    ## for as.DNAbin, snps must be a character (not numeric) matrix... 
    if(class(try(as.DNAbin(snps), silent=TRUE))=="try-error"){
      snps.bin <- as.DNAbin(matrix(as.character(snps), ncol=ncol(snps), byrow=FALSE))
    }else{
      snps.bin <- as.DNAbin(snps)
    }
    if(is.null(snps.gen)){
      snps.gen <- DNAbin2genind(snps.bin, polyThres=0.01)
      snps <- snps.gen@tab
    }else{
      snps <- snps.gen@tab
    }
  }
  
  if(class(x)=="DNAbin"){
    if(is.null(row.names(x))) row.names(x) <- paste("ind", c(1:length(x)), sep=".")
    gen.size <- ncol(x)
    snps.bin <- x
    if(is.null(snps.gen)){
      snps.gen <- DNAbin2genind(snps.bin, polyThres=0.01)
      snps <- snps.gen@tab
    }else{
      snps <- snps.gen@tab
    }
  }
  
  if(class(x)=="genind"){
    if(is.null(row.names(x@tab))) row.names(x@tab) <- paste("ind", c(1:length(x)), sep=".")
    gen.size <- ncol(x@tab) ## CAREFUL--this could be the number of SNPs instead of the GENOME size!
    dat.gen <- x
    ## need snps matrix to be character (nts not binary) for as.DNAbin to be happy
    snps.bin <- as.DNAbin(matrix(as.character(dat.gen@tab), ncol=ncol(dat.gen@tab), byrow=FALSE))
    ## seg.sites ## --> ALTERNATIVE TO USING DNAbin2genind w polyThresh!!
    ## convert from genomes to SNPs (unless separate snps.gen also provided)
    if(is.null(snps.gen)){
      snps.gen <- DNAbin2genind(snps.bin, polyThres=0.01)
      snps <- snps.gen@tab
    }else{
      snps <- snps.gen@tab
    }
  }
  
  ## NOTE: Do we want to use all genomic data for (optional) tree reconstruction
  ## OR just use SNPs matrix (if latter, ensure snps.bin contains SNPs (not genomes)):
  #if(length(snps.bin[[1]]) != ncol(snps)) snps.bin <- as.DNAbin(snps) ## TODO: ADD CHECK FOR CHARACTER MATRIX AND CONVERT IF NEEDED
  
  
  
  ## convert phenotype to factor
  y <- as.factor(y)
  
  ## set n.ind:
  n.ind <- length(y)
  inds <- c(1:n.ind)
  
  if(is.null(tree)){
    ## reconstruct phylogenetic tree if necessary:
    
    ## if no tree.method selected, use maximum-likelihood as default:
    if(is.null(tree.method)) tree.method <- "ml"
    
    ## working with data in DNAbin format...
    ## get distance matrix for terminal node individuals
    D <- dist.dna(snps.bin[1:n.ind,], model="JC69")
    
    ## MAXIMUM-LIKELIHOOD TREE
    if(tree.method=="ml"){
      dna4 <- as.phyDat(snps.bin[1:n.ind,])
      tre.ini <- nj(D)
      fit.ini <- pml(tre.ini, dna4, k=n.ind)
      fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, 
                       optQ = TRUE, optGamma = TRUE)
      ## STORE RESULTS for maximum-likelihood tree reconstruction:
      ml.results <- list(anova(fit.ini, fit), AIC(fit.ini), AIC(fit))
      names(ml.results) <- c("ANOVA for initial vs. optimal fit", "AIC for initial fit", "AIC for optimal fit")
      
      tree <- fit$tree
      tree <- midpoint(ladderize(tree))
      
      ## plot tree
      if(plot.tree==TRUE){
        plot(tree, show.tip=TRUE, edge.width=2)        
        title("Maximum-likelihood tree")
        axisPhylo()
      } # end plot.tree ml
    } # end max-likelihood tree
    
    
    ## NEIGHBOUR-JOINING TREE
    if(tree.method=="nj"){
      tree <- nj(D)
      tree <- midpoint(ladderize(tree))
      ## plot tree
      if(plot.tree==TRUE){
        plot(tree, edge.width=2)
        title("Neighbour-joining tree")
        axisPhylo()
      } # end plot.tree nj
    } # end neighbour-joining tree
    
    
    ## UPGMA TREE
    if(tree.method=="UPGMA"){
      tree <- hclust(D, method="average") 
      tree <- as.phylo(tree)
      tree <- midpoint(ladderize(tree))
      ## plot tree
      if(plot.tree==TRUE){
        plot(tree, main="")
        title("UPGMA tree")
      } # end plot.tree UPGMA
    } # end UPGMA tree    
    
    # end (optional) tree reconstruction (and (optional) plotting of reconstructed tree)
  }else{
    
    ## If user has already submitted a tree as input:
    
    ## Work with a centered phylo tree for consistency and visualisation's sake:
    if(class(tree) != "phylo") tree <- as.phylo(tree)
    #     ## if you can ladderize the tree without causing problems, do so:
    #     ## NOTE: for some reason plot(x, plot=FALSE) causes the current plot to be replaced with an empty plot window... 
    #     if(class(try(p <- plot(ladderize(tree), plot=FALSE), silent=TRUE))!="try-error") tree <- ladderize(tree)
    ## if the tree is not already rooted, root it:
    if(is.rooted(tree)==FALSE) tree <- midpoint(tree)

    
    if(plot.tree==TRUE){
      plot(tree)
      title("Phylogenetic tree (original)")
      axisPhylo()
    } # end plot.tree
  } # end tree...
  
  ## tree's tip.labels must be numeric (for Fitch parsimony step)
  ## NOTE: maybe just do this at that step and then reverse it (keeping user's original labelling)??
  tree.ori <- tree
  tree$tip.label <- c(1:length(tree$tip.label))
  tree$node.label <- c((n.ind+1):(n.ind+tree$Nnode))
  
  
  ## MT.RATE ##
    
  ## set default mt rate
  if(missing(mt.rate)) mt.rate <- NULL
  
  if(sim.by=="branch"){
    ## USE FITCH.N.MTS AS THETA (MT RATE) IN treeWAS...
    if(is.null(mt.rate)){ # mt.rate <- 1  
      ## TODO: ADD CHECK FOR BIALLELIC SNP, CHANGE LEVELS IF NOT BIALLELIC!!
      snps.phyDat <- as.phyDat(as.matrix(snps), type="USER", levels=c(0,1))
      mt.rate <- parsimony(tree, snps.phyDat, method="fitch") ## --> single-number parsimony score (for indian crypto = 1658889)
    }
  } # end sim.by branch
  
  if(sim.by=="locus"){
    if(is.null(mt.rate)) mt.rate <- 1
  } # end sim.by locus
  
  
  
  
  ##################################################################################################################################
  ## 1) Simulate multiple trees to compare your real/ reconstructed original tree to ###############################################
  ##################################################################################################################################
  
  ##################################################################################################################################
  #### A) Sheppard Method: #########################################################################################################
  ###### - Modify the coalescent method you previously developed
  ###### - Generate N comparator datasets using Mt.rate = theta/2 = 1/(sum of branch lengths) --> n.mts = Poisson dist'd w parameter 1
  ##################################################################################################################################
  
  if(method=="Didelot"){
    
    out <- list()
    genomes <- list()
    trees <- list()
    snps.mat <- list()
    
    for(i in 1:n.reps){
      if(is.null(sim.gen.size)) sim.gen.size <- gen.size
      if(sim.by=="branch"){
        if(sim.gen.size != gen.size) mt.rate <- mt.rate*(sim.gen.size/gen.size)
      }
      
      ## SIMULATE A DATASET | your tree ##
      out[[i]] <- tree.sim(n.ind=n.ind, tree=tree, gen.size=sim.gen.size, sim.by=sim.by, theta=mt.rate*2, 
                           biallelic=TRUE, seed=NULL, plot=FALSE, heatmap=FALSE, plot2=FALSE)
      
      ## NOTE: Do we want to set seed=i (for reproducibility) or not (for randomness' sake)? (OR Make this an option for user to decide?)
      
      genomes[[i]] <- out[[i]][[1]]
      trees[[i]] <- out[[i]][[2]]
      
      ## Convert genomes into snps matrices
      dat <- genomes[[i]]
      dat <- do.call("rbind", dat)
      dat.bin <- as.DNAbin(dat)
      dat.gen <- DNAbin2genind(dat.bin, polyThres=0.01)
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
    
  } # end Didelot method
  
  ## check avg n.mts
  #mean(sapply(c(1:length(snps.mat)), function(e) ncol(snps.mat[[e]]))) # /2 # ?
  
  ################################################################################################################################
  #### B) Farhat Method: #########################################################################################################
  ###### - Modify the method implemented in Farhat et al 2013
  ###### - Determine n.mts per site* with the Fitch algorithm
  ############ *(ie. The site-specific n.mts on the tree required to explain your SNP without knowing where on the tree they occurred)
  ###### - Assign this n.mts to new branches at random (with Prob corr. branch length) (repeat Nx to get N trees)
  ##################################################################################################################################
  
  ## TBD: 
  #### - Implement Fitch method of determining n.mts per site
  #### - Implement randomization procedure (w fixed n.mts/site, randomization ~ branch length)
  #### - Figure out how to create trees with this data
  #### - Save data in 2 lists: genomes and trees
  
  if(method=="Fitch"){
    
    #########################################################
    ## get n.mts per site with fitch's parsimony algorith: ##
    #########################################################
    snps.phyDat <- list()
    fitch.n.mts <- list()
    
    ## NOTE--snps currently contains a matrix with 2 columns for every SNP site. 
    #### With biallelic SNPs, we only need to retain the first of these pairs of columns, ie.:    
    toKeep <- seq(1, ncol(snps), 2)
    
    if(length(toKeep)!=0){
      for(i in 1:length(toKeep)){      
        snps.phyDat[[i]] <- as.phyDat(as.matrix(snps[,toKeep[i]]), type="USER", levels=c(0,1))
        fitch.n.mts[[i]] <- parsimony(tree, snps.phyDat[[i]], method="fitch")
        ###########################################################
        ## length(toKeep) 113382 --> system.time = 15 hours! :(  ##
        ###########################################################
        gc()
        
      } # end for loop
    }
    
    fitch.n.mts <- as.vector(unlist(fitch.n.mts))
    
    ## inspect distribution of n.mts:
    #table(fitch.n.mts)
    
    
    #######################################################################################
    ## Randomize this n.mts per site onto different branches with prob ~ branch lenghts: ##
    #######################################################################################
    
    out <- list()
    genomes <- list()
    trees <- list()
    snps.mat <- list()
    
    for(i in 1:n.reps){
      if(is.null(sim.gen.size)) sim.gen.size <- gen.size
      out[[i]] <- fitch.sim(snps, fitch.n.mts, tree, gen.size=sim.gen.size, biallelic=TRUE, seed=NULL)
      
      ## do we want to set seed=i (for reproducibility) or not (for randomness' sake)?
      
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
  
  
  ##################################################################################################################################
  ## 3) Get results:################################################################################################################
  #### Determine the PC-p-values for all SNP loci | null distributions of correlations from simulated data
  #### Synthesize results output: List of all significant SNPs, their names/locations, their p-values for this phenotype
  ##################################################################################################################################
  
  if(method=="Didelot"){
    ## Calculate correlation between simulated SNPs and phenotype, then between real SNPs and phenotype:
    
    ## convert phenotype factor(y) to numeric:
    phen <- as.numeric(y)
    ## for ease of interpretation, if phen has 2 levels, 1 and 2, make these 0 and 1:
    if(length(unique(phen))==2){
      if(length(phen[-c(which(phen==1), which(phen==2))])==0){
        phen <- replace(phen, which(phen==1), 0)
        phen <- replace(phen, which(phen==2), 1)
      }
    }  
    
    ############################################################################################################################
    ## calculate a correlation btw the (simulated) SNPs in matrix[[i]] and the phenotype, keeping only positive correlations: ##
    ############################################################################################################################
    
    ## create empty list to store correlation results for each of the list of simulated SNPs matrices
    corr.sim <- list()
    
    ## for loop to get correlations:
    for(i in 1:length(snps.mat)){
      ## isolate this element of the list of simulated SNPs matrices
      foo <- snps.mat[[i]]
      
      ## get correlation ##      
      if(!is.null(foo)){
        
        ## CORRELATION ##
        if(test=="cor"){
          all.corr <- cor(foo, phen)  # regular correlation... 
        } # end test cor
        
        ## SCORE ##
        if(test=="score"){
          # ~ Correlation "SCORE" = ((nS1P1 + nS0P0) - (nS1P0 + nS0P1) / (n.total)) ## must be calculated for each SNP individually...   
          all.corr <- sapply(c(1:ncol(foo)), function(e) 
            (((length(which(foo[which(phen==1),e]==1)) + length(which(foo[which(phen==0),e]==0)))
              - (length(which(foo[which(phen==1),e]==0)) + length(which(foo[which(phen==0),e]==1))))
              / nrow(foo))) ## hist(all.corr) ## should be symmetric about 0, and in [0,1]...
        } # end test score
        
        ## FISHER'S EXACT TEST ##
        if(test=="fisher"){
          # fisher's exact test? ##
          all.corr <- sapply(c(1:ncol(foo)), function(e) 
            fisher.test(foo[,e], y=phen, alternative="two.sided")$p.value) 
              ## two.sided bc we want to know if inds w the phen have EITHER more 1s or 0s 
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
      corr.dat <- sapply(c(1:ncol(snps)), function(e) cor(snps[,e], phen)) # regular correlation...
    } # end test cor
    
    ## SCORE ##
    if(test=="score"){
      # ~ Correlation "SCORE" = ((nS1P1 + nS0P0) - (nS1P0 + nS0P1) / (n.total)) ## must be calculated for each SNP individually...   
      corr.dat <- sapply(c(1:ncol(snps)), function(e) 
        (((length(which(snps[which(phen==1),e]==1)) + length(which(snps[which(phen==0),e]==0)))
          - (length(which(snps[which(phen==1),e]==0)) + length(which(snps[which(phen==0),e]==1)))) 
         / nrow(snps)))
    } # end test score
    
    ## FISHER'S EXACT TEST ##
    if(test=="fisher"){
      corr.dat <- sapply(c(1:ncol(snps)), 
                         function(e) fisher.test(snps[,e], y=phen, alternative="two.sided")$p.value) 
                            ## two.sided bc we want to know if inds w the phen have EITHER more 1s or 0s
            
      #             ## TEMPORARY: for african crypto whose snps may not all be snps right now...
      #             corr.dat <- list()
      #             for(i in 1:ncol(snps)){
      #               if(length(levels(as.factor(snps[,i])))< 2){
      #                 corr.dat[[i]] <- NA
      #               }else{
      #                 corr.dat[[i]] <- fisher.test(snps[,i], y=phen, alternative="two.sided")$p.value
      #               }
      #             }
      #             corr.dat <- as.vector(unlist(corr.dat))
      #             corr.dat <- corr.dat[!is.na(corr.dat)]    
      
    } # end test fisher
    

    ## USE ABSOLUTE VALUE
    corr.dat <- abs(corr.dat)
    
    
    ## get empirical p-value based on corr.sim ##
    ## INTERPRETATION: you want the threshold above which only 0.001 proportion (ie. 0.1%) of the data lies:
    
    corr.sim <- as.vector(unlist(corr.sim))
    
    ## CORRELATION ##
    if(test=="cor"){
      ## Identify threshold of significance
      thresh <- quantile(corr.sim, probs=1-p.value)
      ## Identify (real) SNPs w correlations > thresh:
      sig.snps <- which(corr.dat > thresh)
    } # end test cor
    
    ## SCORE ##
    if(test=="score"){
      ## Identify threshold of significance
      thresh <- quantile(corr.sim, probs=1-p.value)
      ## Identify (real) SNPs w correlations > thresh:
      sig.snps <- which(corr.dat > thresh)
    } # end test score
    
    ## FISHER'S EXACT TEST ##
    if(test=="fisher"){  
      ## Identify threshold of significance
      thresh <- quantile(corr.sim, probs=p.value)
      ## Identify (real) SNPs w correlations > thresh:
      sig.snps <- which(corr.dat < thresh)
    } # end test fisher

    
    
    ## get list of those correlation values
    sig.corrs <- corr.dat[sig.snps]
    ## get the list of those SNPs (ie. their locus names)
    sig.snps <- dimnames(snps)[[2]][sig.snps]
    ## re-order list of sig.snps and sig.corrs by value of sig.corr
    NWO <- order(sig.corrs, decreasing=FALSE)
    sig.snps <- sig.snps[NWO]
    sig.corrs <- sig.corrs[NWO]
    
    
    gc()
    
    ##
    
  } # end  method Didelot
  
  
  ########################################
  
  
  
  
  ## for Fitch, we want to calculate the site-specific correlations between the SNP at each locus and the trait, getting
  ## the distributions, across all reps, of the correlations for each site. 
  
  if(method=="Fitch"){
    ## Calculate correlation between simulated SNPs and phenotype
    ## convert phenotype factor(y) to numeric:
    phen <- as.numeric(y)
    ## create empty list to store correlation results for each of the list of simulated SNPs matrices
    corr.sim <- list()
    ## and an empty list, each element of which will contain all the SNPs in column i from each snp set made (ie from each rep)
    foo <- list()
    
    ## for loop to get correlations:
    for(i in 1:ncol(snps.mat[[1]])){ ## all SNPs matrices wil have same ncol as they are generated by permutation...
      ## for each site in the SNPs matrix, get the distribution of correlations across all reps:
      ## make foo a set of all the SNP column i's from each element of the list, snps.mat
      foo[[i]] <- sapply(c(1:length(snps.mat)), function(e) snps.mat[[e]][,i])
      ## get the distribution of the correlations between those snp columns and phen, keeping only positive correlations:
      if(!is.null(foo[[i]])){
        all.corr <- cor(foo[[i]], phen)  # regular correlation...
        ## USE ABSOLUTE VALUE??
        corr.sim[[i]] <- abs(all.corr)
      }else{
        corr.sim[[i]] <- NULL
      }
    } # end for loop
    
    
    ## Calculate correlations between real SNPs and phenotype
    corr.dat <- sapply(c(1:ncol(snps)), function(e) cor(snps[,e], phen)) # regular correlation...
    ## USE ABSOLUTE VALUE??
    corr.dat <- abs(corr.dat)
    
    ## get empirical p-value based on corr.sim FOR EACH SNP LOCUS
    ## make a list to store the threshold sig cut-offs for each SNP's correlation w the phenotype
    thresh <- list()
    corr.sim.old <- corr.sim
    corr.sim <- list()
    are.sig.snps <- list()
    sig.snps <- NA
    sig.corrs <- NA
    ## we only need to calculate for every other site:
    toKeep <- seq(1, ncol(snps), 2)
    for(i in toKeep){
      ## INTERPRETATION: you want the threshold above which only 0.001 proportion (ie. 0.1%) of the data lies:
      round <- which(toKeep==i)
      corr.sim[[round]] <- as.vector(unlist(corr.sim.old[[i]]))
      thresh[[round]] <- quantile(corr.sim[[round]], probs=1-0.001)
    } # end for loop
    
    for(i in 1:length(thresh)){
      ## Identify (real) SNPs w correlations > thresh:
      are.sig.snps[[i]] <- corr.dat[i] > thresh[[i]]
      ## get list of those correlation values
      if(are.sig.snps[[i]]==TRUE){
        if(is.na(sig.snps[1])){
          sig.snps[1] <- dimnames(snps)[[2]][toKeep[i]]
        }else{
          sig.snps[(length(sig.snps)+1)] <- dimnames(snps)[[2]][toKeep[i]]
        }
        if(is.na(sig.corrs[1])){
          sig.corrs[1] <- corr.dat[i]
        }else{
          sig.corrs[(length(sig.corrs)+1)] <- corr.dat[i]
        }  
      }
    } # end for loop
    
    ## re-order list of sig.snps and sig.corrs by value of sig.corr  
    #NWO <- order(sig.corrs, decreasing=FALSE)
    ## TO DO: Would be better to re-order by ~ p-value (or the extent to which 
    ## the correlations exceed the threshold for their SNP... )
    ## something like: 
    keepMe <- which(dimnames(snps)[[2]][toKeep] %in% sig.snps)
    thresh <- as.vector(unlist(thresh[keepMe]))
    NWO <- order(sapply(c(1:length(sig.corrs)), function(e) thresh[e]/sig.corrs[e]), decreasing=TRUE)
    sig.snps <- sig.snps[NWO]
    sig.corrs <- sig.corrs[NWO]
    thresh <- thresh[NWO]
    
    
    gc()
    
    ##
    
  } # end method Fitch
  
  
  
  
  
  
  ##################################################################################################################################
  ## 4) (A) Plot the distribution (ie. the frequency of each level (value) of correlation between the simulated SNPs and the phenotype)
  #### Add an annotations to this plot to identify the locations of: 
  ###### (1) the cut-off threshold for significance for this distribution*
  ########## * given some correction for multiple testing (????) -- argument for p < 0.001 = exploratory analysis only (hypothesis-generating)
  ####### (2) the value of the correlations for significant SNPs
  ##################################################################################################################################
  
  if(method=="Didelot"){
    
    
    ## plot.dist ###################################################################################################################
    
    
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
      ## counting correlations==upper limit of each bin as falling within that bin...
      sig.counts <- sapply(c(1:(length(h$breaks)-1)), function(e) length(X[which(X[which(X <= h$breaks[(e+1)])] > h$breaks[e])]))
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
              ## divide up the space on the y-axis into increments (adding 1 to the denomenator s.t y-max not exceeded)
              increment <- max(h$counts)/(sig.counts[i]+1)
              Y[[i]] <- increment*c(1:sig.counts[i])
            }  
          }
        } # end for loop
      }
      Y <- as.vector(unlist(Y))
      
      if(plot.dist==TRUE){
        
        ## plot histogram of correlations btw real SNPs and phenotype: ##
        plot(h, main="Distribution of SNP-trait correlations \n (with significant SNPs indicated if present)")
        ## ADD threshold line in red on x-axis where thresh hits... 
        abline(v=thresh, col="grey", lwd=2, lty=2)
        ## label threshold line(?)
        text(x=thresh, y=(max(h$counts)), labels="significance threshold", col="grey", pos=2, font=4)
        
        ## Only ask to draw arrows if at least 1 significant SNP:
        if(length(sig.snps) > 0){
          ## ADD arrows pointing from each label to position on X-axis:
          arrows(x0=X , y0=(Y-(max(h$counts)/50)) , x1=X , y1=0 , col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs > threshold at their location on the x-axis:
          text(x=X, y=Y, labels=sig.snps, col="red", font=2)
        }else{
          text(x=thresh, y=(max(h$counts)/4)*3, labels="no significant SNPs found", col="red", font=2, pos=2)
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
      ## counting correlations==upper limit of each bin as falling within that bin...
      sig.counts <- sapply(c(1:(length(h$breaks)-1)), function(e) length(X[which(X[which(X <= h$breaks[(e+1)])] > h$breaks[e])]))
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
              ## divide up the space on the y-axis into increments (adding 1 to the denomenator s.t y-max not exceeded)
              increment <- max(h$counts)/(sig.counts[i]+1)
              Y[[i]] <- increment*c(1:sig.counts[i])
            }  
          }
        } # end for loop
      }
      Y <- as.vector(unlist(Y))
      
      if(plot.dist==TRUE){
        
        ## plot histogram of correlations btw real SNPs and phenotype: ##
        plot(h, main="Distribution of SNP-trait correlations \n (with significant SNPs indicated if present)")
        ## ADD threshold line in red on x-axis where thresh hits... 
        abline(v=thresh, col="grey", lwd=2, lty=2)
        ## label threshold line(?)
        text(x=thresh, y=(max(h$counts)), labels="significance threshold", col="grey", pos=2, font=4)
        
        ## Only ask to draw arrows if at least 1 significant SNP:
        if(length(sig.snps) > 0){
          ## ADD arrows pointing from each label to position on X-axis:
          arrows(x0=X , y0=(Y-(max(h$counts)/50)) , x1=X , y1=0 , col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs > threshold at their location on the x-axis:
          text(x=X, y=Y, labels=sig.snps, col="red", font=2, pos=2)
        }else{
          text(x=thresh, y=(max(h$counts)/4)*3, labels="no significant SNPs found", col="red", font=2, pos=2)
        }
      } # end plot.dist
    } # end test score
    
    
    #########################
    ## FISHER'S EXACT TEST ##
    #########################
    if(test=="fisher"){  
      ## plot.dist ##
      
      #################################################################
      ## plot histogram of correlations btw real SNPs and phenotype: ##
      #################################################################
      
      ## get histogram of correlations btw real SNPs and phenotype:
      h <- hist(corr.dat, plot=FALSE)
      ## get X and Y coords for labelling positions of all sig SNPs
      X <- thresh
      if(length(sig.corrs) > 0) X <- sig.corrs
      ## get the number of sig SNPs in each bin of the histogram
      ## counting correlations==upper limit of each bin as falling within that bin...
      sig.counts <- sapply(c(1:(length(h$breaks)-1)), function(e) length(X[which(X[which(X <= h$breaks[(e+1)])] > h$breaks[e])]))
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
              ## divide up the space on the y-axis into increments (adding 1 to the denomenator s.t y-max not exceeded)
              increment <- max(h$counts)/(sig.counts[i]+1)
              Y[[i]] <- increment*c(1:sig.counts[i])
            }  
          }
        } # end for loop
      }
      Y <- as.vector(unlist(Y))
      
      if(plot.dist==TRUE){
        
        ## plot histogram of correlations btw real SNPs and phenotype: ##
        plot(h, main="Distribution of SNP-trait correlations \n (with significant SNPs indicated if present)")
        ## ADD threshold line in red on x-axis where thresh hits... 
        abline(v=thresh, col="grey", lwd=2, lty=2)
        ## label threshold line(?)
        text(x=thresh, y=(max(h$counts)), labels="significance threshold", col="grey", pos=4, font=4)
        
        ## Only ask to draw arrows if at least 1 significant SNP:
        if(length(sig.snps) > 0){
          ## ADD arrows pointing from each label to position on X-axis:
          arrows(x0=X , y0=(Y-(max(h$counts)/50)) , x1=X , y1=0 , col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs > threshold at their location on the x-axis:
          text(x=X, y=Y, labels=sig.snps, col="red", font=2, pos=4)
        }else{
          text(x=thresh, y=(max(h$counts)/4)*3, labels="no significant SNPs found", col="red", font=2, pos=4)
        }
      } # end plot.dist
    } # end test fisher
    

    
    
    
    ## plot.null.dist ###############################################################################################################
    
    ###############################################################
    ## OR--add SNP annotations to histogram of null distribution ##
    ###############################################################
    
    #################
    ## CORRELATION ##
    #################
    if(test=="cor"){
      ## plot.null.dist ##
      
      #############################################################################
      ## plot (null) histogram of correlations btw SIMULATED SNPs and phenotype: ##
      #############################################################################
      ## plot correlations btw simulated SNPs and phenotype:
      h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)
      
      ## get alternate (null dist) label heights:
      X <- thresh
      if(length(sig.corrs) > 0) X <- sig.corrs
      Y.null <- runif(n=length(sig.snps), min=0.00001, max=max(h.null$counts))

      
      if(plot.null.dist==TRUE){
        if(plot.dist==TRUE) par(ask=TRUE)
        
        ## if the true correlation value for SNP i is > max bin, then extend the x-axis of the plot to accommodate annotation:
        if(max(X) > max(h.null$breaks)){
          ## plot histogram of correlations btw real SNPs and phenotype: ##
          ## EXTENDING THE X-AXIS
          plot(h.null, 
               main="Null distribution of correlations 
           \n (with significant SNPs indicated)",
               xlab="Correlation",
               xlim=c(min(h.null$breaks), max(X)+.05))
        }else{
          ## plot histogram of correlations btw real SNPs and phenotype: ##
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
        text(x=thresh, y=(max(h.null$counts)), labels="significance threshold", pos=2, col="grey", font=4)
        if(length(sig.snps) > 0){
          ## ADD arrows pointing from each label to position on X-axis:
          arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) , x1=X , y1=0 , col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs > threshold at their location on the x-axis:
          text(x=X, y=Y.null, labels=sig.snps, col="red", font=2)
        }else{
          text(x=thresh, y=(max(h.null$counts)/4)*3, labels="no significant SNPs found", col="red", font=2, pos=2)
        }
      } # end plot.null.dist
    } # end test cor
    
    
    ###########
    ## SCORE ##
    ###########
    if(test=="score"){
      ## plot.null.dist ##
      
      #############################################################################
      ## plot (null) histogram of correlations btw SIMULATED SNPs and phenotype: ##
      #############################################################################
      ## plot correlations btw simulated SNPs and phenotype:
      h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)
      
      ## get alternate (null dist) label heights:
      X <- thresh
      if(length(sig.corrs) > 0) X <- sig.corrs
      Y.null <- runif(n=length(sig.snps), min=0.00001, max=max(h.null$counts))

      
      if(plot.null.dist==TRUE){    
        if(plot.dist==TRUE) par(ask=TRUE)
        
        ## if the true correlation value for SNP i is > max bin, then extend the x-axis of the plot to accommodate annotation:
        if(max(X) > max(h.null$breaks)){
          ## plot histogram of correlations btw real SNPs and phenotype: ##
          ## EXTENDING THE X-AXIS
          plot(h.null, 
               main="Null distribution of correlation scores 
           \n (with significant SNPs indicated)",
               xlab="(Correlation) Score",
               xlim=c(min(h.null$breaks), max(X)+.05))
        }else{
          ## plot histogram of correlations btw real SNPs and phenotype: ##
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
        text(x=thresh, y=(max(h.null$counts)), labels="significance threshold", pos=2, col="grey", font=4)
        if(length(sig.snps) > 0){
          ## ADD arrows pointing from each label to position on X-axis:
          arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) , x1=X , y1=0 , col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs > threshold at their location on the x-axis:
          text(x=X, y=Y.null, labels=sig.snps, col="red", font=2, pos=2)
        }else{
          text(x=thresh, y=(max(h.null$counts)/4)*3, labels="no significant SNPs found", col="red", font=2, pos=2)
        }
      } # end plot.null.dist
    } # end test score
    
    
    #########################
    ## FISHER'S EXACT TEST ##
    #########################
    if(test=="fisher"){ 
      ## plot.null.dist ##
      
      #############################################################################
      ## plot (null) histogram of correlations btw SIMULATED SNPs and phenotype: ##
      #############################################################################
      ## plot correlations btw simulated SNPs and phenotype:
      h.null <- hist(as.vector(unlist(corr.sim)), plot=FALSE)
      
      ## get alternate (null dist) label heights:
      X <- thresh
      if(length(sig.corrs) > 0) X <- sig.corrs
      Y.null <- runif(n=length(sig.snps), min=0.00001, max=max(h.null$counts))
      
      if(plot.null.dist==TRUE){    
        if(plot.dist==TRUE) par(ask=TRUE)
        
        ## if the true correlation value for SNP i is < min bin, then extend the x-axis of the plot to accommodate annotation:
        if(min(X) < min(h.null$breaks)){
          ## plot histogram of correlations btw real SNPs and phenotype: ##
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
        text(x=thresh, y=(max(h.null$counts)), labels="significance threshold", pos=4, col="grey", font=4)
        ## only ask to draw arrows if sig snps exist
        if(length(sig.snps) > 0){
          ## ADD arrows pointing from each label to position on X-axis:
          arrows(x0=X , y0=(Y.null-(min(Y.null)/50)) , x1=X , y1=0 , col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs > threshold at their location on the x-axis:
          text(x=X, y=Y.null, labels=sig.snps, col="red", font=2, pos=4)
        }else{
          text(x=thresh, y=(max(h.null$counts)/4)*3, labels="no significant SNPs found", col="red", font=2, pos=4)
        }
      } # end plot.null.dist
    } # end test fisher
    
    

    
    
    } # end Didelot method
  
  
  ##################################################################################################################################
  #### B) Plot the distributions* (ONLY for the SIGNIFICANT SNPs) of the correlations btw the phen and the randomized SNPs
  ########## *as these will be specific to each site...?
  #### Add the 2 annotations as above for Didelot method... 
  ##################################################################################################################################
  
  if(method=="Fitch"){
    
    ## plot Fitch ##
    
    ## Only call to plot if there are some significant SNPs:
    if(length(sig.snps) > 0){
      
      ## we only need to calculate for every other site:
      toKeep <- seq(1, ncol(snps), 2)
      ## get indices of sig snps
      snpsToKeep <- which(dimnames(snps)[[2]][toKeep] %in% sig.snps)
      
      ## for each of the snps in snpsToKeep, plot the null distribution of the correlations btw that SNP and the phen
      if(length(snpsToKeep) > 1) par(ask=TRUE)
      for(i in snpsToKeep[NWO]){
        round <- which(snpsToKeep[NWO]==i)
        ######################################################################
        ## plot histogram of correlations btw simulated SNPs and phenotype: ##
        ######################################################################
        ## get histogram of correlations btw simulated SNPs and phenotype:
        h <- hist(corr.sim[[i]], plot=FALSE)
        ## get X and Y coords for labelling position of sig SNP
        ## get position on the x-axis of the labels
        X <- sig.corrs[round]
        ## get height of labels
        Y <- max(h$counts)/4
        
        if(plot.Fitch==TRUE){
          
          ## if the true correlation value for SNP i is > max bin, then extend the x-axis of the plot to accommodate annotation:
          if(X > max(h$breaks)){
            ## plot histogram of correlations btw real SNPs and phenotype: ##
            ## EXTENDING THE X-AXIS
            plot(h, 
                 main=paste("Distribution of SNP-trait correlations for SNP", 
                            dimnames(snps)[[2]][toKeep][i], 
                            "\n (with true correlation indicated)"),
                 xlab="Correlation",
                 xlim=c(min(h$breaks), X+.05))
          }else{
            ## plot histogram of correlations btw real SNPs and phenotype: ##
            ## WITHOUT EXTENDING THE X-AXIS
            plot(h, 
                 main=paste("Distribution of SNP-trait correlations for SNP", 
                            dimnames(snps)[[2]][toKeep][i], 
                            "\n (with true correlation indicated)"))           
          }
          
          ## ADD threshold line in red on x-axis where thresh hits... 
          abline(v=thresh[[round]], col="grey", lwd=2, lty=2)
          ## label threshold line(?)
          text(x=thresh[[round]], y=(max(h$counts)), labels="significance threshold", col="grey", pos=2, font=4)
          ## ADD arrows pointing from each label to position on X-axis:
          arrows(x0=X , y0=Y , x1=X , y1=0 , col="blue", length=0.1, lwd=1)
          ## add annotation text labelling SNPs > threshold at their location on the x-axis:
          text(x=X, y=Y, labels=sig.snps[round], col="red", font=2)
        } # end plot.Fitch
      } # end for loop
    }
    
  } # end Fitch/ modified-Farhat method
  
  par(ask=FALSE)
  
  ##################################################################################################################################
  ## 5) Return results list (ie. "call" (parameters input), sig SNPs found), print plots (if plots=TRUE) ###########################
  ##################################################################################################################################
  
  if(length(sig.snps)==0) sig.snps <- sig.corrs <- NULL
  
  #   ## CORRELATION ##
  #   if(test=="cor"){
  #     
  #   } # end test cor
  #   
  #   ## SCORE ##
  #   if(test=="score"){
  #     
  #   } # end test score
  #   
  #   ## FISHER'S EXACT TEST ##
  #   if(test=="fisher"){  
  #     
  #   } # end test fisher
  
  
  ## Get counts for n.sig.snps in each cell of the contingency table: 
  snps.toKeep <- snps[,which(dimnames(snps)[[2]] %in% sig.snps)]
  S1P1 <- sapply(c(1:ncol(snps.toKeep)), function(e) length(which(snps[which(phen==1),e]==1)))
  S0P0 <- sapply(c(1:ncol(snps.toKeep)), function(e) length(which(snps[which(phen==0),e]==0)))
  S1P0 <- sapply(c(1:ncol(snps.toKeep)), function(e) length(which(snps[which(phen==0),e]==1)))
  S0P1 <- sapply(c(1:ncol(snps.toKeep)), function(e) length(which(snps[which(phen==1),e]==0)))
  
  ## TBD: calculate p-values for sig.snps
  p.vals <- rep(0, length(sig.snps))
  
  ## make a data.frame containing all relevant output for sig.snps
  if(length(sig.snps) > 0){
    df <- data.frame(sig.snps, p.vals, sig.corrs, S1P1, S0P0, S1P0, S0P1)
    names(df) <- c("SNP.locus", "p.value", "Test.statistic", "S1P1", "S0P0", "S1P0", "S0P1")
  }else{
    df <- "No significant SNPs found."
  }

  ## TO DO: 
  ## ADD MANHATTAN PLOT
  ## CALCULATE and ADD p-values FOR EACH SIG SNP...
  ## ALSO RETURN ALL 4 OF THE NUMBERS IN ALL.CORR FOR EACH SIG SNP (ie S0P0, S1P1, S0P1, S1P0)
  
  if(method=="Didelot"){
    results <- list()
    results[[1]] <- h
    results[[2]] <- h.null
    results[[3]] <- summary(corr.sim)
    results[[4]] <- thresh
    results[[5]] <- df
    

    names(results) <- c("Histogram of true SNP-trait correlations",
                        "Histogram of null distribution of simulated SNP-trait correlations",
                        "Summary of null distribution",
                        "Significance threshold",
                        "Significant SNPs")

    
  } # end method Didelot
  ##
  
  
  
  if(method=="Fitch"){
    results <- list()
    results[[1]] <- sig.snps
    results[[2]] <- sig.corrs
    results[[3]] <- thresh
    #   results[[4]] <- h  
    names(results) <- c("SNPs significantly associated with trait (lowest-highest excess of corr' > thresh)",
                        "Correlations between (significant) SNPs and phenotype/trait",
                        "Significance threshold")
    #                       ,"histogram of true SNP-trait correlations")
    
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
