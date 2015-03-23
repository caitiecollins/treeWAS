

##############
## tree.sim ##
##############

#'tree.sim
#': simulate genotypes based on a given phylogenetic tree.
#'
#'This function simulates genetic data given an inputed phylogenetic tree structure.  
#' 
#' @param n.ind The number of individuals in the tree. (Defaults to length(tree$tip.label).)
#' @param tree The phylogenetic tree to be used as the basis for the genetic relationship between
#'  the genotypes of the individuals in your simulated dataset.
#' @param sim.by The method by which the simulated data is generated (when \code{method = "Didelot"}).
#' \itemize{
#'  \item If \code{sim.by = "locus"}, then \code{mt.rate} denotes the expected number of mutations to occur at each locus.
#'  \item If \code{sim.by = "branch"}, then \code{mt.rate} denotes the expected number of mutations to occur in the whole genome.
#'  }
#' @param theta An optional integer specifying the mutation rate to be used when simulating comparator datasets.
#'  If \code{mt.rate = NULL}, when \code{sim.by = "locus"}, \code{mt.rate} will be set to 1*2; 
#'  and when \code{sim.by = "branch"}, \code{mt.rate} will be set to \code{gen.size}*2.
#' @param biallelic A logical indicating whether SNPs should be biallelic or otherwise. (TRUE by default.)
#' @param seed An optional number to be used as the seed of pseudo-random number generation processes that occur within the \code{tree.sim} 
#'  function. Identical pseudo-random processes generated with the same seed will have identical results. 
#' @param plot A logical indicating whether to plot the phylogenetic tree.
#' @param heatmap A logical indicating whether to plot a heatmap of the genetic distances between the individuals in your tree and simulated dataset.
#' @param plot2 An optional character specifying the method of phylogenetic reconstruction run to confirm that a 
#'  phylogenetic tree which is similar to the inputed phylogeny can be generated based on the simulated genomes alone. 
#' 
#' @author Caitlin Collins \email{caitlin.collins12@@imperial.ac.uk}
#' @export
#' @examples
#' \dontrun{
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
#' ## run tree.sim 
#' result <- tree.sim(n.ind=length(tree$tip.label), tree=tree, 
#'            gen.size=1000, sim.by="locus", theta=NULL, biallelic=TRUE, seed=NULL, plot=TRUE,
#'            heatmap=FALSE, plot2="ml")
#' 
#' }
#' @import adegenet ape phangorn 


# n.ind <- 10 # n.genomes you want to end up with
# gen.size <- 1000000 # bases
# theta <- gen.size*2 # (if sim.by=="branch")# OR # 1*2 # (if sim.by=="locus") 
# biallelic <- TRUE # if TRUE, select ONLY complementary nt; if FALSE, select from 3 alternatives (ie. A/C/G/T-current nt)
# seed <- 1 # allow user to control randomization to get reproducible results.
# sim.by <- c("branch", "locus")

#tree <- tree.indian
#out <- tree.sim(n.ind=length(tree$tip.label), tree=tree, gen.size=10000, theta=100*2, biallelic=FALSE, seed=2, plot=TRUE, heatmap=FALSE, plot2="ml")

tree.sim <- function(n.ind=length(tree$tip.label), tree=tree, gen.size=1000, 
                     sim.by="locus", theta=NULL, biallelic=TRUE, seed=NULL, plot=TRUE,
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

  ##################################
  ## REGULAR simulation procedure ##
  ##################################
  
  snps.loci <- list()
  
  ## draw n.mts per branch and their loci 
  for(i in rev(1:length(tree$edge.length))){  
    ## get branch length i 
    L[i] <- rev(tree$edge.length)[i]
    
    ## get lambda, the mean of the poisson dist 
    ## from which we draw the n.mts to occur on branch i (proportional to branch length i out of the total sum of all branch lengths in the tree)
    lambda[i] <- ((theta/2) * (L[i] / time.total))
    
    ## get n.mts from poisson dist: 
    n.mts[i] <- rpois(n=1, lambda[i])
    
    ## draw n.mts locations in the genome (allowing for replacement)
    snps.loci[[i]] <- sample(c(1:gen.size), n.mts[i], replace=TRUE)  
  } # end for loop

  } # end sim.by branch
  
  
  ##############
  ## BY LOCUS ## ## NOTE: ## for sim.by LOCUS, MT.RATE ~= the EXPECTED n.mts AT EACH SITE (eg. 1(*2)) !
  ##############             #######################################################################
  if(sim.by=="locus"){
    
    #####################################
    ## SNP-BY-SNP simulation procedure ##
    #####################################
    
    ## draw the number of mutations to occur at each site:
    n.mts <- rpois(n=gen.size, lambda=(theta/2))
    
    ## for each site, draw the branches to which you will assign the mts for this site (~ branch length):
    snps.loci <- sapply(c(1:length(n.mts)), function(e) sample(c(1:length(tree$edge.length)), n.mts[e], replace=TRUE, prob=tree$edge.length))
    ## rearrange snps.loci s.t it becomes a list of length tree$edge.length, 
    ## each element of which contains the locations of the mutations that will occur on that branch
    snps.loci <- sapply(c(1:length(tree$edge.length)), function(f) seq_along(snps.loci)[sapply(snps.loci, function(e) f %in% e)])

  } # end sim.by locus
  
  ## mini fn testing for output==integer(0) ##
  .is.integer0 <- function(x){
    is.integer(x) && length(x) == 0L
  } # end .is.integer0()
  
  
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
  for(i in 1:length(x)){
    ## for all genomes other than root, we mutate the 
    ## genome of the node preceding it, according to snps.loci.
    ## Draw new nts for each locus selected for mutation: 
    if(!.is.integer0(snps.loci.unique[[i]])){
      if(biallelic==FALSE){
        new.nts[[i]] <- sapply(c(1:length(snps.loci.unique[[i]])), function(e) 
          sample(c("a", "c", "g", "t")[-which(c("a", "c", "g", "t") 
                                              %in% genomes[[tree$edge[x[i],1]]]
                                              [snps.loci.unique[[i]][e]])], 1))
      }else{
        new.nts[[i]] <- sapply(c(1:length(snps.loci.unique[[i]])), function(e) 
          selectBiallelicSNP(c("a", "c", "g", "t")[which(c("a", "c", "g", "t") 
                                                         %in% genomes[[tree$edge[x[i],1]]]
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
      temp <- genomes[[tree$edge[x[i],1]]]
      temp[snps.loci.unique[[i]]] <- new.nts[[i]]    
      genomes[[tree$edge[x[i],2]]] <- temp 
      
    }else{
      ## if no mts occur on branch, set genotype of 
      ## downstream individual to be equal to ancestor's
      genomes[[tree$edge[x[i],2]]] <- genomes[[tree$edge[x[i],1]]]    
    }  
  } # end of for loop selecting new nts at mutator loci
  
  
  
  
  
  dna <- as.DNAbin(genomes)
  names(dna) <- c(1:length(genomes))
  
  
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
  
  if(plot==TRUE){
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
      edgelabels(text=sapply(c(1:length(snps.loci)), function(e) length(snps.loci[[e]])), col="red", frame="none", cex=1.1, adj=c(1,-0.5))
    } # end sim.by locus
    
  }
  
  if(plot=="simple"){
    plot(tree, show.tip=TRUE, edge.width=2, cex=0.6, adj=.25)
    title("Coalescent tree") 
  }
  
  D <- dist.dna(dna[1:n.ind], model="JC69")
  
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
    
    ## NOTE--you may want to store these in a results.ml list and return it with your results instead of printing
    ## OR at least print a message (eg. "Printing maximum-likelihood calculations...") before printing these numbers... 
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
  
  genomes <- genomes[1:n.ind]
  
  out <- list(genomes, tree)
  return(out)
  
} # end coalescent.sim





#########################################################################
## selectBiallelicSNP: fn that returns the alternative nt| the nt inut ##
########################################################################

## NOTE-- while this is not inherently the definition of a biallelic SNP, it is currently suiting my purposes
#### by fulfilling the function of ensuring that sites always revert back and for between one state and ONE other,
#### hence never creating any triallelic sites or tetralellic sites --> binary encoding guaranteed to work fine. 

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
##########################################################################