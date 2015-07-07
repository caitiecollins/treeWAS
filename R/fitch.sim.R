##############################################
##############################################
## fitch.sim: simulation/ randomization fn ###
##############################################
##############################################


fitch.sim <- function(snps, fitch.n.mts, tree, gen.size, biallelic=TRUE, seed=NULL){
  
  require(adegenet)
  require(ape)
  require(phangorn)
  
  
  ## NOTE--Seed worry: could it cause the n.mts mts to be pseudorandomly assigned to the same branches for each site...?)
  if(!is.null(seed)) set.seed(seed) 
  
  ## simulate genotype for root individual:
  gen.root <- sample(c("a", "c", "g", "t"), gen.size, replace=TRUE)
  ## get the total time (ie. branch length from root to tip)
  # CHECK FOR ULTRAMETRICITY (ie. that all branch.lengths are equal):
  if(is.ultrametric(tree)==TRUE){
    ## if we have a real ultrametric!!!!! tree (eg. the inputted one):
    branch.lengths <- distRoot(tree)
  }else{
    ## enfore ultrametricity (or find alternative...)
    tree <- compute.brtime(tree)
    ## if we have a real ultrametric!!!!! tree (eg. the inputted one):
    branch.lengths <- distRoot(tree)
  }
  
  ## tree's tip.labels must be numeric (for Fitch parsimony step)
  tree$tip.label <- c(1:length(tree$tip.label))
  
  
  time.total <- branch.lengths[1]
  ## make dummy variables in which to store the resulting n.mts variables:
  L <- lambda <- n.mts <- new.nt <- NA
  snps.loci <- list()
  

## generate mts for each of n.ind genomes: 
## NOTE: We're working from the root down,  so from the bottom row of 
## tree$edge up to the top row --> more Mts on longer branches... 

# if(method=="fitch"){
  
  #####################################################
  ## FITCH pseudo-simulation randomization procedure ##
  #####################################################
  
  ## n.mts and sites of those mts are fixed/ submitted as input to the fn  
  ## BUT, the branches that they fall on are selected w Prob ~ branch lengths
  ## ie. use tree$edge.length as a probability distribution for EACH mt site
  ## so: Sample, from the set of all edges (from 1 to nrow(tree$edge)), the total n.mts to occur at a given site,
  ## with the probability of choosing a given branch for any or all (therefore sample w replacement)
  ## of those mts being proportional to its relative branch length... 
  
  ##################################
  ## required INPUT from treeWAS: ##
  ##################################
  n.mts <- fitch.n.mts ## from treeWAS
  #   snps <- snps
  #   tree <- tree
  ###################################
  
  toKeep <- seq(1, ncol(snps), 2)
  #snps.loci <- colnames(snps)[toKeep]
  
  branches <- list()
  snps.loci <- list()
  for(i in 1:length(n.mts)){
    ## branches i will contain the branches on which the mts at site colnames(snps)[toKeep][i] will occur
    branches[[i]] <- sample(c(1:nrow(tree$edge)), n.mts[i], replace=TRUE, prob=(tree$edge.length)/sum(tree$edge.length))
    snps.loci[[i]] <- rep(colnames(snps)[toKeep][i], n.mts[i])
  }
  ## collapse lists:
  branches <- as.vector(unlist(branches))
  snps.loci <- as.vector(unlist(snps.loci))
  ## re-order s.t branches listed 1 to n.branches and mt sites sorted to match...
  snps.loci <- snps.loci[order(branches, decreasing=FALSE)]
  branches <- sort(branches, decreasing=FALSE)
  
  ## rebuild output in same format as standard simulation approach output...
  ## want to end up with:
  #### - n.mts on branches i (ie branches 1 to n.branches)
  #### - the sites of the mts occuring on branches i
  N.MTS <- list()
  SNPS.LOCI <- list()  
  for(i in 1:nrow(tree$edge)){
    if(!any(names(table(branches))==i)){
      N.MTS[[i]] <- 0
    }else{
      N.MTS[[i]] <- table(branches)[which(names(table(branches))==i)]
      names(N.MTS[[i]]) <- NULL
    }
    if(!any(branches==i)){
      SNPS.LOCI[[i]] <- sample(c(1:100),0)
    }else{
      ## check:
      #all(colnames(snps)[which(colnames(snps) %in% snps.loci[which(branches==i)])]==snps.loci[which(branches==i)])
      ## check 2: --> THUS, using which() on only the ODD colnames of SNPs --> INDICES with the same numbers as are in the loci NAMES! 
      #which(colnames(snps)[seq(1, ncol(snps), 2)] %in% snps.loci[which(branches==i)])
      SNPS.LOCI[[i]] <- which(colnames(snps)[seq(1, ncol(snps), 2)] %in% snps.loci[which(branches==i)]) 
    }
  } # end for loop
  
  snps.loci <- n.mts <- NULL
  ## over-write preliminary lists with final output:
  n.mts <- N.MTS
  snps.loci <- SNPS.LOCI
  
## check:  
# sapply(c(1:ncol(snps)), function(e) length(which(as.vector(unlist(snps.loci))==e)))
# fitch.n.mts

  # end of FITCH procedure of pseudo-simulation randomization w fixed n.mts/site. #
  ################################################################################# 
  
# }else{
#   
#   ##################################
#   ## REGULAR simulation procedure ##
#   ##################################
#   
#   ## draw n.mts per branch and their loci 
#   for(i in rev(1:length(tree$edge.length))){  
#     ## get branch length i 
#     L[i] <- rev(tree$edge.length)[i]
#     
#     ## get lambda, the mean of the poisson dist 
#     ## from which we draw the n.mts to occur on branch i
#     lambda[i] <- (theta/2) * L[i]
#     
#     ## get n.mts from poisson dist: 
#     n.mts[i] <- rpois(n=1, lambda[i])
#     
#     ## draw n.mts locations in the genome (allowing for replacement)
#     snps.loci[[i]] <- sample(c(1:gen.size), n.mts[i], replace=TRUE)  
#   } # end for loop
# }


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

#print(length(genomes))

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


# 
# 
# 
# 
# dna <- as.DNAbin(genomes)
# names(dna) <- c(1:length(genomes))
# 
# 
# if(heatmap==TRUE){
#   ## get a distance matrix between the genomes
#   D <- dist.dna(dna, model="JC69")
#   
#   mat <- t(as.matrix(D))
#   mat <- mat[,ncol(mat):1]
#   par(mar=c(1,5,5,1))
#   image(x=1:ncol(mat), y=1:ncol(mat), mat, 
#         col=rev(heat.colors(100)), 
#         xaxt="n", yaxt="n", xlab="", ylab="")
#   axis(side=2, at=c(1:ncol(mat)), 
#        lab=rev(names(dna)), las=2, cex.axis=1)
#   axis(side=3, at=c(1:ncol(mat)), 
#        lab=names(dna), las=1, cex.axis=1)
#   ## return margin parameter to default:
#   par(mar=c(5,4,4,2)+0.1)
# }
# 
# if(plot==TRUE){
#   plot(tree, show.tip=FALSE, edge.width=2)
#   title("Coalescent tree")
#   tiplabels(text=tree$tip.label, cex=1, adj=-.5) 
#   nodelabels(text=rev(unique(tree$edge[,1])), cex=0.75)
#   edgelabels(text=paste("e", c(1:nrow(tree$edge)), sep="."), cex=0.66)
#   edgelabels(text=rev(n.mts), col="red", frame="none", cex=1.1, adj=c(1,-0.5))
# }
# 
# if(plot=="simple"){
#   plot(tree, show.tip=TRUE, edge.width=2, cex=0.6, adj=.25)
#   title("Coalescent tree") 
# }
# 
# D <- dist.dna(dna[1:n.ind], model="JC69")
# 
# if(plot2=="nj"){
#   tree1 <- nj(D)
#   #tree1 <- root(tree1, node=(n.ind+5)) 
#   tree1 <- ladderize(tree1)
#   plot(tree1, edge.width=2)
#   title("Neighbour-joining tree")
#   axisPhylo()
# }
# 
# if(plot2=="UPGMA"){
#   tree2 <- hclust(D, method="average") 
#   plot(tree2, main="")
#   title("UPGMA tree")
# }
# 
# if(plot2=="ml"){
#   dna4 <- as.phyDat(dna[1:n.ind])
#   tre.ini <- nj(dist.dna(dna[1:n.ind], model="JC69"))
#   fit.ini <- pml(tre.ini, dna4, k=n.ind)
#   fit <- optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, 
#                    optQ = TRUE, optGamma = TRUE)
#   
#   ## NOTE--you may want to store these in a results.ml list and return it with your results instead of printing
#   ## OR at least print a message (eg. "Printing maximum-likelihood calculations...") before printing these numbers... 
#   anova(fit.ini, fit)
#   AIC(fit.ini)
#   AIC(fit) 
#   
#   tree3 <- fit$tree
#   tree3 <- root(tree3, node=(n.ind+5))
#   plot(tree3, show.tip=TRUE, edge.width=2)
#   title("Maximum-likelihood tree")
#   axisPhylo()
# }
# 
# par(ask=FALSE)

genomes <- genomes[1:n.ind]

out <- list(genomes, tree)
# out <- genomes
return(out)

} # end fitch.sim































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