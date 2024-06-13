
## WARNING:
## SOME LINES OF SNP.SIM.Q.R HAVE FALLEN OUT OF DATE w SNP.SIM.R.
## BEFORE USING SNP.SIM.Q, UPDATE (at least!):
## sapply, add n.mts > l.edge check,
## convert snps to logical, for loop (selectBiallelicSNP --> !l[[i]])

###############
## snp.sim.Q ##
###############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Aternative SNPs simulation fn.
#'
#' Currently under development. Please use the regular snp.sim function to simulate genetic data.
#'
#' @param n.snps An integer specifying the number of snps columns to be simulated.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#'
#' @examples
#' ## Example ##
#'
#' @import adegenet
#' @rawNamespace import(ape, except = zoom)
#' @importFrom Hmisc all.is.numeric
#' @importFrom phangorn midpoint
#'
#' @export

########################################################################
#  @useDynLib phangorn, .registration = TRUE

snp.sim.Q <- function(n.snps = 10,
                      tree = coalescent.tree.sim(100),
                      phen.reconstruction, ## provide phen.rec and use this to set s=n.phen.subs... 
                      s = 15, ## n.subs for correlated snps.assoc
                      af = 10, ## association factor (relative odds of assoc vs. non-assoc in Q mat)
                      snp.root = NULL,
                      # assoc.prob = 100,
                      heatmap = FALSE,
                      reconstruct = FALSE,
                      dist.dna.model = "JC69",
                      grp.min = 0.25,
                      row.names = NULL,
                      set=3,
                      seed=1){
  
  # require(adegenet)
  # require(ape)
  
  temp <- snps <- snps.assoc <- tree.reconstructed <- sets <- phen <- phen.nodes <- NULL
  
  n.snps.assoc <- n.snps
  
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
  
  ## check phen.rec:
  if(!is.null(phen.reconstruction)){
    if(length(unique(phen.reconstruction[!is.na(phen.reconstruction)])) == 2){
      phen.reconstruction <- as.logical(as.numeric(phen.reconstruction))
      
      ## Infer n.phen.subs:
      n.phen.subs <- length(which(phen.reconstruction[tree$edge[,1]] != phen.reconstruction[tree$edge[,2]]))
    }else{
      phen.reconstruction <- NULL
      n.phen.subs <- NULL
      warning("phen.reconstruction must be binary. Ignorning (sorry).")
    }
  }else{
    phen.reconstruction <- NULL
    n.phen.subs <- NULL
  }
  
  
  ## For n.subs = n or = dist approaches:
  gen.root <- sample(c(TRUE, FALSE), gen.size, replace=TRUE)
  if(is.null(phen.reconstruction)){
    phen.root <- sample(c(TRUE, FALSE), gen.size, replace=TRUE)
  }else{
    phen.root <- phen.reconstruction[min(tree$edge[,1])]
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
    
    
  ######################################################################################################################################################################
  ######################################################################################################################################################################
  
  
  #########################
  ## SIM ASSOCIATED SNPS ##
  #########################
  
  ## Need to treat ASSOCIATED SNPs differently:
  ## (non.assoc.snps do NOT need to pass the "while" check;
  ## they just need to match phen.loci at this point.)
  
  snps.assoc.nodes <- phen.nodes <- NULL
  
  if(n.snps.assoc != 0){
    
    ## (*) Now assuming phen.rec provided as input
    ## If s provided, s overrides n.phen.subs; 
    ## else s=n.phen.subs as inferred from phen.reconstruction: 
    if(is.null(s)){
      if(!is.null(n.phen.subs)){
        s <- n.phen.subs
      }else{
        s <- 15 
        warning("Setting s (n.subs for snps.assoc in Q rate matrix) as 15.
                User must provide phen.reconstruction (s=n.phen.subs) or set s argument to specify manually.")
      }
    }
    ## Scale s to account for total edge.length:
    s <- s/sum(tree$edge.length)
    
    ## Set af = the association factor:
    ## *af = the relative odds of var1 changing state towards being in-phase w var2 vs. 
    ## changing to the opposite state, out-of-step w var2, over a given branch length...  
    ## (Q is actually an instantaneous transition rate matrix, but we will account for branch lengths later.)
    if(is.null(af)){
      af <- 10 
      warning("Setting association factor af = 10.")
    }
    
    
    ## Get Q, the dependent transition rate/prob matrix:
    Q.mat <- matrix(c(NA,     1*s, 1*s, 0,
                      1*af*s, NA,  0,   1*af*s,
                      1*af*s, 0,   NA,  1*af*s,
                      0,      1*s, 1*s, NA),
                    nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
    
    diag(Q.mat) <- sapply(c(1:nrow(Q.mat)), function(e) -sum(Q.mat[e, c(1:ncol(Q.mat))[-e]]))
    Q <- Q.mat
    
    
    
    ## get snps.loci for all ASSOCIATED snps, conditional on phen, Q: ##
    snps.assoc.nodes <- list()
    N.OVERLAP <- list()
    
    Qt <- list()
    
    ## get snps.loci for FIRST ASSSOCIATED snp WITH PHEN.loci:
    i <- 1
    edges <- tree$edge
    root.nt <- gen.root[snps.assoc[i]]
    
    
    ## Need to update tree$edge.length probs to mirror sampling without replacement...
    ## BUT only remove branches from the calculation once they have had a sub on them...?
    # N.SUBS.TOTAL <- rpois(1, n.phen.subs)  ############
    
    ###############
    ## FOR LOOP: ## 
    ###############
    for(i in 1:n.snps.assoc){
      ## get nt for root at this locus:
      root.nt <- gen.root[snps.assoc[i]]
      
      snp.node <- as.list(rep(NA, length(unique(as.vector(edges)))))
      
      snp.node[[edges[x[1], 1]]] <- root.nt
      
      ## go from last to first edge in edges:
      for(e in x){
        
        probs <- NULL
        
        ####################################################
        ## get conditional probs for each edge w matexpo! ##
        ####################################################
        ## (run within code...)
        Qt[[e]] <- matexpo(Q*tree$edge.length[e])
        
        P <- Qt[[e]]
        rownames(P) <- rownames(Qt[[e]]) <- rownames(Q)
        colnames(P) <- colnames(Qt[[e]]) <- colnames(Q)
        
        if(snp.node[[edges[e, 1]]] == FALSE & phen.reconstruction[[edges[e, 1]]] == FALSE) probs <- P[1,]
        if(snp.node[[edges[e, 1]]] == FALSE & phen.reconstruction[[edges[e, 1]]] == TRUE) probs <- P[2,]
        if(snp.node[[edges[e, 1]]] == TRUE & phen.reconstruction[[edges[e, 1]]] == FALSE) probs <- P[3,]
        if(snp.node[[edges[e, 1]]] == TRUE & phen.reconstruction[[edges[e, 1]]] == TRUE) probs <- P[4,]
        
        ## Now we KNOW, A PRIORI, the phen.node for the DESCENDANT!
        probs.mod <- replace(probs, which(!as.logical(as.numeric(keepLastN(colnames(Q), 1))) == phen.reconstruction[[edges[e, 2]]]), 0)
        SP.dec <- sample(colnames(Q), 1, prob = probs.mod)
        
        S.dec <- as.logical(as.numeric(keepFirstN(SP.dec, 1)))
        names(S.dec) <- NULL
        
        snp.node[[edges[e, 2]]] <- S.dec
        
      } # end for (e) loop
      
      ## STORE SNPS.ASSOC (FOR ALL NODES):
      snps.assoc.nodes[[i]] <- as.vector(unlist(snp.node))
      
      ## Get proportion overlap btw phen and snps.assoc.i:
      N.overlap <- length(which(phen.reconstruction[1:n.ind] == snps.assoc.nodes[[i]][1:n.ind]))
      N.overlap <- max(N.overlap, (n.ind-N.overlap))
      N.OVERLAP[[i]] <- N.overlap
      
    } # end for loop
    
    ## Bind SNPs.ASSOC into matrix:
    snps.assoc.nodes <- do.call("cbind", snps.assoc.nodes)
    
    
    ###################################################
    ## TEMP -- COMPARE PHEN & ALL SNPS.ASSOC w PLOT: ##
    ###################################################
    par(mfrow=c(2,6))
    plot_phen(tree, phen.nodes=phen.reconstruction, main.title="phen")
    if(n.snps.assoc >= 10){
      ## just plot first 10:
      for(i in 1:5){
        plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE,
                  main.title=paste("snp.assoc", i, sep=" "))
        title(N.OVERLAP[[i]], line=0, font.main=1)
      }
      plot_phen(tree, phen.nodes=phen.reconstruction, main.title="phen")
      for(i in 6:10){
        plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE,
                  main.title=paste("snp.assoc", i, sep=" "))
        title(N.OVERLAP[[i]], line=0, font.main=1)
      }
    }else{
      ## plot up to 5:
      if(n.snps.assoc <= 5){
        par(mfrow=c(1,(n.snps.assoc+1)))
        for(i in 1:n.snps.assoc){
          plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE,
                    main.title=paste("snp.assoc", i, sep=" "))
          title(N.OVERLAP[[i]], line=0, font.main=1)
        }
      }else{
        ## plot btw 5 and 10:
        for(i in 1:5){
          plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE,
                    main.title=paste("snp.assoc", i, sep=" "))
          title(N.OVERLAP[[i]], line=0, font.main=1)
        }
        plot_phen(tree, phen.nodes=phen.reconstruction, main.title="phen")
        for(i in 6:n.snps.assoc){
          plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE,
                    main.title=paste("snp.assoc", i, sep=" "))
          title(N.OVERLAP[[i]], line=0, font.main=1)
        }
      }
    }
    
    
    par(mfrow=c(1,1)) # end temp panel plot
    
    gc()
    
  } # end of snps.assoc generation
  
  
  
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  
  
  ###########################################
  ## GET COMPLETE SNPS MATRIX ("genomes"): ##
  ###########################################
  
  if(!is.null(temp)){
    ## Create snps matrix:
    snps <- temp
    rm(temp)
    
    ## Attach snps.assoc loci to last column:
    if(!is.null(snps.assoc.nodes)){
      snps <- cbind(snps[,c(1:(ncol(snps)-(n.snps.assoc)))], snps.assoc.nodes[1:n.ind,])
    }
    ## keep only rows containing terminal individuals:
    snps <- snps[1:n.ind, ]
  }else{
    
    ## If only snps.assoc simulated (no non-assoc snps):
    snps <- snps.assoc.nodes[1:n.ind,]
  }
  
  
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
  
  ## (NO LONGER) CONVERT TO NUMERIC: ##
  ## Convert from logical to binary SNPs (for terminal nodes only):
  # snps <- replace(snps, which(snps == TRUE), 1)
  
  ## Reassort snps.assoc to new columns:
  if(!is.null(snps.assoc)){
    
    if(!is.null(temp)){
      rm(temp)
    
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
    }
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
  
  ##################
  ## get RESULTS: ##
  ##################
  out <- list(snps, snps.assoc, snps.assoc.nodes, tree.reconstructed, phen, phen.nodes)
  names(out) <- c("snps", "snps.assoc", "snps.assoc.nodes", "tree.reconstructed",  "phen", "phen.nodes")
  
  return(out)
  
} # end snp.sim.Q




