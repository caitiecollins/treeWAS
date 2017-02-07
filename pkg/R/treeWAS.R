



###################
## print.treeWAS ##
###################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Print \code{treeWAS} output.
#'
#' Print the results of \code{treeWAS}, excluding longer data elements within the output.
#'
#' @param x The output returned by \code{treeWAS}.
#' @param sort.by.p A logical indicating whether to sort the results by decreasing p-value (\code{TRUE}) or by locus (\code{FALSE}, the default).
#'
#' @examples
#' ## Example ##
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @import adegenet ape phangorn

########################################################################

print.treeWAS <- function(x, sort.by.p = FALSE){
  cat("\t#################### \n")
  cat("\t## treeWAS output ## \n")
  cat("\t#################### \n")
  cat("\t \n")

  cat("\t#################### \n")
  cat("\t## All findings:  ## \n")
  cat("\t#################### \n")
  cat("Number of significant loci: ")
  print(length(x$treeWAS.combined$treeWAS.combined))

  res <- x$treeWAS.combined$treeWAS.combined
  if(all.is.numeric(res)){
    res <- as.character(sort(as.numeric(res, decreasing=FALSE)))
  }else{
    if(all.is.numeric(removeLastN(res, 2))){
      ord <- match(sort(as.numeric(removeLastN(res, 2), decreasing=FALSE)), removeLastN(res, 2))
      res <- res[ord]
    }
  }
  cat("Significant loci: \n")
  print(res)

  cat("\t \n")

  cat("\t######################## \n")
  cat("\t## Findings by test:  ## \n")
  cat("\t######################## \n")

  N <- c(2:(length(x)-1))

  for(i in N){
    test <- names(x)[i]

    cat("\t", paste(c("######", rep("#", nchar(test)), "######"), collapse=""), " \n")
    cat("\t ## ", test, "test ## \n")
    cat("\t", paste(c("######", rep("#", nchar(test)), "######"), collapse=""), " \n")

    res <- x[[i]]$sig.snps

    cat("Number of significant loci: \n")
    if(is.vector(res)){
      print(0)

      cat("Significance threshold: \n")
      print(x[[i]]$sig.thresh)

    }else{

      print(length(res$SNP.locus))

      cat("Significance threshold: \n")
      print(x[[i]]$sig.thresh)


      cat("Significant loci: \n")
      if(sort.by.p == TRUE){
        print(res)
      }else{
        sl <- sort(res$SNP.locus, decreasing=FALSE)
        ord <- match(sl, res$SNP.locus)
        print(res[ord, ])
      }

    }

  } # end for loop

} # end print.treeWAS

#############
## treeWAS ##
#############


###################################################################################################################################

###################
## DOCUMENTATION ##
###################

#' Phylogenetic tree-based GWAS for microbes.
#'
#' This function implements a phylogenetic approach to genome-wide association studies (GWAS) designed for use in bacteria and viruses.
#' The \code{treeWAS} approach allows for the identification of significant asociations between genotype and phenotype, while accounting for the
#' confounding effects of clonal population structure, stratification (overlap between the population structure and phenotypic distribution),
#' and homologous recombination.
#'
#'
#' @param snps A matrix containing binary genetic data, with individuals in the rows and genetic loci in the columns and both rows and columns labelled.
#' @param phen A vector containing the phenotypic state of each individual, whose length is equal to the number of rows in \code{snps} and which is named with the same set of labels.
#'              The phenotype can be either binary (character or numeric) or continuous (numeric).
#' @param tree A \code{phylo} object containing the phylogenetic tree; or, a character string, one of \code{"nj"}, \code{"ml"}, or \code{"UPGMA"} (the default),
#'                specifying the method of phylogenetic reconstruction.
#' @param n.subs A numeric vector containing the homoplasy distribution (if known, see details), or NULL (the default).
#' @param sim.n.snps An integer specifying the number of loci to be simulated for estimating the null distribution
#'                      (by default \code{10*ncol(snps)}). If memory errors arise during the analyis of a large dataset,
#'                      it may be necesary to reduce \code{sim.n.snps} from a multiple of 10 to, for example, 5x the number of loci.
#' @param test A character string or vector containing one or more of the following available tests of association:
#'              \code{"terminal"}, \code{"simultaneous"}, \code{"subsequent"}, \code{"cor"}, \code{"fisher"}. By default, the first three tests are run (see details).
#'
#' @param snps.reconstruction Either a character string specifying \code{"parsimony"} (the default) or \code{"ML"} (maximum likelihood)
#'                              for the ancestral state reconstruction of the genetic dataset,
#'                              or a matrix containing this reconstruction if it has been performed elsewhere.
#'
#' @param snps.sim.reconstruction A character string specifying \code{"parsimony"} (the default) or \code{"ML"} (maximum likelihood)
#'                                  for the ancestral state reconstruction of the simulated null genetic dataset.
#'
#' @param phen.reconstruction Either a character string specifying \code{"parsimony"} (the default) or \code{"ML"} (maximum likelihood)
#'                              for the ancestral state reconstruction of the phenotypic variable,
#'                              or a vector containing this reconstruction if it has been performed elsewhere.
#' @param na.rm A logical indicating whether columns in \code{snps} containing more than 50\% \code{NA}s should be removed at the outset (TRUE, the default) or not (FALSE).
#'
#' @param p.value A number specifying the base p-value to be set the threshold of significance (by default, \code{0.01}).
#' @param p.value.correct A character string, either \code{"bonf"} (the default) or \code{"fdr"}, specifying whether correction for multiple testing
#'                          should be performed by Bonferonni correction (recommended) or the False Discovery Rate.
#' @param p.value.by A character string specifying how the upper tail of the p-value distribution is to be identified.
#'                      Either \code{"count"} (the default, recommended) for a simple count-based approach or \code{"density"} for a kernel-density based approximation.
#' @param dist.dna.model A character string specifying the type of model to use in reconstructing the phylogenetic tree for
#'                          calculating the genetic distance between individual genomes, only used if \code{tree} is a character string (see ?dist.dna).
#' @param plot.tree A logical indicating whether to generate a plot of the phylogenetic tree (\code{TRUE}) or not (\code{FALSE}, the default).
#' @param plot.manhattan A logical indicating whether to generate a manhattan plot for each association score (\code{TRUE}, the default) or not (\code{FALSE}).
#' @param plot.null.dist A logical indicating whether to plot the null distribution of association score statistics (\code{TRUE}, the default) or not (\code{FALSE}).
#' @param plot.dist A logical indicating whether to plot the true distribution of association score statistics (\code{TRUE}) or not (\code{FALSE}, the default).
#' @param snps.assoc An optional character string or vector specifying known associated loci to be demarked in results plots (e.g., from previous studies or if data is simulated); else NULL.
#' @param filename.plot An optional character string denoting the file location for saving any plots produced; else \code{NULL}.
#' @param seed An optional integer to control the pseudo-randomisation process and allow for identical repeat runs of the function; else \code{NULL}.
#'
###########################################################
#'
#' @details
#' \strong{Data Cleaning}
#'
#' The genetic data matrix, phenotype, and phylogenetic tree (terminal nodes)
#' should all be labelled with corresponding sets of names for the individuals.
#' The order of individuals does not matter, as they will be rearranged to match within the function.
#'
#' Any individual that is not present in one or more of the genetic data matrix,
#' the phenotypic variable, and/or the phylogenetic tree must be removed.
#'
#' If, in the genetic data matrix, redundant columns are present for binary loci
#' (ie. Denoting the state of the second allele as the inverse of the previous colummn), these should be removed.
#' For loci with more than two alleles, all columns should be retained. The removal of redundant binary columns
#' can be done by hand; but, the function \code{get.binary.snps} may also be used. This function requires
#' column names to have a two-character suffix and unique locus identifiers. The function expects the suffixes
#' ".a", ".c", ".g", ".t" (e.g., "Locus_123243.a", "Locus_123243.g"), though alternative two-character suffixes can be
#' used (e.g., "Locus_123243_1", "Locus_123243_2") by setting the argument \code{force = TRUE}.
#' Please also be careful not to accidentally remove any purposeful duplications with repeated names;
#' for example, if you have deliberately duplicated columns without subsequently generating unique column names
#' (e.g., by expanding unique columns according to an index returned by ClonalFrameML).
#'
#' Missing data is permitted (denoted by NA values only) in the genetic data matrix,
#' but if more than 50% of any column is composed of NAs,
#' we recommend that this column be removed from the dataset.
#' Note that the removal of majority-missing columns will be performed automatically within \code{treeWAS}.
#' If, for some reason, you do not wish this to be the case, set the \code{na.rm} argument to \code{FALSE}.
#'
#'
###########################################################
#' \strong{Homoplasy Distribution}
#'
#' The homoplasy distribution contains the number of substitutions per site.
#'
#' If this information is not know, it will be reconstructed within \code{treeWAS} using Fitch's parsimony.
#'
#' If this information is known (i.e., it has been estimated elsewhere through a parsimonious reconstruction),
#' it can be provided in the \code{n.subs} argument.
#' It must be in the form of a \emph{named} vector (or table), or a vector in which the \emph{i}'th element
#' contains the number of \emph{loci} that have been estimated to undergo \emph{i} substitutions on the tree.
#' The vector must be of length \emph{max n.subs}, and "empty" indices must contain zeros.
#' For example: the vector \code{n.subs = c(1833, 642, 17, 6, 1, 0, 0, 1)}, could be used to define the homoplasy distribution for a dataset with 2500 loci,
#' where the maximum number of substitutions to be undergone on the tree by any locus is 8, and no loci undergo either 6 or 7 substitutions.
#'
#'
###########################################################
#' \strong{Ancestral State Reconstrution}
#'
#' If ancestral state reconstruction has been performed outside of \code{treeWAS} for the \code{snps} and/or \code{phen} variable,
#' these reconstructions can be submitted in the form of a matrix and a vector, respectively, to the
#' \code{snps.reconstruction} and \code{phen.reconstruction} arguments. Please note the formatting requirements.
#'
#' If provided by the user, \code{snps.reconstruction} should contain \code{snps} in its first \code{nrow(snps)} rows,
#' and then have the reconstructed states in rows \code{nrow(snps)+1} to \code{nrow(snps)+tree$Nnode}
#'
#' If provided by the user, \code{phen.reconstruction} should contain \code{phen} in its first \code{length(phen)} elements,
#' and then have the reconstructed states in elements \code{length(phen)+1} to \code{length(phen)+tree$Nnode}
#'
#' If created externally, the \code{snps.reconstruction} must have been generated through either parsimony or ML, and the
#' \code{snps.sim.reconstruction} argument should be set to match (as either \code{"parsimony"} or \code{"ML"}), so that a direct comparison can be made.
#' At least for small datasets, it may be worth (re-)running this reconstruction within \code{treeWAS} instead, in case any inconsistencies exist between the
#' external and internal methods of reconstruction.
#'
#'
###########################################################
#' \strong{Tests of Association}
#'
#' \emph{Notation used in the equations below:}
#'            \deqn{G = Genotypic state...}
#'            \deqn{P = Phenotypic state...}
#'            \deqn{t = ... at terminal nodes}
#'            \deqn{a = ... at ancestral nodes}
#'            \deqn{d = ... at descendant nodes}
#'            \deqn{Nterm = Number of terminal nodes}
#'
#' \describe{
#' \item{\code{terminal}}{The \code{terminal} test solves the following equation, for each genetic locus, at the terminal nodes of the tree only:
#'                        \deqn{Terminal = | (1/Nterm)*(Pt*Gt - (1 - Pt)*Gt - Pt*(1 - Gt) + (1 - Pt)*(1 - Gt)) |}
#'                        The \code{terminal} test is a sample-wide test of association that
#'                        seeks to identify broad patterns of correlation between genetic loci and the phenotype,
#'                        without relying on inferences drawn from reconstructions of the ancestral states.}
#'
#' \item{\code{simultaneous}}{The \code{simultaneous} test solves the following equation, for each genetic locus, across each branch in the tree:
#'                            \deqn{Simultaneous = | (Pa - Pd)*(Ga - Gd) |}
#'                            This allows for the identification of simultaneous substitutions (i.e., substitutions occuring in
#'                            both the genetic locus and phenotypic variable on the same branch of the phylogenetic tree,
#'                            or parallel change on the same branch if non-binary data). Simultaneous substitutions are an
#'                            indicator of a deterministic relationship between genotype and phenotype. Moreover, because
#'                            this score measures \emph{only} simultaneous change, and is not negatively impacted by
#'                            the lack of association on other branches, it may be able to detect associations occurring
#'                            within some clades but not others and, therefore, to identify
#'                            loci giving rise to the phenotype through complementary pathways.}
#'
#' \item{\code{subsequent}}{The \code{subsequent} test solves the following equation, for each genetic locus, across each branch in the tree:
#'                            \deqn{Subsequent = | 4/3(Pa*Ga) + 2/3(Pa*Gd) + 2/3(Pd*Ga) + 4/3(Pd*Gd) - Pa - Pd - Ga - Gd + 1 |}
#'                            Calculating this metric across all branches of the tree allows us to measure in what
#'                            proportion of tree branches we expect the genotype and phenotype to be in the same state.}
#'
#'    }
#'
#'
###########################################################
#'
#' @return A list is returned as the output of \code{treeWAS} containing the data used within it,
#' and results including information about significant loci identified, if any.
#'
#'
#'
###########################################################
#'
#' @examples
#'
#' \dontrun{
#' ## load example homoplasy distribution
#' data(dist_0)
#' str(dist_0)
#'
#'
#' ## simulate a tree, phenotype, and genetic
#' ## data matrix with 10 associated loci:
#' dat <- coalescent.sim(n.ind = 100,
#'                         n.snps = 1000,
#'                         n.subs = dist_0,
#'                         n.snps.assoc = 10,
#'                         assoc.prob = 90,
#'                         n.phen.subs = 15,
#'                         phen = NULL,
#'                         plot = TRUE,
#'                         heatmap = FALSE,
#'                         reconstruct = FALSE,
#'                         dist.dna.model = "JC69",
#'                         grp.min = 0.25,
#'                         row.names = NULL,
#'                         coaltree = TRUE,
#'                         s = NULL,
#'                         af = NULL,
#'                         filename = NULL,
#'                         set = 1,
#'                         seed = 1)
#'
#' ## isolate elements of output:
#' snps <- dat$snps
#' phen <- dat$phen
#' snps.assoc <- dat$snps.assoc
#' tree <- dat$tree
#'
#' ## run treeWAS:
#' out <- treeWAS(snps = snps,
#'                 phen = phen,
#'                 tree = tree,
#'                 n.subs = dist_0,
#'                 sim.n.snps = ncol(snps)*10,
#'                 test = c("terminal", "simultaneous", "subsequent"),
#'                 snps.reconstruction = "parsimony",
#'                 snps.sim.reconstruction = "parsimony",
#'                 phen.reconstruction = "parsimony",
#'                 p.value = 0.01,
#'                 p.value.correct = "bonf",
#'                 p.value.by = "count",
#'                 dist.dna.model = NULL,
#'                 plot.tree = FALSE,
#'                 plot.manhattan = TRUE,
#'                 plot.null.dist = TRUE,
#'                 plot.dist = FALSE,
#'                 snps.assoc = NULL,
#'                 filename.plot = NULL,
#'                 seed = NULL)
#'
#' ## examine output:
#' print(out)
#'
#' }
#'
#'
###########################################################
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @import adegenet ape phangorn
#' @importFrom Hmisc all.is.numeric
#'
#' @export

###################################################################################################################################







# \deqn{Terminal = |\sum_{i = 1}^{n_{term}}\frac{1}{n_{term}} (p_{i}^{des}s_{i}^{des}\ -\ (1 - p_{i}^{des})s_{i}^{des}\ -\ p_{i}^{des}(1 - s_{i}^{des})\ +\ (1 - p_{i}^{des})(1 - s_{i}^{des}))\ |}

#################
## PARAMETERS: ##
#################

# n.subs <- NULL
# dist.dna.model <- "JC69"
# plot.tree <- FALSE
# test <- c("terminal", "simultaneous", "subsequent")
# p.value <- 0.001
# p.value.correct <- "fdr"
# p.value.by <- "count"
# sim.n.snps <- ncol(snps)*10
# n.reps <- 1
# plot.manhattan <- TRUE
# plot.null.dist <- TRUE
# plot.dist <- FALSE
# snps.reconstruction <- "parsimony"
# snps.sim.reconstruction <- "parsimony"
# phen.reconstruction <- "parsimony"


##############
## EXAMPLE: ##
##############

# foo  <- treeWAS(snps,
#                 phen,
#                 tree = tree,
#                 n.subs = NULL,
#                 sim.n.snps = ncol(snps)*10,
#                 test = c("terminal", "simultaneous", "subsequent"),
#                 p.value = 0.01,
#                 p.value.correct = "bonf",
#                 p.value.by = "count",
#                 dist.dna.model = "JC69",
#                 plot.tree = FALSE,
#                 plot.manhattan = TRUE,
#                 plot.null.dist = TRUE,
#                 plot.dist = FALSE,
#                 snps.assoc = NULL, # for (manhattan) plot
#                 snps.reconstruction = "parsimony",
#                 phen.reconstruction = "parsimony",
#                 filename.plot = NULL,
#                 seed = NULL)



treeWAS <- function(snps,
                    phen,
                    tree = c("UPGMA", "nj", "ml"),
                    n.subs = NULL,
                    sim.n.snps = ncol(snps)*10,
                    test = c("terminal", "simultaneous", "subsequent"),
                    snps.reconstruction = "parsimony",
                    snps.sim.reconstruction = "parsimony",
                    phen.reconstruction = "parsimony",
                    na.rm = TRUE,
                    p.value = 0.01,
                    p.value.correct = c("bonf", "fdr", FALSE),
                    p.value.by = c("count", "density"),
                    dist.dna.model = "JC69",
                    plot.tree = FALSE,
                    plot.manhattan = TRUE,
                    plot.null.dist = TRUE,
                    plot.dist = FALSE,
                    snps.assoc = NULL, # for (manhattan) plot
                    filename.plot = NULL,
                    seed = NULL){

  ###################
  ## LOAD PACKAGES ##
  ###################
  # require(adegenet)
  # require(phangorn)
  # require(ape)
  # # require(ade4) #?
  # require(Hmisc) # all.is.numeric

  snps.sim <- snps.rec <- snps.REC <- snps.sim.rec <- snps.sim.REC <- phen.rec <- phen.REC <- NULL

  #####################################################################
  ## 0) HANDLE INPUT DATA #############################################
  #####################################################################

  #####################
  ## HANDLE TEST ARG ##
  #####################

  ## Allow partial matching of argument names:
  test <- match.arg(arg = test,
                    choices = c("terminal", "simultaneous", "subsequent", "cor", "fisher"),
                    several.ok = TRUE)

  ########################
  ## HANDLE SNPS & PHEN ##
  ########################
  if(!is.matrix(snps)) snps <- as.matrix(snps)
  # x <- snps
  n.snps <- ncol(snps)

  ## convert phenotype to factor
  phen <- as.factor(phen)
  # y <- phen

  ## set n.ind:
  n.ind <- length(phen)
  inds <- c(1:n.ind)

  #################
  ## HANDLE TREE ##
  #################

  ## RECONSTRUCTED TREE ##

  if(class(tree) == "character"){

    tree <- tolower(tree)

    if(!any(c("upgma", "nj", "ml") %in% tree)){
      warning("If tree is not a phylo object,
              it should be one of 'UPGMA', 'NJ', 'ML',
              specifying which method is to be used to
              reconstruct the phylogenetic tree from the snps provided.
              Choosing 'UPGMA' by default.")
      tree <- "upgma"
    }
    tree <- tree.reconstruct(snps,
                             method = tree,
                             dist.dna.model = dist.dna.model,
                             plot = plot.tree)

    ## HANDLE TREE: ##
    ## Always work with trees in "pruningwise" order:
    tree <- reorder.phylo(tree, order="pruningwise")
    ## Trees must be rooted:
    if(!is.rooted(tree)) tree <- midpoint(tree)

  }else{

    ## USER-PROVIDED TREE ##

    ## If user has already submitted a tree as input:
    ## Work with a centered phylo tree for
    ## consistency and visualisation's sake:
    if(class(tree) != "phylo") tree <- as.phylo(tree)

    ## HANDLE TREE: ##
    ## Always work with trees in "pruningwise" order:
    tree <- reorder.phylo(tree, order="pruningwise")
    ## Trees must be rooted:
    if(!is.rooted(tree)) tree <- midpoint(tree)


    if(plot.tree==TRUE){
      plot(tree)
      title("Phylogenetic tree (original)")
      axisPhylo()
    } # end plot.tree

  }# end tree...


  ####################################################################
  ## Check for COALESCENT or RTREE-TYPE ORDERING before SIMULATING: ##
  ####################################################################

  # if(coaltree == FALSE){
  ## Simulation should start from the lowest internal node index (ie n.terminal+1):
  if(unique(tree$edge[,1])[1] == (tree$Nnode+2)){
    ## Simulate from top:bottom?
    x <- 1:nrow(tree$edge)
  }else{
    ## Extra check:
    if(unique(tree$edge[,1])[length(unique(tree$edge[,1]))] == (tree$Nnode+2)){
      ## Simulate from bottom:top?
      x <- rev(c(1:nrow(tree$edge)))
    }else{
      stop("This simulation procedure expects to find the root node/first internal node
           (ie. n.terminal+1) in either the FIRST or LAST row of tree$edge[,1],
           once the tree has been reordered to be in 'pruningwise' configuration.
           This is NOT the case with your tree. Please check.")
    }
    }
  ####################################################################




  ####################################################################
  ########################
  ## RUN CHECKS ON DATA ##
  ########################
  snps.ori <- snps
  snps.reconstruction.ori <- snps.reconstruction
  phen.ori <- phen
  tree.ori <- tree
  ####################################################################
  ## CHECK TO ENSURE ALL CONTAIN SAME SET OF LABELS:
  ## Check that labels are present:
  if(is.null(tree$tip.label)) stop("Trees must have tip.labels corresponding to rownames(snps).")
  if(is.null(rownames(snps))) stop("SNPs must have rownames corresponding to tree$tip.label.")
  if(is.null(names(phen))) stop("Phen must have names corresponding to tree$tip.label.")

  ## Cross-check labels with each other:
  if(!all(tree$tip.label %in% rownames(snps))) stop("Some elements of tree$tip.label
                                                    are absent from rownames(snps).")
  if(!all(rownames(snps) %in% tree$tip.label)) stop("Some elements of rownames(snps)
                                                    are absent from tree$tip.label.")
  if(!all(tree$tip.label %in% names(phen))) stop("Some elements of tree$tip.label
                                                    are absent from names(phen).")
  if(!all(names(phen) %in% tree$tip.label)) stop("Some elements of names(phen)
                                                    are absent from tree$tip.label.")
  if(!all(names(phen) %in% rownames(snps))) stop("Some elements of names(phen)
                                                    are absent from rownames(snps).")
  if(!all(rownames(snps) %in% names(phen))) stop("Some elements of rownames(snps)
                                                    are absent from names(phen).")
  ####################################################################
  ################
  ## CHECK PHEN ##
  ################
  ####################################################################
  ## CHECK IF ANY PHEN MISSING:
  if(any(is.na(phen))){
    toRemove <- names(which(is.na(phen)))
    warning(c("The phenotypic variable for individual(s) ", toRemove, " is missing.
              Removing these individuals from the analysis."))
    ## remove individuals from phen:
    phen <- phen[-which(names(phen) %in% toRemove)]
    ## remove individuals from snps:
    snps <- snps[-which(rownames(snps) %in% toRemove), ]
    ## remove individuals from tree:
    tree <- drop.tip(tree, tip = toRemove)
  }

  ####################################################################
  #################
  ## CHECK SNPS: ##
  #################
  ####################################################################
  ## CHECK IF ALL POLYMORPHIC: (not necessary...)
  #   tab <- sapply(c(1:ncol(snps)), function(e) min(table(snps[,e]))/nrow(snps))
  #
  #   ## TO DO (?): COULD REPLACE 0.01 w a POLYTHRESH ARGUMENT
  #   if(length(which(tab < 0.01)) > 0){
  #     toRemove <- which(tab < 0.01)
  #     if(length(toRemove) > 0){
  #       snps <- snps[, -toRemove]
  #     }
  #     ## (+ snps.reconstruction)
  #     if(is.matrix(snps.reconstruction)){
  #       snps.reconstruction <- snps.reconstruction[, -toRemove]
  #     }
  #   }

  ####################################################################
  ## CHECK IF ALL CONTAIN MINORITY OF NAs:
  if(is.null(na.rm)) na.rm <- TRUE
  if(na.rm != FALSE){
  NA.tab <- sapply(c(1:ncol(snps)), function(e) length(which(is.na(snps[,e]))))

  ## What proportion of individuals w NAs is acceptable?
  toRemove <- which(NA.tab > floor(nrow(snps)/2)) # Remove columns w > 50% of NAs
  if(length(toRemove) > 0){
    snps <- snps[, -toRemove]

    ## (+ snps.reconstruction)
    if(is.matrix(snps.reconstruction)){
      snps.reconstruction <- snps.reconstruction[, -toRemove]
    }
  }
  }
  ####################################################################
  ## CHECK FOR NON-BINARY SNPS (?): ## DO THIS OUTSIDE OF THE treeWAS PIPELINE -- TOO MANY OPPORTUNITIES FOR ERRORS IF DONE AUTOMATICALLY.
  ## (This will only work if all snps colnames end in one of .a/.c/.g/.t)
  #   snps <- get.binary.snps(snps)
  #   ## (+ snps.reconstruction)
  #   if(is.matrix(snps.reconstruction)){
  #     snps.reconstruction <- get.binary.snps(snps.reconstruction)
  #   }
  ####################################################################

  ## REORDER SNPS TO MATCH TREE$TIP.LABEL
  if(!identical(rownames(snps), tree$tip.label)){
    ord <- match(tree$tip.label, rownames(snps))
    snps <- snps[ord,]
    ## check:
    if(!identical(rownames(snps), tree$tip.label)){
      stop("Unable to rearrange snps such that rownames(snps)
            match content and order of tree$tip.label.
            Please do this by hand.")
    }
  }

  ## REORDER PHEN TO MATCH TREE$TIP.LABEL
  if(!identical(names(phen), tree$tip.label)){
    ord <- match(tree$tip.label, names(phen))
    phen <- phen[ord]
    ## check:
    if(!identical(names(phen), tree$tip.label)){
      stop("Unable to rearrange phen such that names(phen)
           match content and order of tree$tip.label.
           Please do this by hand.")
    }
    }

  ####################################################################





  ###################
  ## HANDLE N.SUBS ##
  ###################

  if(is.null(n.subs)){

    ## if n.subs is a vector (ie. distribution) ##
    ## we use this distribution directly (but in proportion with the number of sites)
    ## to specify the n.subs per site. (Handled within snp.sim fn.)

    ## if n.subs is NULL ##
    ## we compute the distribution of the n.subs-per-site
    ## using the Fitch parsimony score calculation fns from phangorn.


    ###########
    ## TO DO ##   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
    ###########

    ## if either test 2 or test 3 will be run with parsimonious/user-provided reconstruction,
    ## get n.subs from this to avoid duplication..?
    #     if(any(c("simultaneous", "subsequent") %in% test) & snps.reconstruction != "ml"){
    #
    #       ## run get.ancestral.pars
    #       snps.pars <- get.ancestral.pars(var=snps, tree=tree)
    #
    #       ## get elements of output
    #       snps.rec <- snps.pars$var.rec
    #       snps.subs.edges <- snps.pars$subs.edges
    #
    #       ## CHECK--Compare costs:
    #       cost1 <- get.fitch.n.mts(snps=snps, tree=tree)
    #       cost2 <- sapply(c(1:length(snps.subs.edges)), function(e) length(snps.subs.edges[[e]][["total"]]))
    #       table(cost1)
    #       table(cost2) ## longer tail...
    #
    #     }else{

    ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

    ## get parsimomy cost for each SNP locus using fitch:
    n.subs <- get.fitch.n.mts(snps=snps, tree=tree)
    n.subs <- table(n.subs)
    ## handle n.subs "levels" with 0 SNP loci at those levels:
    noms <- as.numeric(names(n.subs))
    temp <- rep(0, max(noms))
    for(i in 1:max(noms)){
      if(i %in% noms) temp[i] <- n.subs[which(noms==i)]
    }
    n.subs <- temp
    # names(n.subs) <- 1:length(n.subs)
    ## check?
    # sum(n.subs) == ncol(snps) ## should be TRUE
    # }
  }

  ## check:
  # barplot(n.subs, col=transp("blue", 0.5), names=c(1:length(n.subs)))
  # title("Homoplasy distribution (snps)")

  #####################################################
  ## 1) Simulate multiple snps datasets to compare your
  ## real correlations w phen to  #####################
  #####################################################

  ## check sim.n.snps:
  if(is.null(sim.n.snps)){
    sim.n.snps <- ncol(snps)*10
  }else{
    ## Update sim.n.snps to match the REAL ncol(snps) after checks...
    if(n.snps != ncol(snps)){
      if(sim.n.snps == n.snps){
        sim.n.snps.ori <- sim.n.snps
        sim.n.snps <- ncol(snps)
        warning(paste("Note: Updating sim.n.snps to match the number of real loci after data cleaning.
                  Input:", sim.n.snps.ori, " -->
                      Updated:", sim.n.snps))
      }else{
        Nx <- c(2:10)
        if(any(n.snps*Nx %in% sim.n.snps)){
          sim.n.snps.ori <- sim.n.snps
          Nx <- Nx[which(Nx == (sim.n.snps/n.snps))]
          sim.n.snps <- ncol(snps)*Nx
          warning(paste("Note: Updating sim.n.snps to match ", Nx, "x the number of real loci after data cleaning.
                  Input: ", sim.n.snps.ori, " -->
                  Updated: ", sim.n.snps, sep=""))
          ## TO DO: ##
          ## If the user, for whatever reason, wanted to simulate the number they requested
          ## (ie. matching the input n.snps (coincidentally?)), could either
          ## (A) Add an argument like upate.sim.n.snps = FALSE,
          ## (B) Tell them to do their own data cleaning exactly as I do but outside of/before running the treeWAS function.
          ## (C) Fudge--tell them to add 1, eg. 10,000 --> 10,001?
        }
      }
    }
  }
  out <- genomes <- snps.mat <- list()

  ## SIMULATE A DATASET | your tree ##
  if(!is.null(seed)) set.seed(seed)
  out[[1]] <- snp.sim(n.snps = sim.n.snps,
                      n.subs = n.subs,
                      n.snps.assoc = 0,
                      assoc.prob = 100,
                      tree = tree,
                      phen.loci = NULL,
                      heatmap = FALSE,
                      reconstruct = FALSE,
                      dist.dna.model = dist.dna.model,
                      row.names = rownames(snps),
                      seed = seed)

  genomes[[1]] <- out[[1]][[1]]

  ## Modify genomes/snps matrices
  if(!is.null(genomes[[1]])){
    snps.mat[[1]] <- genomes[[1]]
  }else{
    snps.mat[[1]] <- NULL
  }

  gc()


  print("treeWAS snps sim done.")


  ################################################################
  ## 3) Get results:##############################################
  #### Determine the phylogenetially correct p-values for SNPs | #
  ##   null distributions of correlations from simulated data ####
  #### Synthesize results output: List of all significant SNPs, ##
  ##   their names/locations, their p-values for this phenotype ##
  ################################################################

  ##################################################
  ## RUN CHECKS ONCE BEFORE get.sig.snps FOR LOOP ##
  ##################################################

  ## NOTE: These checks are repeated within the get.sig.snps fn
  ## as an extra layer of safety/ in case users want to use it alone,
  ## but it is more economical to run them once outside of the for loop..

  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

  snps.unique <- snps.index <- snps.sim.unique <- snps.sim.index <- NULL

  #################
  ## Handle snps ##
  #################
  ## Check snps column names
  if(is.null(colnames(snps))) colnames(snps) <- c(1:ncol(snps))

  ################################
  ## Handle snps.sim --> matrix ##
  ################################
  snps.sim <- snps.mat

  ## Handle matrix/list input:
  if(class(snps.sim) == "list"){
    ## If list of length 1...
    if(length(snps.sim) == 1){
      ## keep matrix:
      snps.sim <- snps.sim[[1]]
    }else{
      ## If list of multiple matrices...
      ## merge all elements into one big matrix
      ## by pasting columns together:
      snps.sim <- do.call("cbind", snps.sim)
    }
  }

  ## check n.subs:
  # n.subs.sim <- get.fitch.n.mts(snps.sim, tree=tree)
  # n.subs.sim2 <- table(n.subs.sim)
  # barplot(n.subs.sim2, col=transp("blue", 0.5), names=c(1:length(n.subs.sim2)))
  # title("Homoplasy distribution \n (snps.sim)")

  #################
  ## Handle phen ##
  #################
  ## convert phenotype to numeric:
  ## NOTE--this is also necessary for returning results in step (5)!
  phen.ori <- phen
  if(!is.numeric(phen)) phen <- as.numeric(phen)
  ## for ease of interpretation,
  ## if phen has 2 levels, 1 and 2,
  ## make these 0 and 1:
  if(length(unique(phen))!=2){
    stop("This function is only designed for phenotypes with two levels.")
  }else{
    if(length(phen[-c(which(phen==1), which(phen==2))])==0){
      phen <- replace(phen, which(phen==1), 0)
      phen <- replace(phen, which(phen==2), 1)
    }
  }
  ## ensure ind names not lost
  names(phen) <- names(phen.ori)

  ##############################################################################################
  ## Reconstruct ancestral SNPs & phen by parsimony/ML (for tests simultaneous & subsequent) ##
  ##############################################################################################

  ## Ensure we are only reconstructing ancestral states ONCE here, to be used in MULTIPLE tests later.
  snps.REC <- snps.sim.REC <- phen.REC <- NULL

  if(any(c("simultaneous", "subsequent") %in% test)){

    #######################
    ## Reconstruct SNPs: ##
    #######################

    ## By PARSIMONY or ML: ##
    ############################
    ## Reconstruct REAL SNPs: ##
    ############################
    ## If user-provided reconsruction:
    if(is.matrix(snps.reconstruction)){
      ## CHECK:
      if(nrow(snps.reconstruction) != (nrow(snps)+tree$Nnode)){
        warning("The number of rows in the provided snps.reconstruction is not equal to the
                                                            total number of nodes in the tree. Performing a new parsimonious reconstruction instead.")
        snps.reconstruction <- "parsimony"
      }
      if(ncol(snps.reconstruction) != ncol(snps)){
        warning("The number of columns in the provided snps.reconstruction is not equal to the number of
                                                            columns in the snps matrix. Performing a new parsimonious reconstruction instead.")
        snps.reconstruction <- "parsimony"
      }
    }

    ## If not user-provided or checks failed, reconstruct ancestral states:
    if(is.matrix(snps.reconstruction)){
      snps.rec <- snps.reconstruction
    }else{
    # system.time( # 274
      snps.REC <- asr(var = snps, tree = tree, type = snps.reconstruction)
    # )
    snps.rec <- snps.REC$var.rec
    }

    #################################
    ## Reconstruct SIMULATED SNPs: ##
    #################################
    # system.time(
      snps.sim.REC <- asr(var = snps.sim, tree = tree, type = snps.sim.reconstruction)
    # )
    snps.sim.rec <- snps.sim.REC$var.rec


    #######################
    ## Reconstruct phen: ##
    #######################
    ## If user-provided reconsruction:
    if(length(phen.reconstruction) > 1){
      ## CHECK:
      if(length(phen.reconstruction) != (length(phen)+(tree$Nnode))){
        warning("The number of individuals in the provided phen.reconstruction is not equal to the
                total number of nodes in the tree. Performing a new parsimonious reconstruction instead.")
        phen.reconstruction <- "parsimony"
      }
    }

    ## If not user-provided or checks failed, reconstruct ancestral states:
    if(length(phen.reconstruction) > 1){
      phen.rec <- phen.reconstruction
    }else{
      ## By PARSIMONY or ML: ##
      phen.REC <- asr(var = phen, tree = tree, type = phen.reconstruction)
      phen.rec <- phen.REC$var.rec
    }

  } # end reconstruction for tests 2 & 3


  ## !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ###


  ###########################
  ## GET UNIQUE SNPS(.SIM) ##
  ###########################

  ## Get UNIQUE snps + index
  snps.complete <- snps
  temp <- get.unique.matrix(snps, MARGIN=2)
  snps.unique <- temp$unique.data
  snps.index <- temp$index

  ## Get UNIQUE snps.sim + index
  snps.sim.complete <- snps.sim
  temp <- get.unique.matrix(snps.sim, MARGIN=2)
  snps.sim.unique <- temp$unique.data
  snps.sim.index <- temp$index

  ## Get UNIQUE snps.reconstruction
  snps.rec.complete <- snps.rec
  temp <- get.unique.matrix(snps.rec, MARGIN=2)
  snps.rec <- temp$unique.data
  snps.rec.index <- temp$index
  if(!identical(snps.rec.index, snps.index)){
    warning("Careful-- snps and snps.rec should have the same index when reduced
            to their unique forms!") ## SHOULD THIS BE A "STOP" INSTEAD? OR IS THIS ERROR NOT FATAL OR NOT POSSIBLE????
  }

  ## Get UNIQUE snps.sim.reconstruction
  snps.sim.rec.complete <- snps.sim.rec
  temp <- get.unique.matrix(snps.sim.rec, MARGIN=2)
  snps.sim.rec <- temp$unique.data
  snps.sim.rec.index <- temp$index
  if(!identical(snps.sim.rec.index, snps.sim.index)){
    warning("Careful-- snps.sim and snps.sim.rec should have the same index when reduced
            to their unique forms!") ## SHOULD THIS BE A "STOP" INSTEAD? OR IS THIS ERROR NOT FATAL OR NOT POSSIBLE????
  }

  print("reconstructions done")

  #######################
  ## identify sig.snps ##
  #######################
  ## Note: UNIQUE snps & snps.sim are identified WITHIN the get.sig.snps fn
  ## to reduce computational time, but results are identified on the basis of all
  ## ORIGINAL snps & snps.sim columns inputted.

  sig.list <- list()

  # test <- c("terminal", "simultaneous", "subsequent")
  TEST <- as.list(test)

  ## Run get.sig.snps fn once for each association test:
  system.time( # 100 - 164 (why such a difference?)
    for(i in 1:length(TEST)){
      sig.list[[i]] <- get.sig.snps(snps = snps,
                                    snps.unique = snps.unique,
                                    snps.index = snps.index,
                                    snps.sim = snps.sim,
                                    snps.sim.unique = snps.sim.unique,
                                    snps.sim.index = snps.sim.index,
                                    phen = phen,
                                    tree = tree,
                                    test = TEST[[i]],
                                    n.tests = length(TEST),
                                    p.value = p.value,
                                    p.value.correct = p.value.correct,
                                    p.value.by = p.value.by,
                                    snps.reconstruction = snps.rec,
                                    snps.sim.reconstruction = snps.sim.rec,
                                    phen.reconstruction = phen.rec)
    }
  )

  names(sig.list) <- test
  # str(sig.list)

  print("get sig snps done.")

  ## DOUBLE CHECKING ##
  #   str(sig.list[[i]])
  #   sig.list[[i]]$sig.snps
  #   sig.list[[i]]$sig.corrs
  #   ## plot
  # hist(sig.list[[i]][[1]]$corr.sim)
  # hist(sig.list[[i]][[1]]$corr.dat)

  # sig.list[[2]][[1]]$corr.dat[snps.assoc]
  # sig.list[[3]][[1]]$corr.dat[snps.assoc]

  ## BUG CHECKING ##
  ## get.sig.snps
  #     snps <-  snps
  #     snps.sim <- snps.sim
  #     phen <- phen
  #     tree <- tree
  #     test <- "simultaneous"
  #     p.value <- p.value
  #     p.value.correct <- p.value.correct
  #     p.value.by <- p.value.by
  #     snps.reconstruction <- snps.rec
  #     snps.sim.reconstruction <- snps.sim.rec
  #     phen.reconstruction <- phen.rec

  #################
  ## GET RESULTS ##
  #################

  ## set margins for plotting:
  par.mar.ori <- par()$mar
  par(mar=c(5, 2, 4, 1)+0.1)

  RES <- list()

  ## get results for each test run:
  for(i in 1:length(sig.list)){

    #############################################
    ## isolate elements of get.sig.snps output ##
    #############################################

    corr.dat <- sig.list[[i]]$corr.dat
    corr.sim <- sig.list[[i]]$corr.sim
    p.vals <- sig.list[[i]]$p.vals
    sig.snps.names <- sig.list[[i]]$sig.snps.names
    sig.snps <- sig.list[[i]]$sig.snps
    sig.corrs <- sig.list[[i]]$sig.corrs
    sig.p.vals <- sig.list[[i]]$sig.p.vals
    min.p <- sig.list[[i]]$min.p
    sig.thresh <- sig.list[[i]]$sig.thresh

    ########################################

    ##########################
    ## NEW: MANHATTAN PLOT! ##
    ##########################
    if(plot.manhattan == TRUE){

      manhattan.plot(p.vals = corr.dat,
                     col = "funky",
                     transp = 0.75,
                     sig.thresh = sig.thresh,
                     thresh.col="red",
                     snps.assoc = NULL,
                     snps.assoc.col = "red",
                     jitter.amount = 0.00001,
                     min.p = NULL,
                     log10=FALSE,
                     ylab=paste(TEST[[i]], "score", sep=" "))
    } # end plot manhattan


    ##################################
    ## 4) (A) Plot the distribution ##
    ##################################

    ## Generate one histogram per test:
    plot.sig.snps(corr.dat = corr.dat,
                  corr.sim = corr.sim,
                  corr.sim.subset = NULL,
                  sig.corrs = corr.dat[sig.snps],
                  sig.snps = sig.snps.names,
                  sig.thresh = sig.thresh,
                  test = TEST[[i]],
                  sig.snps.col = "black",
                  hist.col = rgb(0,0,1,0.5),
                  hist.subset.col = rgb(1,0,0,0.5),
                  thresh.col = "red",
                  snps.assoc = NULL,
                  snps.assoc.col = "blue",
                  bg = "lightgrey",
                  grid = TRUE,
                  freq = FALSE,
                  plot.null.dist = TRUE,
                  plot.dist = FALSE)


    ########################################
    ## 5) Return results list ##############
    ########################################

    if(length(sig.snps)==0) sig.snps <- sig.corrs <- NULL

    ###########
    ## make a data.frame containing all relevant output for sig.snps
    if(length(sig.snps) > 0){

      ## Get counts for n.sig.snps in each cell of the contingency table:
      #     toKeep <- sapply(c(1:length(sig.snps)),
      #                      function(e)
      #                        which(dimnames(snps)[[2]] == sig.snps))
      toKeep <- sig.snps
      snps.toKeep <- snps[,toKeep]

      ##
      if(length(toKeep) > 1){
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
      }else{
        ## if only ONE sig snp (haploid) identified:
        S1P1 <- length(which(snps.toKeep[which(phen==1)]==1))
        S0P0 <- length(which(snps.toKeep[which(phen==0)]==0))
        S1P0 <- length(which(snps.toKeep[which(phen==0)]==1))
        S0P1 <- length(which(snps.toKeep[which(phen==1)]==0))

      }
      df <- data.frame(sig.snps,
                       sig.p.vals,
                       sig.corrs,
                       S1P1, S0P0, S1P0, S0P1)
      names(df) <- c("SNP.locus",
                     "p.value",
                     "Test.statistic",
                     "S1P1", "S0P0", "S1P0", "S0P1")

      ## NOTE: Could return sig.snps.names somewhere here
      ## in addition to sig.snps loci ####    ####    ####    ####

    }else{
      df <- "No significant SNPs found."
    }

    ## 0 p.vals
    #   min.p <- paste("p-values listed as 0 are <",
    #                  1/length(corr.sim), sep=" ")
    min.p <- 1/length(corr.sim)
    names(min.p) <- c("p-values listed as 0 are less than:")

    ## TO DO:
    ## ADD MANHATTAN PLOT

    results <- list()
    results[[1]] <- corr.dat
    results[[2]] <- corr.sim
    results[[3]] <- p.vals
    results[[4]] <- sig.thresh
    results[[5]] <- df
    results[[6]] <- min.p

    names(results) <- c("corr.dat",
                        "corr.sim",
                        "p.vals",
                        "sig.thresh",
                        "sig.snps",
                        "min.p.value")

    RES[[i]] <- results
  } # end for loop


  ## return plot margins to their original state:
  par(mar=par.mar.ori)


  ################################
  ## NEW! Get COMBINED results: ##
  ################################
  ## get uniques(sig.snps) for terminal, simultaneous, subsequent results combined
  ## for best thresh only (ie. pval.0.01.bonf.count.10.x.n.snps)
  SNP.loci <- vector("list", length=length(test))
  names(SNP.loci) <- test
  for(t in 1:length(test)){
    temp <- RES[[t]]$sig.snps
    if(is.vector(temp)){
      SNP.loci[[t]] <- NA
    }else{
      SNP.loci[[t]] <- rownames(temp) # temp$SNP.locus
    }
  }
  ## store 3 tests in separate list:
  SNP.loci.ori <- SNP.loci

  ## make list of length 2 (all, separately):
  SNP.loci <- vector("list", length=2)
  names(SNP.loci) <- c("treeWAS.combined", "treeWAS")
  if(length(as.vector(unlist(SNP.loci.ori))) > 0 & !all(is.na(as.vector(unlist(SNP.loci.ori))))){
    SNP.loci[[1]] <- sort(unique(as.vector(unlist(SNP.loci.ori))), decreasing=FALSE)
  }else{
    SNP.loci[[1]] <- NA
  }
  SNP.loci[[2]] <- SNP.loci.ori


  ## Assign treeWAS.combined to FIRST element of RES:
  RES.ORI <- RES
  RES <- list()
  RES[[1]] <- SNP.loci
  RES[(2:(length(RES.ORI)+1))] <- RES.ORI


  ## Also return data (in the form we were working with):
  dat <- list("snps" = snps,
              "snps.reconstruction" = snps.rec,
              "snps.sim" = snps.sim,
              "snps.sim.reconstruction" = snps.sim.rec,
              "phen" = phen,
              "phen.reconstruction" = phen.rec,
              "tree" = tree,
              "n.subs" = n.subs)

  ## Assign data to last element of RES:
  RES[[(length(RES)+1)]] <- dat

  ## name elements of RES:
  names(RES) <- c("treeWAS.combined", test, "dat")
  results <- RES

  class(results) <- "treeWAS"

  return(results)

} # end treeWAS



# corr.dat <- foo$terminal$corr.dat
# corr.sim <- foo$terminal$corr.sim
# hist(corr.sim, col=transp("red", 0.5), freq=F, xlim=c(0, ceiling(max(c(corr.dat, corr.sim)))))
# hist(corr.dat, col=transp("blue", 0.5), freq=F, add=T, xlim=c(0, ceiling(max(c(corr.dat, corr.sim)))))
#
# corr.dat <- foo$simultaneous$corr.dat
# corr.sim <- foo$simultaneous$corr.sim
# hist(corr.sim, col=transp("red", 0.5), freq=F, xlim=c(0, ceiling(max(c(corr.dat, corr.sim)))))
# hist(corr.dat, col=transp("blue", 0.5), freq=F, add=T, xlim=c(0, ceiling(max(c(corr.dat, corr.sim)))))
#
# corr.dat <- foo$subsequent$corr.dat
# corr.sim <- foo$subsequent$corr.sim
# hist(corr.sim, col=transp("red", 0.5), freq=F, xlim=c(0, ceiling(max(c(corr.dat, corr.sim)))))
# hist(corr.dat, col=transp("blue", 0.5), freq=F, add=T, xlim=c(0, ceiling(max(c(corr.dat, corr.sim)))))


##############################################################################################
## legend (temporary?)
## Not necessarily needed--usually only a couple UNIQUE thresholds...
#     par(mfrow=c(1,2))
#     par(oma = c(5,4,0,0) + 0.1)
#     par(mar = c(0,0,1,1) + 0.1)
#     ## column 1:
#     midpoints1 <- barplot(rep(10, length(THRESH)/2),
#                         col = seasun(length(THRESH))[1:(length(THRESH)/2)],
#                         horiz=TRUE)
#     ##overlay names:
#     text(3, midpoints1, labels=names(THRESH)[1:(length(THRESH)/2)], cex=0.75, adj=0.3)
#
#     ## column 2:
#     midpoints2 <- barplot(rep(10, length(THRESH)/2),
#                           col = seasun(length(THRESH))[((length(THRESH)/2)+1):length(THRESH)],
#                           horiz=TRUE)
#     ##overlay names:
#     text(3, midpoints2, labels=names(THRESH)[((length(THRESH)/2)+1):length(THRESH)], cex=0.75, adj=0.3)
#
#     ## return to original par settings:
#     par(mfrow=c(1,1))
#
#     ## add title
#     title("Legend: Significance Thresholds")
#
#     par(oma=c(0,0,0,0))
#     par(mar=c(5,4,4,2)+0.1)
##############################################################################################



# ## only 2 unique sets of p.vals for each test
# ## (for n.snps & 10x n.snps):
# pv <- list()
# summ <- list()
# for(j in 1:length(results)){
#   pv[[j]] <- list()
#   summ[[j]] <- list()
# for(i in 1:length(results[[j]])){
#   pv[[j]][[i]] <- results[[j]][[i]]$p.vals
#   summ[[j]][[i]] <- summary(pv[[j]][[i]])
# }
# }
#
# length(unique(summ[[2]]))
# length(unique(summ[[2]][seq(2, length(summ[[3]]), 2)]))
# length(unique(summ[[2]][seq(2, length(summ[[3]]), 2)]))




###############
## CHECK!!!! ##
###############
## w assoc.prob == 100, still getting low/0(!) scores for "snps.assoc"
## get original phen for all terminal + internal nodes and edges...

## checking simultaneous score:

# str(foo)
# phen.str <- foo$phen.plot.col
# phen1 <- foo$phen
# phen2 <- phen.str$all.nodes
#
# head(phen1, 20)
# head(phen2, 20)
#
# phen2[which(phen2 == "blue")] <- "A"
# phen2[which(phen2 == "red")] <- "B"
# phen2 <- as.factor(phen2)
# names(phen2) <- paste("ind", 1:length(phen2), sep=".")
#
# phen.edges <- phen.str$edges
# phen.edges[which(phen.edges == "blue")] <- "A"
# phen.edges[which(phen.edges == "red")] <- "B"
# # phen.edges[which(phen.edges == "green")] <- "C"
# phen.edges <- as.factor(phen.edges)
# names(phen.edges) <- paste("ind", 1:length(phen.edges), sep=".")
# head(phen.edges)
#
# ## which edges should/do contain phen subs?
# which(phen.edges == "green")
# ## run relevant code in reconstruct to get phen.subs.edges list (to see which edges are identified):
# phen.subs.edges$total
#
# snps.diffs <- list()
# snps.assoc.index <- index[snps.assoc] ## GAK! -- all snps.assoc have the same index (915)! set.seed problem, for a start...
# for(i in snps.assoc.index){
#   snps.diffs[[i]] <- get.branch.diffs(var = snps.rec[,i],
#                                       edges = edges)
# }
# which(abs(snps.diffs[[i]]) == 1)
# ## ALSO--important to NOTE that a lot of the reconstructed edges are only off by one...
# ## (thus should the subsequent score not be doing much better????)
# length(which(which(abs(snps.diffs[[i]]) == 1) %in% phen.subs.edges$total)) ## SO Shouldn't the max corr.dat score2 be 11????

###############
## CHECK!!!! ##
###############
## The FPR for Score 2 (simultaneous) should NOT be that high-- something wrong with the threshold selection??
#
# snps.assoc.ori.ori <- snps.assoc
#
# sim.dat <- sim.dat.ori <- results$dat$simultaneous
# corr.dat <- corr.dat.ori <- sim.dat$corr.dat
# corr.sim <- corr.sim.ori <- sim.dat$corr.sim
# p.vals <- p.vals.ori <- sim.dat$p.vals
#
# sim.res <- sim.res.ori <- results$res$simultaneous
# str(sim.res[[1]])
# str(sim.res[[32]])
# thresh.ori <- thresh <- sim.res[[1]]$sig.thresh
# thresh.ori <- thresh <- sim.res[[32]]$sig.thresh
#
# table(corr.sim)
#
# hist(corr.sim)
# lines(density(corr.sim), col="red", lwd=2)
#
# str(density(corr.sim))
#
# ## with first 10000 only??
# corr.sim.ori <- corr.sim
# corr.sim <- corr.sim.ori[1:10000]
# table(corr.sim)



# ###
#
# t.corr.sim <- results$dat$terminal$corr.sim
#
# ########################
# ## get density curve: ##
# ########################
#
# ## SIMULTANEOUS SCORE: ##
# dat <- corr.sim[1:10000]
# d <- density(dat)
# # from=0 may be necessary for aligning polygon w hsit alon x-axis, BUT may cause problems for polygon drawing (try lines instead?)
# xmax <- max(hist(dat)$mids)+min(hist(dat)$mids)
# ymax <- ceiling(max(d$y))
# hist(dat, freq=F, xlim=c(0,xmax), ylim=c(0,ymax))
# # lines(d, col="red", lwd=2, xlim=c(0,1), ylim=c(0,ymax))
# polygon(d, col=transp("red", 0.25), border="red", lwd=2, xlim=c(0,xmax), ylim=c(0,ymax))
#
# ## TERMINAL SCORE ##
# dat <- t.corr.sim
# d <- density(dat, from=0)
# ymax <- ceiling(max(d$y))
# hist(dat, freq=F, xlim=c(0,1), ylim=c(0,ymax))
# # lines(d, col="red", lwd=2, xlim=c(0,1), ylim=c(0,ymax))
# polygon(d, col=transp("red", 0.25), border="red", lwd=2, xlim=c(0,1), ylim=c(0,ymax))


###




###
