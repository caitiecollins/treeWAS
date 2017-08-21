

###################
## TO DO (2017): ##
###################


## dependencies problems?ssh caitlin@131.251.130.191

## Hmisc new version requiring a non-existant version of pkg "survival"..
## -->
## Soln:
# library(devtools)
# install_version("Hmisc", version = "3.9-1", repos = "http://cran.us.r-project.org")



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
#' @import adegenet ape

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
  l <- length(x$treeWAS.combined$treeWAS.combined)
  if(any(is.na(x$treeWAS.combined$treeWAS.combined))){
    l <- l - length(which(is.na(x$treeWAS.combined$treeWAS.combined)))
  }
  print(l)

  res <- x$treeWAS.combined$treeWAS.combined
  if(all.is.numeric(res)){
    res <- as.character(sort(as.numeric(res, decreasing=FALSE)))
  }else{
    if(all.is.numeric(removeLastN(res, 2))){
      ord <- match(sort(as.numeric(removeLastN(res, 2), decreasing=FALSE)), removeLastN(res, 2))
      res <- res[ord]
    }
  }
  if(l > 0){
    cat("Significant loci: \n")
    print(res)
  }

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

      print(nrow(res))

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
#' This function implements a phylogenetic approach to genome-wide association studies
#'(GWAS) designed for use in bacteria and viruses.
#' The \code{treeWAS} approach allows for the identification of
#' significant asociations between genotype and phenotype, while accounting for the
#' confounding effects of clonal population structure,
#' stratification (overlap between the population structure and phenotypic distribution),
#' and homologous recombination.
#'
#'
#' @param snps A matrix containing binary genetic data, with individuals in the rows
#'                and genetic loci in the columns and both rows and columns labelled.
#' @param phen A vector containing the phenotypic state of each individual, whose length is equal to the number of rows in
#'              \code{snps} and which is named with the same set of labels.
#'              The phenotype can be either binary (character or numeric) or continuous (numeric).
#' @param tree A \code{phylo} object containing the phylogenetic tree; or, a character string,
#'                one of \code{"NJ"}, \code{"BIONJ"} (the default), \code{"ML"}, or \code{"UPGMA"},
#'                or, if NAs are present in the distance matrix, one of: \code{"NJ*"} or \code{"BIONJ*"},
#'                specifying the method of phylogenetic reconstruction.
#' @param n.subs A numeric vector containing the homoplasy distribution (if known, see details), or NULL (the default).
#' @param n.snps.sim An integer specifying the number of loci to be simulated for estimating the null distribution
#'                      (by default \code{10*ncol(snps)}). If memory errors arise during the analyis of a large dataset,
#'                      it may be necesary to reduce \code{n.snps.sim} from a multiple of 10 to, for example, 5x the number of loci.
#' @param chunk.size An integer indicating the number of \code{snps} loci to be analysed at one time. This provides a solution for
#'                    machines with insufficient memory to analyse the dataset at hand.
#'                    Note that smaller values of \code{chunk.size} will increase the computational time required
#'                    (e.g., for \code{chunk.size = ncol(snps)/2}, treeWAS will take twice as long to complete).
#' @param test A character string or vector containing one or more of the following available tests of association:
#'              \code{"terminal"}, \code{"simultaneous"}, \code{"subsequent"}, \code{"cor"}, \code{"fisher"}.
#'              By default, the first three tests are run (see details).
#'
#' @param snps.reconstruction Either a character string specifying \code{"parsimony"} (the default) or \code{"ML"} (maximum likelihood)
#'                              for the ancestral state reconstruction of the genetic dataset,
#'                              or a matrix containing this reconstruction if it has been performed elsewhere
#'                              \emph{and} you provide the tree.
#'
#' @param snps.sim.reconstruction A character string specifying \code{"parsimony"} (the default) or \code{"ML"} (maximum likelihood)
#'                                  for the ancestral state reconstruction of the simulated null genetic dataset.
#'
#' @param phen.reconstruction Either a character string specifying \code{"parsimony"} (the default) or \code{"ML"} (maximum likelihood)
#'                              for the ancestral state reconstruction of the phenotypic variable,
#'                              or a vector containing this reconstruction if it has been performed elsewhere.
#' @param phen.type An optional character string specifying whether the ancestral state reconstruction
#'                  of the phenotypic variable, if performed via ML, should treat the phenotype
#'                  as either \code{"discrete"} or \code{"continuous"}. By default,
#'                  \code{phen.type} is \code{NULL}, in which case ML reconstructions will be "continuous" for any
#'                  non-binary phenotypes.
#' @param na.rm A logical indicating whether columns in \code{snps} containing more than 75\% \code{NA}s
#'                should be removed at the outset (TRUE, the default) or not (FALSE).
#'
#' @param p.value A number specifying the base p-value to be set the threshold of significance (by default, \code{0.01}).
#' @param p.value.correct A character string, either \code{"bonf"} (the default) or \code{"fdr"},
#'                          specifying whether correction for multiple testing
#'                          should be performed by Bonferonni correction (recommended) or the False Discovery Rate.
#' @param p.value.by A character string specifying how the upper tail of the p-value distribution is to be identified.
#'                      Either \code{"count"} (the default, recommended) for a simple count-based approach or
#'                      \code{"density"} for a kernel-density based approximation.
#' @param dist.dna.model A character string specifying the type of model to use in reconstructing the phylogenetic tree for
#'                          calculating the genetic distance between individual genomes,
#'                          only used if \code{tree} is a character string (see ?dist.dna).
#' @param plot.tree A logical indicating whether to generate a plot of the phylogenetic tree
#'                    (\code{TRUE}, the default) or not (\code{FALSE}).
#' @param plot.manhattan A logical indicating whether to generate a manhattan plot for each association score
#'                         (\code{TRUE}, the default) or not (\code{FALSE}).
#' @param plot.null.dist A logical indicating whether to plot the null distribution of association score statistics
#'                        (\code{TRUE}, the default) or not (\code{FALSE}).
#' @param plot.dist A logical indicating whether to plot the true distribution of association score statistics
#'                    (\code{TRUE}) or not (\code{FALSE}, the default).
#' @param snps.assoc An optional character string or vector specifying known associated loci to be demarked in
#'                      results plots (e.g., from previous studies or if data is simulated); else NULL.
#' @param filename.plot An optional character string denoting the file location for
#'                        saving any plots produced; else \code{NULL}.
#' @param seed An optional integer to control the pseudo-randomisation process and allow
#'                for identical repeat runs of the function; else \code{NULL}.
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
#' Missing data is permitted (denoted by NA values only) in the genetic data matrix.
#' However, any row or column that is entirely missing will be automatically removed within \code{treeWAS}.
#' In addition, any column that is more than 75% NAs will be removed by default
#' (though, if you do not wish this to be the case, you can set \code{na.rm} to \code{FALSE}).
#'
#' The phylogenetic tree, if provided by the user, should contain only the terminal nodes corresponding to the individuals under analysis.
#' Any additional individuals, including the outgroup, if not under analysis should be removed prior to running \code{treeWAS}.
#' The tree can be either rooted or unrooted.
#'
#'
###########################################################
#' \strong{Homoplasy Distribution}
#'
#' The homoplasy distribution contains the number of substitutions per site.
#'
#' If this information is not known, it will be reconstructed within \code{treeWAS} using Fitch's parsimony.
#'
#' If this information is known (i.e., it has been estimated elsewhere through a parsimonious reconstruction),
#' it can be provided in the \code{n.subs} argument.
#' It must be in the form of a \emph{named} vector (or table), or a vector in which the \emph{i}'th element
#' contains the number of \emph{loci} that have been estimated to undergo \emph{i} substitutions on the tree.
#' The vector must be of length \emph{max n.subs}, and "empty" indices must contain zeros.
#' For example: the vector \code{n.subs = c(1833, 642, 17, 6, 1, 0, 0, 1)},
#' could be used to define the homoplasy distribution for a dataset with 2500 loci,
#' where the maximum number of substitutions to be undergone on the tree
#' by any locus is 8, and no loci undergo either 6 or 7 substitutions.
#'
#'
###########################################################
#' \strong{Ancestral State Reconstrution}
#'
#' If ancestral state reconstruction has been performed
#' outside of \code{treeWAS} for the \code{snps} and/or \code{phen} variable,
#' these reconstructions can be submitted in the form of a matrix and a vector, respectively, to the
#' \code{snps.reconstruction} and \code{phen.reconstruction} arguments.
#' Please note the formatting requirements.
#'
#' If provided by the user, \code{snps.reconstruction} should contain \code{snps} in its first \code{nrow(snps)} rows,
#' and then have the reconstructed states in rows \code{nrow(snps)+1} to \code{max(tree$edge)}
#'
#' If provided by the user, \code{phen.reconstruction} should contain \code{phen} in its first \code{length(phen)} elements,
#' and then have the reconstructed states in elements \code{length(phen)+1} to \code{max(tree$edge)}
#'
#' If created externally, the \code{snps.reconstruction} must have
#' been generated through either parsimony or ML, and the
#' \code{snps.sim.reconstruction} argument should be set to match
#' (as either \code{"parsimony"} or \code{"ML"}), so that a direct comparison can be made.
#' At least for small datasets, it may be worth (re-)running this
#' reconstruction within \code{treeWAS} instead, in case any inconsistencies exist between the
#' external and internal methods of reconstruction.
#'
#' In addition, if either \code{snps.reconstruction} or \code{phen.reconstruction} is being provided as input,
#' the user must ensure that \code{tree$node.label} contains labels for all internal nodes and
#' that this same set of names is used to label the rows or indices correspondng
#' to the internal nodes in \code{snps.reconstruction} and/or \code{phen.reconstruction}.
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
#' @return The output of `treeWAS` contains the set of significant loci identified
#' as well as all relevant information used by or generated within \code{treeWAS}.
#' To examine the significant findings alone, we recommend using the \code{print} function.
#'
#' The \code{treeWAS} function returns a list object, which takes on the following general structure:
#'
#' \describe{
#'     \item{\strong{$treeWAS.combined}}{
#'     The first element is a list of length two containing the identities of significant findings:
#'
#'     \itemize{
#'           \item{\code{$treeWAS.combined}
#'           The pooled set of significant loci identified by any association score.}
#'
#'           \item{\code{$treeWAS}
#'           A list with the sets of significant loci identified by each association score individually.}
#'           }
#'           }
#'
#'
#'     \item{\strong{$[SCORE]}}{
#'     There are list elements for each association score, containing the original score values
#'      for each locus and additional information for significant loci.
#'      By default, there will be three such \code{$[SCORE]}-type elements called
#'      \code{$terminal}, \code{$simultaneous}, and \code{$subsequent},
#'      each of which will have the following elements:
#'
#'      \itemize{
#'           \item{\code{$corr.dat}
#'           The association score values for loci in the empirical genetic dataset.}
#'
#'           \item{\code{$corr.sim}
#'           The association score values for loci in the simulated genetic dataset.}
#'
#'           \item{\code{$p.vals}
#'           The p-values associated with the loci in the empirical genetic dataset for this association score.}
#'
#'           \item{\code{$sig.thresh}
#'           The significance threshold for this association score.}
#'
#'           \item{\code{$sig.snps}
#'           A data frame describing the genetic loci identified as significant.
#'           The last four columns will only be present if the data is binary, in which case they will
#'           contain the cell counts of a 2x2 table of genotypic and phenotypic states for each significant locus.
#'           \itemize{
#'           \item{\code{row.names}: The column names of significant loci.}
#'           \item{\code{$SNP.locus}: The column positions of significant loci in \code{dat$snps} (see below).}
#'           \item{\code{$p.value}: The p-values for significant loci.}
#'           \item{\code{$score}: The association score values for significant loci.}
#'           \item{\code{$G1P1}: N.individuals with genotype = 1 and phenotype = 1 at this locus.}
#'           \item{\code{$G0P0}: N.individuals with genotype = 0 and phenotype = 0 at this locus.}
#'           \item{\code{$G1P0}: N.individuals with genotype = 1 and phenotype = 0 at this locus.}
#'           \item{\code{$G0P1}: N.individuals with genotype = 0 and phenotype = 1 at this locus.}
#'           }
#'           }
#'
#'           \item{\code{$min.p.value}
#'           The minimum p-value. P-values listed as zero can only truly be defined as below this value.}
#'           }
#'           }
#'
#'
#'     \item{\strong{$dat}}{
#'     The final element contains all of the data either used by or generated within \code{treeWAS}.
#'     Objects that were provided as inputs to the \code{treeWAS} function will be returned here
#'     in the form in which they were analysed (i.e., after data cleaning within \code{treeWAS}).
#'
#'     \itemize{
#'           \item{\code{$snps}
#'           The empirical genetic data matrix.}
#'
#'           \item{\code{$snps.reconstruction}
#'           The ancestral state reconstruction of the empirical genetic data matrix.}
#'
#'           \item{\code{$snps.sim}
#'           The simulated genetic data matrix.}
#'
#'           \item{\code{$snps.sim.reconstruction}
#'           The ancestral state reconstruction of the simulated genetic data matrix.}
#'
#'           \item{\code{$phen}
#'           The phenotypic variable.}
#'
#'           \item{\code{$phen.reconstruction}
#'           The ancestral state reconstruction of the phenotype.}
#'
#'           \item{\code{$tree}
#'           The phylogenetic tree.}
#'
#'           \item{\code{$n.subs}
#'           The homoplasy distribution. Each element represents a number of substitutions (from 1 to \code{length(n.subs)})
#'           and contains the number of loci that have been inferred to undergo that many substitutions.}
#'           }
#'           }
#'
#'           }
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
#'                 n.snps.sim = ncol(snps)*10,
#'                 test = c("terminal", "simultaneous", "subsequent"),
#'                 snps.reconstruction = "parsimony",
#'                 snps.sim.reconstruction = "parsimony",
#'                 phen.reconstruction = "parsimony",
#'                 phen.type = NULL,
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
#' @importFrom phangorn midpoint
#' @importFrom scales rescale
#'
#' @export

###################################################################################################################################
# @useDynLib phangorn, .registration = TRUE


treeWAS <- function(snps,
                    phen,
                    tree = c("BIONJ", "NJ", "UPGMA", "ML", "BIONJ*", "NJ*"),
                    n.subs = NULL,
                    n.snps.sim = ncol(snps)*10,
                    chunk.size = ncol(snps),
                    test = c("terminal", "simultaneous", "subsequent"),
                    snps.reconstruction = "parsimony",
                    snps.sim.reconstruction = "parsimony",
                    phen.reconstruction = "parsimony",
                    phen.type = NULL,
                    na.rm = TRUE,
                    p.value = 0.01,
                    p.value.correct = c("bonf", "fdr", FALSE),
                    p.value.by = c("count", "density"),
                    dist.dna.model = "JC69",
                    plot.tree = TRUE,
                    plot.manhattan = TRUE,
                    plot.null.dist = TRUE,
                    plot.dist = FALSE,
                    snps.assoc = NULL, # for (manhattan) plot
                    filename.plot = NULL,
                    seed = NULL){

  snps.sim <- snps.rec <- snps.REC <- snps.sim.rec <- snps.sim.REC <- phen.rec <- phen.REC <- NULL

  #####################################################################
  ## 0) HANDLE INPUT DATA #############################################
  #####################################################################

  #################
  ## HANDLE ARGS ##
  #################

  ## Allow partial matching of argument names:

  ## TEST ##
  test <- tolower(test)
  test <- match.arg(arg = test,
                    choices = c("terminal", "simultaneous", "subsequent", "cor", "fisher"),
                    several.ok = TRUE)

  ## TREE ##
  if(is.character(tree)){
    tree <- tolower(tree)
    if(class(try(match.arg(arg = tree,
                           choices = c("bionj", "nj", "upgma", "ml", "nj*", "bionj*"),
                           several.ok = FALSE), silent=TRUE)) == "try-error"){
      tree <- "bionj"
      cat("If tree is not a phylo object, please specify one of the following reconstruction methods:
          'UPGMA', 'NJ', 'BIONJ', ML', 'NJ*', 'BIONJ*'. Choosing 'BIONJ' by default.\n")
    }else{
      tree <- match.arg(arg = tree,
                        choices =  c("bionj", "nj", "upgma", "ml", "nj*", "bionj*"),
                        several.ok = FALSE)
    }
  }


  ## P-VALUE CORRECT ##
  p.value.correct <- tolower(p.value.correct)
  p.value.correct <- match.arg(arg = p.value.correct,
                               choices = c("bonf", "fdr", "false"),
                               several.ok = FALSE)
  if(p.value.correct == "false") p.value.correct <- FALSE


  ## P-VALUE BY ##
  p.value.by <- tolower(p.value.by)
  p.value.by <- match.arg(arg = p.value.by,
                          choices = c("count", "density"),
                          several.ok = FALSE)

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

    ## Check Rows first: (mandatory) ##
    NA.tab <- sapply(c(1:nrow(snps)), function(e) length(which(is.na(snps[e,]))))
    ## Remove ENTIRELY missing rows...
    toRemove <- rownames(snps)[which(NA.tab == ncol(snps))]
    ## If there are any...
    if(length(toRemove) > 0){
      ## print notice:
      cat("Removing", length(toRemove), "missing individual(s); only NAs in row(s).\n", sep=" ")
      ## remove individuals from snps:
      snps.ini <- snps
      snps <- snps[-which(rownames(snps) %in% toRemove), ]
      ## remove individuals from phen:
      phen <- phen[-which(names(phen) %in% toRemove)]
    } # end row NA check 1

    ## Ensure snps.rec is NOT user-provided if tree not provided:
    if(class(snps.reconstruction) != "character") snps.reconstruction <- "parsimony"

    ## Reconstruct tree: ##
    ## NOTE: Should really request WHOLE GENOMES here (add separate argument (seqs)) <------------- ##### TO DO #####
    tree <- tree.reconstruct(snps,
                             method = tree,
                             dist.dna.model = dist.dna.model,
                             plot = FALSE)

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

  } # end tree...


  ####################################################################
  ############################
  ## Get Anc-Des EDGE ORDER ##
  ############################
  ## Get sequence from lowest ("root", Nterm+1) to highest ancestral node:
  ix <- c(min(tree$edge[,1]):max(tree$edge[,1]))
  ## Get for loop index of rows in tree$edge[,1], in pairs, from lowest to highest:
  x <- as.vector(unlist(sapply(c(1:length(ix)), function(e) which(tree$edge[,1] == ix[e]))))
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
  ##############################################
  ## ASSIGN NODE LABELS to TREE (& SNPS.REC): ##
  ##############################################
  if(is.null(tree$node.label)){
    tree$node.label <- paste("NODE", c((length(tree$tip.label)+1):max(tree$edge)), sep="_")
  }

  ## get index for terminal nodes:
  ixt <- c(1:length(tree$tip.label))
  ## get index for internal nodes:
  ixi <- c((nrow(snps)+1):max(tree$edge[,2]))

  ## Check phen.rec names:
  if(length(phen.reconstruction) == max(tree$edge[,2])){
    if(!identical(names(phen.reconstruction)[ixi], tree$node.label)){
      ## rearrange internal nodes if possible:
      ord <- match(names(phen.reconstruction)[ixi], tree$node.label)
      if(length(which(is.na(ord))) == 0){
        phen.reconstruction[ixi] <- phen.reconstruction[ixi][ord]
        names(phen.reconstruction)[ixi] <- names(phen.reconstruction[ixi])[ord]
      }else{
        if(identical(names(phen.reconstruction)[ixt], tree$tip.label)){
          names(phen.reconstruction)[ixi] <- tree$node.label
          cat("Assuming phen.reconstruction[Nterminal+1:Ntotal] correspond to tree$node.label,
              although labels do not match.\n")
        }else{
          warning("The names of phen.reconstruction[Nterminal+1:Ntotal] do not correspond to tree$node.label.
                  Reconstructing phen internally instead.\n")
          phen.reconstruction <- "parsimony"
        }
      }
    }
  }

  ## Check snps.rec rownames:
  if(is.matrix(snps.reconstruction)){
    if(!identical(rownames(snps.reconstruction)[ixi], tree$node.label)){
      ## rearrange internal nodes if possible:
      ord <- match(rownames(snps.reconstruction)[ixi], tree$node.label)
      if(length(which(is.na(ord))) == 0){
        snps.reconstruction[ixi,] <- snps.reconstruction[ixi[ord],]
        rownames(snps.reconstruction)[ixi] <- rownames(snps.reconstruction)[ixi][ord]
      }else{
        if(identical(rownames(snps.reconstruction)[ixt], tree$tip.label)){
          rownames(snps.reconstruction)[ixi] <- tree$node.label
          cat("Assuming snps.reconstruction[Nterminal+1:Ntotal,] correspond to tree$node.label,
              although labels do not match.\n")
        }else{
          warning("The names of snps.reconstruction[Nterminal+1:Ntotal,] do not correspond to tree$node.label.
                  Reconstructing snps internally instead.\n")
          snps.reconstruction <- "parsimony"
        }
      }
    }
  }


  ####################################################################
  ################
  ## CHECK PHEN ##
  ################
  ####################################################################
  ## CHECK IF ANY PHEN MISSING:
  if(any(is.na(phen))){
    toRemove <- names(which(is.na(phen)))
    warning(c("The phenotypic variable for individual(s) ", toRemove, " is missing.
              Removing these individuals from the analysis.\n"))
    ## remove individuals from phen:
    phen <- phen[-which(names(phen) %in% toRemove)]
    ## remove individuals from snps:
    snps.ini <- snps # need snps.ini below
    snps <- snps[-which(rownames(snps) %in% toRemove), ]

    ## remove individuals from tree:
    tree <- drop.tip(tree, tip = toRemove)

    ## (+ snps.reconstruction)
    if(is.matrix(snps.reconstruction)){
      toKeep <- which(rownames(snps.reconstruction) %in% c(rownames(snps), tree$node.label))
      snps.reconstruction <- snps.reconstruction[toKeep, ]
    }
  }

  ####################################################################
  #################
  ## CHECK SNPS: ##
  #################
  ####################################################################
  ## CHECK IF BINARY:
  if(length(unique(as.vector(unlist(snps[!is.na(snps)])))) != 2){
    stop("snps must be a binary matrix")
  }else{
    ## Convert to numeric, binary:
    if(!is.numeric(snps)){
      na.before <- length(which(is.na(snps)))
      if(!all.is.numeric(snps)){
        r.noms <- rownames(snps)
        c.noms <- colnames(snps)
        snps <- matrix(as.numeric(as.factor(snps))-1, nrow=nrow(snps), ncol=ncol(snps))
        rownames(snps) <- r.noms
        colnames(snps) <- c.noms
      }else{
        r.noms <- rownames(snps)
        c.noms <- colnames(snps)
        snps <- matrix(as.numeric(as.character(snps)), nrow=nrow(snps), ncol=ncol(snps))
        rownames(snps) <- r.noms
        colnames(snps) <- c.noms
      }
      na.after <- length(which(is.na(snps)))
      if(na.after > na.before){
        stop("NAs created in converting snps to numeric.")
      }
    }
  }
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
  #######################################################
  ## CHECK IF ANY INDIVIDUALS ARE 100% (or 75%??) NAS: ##
  #######################################################
  ## Rows: (mandatory)
  NA.tab <- sapply(c(1:nrow(snps)), function(e) length(which(is.na(snps[e,]))))

  ## Remove any ENTIRELY missing rows...
  toRemove <- rownames(snps)[which(NA.tab == ncol(snps))]

  ## Remove any NON-POLYMORPHIC rows? (i.e., entirely 1 or 0)
  # rs <- rowSums(snps, na.rm=TRUE)
  # # toRemove2 <- rownames(snps)[which(rs %in% c(0, ncol(snps)))]
  # l <- rep(ncol(snps), nrow(snps)) - NA.tab
  # toRemove2 <- rownames(snps)[which(rs == 0 | rs == l)]
  # toRemove <- c(toRemove, toRemove2)

  if(length(toRemove) > 0){
    ## print notice:
    cat("Removing", length(toRemove), "missing individual(s); only NAs in row(s).\n", sep=" ")
    ## remove individuals from snps:
    snps.ini <- snps # need snps.ini below
    snps <- snps[-which(rownames(snps) %in% toRemove), ]

    ## remove individuals from phen:
    phen <- phen[-which(names(phen) %in% toRemove)]
    ## remove individuals from tree:
    tree <- drop.tip(tree, tip = toRemove)

    ## (+ snps.reconstruction)
    if(is.matrix(snps.reconstruction)){
      toKeep <- which(rownames(snps.reconstruction) %in% c(rownames(snps), tree$node.label))
      snps.reconstruction <- snps.reconstruction[toKeep, ]
    }
    rm(snps.ini)
  }

  ####################################################################
  ##########################################
  ## CHECK IF ALL LOCI CONTAIN < 75% NAs: ##
  ##########################################
  ## Columns:
  if(is.null(na.rm)) na.rm <- TRUE

  ## Get number of NAs per column:
  NA.tab <- sapply(c(1:ncol(snps)), function(e) length(which(is.na(snps[,e]))))

  ## Remove any ENTIRELY missing columns
  toRemove <- colnames(snps)[which(NA.tab == nrow(snps))]

  ## Remove any NON-POLYMORPHIC columns? (i.e., entirely 1 or 0, -NAs)
  # cs <- colSums(snps, na.rm=TRUE)
  # l <- rep(nrow(snps), ncol(snps)) - NA.tab
  # toRemove2 <- colnames(snps)[which(cs == 0 | cs == l)]
  # toRemove <- c(toRemove, toRemove2)

  if(length(toRemove) > 0){
    snps <- snps[, -which(colnames(snps) %in% toRemove)]

    ## (+ snps.reconstruction)
    if(is.matrix(snps.reconstruction)){
      snps.reconstruction <- snps.reconstruction[, -which(colnames(snps.reconstruction) %in% toRemove)]
    }
    ## + update NA.tab:
    NA.tab <- NA.tab[-which(NA.tab == nrow(snps))]
  }

  ## Unless user declined, remove columns missing >= 75%
  if(na.rm != FALSE){
    ## What proportion of individuals w NAs is acceptable?
    toRemove <- colnames(snps)[which(NA.tab > floor(nrow(snps)*(3/4)))] # Remove columns w > 75% of NAs
    if(length(toRemove) > 0){
      snps <- snps[, -which(colnames(snps) %in% toRemove)]

      ## (+ snps.reconstruction)
      if(is.matrix(snps.reconstruction)){
        snps.reconstruction <- snps.reconstruction[, -which(colnames(snps.reconstruction) %in% toRemove)]
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
  ##########################################
  ## REORDER SNPS TO MATCH TREE$TIP.LABEL ##
  ##########################################
  if(!identical(as.character(rownames(snps)), as.character(tree$tip.label))){
    ord <- match(tree$tip.label, rownames(snps))
    snps <- snps[ord,]
    ## check:
    if(!identical(as.character(rownames(snps)), as.character(tree$tip.label))){
      stop("Unable to rearrange snps such that rownames(snps)
           match content and order of tree$tip.label.
           Please check that these match.")
    }
    }

  ## (+ snps.reconstruction)
  if(is.matrix(snps.reconstruction)){
    if(!identical(as.character(rownames(snps.reconstruction)[1:nrow(snps)]), as.character(tree$tip.label))){
      ord <- match(tree$tip.label, rownames(snps.reconstruction)[1:nrow(snps)])
      snps.reconstruction <- snps.reconstruction[c(ord, (length(ord)+1):nrow(snps.reconstruction)),]
      ## check:
      if(!identical(as.character(rownames(snps.reconstruction)[1:nrow(snps)]), as.character(tree$tip.label))){
        stop("Unable to rearrange snps.reconstruction such that
              rownames(snps.reconstruction)[1:nrow(snps)]
              match content and order of tree$tip.label.
              Please check that these match.")
      }
    }

  }

  ## REORDER PHEN TO MATCH TREE$TIP.LABEL
  if(!identical(as.character(names(phen)), as.character(tree$tip.label))){
    ord <- match(tree$tip.label, names(phen))
    phen <- phen[ord]
    ## check:
    if(!identical(as.character(names(phen)), as.character(tree$tip.label))){
      stop("Unable to rearrange phen such that names(phen)
           match content and order of tree$tip.label.
           Please check that these match.")
    }
    }

  ####################################################################
  ############################################################
  ## CHECK TREE: (set NEGATIVE branch lengths to zero (??)) ##
  ############################################################
  toChange <- which(tree$edge.length < 0)
  if(length(toChange) > 0){
    tree$edge.length[toChange] <- 0
    cat("Setting", length(toChange), "negative branch lengths to zero.\n", sep=" ")
  }

  ####################################################################
  #################
  ## Handle phen ##
  #################
  phen.rec.method <- levs <- NULL
  ## convert phenotype to numeric:
  ## NOTE--this is also necessary for returning results in step (5)!
  phen.ori <- phen

  ## Convert to numeric (required for assoc tests):
  na.before <- length(which(is.na(phen)))

  ## CHECK for discrete vs. continuous:
  ## NB: can only be binary or continuous at this point...
  levs <- unique(as.vector(unlist(phen)))
  n.levs <- length(levs[!is.na(levs)])
  ## BINARY: ##
  if(n.levs == 2){
    ## Convert phen to numeric:
    if(!is.numeric(phen)){
      if(all.is.numeric(phen)){
        phen <- as.numeric(as.character(phen))
      }else{
        phen <- as.numeric(as.factor(phen))
      }
    }
    ## Set phen.rec.method: ##
    phen.rec.method <- "discrete"
    if(!is.null(phen.type)){
      if(phen.type == "continuous"){
        phen.rec.method <- "continuous"
        warning("phen is binary. Are you sure phen.type is 'continuous'?")
      }
    } # end phen.type (binary)
  }else{
    ## DISCRETE or CONTINUOUS: ##
    ## Convert phen to numeric:
    if(!is.numeric(phen)){
      if(all.is.numeric(phen)){
        phen <- as.numeric(as.character(phen))
      }else{
        stop("phen has more than 2 levels but is not numeric (and therefore neither binary nor continuous).")
      }
    }
    ## Set phen.rec.method: ##
    phen.rec.method <- "continuous"
    ## Get proportion unique:
    prop.u <- length(unique(phen))/length(phen)
    if(!is.null(phen.type)){
      if(phen.type == "discrete"){
        phen.rec.method <- "discrete"
        if(prop.u > 0.5){
          cat("Performing *discrete* reconstruction, although phen is ", round(prop.u, 2)*100, "% unique: Are you sure phen.type is 'discrete'?", sep="")
        }
      }
    } # end phen.type (discrete/continuous)
  }
  ## ensure ind names not lost
  names(phen) <- names(phen.ori)

  ## Check that no errors occurred in conversion:
  na.after <- length(which(is.na(phen)))
  if(na.after > na.before){
    stop("NAs created while converting phen to numeric.")
  }
  ####################################################################
  #######################
  ## Reconstruct phen: ##
  #######################
  phen.rec <- NULL
  if(any(c("simultaneous", "subsequent") %in% test)){
    ## If user-provided reconsruction:
    if(length(phen.reconstruction) > 1){
      ## CHECK:
      if(length(phen.reconstruction) != max(tree$edge[,2])){
        warning("The number of individuals in the provided phen.reconstruction is not equal to the
                total number of nodes in the tree. Performing a new reconstruction instead.\n")
        if(phen.rec.method == "discrete") phen.reconstruction <- "parsimony"
        if(phen.rec.method == "continuous") phen.reconstruction <- "ml"
      }else{
        ## If phen provided and correct length:
        phen.reconstruction.ori <- phen.reconstruction
        levs <- unique(as.vector(unlist(phen.reconstruction)))
        n.levs <- length(levs[!is.na(levs)])
        ## Convert phen.reconstruction to numeric:
        if(!is.numeric(phen.reconstruction)){
          ## convert to numeric if possible:
          if(all.is.numeric(phen.reconstruction)){
            phen.reconstruction <- as.numeric(as.character(phen.reconstruction))
          }else{
            ## if binary, convert from factor:
            if(n.levs == 2){
              phen.reconstruction <- as.numeric(as.factor(phen.reconstruction))
            }else{
              ## if neither binary nor numeric, warning + reconstruct:
              warning("phen.reconstruction has more than 2 levels but is not numeric. Reconstructing phen internally instead.\n")
              if(phen.rec.method == "discrete") phen.reconstruction <- "parsimony"
              if(phen.rec.method == "continuous") phen.reconstruction <- "ml"
            }
          }
        }
        if(length(phen.reconstruction) > 1) names(phen.reconstruction) <- names(phen.reconstruction.ori)
      }
    }

    ## If not user-provided or checks failed, reconstruct ancestral states:
    if(length(phen.reconstruction) > 1){
      phen.rec <- phen.reconstruction
    }else{
      ## By PARSIMONY or ML: ##
      phen.rec <- asr(var = phen, tree = tree, type = phen.reconstruction, method = phen.rec.method)
    }
  }
  ####################################################################


  ###################
  ## filename.plot ##
  ###################
  if(!is.null(filename.plot)){
    ## save whatever plots are drawn before dev.off() is called to filename.plot:
    pdf(file = filename.plot, width = 7, height = 11)
  }


  ###############
  ## PLOT TREE ##
  ###############

  if(plot.tree==TRUE){
    ## get plot margins:
    mar.ori <- par()$mar

    ## get tip.col:
    leafCol <- "black"
    if(all.is.numeric(phen)){
      var <- as.numeric(as.character(phen))
    }else{
      var <- as.character(phen)
    }
    levs <- unique(var)
    if(length(levs) == 2){
      ## binary:
      myCol <- c("red", "blue")
      leafCol <- var
      ## for loop
      for(i in 1:length(levs)){
        leafCol <- replace(leafCol, which(leafCol == levs[i]), myCol[i])
      } # end for loop
    }else{
      if(is.numeric(var)){
        ## numeric:
        # myCol <- seasun(length(levs))
        myCol <- num2col(var, col.pal = seasun)
        leafCol <- myCol
      }else{
        ## categorical...
        myCol <- funky(length(levs))
        leafCol <- var
        ## for loop
        for(i in 1:length(levs)){
          leafCol <- replace(leafCol, which(leafCol == levs[i]), myCol[i])
        } # end for loop
      }
    }

    ## PLOT TREE:
    plot(tree, show.tip=T, tip.col=leafCol, align.tip.label=TRUE, cex=0.5)
    title("Phylogenetic tree")
    axisPhylo()

    ## node labels (optional):
    # if(!is.null(tree$node.label)){
    #   nodeLabs <- tree$node.label
    # }else{
    #   nodeLabs <- c((length(tree$tip.label)+1):max(tree$edge))
    # }
    # nodelabels(nodeLabs, cex=0.6, font=2, col="blue", frame="none")

    ## reset plot margins:
    par(mar=mar.ori)
  } # end plot.tree


  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

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
    names(n.subs) <- 1:length(n.subs)
    ## check?
    # sum(n.subs) == ncol(snps) ## should be TRUE
    # }

  }else{

    ## INPUT n.subs ##

    ## Assign names if null:
    if(is.null(names(n.subs))) names(n.subs) <- 1:length(n.subs)
  }

  ## check:
  # barplot(n.subs, col=transp("blue", 0.5), names=c(1:length(n.subs)))
  # title("Homoplasy distribution \n(treeWAS Fitch)")


  ######   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ############################################################   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ####### BEFORE RUNNING CHUNK-BY-CHUNK TREEWAS ...    #######   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ############################################################   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ######   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

  #######################
  ## Handle chunk.size ##
  #######################

  CHUNKS <- list()

  if(is.null(chunk.size)) chunk.size <- ncol(snps)
  if(chunk.size > ncol(snps)) chunk.size <- ncol(snps)
  ## Put all loci in one chunk segment:
  if(chunk.size == ncol(snps)){
    CHUNKS[[1]] <- 1:ncol(snps)
  }else{
    ## get chunks:
    chunks <- seq(1, ncol(snps), chunk.size)
    ## handle last chunk:
    if(chunks[length(chunks)] < ncol(snps)){
      # chunks <- c(chunks, ncol(snps))
      ## If last chunk too small (<= 10%(chunk.size)?), append to penultimate chunk:
      if((ncol(snps) - chunks[length(chunks)]) <= chunk.size*0.1){
        chunks <- chunks[1:(length(chunks)-1)]
        chunks <- c(chunks, (ncol(snps)+1))
      }else{
        chunks <- c(chunks, (ncol(snps)+1))
      }
    }else{
      chunks[length(chunks)] <- chunks[length(chunks)]+1
    }
    ## Make list of chunk segments:
    for(i in 1:(length(chunks)-1)){
      CHUNKS[[i]] <- chunks[i]:(chunks[(i+1)]-1)
    }
  }
  ## check?
  # for(i in 1:length(CHUNKS)) print(range(CHUNKS[[i]]))


  #########################################
  ## BEFORE RUNNING CHUNK-BY-CHUNK LOOP: ##
  #########################################
  ## Store "master" objects, by chunk: ##
  ## - snps
  ## - snps.rec (if present)
  ## - n.subs --> CONVERT s.t. sum(n.subs) == ncol(snps) --> AND SEGREGATE into COLUMNS for each CHUNK

  ## And make lists to store created objects: ##
  ## - corr.dat (x n.tests)
  ## - corr.sim (x n.tests)
  ## - snps.sim
  ## - snps.rec
  ## - snps.sim.rec
  CORR.DAT <- CORR.SIM <- SNPS.SIM <- SNPS.SIM.REC <- N.SNPS.SIM <- SNPS <- SNPS.REC <- N.SUBS <- list()

  #############################
  ## GET VALUES FOR CHUNK(S) ##
  #############################

  ##  One chunk only:
  if(length(CHUNKS) == 1){

    N.SUBS[[1]] <- n.subs
    SNPS[[1]] <- snps
    N.SNPS.SIM[[1]] <- n.snps.sim
    if(length(snps.reconstruction) > 1){
      SNPS.REC[[1]] <- snps.reconstruction
    }else{
      SNPS.REC[[1]] <- snps.reconstruction
    }

  }else{ # end no CHUNKS (one only)

    ## if working w multiple CHUNKS

    #################################################
    ## CONVERT n.subs st sum(n.subs) == ncol(snps) ##
    #################################################
    dist <- n.subs
    gen.size <- ncol(snps)

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
    n.subs <- dist.prop

    ## Remove unnecessary objects...
    rm(dist)
    rm(dist.prop)
    rm(dist.sum)

    ## Assign n.subs for each chunk:
    n.subs.ori <- n.subs
    c.lims <- N.SUBS <- list()
    for(i in 1:length(CHUNKS)) c.lims[[i]] <- range(CHUNKS[[i]])

    ## FOR LOOP to get n.subs segment for each chunk:
    for(i in 1:length(CHUNKS)){
      ## get chunk limits:
      N <- c.lims[[i]][2] - (c.lims[[i]][1]-1)
      start <- 1
      end <- 1
      n.chunk <- n.subs[start:end]
      ## get levels of n.subs for this chunk:
      counter <- 0
      while(sum(n.chunk) < N){
        counter <- counter+1
        # print(counter)
        n.chunk <- n.subs[start:(end+counter)]
      } # end while loop
      ## cut extra from last level:
      if(sum(n.chunk) > N){
        n.chunk[length(n.chunk)] <- n.chunk[length(n.chunk)] - (sum(n.chunk) - N)
      }
      ## upadte n.subs (remaining):
      toRemove <- which(names(n.subs) %in% names(n.chunk))
      n.subs[toRemove] <- n.subs[toRemove]  - n.chunk
      N.SUBS[[i]] <- n.chunk
    } # end for loop

    ## get SNPS, SNPS.REC, N.SNPS.SIM ##
    fac <- n.snps.sim/ncol(snps)
    for(i in 1:length(CHUNKS)){
      SNPS[[i]] <- snps[, CHUNKS[[i]]]
      N.SNPS.SIM[[i]] <- ncol(SNPS[[i]])*fac
      if(length(snps.reconstruction) > 1){
        SNPS.REC[[i]] <- snps.reconstruction[, CHUNKS[[i]]]
      }else{
        SNPS.REC[[i]] <- snps.reconstruction
      }
    } # end for loop

  } # end multiple CHUNKS



  ######   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ############################################################   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ###### FOR LOOP for CHUNK-BY-CHUNK TREEWAS STARTS HERE #####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ############################################################   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ######   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

  for(i in 1:length(CHUNKS)){

    snps <- SNPS[[i]]
    n.snps.sim <- N.SNPS.SIM[[i]]
    snps.reconstruction <- SNPS.REC[[i]]
    n.subs <- N.SUBS[[i]]

    #####################################################
    ## 1) Simulate null genetic dataset to compare your #
    ## real correlations w phen to  #####################
    #####################################################

    ## check n.snps.sim:
    if(is.null(n.snps.sim)){
      n.snps.sim <- ncol(snps)*10
    }else{
      ## Update n.snps.sim to match the REAL ncol(snps) after checks...
      if(n.snps != ncol(snps)){
        if(n.snps.sim == n.snps){
          n.snps.sim.ori <- n.snps.sim
          n.snps.sim <- ncol(snps)

          ## Print update notice:
          cat("Note: Updating n.snps.sim to match the number of real loci after data cleaning.
          Input:", n.snps.sim.ori, " -->
          Updated:", n.snps.sim, "\n")

        }else{
          Nx <- c(2:200)
          if(any(n.snps*Nx == n.snps.sim)){
            n.snps.sim.ori <- n.snps.sim
            Nx <- Nx[which(Nx == (n.snps.sim/n.snps))]
            n.snps.sim <- ncol(snps)*Nx

            ## Print update notice:
            cat("Note: Updating n.snps.sim to match ", Nx, "x the number of real loci after data cleaning.
            Input: ", n.snps.sim.ori, " -->
            Updated: ", n.snps.sim, "\n", sep="")

            ## TO DO: ##
            ## If the user, for whatever reason, wanted to simulate the number they requested
            ## (ie. matching the input n.snps), could either
            ## (A) Add an argument like upate.n.snps.sim = FALSE,
            ## (B) Tell them to do their own data cleaning exactly as I do but outside of/before running the treeWAS function.
            ## (C) Fudge--tell them to add 1, eg. 10,000 --> 10,001?
          }
        }
      }
  }
    out <- genomes <- snps.mat <- list()

    ## SIMULATE A DATASET | your tree ##
    if(!is.null(seed)) set.seed(seed)
    out[[1]] <- snp.sim(n.snps = n.snps.sim,
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

    # gc()


    print(paste("treeWAS snps sim done @", Sys.time()))


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

    ## Handle phen
    ## (moved to initial checks)


    ###########################
    ## GET UNIQUE SNPS(.SIM) ##
    ###########################

    ## Get UNIQUE snps + index
    snps.complete <- snps
    temp <- get.unique.matrix(snps, MARGIN=2)
    snps <- temp$unique.data
    snps.index <- temp$index

    ## Get UNIQUE snps.sim + index
    snps.sim.complete <- snps.sim
    temp <- get.unique.matrix(snps.sim, MARGIN=2)
    snps.sim <- temp$unique.data
    snps.sim.index <- temp$index

    ##############################################################################################
    ## Reconstruct ancestral SNPs & phen by parsimony/ML (for tests simultaneous & subsequent)  ##
    ##############################################################################################

    ## Ensure we are only reconstructing ancestral states ONCE here, to be used in MULTIPLE tests later.
    snps.REC <- snps.sim.REC <- NULL

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
        if(nrow(snps.reconstruction) != max(tree$edge[,2])){
          warning("The number of rows in the provided snps.reconstruction is not equal to the
                  total number of nodes in the tree. Performing a new parsimonious reconstruction instead.\n")
          snps.reconstruction <- snps.sim.reconstruction <- "parsimony"
        }
        if(ncol(snps.reconstruction) != ncol(snps)){
          warning("The number of columns in the provided snps.reconstruction is not equal to the number of
                  columns in the snps matrix. Performing a new parsimonious reconstruction instead.\n")
          snps.reconstruction <- snps.sim.reconstruction <- "parsimony"
        }
      }

      ## If not user-provided or checks failed, reconstruct ancestral states:
      if(is.matrix(snps.reconstruction)){
        snps.rec <- snps.reconstruction
        ## Get UNIQUE snps.reconstruction
        snps.rec.complete <- snps.rec
        temp <- get.unique.matrix(snps.rec, MARGIN=2)
        snps.rec <- temp$unique.data
        snps.rec.index <- temp$index
        if(!identical(snps.rec.index, snps.index)){
          warning("Careful-- snps and snps.rec should have the same index when reduced
              to their unique forms.\n") ## SHOULD THIS BE A "STOP" INSTEAD? OR IS THIS ERROR NOT FATAL OR NOT POSSIBLE????
        }
      }else{
        # system.time( # 274
        snps.rec <- asr(var = snps, tree = tree, type = snps.reconstruction, unique.cols = TRUE)
        # )
      }

      #################################
      ## Reconstruct SIMULATED SNPs: ##
      #################################
      # system.time(
      snps.sim.rec <- asr(var = snps.sim, tree = tree, type = snps.sim.reconstruction, unique.cols = TRUE)
      # )

        } # end reconstruction for tests 2 & 3


    ## !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ###

    # ## Get UNIQUE snps.sim.reconstruction
    # snps.sim.rec.complete <- snps.sim.rec
    # temp <- get.unique.matrix(snps.sim.rec, MARGIN=2)
    # snps.sim.rec <- temp$unique.data
    # snps.sim.rec.index <- temp$index
    # if(!identical(snps.sim.rec.index, snps.sim.index)){
    #   warning("Careful-- snps.sim and snps.sim.rec should have the same index when reduced
    #           to their unique forms!") ## SHOULD THIS BE A "STOP" INSTEAD? OR IS THIS ERROR NOT FATAL OR NOT POSSIBLE????
    # }

    print(paste("Reconstructions completed @", Sys.time()))

    # gc()

    #######################
    ## identify sig.snps ##
    #######################
    ## Note: UNIQUE snps & snps.sim are identified WITHIN the get.sig.snps fn
    ## to reduce computational time, but results are identified on the basis of all
    ## ORIGINAL snps & snps.sim columns inputted.

    # test <- c("terminal", "simultaneous", "subsequent")
    TEST <- as.list(test)
    CORR.DAT[[i]] <- CORR.SIM[[i]] <- vector(mode="list", length=length(test))
    # names(CORR.DAT[[i]]) <- names(CORR.SIM[[i]]) <- test

    ## Run get.assoc.scores for each assoc.test
    ## Then run get.sig.snps (ONLY once ALL of corr.sim, corr.dat have been generated w get.assoc.scores (IFF running chunk-by-chunk!))
    # system.time(
      for(t in 1:length(TEST)){
        assoc.scores <- get.assoc.scores(snps = snps,
                                         snps.sim,
                                         phen = phen,
                                         tree = tree,
                                         test = TEST[[t]],
                                         snps.reconstruction = snps.rec,
                                         snps.sim.reconstruction = snps.sim.rec,
                                         phen.reconstruction = phen.rec,
                                         unique.cols = TRUE)
        ## EXPAND CORR.DAT & CORR.SIM:
        cd <- assoc.scores$corr.dat
        cs <- assoc.scores$corr.sim
        cd <- cd[snps.index]
        cs <- cs[snps.sim.index]
        names(cd) <- colnames(snps.complete)
        names(cs) <- colnames(snps.sim.complete)

        ## STORE DATA FOR EACH TEST:
        CORR.DAT[[i]][[t]] <- cd
        CORR.SIM[[i]][[t]] <- cs
        rm(assoc.scores)
      } # end for (t) loop
    # )

    ## EXPAND DATA BEFORE STORING FOR THIS CHUNK:
    snps.sim <- snps.sim.complete
    snps.rec <- snps.rec[,snps.index]
    colnames(snps.rec) <- colnames(snps.complete)
    snps.sim.rec <- snps.sim.rec[,snps.sim.index]
    colnames(snps.sim.rec) <- colnames(snps.sim.complete)

    ## STORE DATA CREATED FOR THIS CHUNK:
    SNPS.SIM[[i]] <- snps.sim
    SNPS.REC[[i]] <- snps.rec
    SNPS.SIM.REC[[i]] <- snps.sim.rec

    } # end for (i) loop (CHUNKS)
  ######   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ############################################################   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ####### END LOOP/OPTIONAL CHUNK-BY-CHUNK TREEWAS HERE ######   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ############################################################   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  ######   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###

  ############################
  ## (COMBINE CHUNKS FIRST) ##
  ############################

  ## Get corr.dat & corr.sim:
  noms <- names(CORR.DAT[[1]])
  CORR.DAT <- sapply(c(1:length(CORR.DAT[[1]])), function(e) unlist(sapply(CORR.DAT, "[", e)), simplify=FALSE)
  names(CORR.DAT) <- test
  CORR.SIM <- sapply(c(1:length(CORR.SIM[[1]])), function(e) unlist(sapply(CORR.SIM, "[", e)), simplify=FALSE)
  names(CORR.SIM) <- test

  ## Get snps, snps.sim, snps.rec, snps.sim.rec:
  snps <- do.call(cbind, SNPS)
  rm(SNPS)
  snps.sim <- do.call(cbind, SNPS.SIM)
  rm(SNPS.SIM)
  snps.rec <- do.call(cbind, SNPS.REC)
  rm(SNPS.REC)
  snps.sim.rec <- do.call(cbind, SNPS.SIM.REC)
  rm(SNPS.SIM.REC)

  #######################################
  ## Get sig.snps for ALL assoc.scores ##
  #######################################
  sig.list <- list()

  # system.time(
    for(t in 1:length(TEST)){
      sig.list[[t]] <- get.sig.snps(corr.dat = CORR.DAT[[t]],
                                    corr.sim = CORR.SIM[[t]],
                                    snps.names = colnames(snps),
                                    test = TEST[[t]],
                                    n.tests = length(TEST),
                                    p.value = p.value,
                                    p.value.correct = p.value.correct,
                                    p.value.by = p.value.by)
    }
  # )

  names(sig.list) <- test

  ## Remove unnecessary objects (already stored in sig.list):
  rm(CORR.DAT)
  rm(CORR.SIM)


  print(paste("ID of significant loci completed @", Sys.time()))


  #################
  ## GET RESULTS ##
  #################

  ## set margins for plotting:
  par.mar.ori <- par()$mar
  # par(mar=c(5, 2, 4, 1)+0.1)
  par(mar=c(5, 4, 4, 1)+0.1)

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

    ############################
    ## 4) (A) MANHATTAN PLOT! ##
    ############################
    if(plot.manhattan == TRUE){

      manhattan.plot(p.vals = abs(corr.dat),
                     col = "funky",
                     transp = 0.25,
                     sig.thresh = sig.thresh,
                     thresh.col="red",
                     snps.assoc = NULL,
                     snps.assoc.col = "red",
                     jitter.amount = 0.00001,
                     min.p = NULL,
                     log10=FALSE,
                     ylab=paste(TEST[[i]], "score", sep=" "))
      ## Add subtitle:
      title(paste("\n \n(", TEST[[i]], "score)"), cex.main=0.9)
    } # end plot manhattan


    ####################################
    ## 4) (B,C) Plot the distribution ##
    ####################################

    if(plot.null.dist == TRUE){
      ## Generate one histogram per test:
      plot.sig.snps(corr.dat = abs(corr.dat),
                    corr.sim = abs(corr.sim),
                    corr.sim.subset = NULL,
                    sig.corrs = abs(corr.dat[sig.snps]),
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
    } # end plot.null.dist

    if(plot.dist == TRUE){
      ## Generate one histogram per test:
      plot.sig.snps(corr.dat = abs(corr.dat),
                    corr.sim = abs(corr.sim),
                    corr.sim.subset = NULL,
                    sig.corrs = abs(corr.dat[sig.snps]),
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
                    plot.null.dist = FALSE,
                    plot.dist = TRUE)
    } # end plot.dist


    ########################################
    ## 5) Return results list ##############
    ########################################
    phen.curr <- phen # store current phen for output
    if(length(sig.snps)==0) sig.snps <- sig.corrs <- NULL

    ###########
    ## Make a data.frame containing all relevant output for sig.snps
    if(length(sig.snps) > 0){

      ## Get counts for n.sig.snps in each cell of the contingency table:
      toKeep <- sig.snps
      snps.toKeep <- snps[,toKeep]

      ## Only get G1P1 etc. if binary phen...
      G1P1 <- G0P0 <- G1P0 <- G0P1 <- NA
      levs <- unique(as.vector(unlist(phen)))
      levs <- levs[!is.na(levs)]
      n.levs <- length(levs)

      if(n.levs == 2){

        ## store ind names:
        noms <- names(phen)
        ## If binary, convert phen to 0/1:
        phen <- as.numeric(as.factor(phen))
        phen <- rescale(phen, to=c(0,1)) # require(scales)
        # if(length(phen[-c(which(phen==1), which(phen==2))])==0){
        #   phen <- replace(phen, which(phen==1), 0)
        #   phen <- replace(phen, which(phen==2), 1)
        # }
        ## ensure ind names not lost:
        names(phen) <- noms

        if(length(toKeep) > 1){
          G1P1 <- sapply(c(1:ncol(snps.toKeep)),
                         function(e)
                           length(which(snps.toKeep[which(phen==1),e]==1)))
          G0P0 <- sapply(c(1:ncol(snps.toKeep)),
                         function(e)
                           length(which(snps.toKeep[which(phen==0),e]==0)))
          G1P0 <- sapply(c(1:ncol(snps.toKeep)),
                         function(e)
                           length(which(snps.toKeep[which(phen==0),e]==1)))
          G0P1 <- sapply(c(1:ncol(snps.toKeep)),
                         function(e)
                           length(which(snps.toKeep[which(phen==1),e]==0)))
        }else{
          ## if only ONE sig snp (haploid) identified:
          G1P1 <- length(which(snps.toKeep[which(phen==1)]==1))
          G0P0 <- length(which(snps.toKeep[which(phen==0)]==0))
          G1P0 <- length(which(snps.toKeep[which(phen==0)]==1))
          G0P1 <- length(which(snps.toKeep[which(phen==1)]==0))

        }
        df <- data.frame(sig.snps,
                         sig.p.vals,
                         sig.corrs,
                         G1P1, G0P0, G1P0, G0P1)
        names(df) <- c("SNP.locus",
                       "p.value",
                       "score",
                       "G1P1", "G0P0", "G1P0", "G0P1")
      }else{ # end G1P1 etc.

      ## If G1P1 etc. not valid (non-binary phen),
      ## return df without these variables..
      df <- data.frame(sig.snps,
                       sig.p.vals,
                       sig.corrs)
      names(df) <- c("SNP.locus",
                     "p.value",
                     "score")

      ## NOTE: Could return sig.snps.names as column rather than rownames..? ####   (??)   ####
      }

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

  #######################
  ## end filename.plot ##
  #######################
  if(!is.null(filename.plot)){
    dev.off()
  }


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

  phen <- phen.curr

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








# ##############################################
# ## ASSIGN NODE LABELS to TREE (& SNPS.REC): ##
# ##############################################
# ## If NO snps.rec AND tree has NO nodelabs, ASSIGN nodelabs to tree:
# if(!is.matrix(snps.reconstruction)){
#   if(is.null(tree$node.label)){
#     tree$node.label <- paste("NODE", c((length(tree$tip.label)+1):max(tree$edge)), sep="_")
#   }
# }else{
#   ## If snps.rec is MATRIX..##
#   ## get index for terminal nodes:
#   ixt <- c(1:length(tree$tip.label))
#   ## get index for internal nodes:
#   ixi <- c((nrow(snps)+1):nrow(snps.reconstruction))
#   ## Unless tree has node.labs AND they match snps.rec names...
#   if(is.null(tree$node.label) | !identical(tree$node.label, rownames(snps)[ixi])){
#     ## If snps.rec[ixt,] match tree$tip.labs, set matching NODElabs (one or both):
#     if(identical(rownames(snps.reconstruction)[ixt], tree$tip.label)){
#       ## If snps.rec HAS labs, BUT tree does NOT:
#       if(is.null(tree$node.label)){
#         ## If snps.rec labs NOT numeric, assign to nodelabs as is
#         if(!all.is.numeric(rownames(snps.reconstruction)[ixi])){
#           tree$node.label <- rownames(snps.reconstruction)[ixi]
#         }else{
#           ## If snps.rec is NUMERIC, paste NODE_ and assign to both snps.rec and nodelabs:
#           tree$node.label <- rownames(snps.reconstruction)[ixi] <- paste("NODE", rownames(snps.reconstruction)[ixi], sep="_")
#         }
#         ## + print notice:
#         cat("tree$node.label is NULL. Assuming identical to rownames(snps.reconstruction)[Nterminal+1:Ntotal].\n")
#
#       }else{
#         ## If nodelabs & snps.rec rowlabs are NOT identical (while both are PRESENT and TERMINAL labs MATCH)...
#         ## Unless they MATCH in a different order...
#         if(length(which(is.na(match(tree$node.label, rownames(snps.reconstruction)[ixi])))) > 0){
#           ## see if you can remove NODE_ & achieve match, otherwise,
#           ## assign nodelabs to snps.rec rows & print warning (snps.rec order may not match nodelabs order):
#           nodeNs <- removeFirstN(tree$node.label, 5)
#           ord <- match(nodeNs,rownames(snps.reconstruction)[ixi])
#           if(all.is.numeric(nodeNs) & length(which(is.na(ord)))==0){
#             rownames(snps.reconstruction)[ixi] <- nodeNs[ord]
#           }else{
#             rownames(snps.reconstruction)[ixi] <- tree$node.label
#             warning("Rownames of snps.reconstruction[Nterminal+1:Ntotal,] do not match tree$node.label.
#                   Cannot be sure that reconstruction's rows correspond to tree's internal nodes.\n")
#           }
#         }
#       }
#
#
#     }else{
#       ## If snps.rec[ixt] and tiplabs do NOT match (while nodelabs is NULL or does not match snps.rec[ixi])
#       ## assign labs to nodelabs and snps.rec[ixi]
#       ## + print WARNING (snps.rec order may not match nodelabs order):
#       tree$node.label <- rownames(snps.reconstruction)[ixi] <- paste("NODE", ixi, sep="_")
#       warning("Rownames of snps.reconstruction[Nterminal+1:Ntotal,] do not match tree$node.label.
#                 Cannot be sure that reconstruction's rows correspond to tree's internal nodes.\n")
#     }
#   }
# } # end pipeline to check/assign matching snps.rec rownames and tree$node.labels
##############################################################################################


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
