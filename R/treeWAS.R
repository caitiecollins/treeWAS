
### TO DO:
## - Update memreq estimate to be on the basis of UNIQUE snps columns rather than TOTAL ncol...
## - Change snps, SNPS, snps.sim, SNPS.SIM, snps.rec, SNPS.REC to logicals instead of numerics



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
#' @param sort.by.p A logical indicating whether to sort the results by decreasing p-value (\code{TRUE})
#'                  or by locus (\code{FALSE}, the default).
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export
#' @import adegenet

########################################################################
# #' @rdname treeWAS
# #' @method treeWAS print
# #' @S3method treeWAS print
# https://stackoverflow.com/questions/7198758/roxygen2-how-to-properly-document-s3-methods

print.treeWAS <- function(x, sort.by.p = FALSE, digits = 3){
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

  ## only sort/print if any sig loci present:
  if(l > 0){
    if(all.is.numeric(res)){
      res <- as.character(sort(as.numeric(res, decreasing=FALSE)))
    }else{
      ## get corresponding loci:
      locus <- sapply(c(1:length(res)), function(e) which(colnames(x$dat$snps) == res[e]))
      ## reorder by locus order:
      ord <- match(locus, sort(locus, decreasing=F))
      res <- res[ord]
    }

    cat("Significant loci: \n")
    print(res)
  }

  cat("\t \n")

  cat("\t######################## \n")
  cat("\t## Findings by test:  ## \n")
  cat("\t######################## \n")

  if("pairwise" %in% names(x)){
    ## categorical (+pairwise)
    N <- c(2:(length(x)-2))
  }else{
    N <- c(2:(length(x)-1))
  }

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
        print(signif(res, digits = digits))
      }else{
        sl <- sort(res$SNP.locus, decreasing=FALSE)
        ord <- match(sl, res$SNP.locus)
        print(signif(res[ord, ], digits = digits))
      }
    }
    cat("\t \n")

  } # end for (i) loop




  ## Print pairwise results for categorical phens:
  if("pairwise" %in% names(x)){

    cat("\t############################################ \n")
    cat("\t## Pairwise scores for significant snps:  ## \n")
    cat("\t############################################ \n")

    PW <- x$pairwise$pairwise.combined
    PW.res <- x$pairwise$pairwise.tests

    ## No sig snps:
    if(is.null(PW)){

      cat("No snps achieved overall significance, so no pairwise tests performed. \n")

    }else{
      ## PW tests on sig.snps:
      if(is.list(PW)){
        for(p in 1:length(PW)){

          pair <- names(PW)[p]

          cat("\t", paste(c("########", rep("#", nchar(pair))), collapse=""), " \n")
          cat("\t ## ", pair, " ## \n")
          cat("\t", paste(c("########", rep("#", nchar(pair))), collapse=""), " \n")


          cat("Significance thresholds: \n")
          sig.thresh <- NA
          for(t in 1:length(PW.res[[p]])){
            sig.thresh[t] <- PW.res[[p]][[t]]$sig.thresh
          } # end for (t) loop
          names(sig.thresh) <- names(PW.res[[p]])
          st.df <- data.frame(sig.thresh)
          print(st.df)


          print(signif(PW[[p]], digits = digits))

          cat("\t \n")

        } # end for (p) loop
      }
    }

  } # end pairwise

} # end print.treeWAS





###################
## write.treeWAS ##
###################

########################################################################

###################
## DOCUMENTATION ##
###################

#' Write \code{treeWAS} output to a CSV file.
#'
#' Save the results of \code{treeWAS} to a CSV file as a summary table of significant findings and scores
#' (excluding longer data elements within the output).
#' .
#'
#' @param x The output returned by \code{treeWAS}.
#' @param filename A character string containing the path and filename to which the .csv file will be saved;
#'                 by default, \code{filename = "./treeWAS_results"} and so
#'                 would be saved to the current working directory.
#'
#' @examples
#' ## Example ##
#' \dontrun{
#' ## Load data:
#' data(snps)
#' data(phen)
#' data(tree)
#'
#' ## Run treeWAS:
#' out <- treeWAS(snps, phen, tree, seed = 1)
#'
#' ## Save results to home directory:
#' write.treeWAS(x = out, filename = "~/treeWAS_results")
#' }
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @export

########################################################################

write.treeWAS <- function(x, filename="./treeWAS_results"){

  ## get sig locus names:
  name <- x$treeWAS.combined$treeWAS.combined

  ## If NO sig findings, print dummy (NULL) csv:
  if(length(name[!is.na(name)]) == 0){

    ## Make dummy variables:
    name <- locus <- p.value.1 <- p.value.2 <- p.value.3 <-
      score.1 <- score.2 <- score.3 <- G1P1 <- G0P0 <- G1P0 <- G0P1 <- rep(NA, 1)

    ## Make table:
    tab <- data.frame(name, locus, score.1, score.2, score.3,
                      p.value.1, p.value.2, p.value.3, G1P1, G0P0, G1P0, G0P1)

    ## If SIG loci found, print csv:
  }else{
    ## get corresponding loci:
    locus <- sapply(c(1:length(name)), function(e) which(colnames(x$dat$snps) == name[e]))
    ## reorder by locus order:
    ord <- match(locus, sort(locus, decreasing=F))
    locus <- locus[ord]
    name <- name[ord]

    ## Make dummy variables:
    p.value.1 <- p.value.2 <- p.value.3 <-
      score.1 <- score.2 <- score.3 <- G1P1 <- G0P0 <- G1P0 <- G0P1 <- rep(NA, length(name))

    ## Make table:
    tab <- data.frame(name, locus, score.1, score.2, score.3,
                      p.value.1, p.value.2, p.value.3, G1P1, G0P0, G1P0, G0P1)

    ## Add relevant values for each test: ##
    ## Terminal:
    if(length(x$terminal$sig.snps) > 1){
      df <- x$terminal$sig.snps
      toKeep <- match(df$SNP.locus, locus)
      tab$p.value.1[toKeep] <- df$p.value
      tab$score.1[toKeep] <- df$score
      if("G1P1" %in% names(df)){
        tab$G1P1[toKeep] <- df$G1P1
        tab$G0P0[toKeep] <- df$G0P0
        tab$G1P0[toKeep] <- df$G1P0
        tab$G0P1[toKeep] <- df$G0P1
      }
    }

    ## Simultaneous:
    if(length(x$simultaneous$sig.snps) > 1){
      df <- x$simultaneous$sig.snps
      toKeep <- match(df$SNP.locus, locus)
      tab$p.value.2[toKeep] <- df$p.value
      tab$score.2[toKeep] <- df$score
      if("G1P1" %in% names(df)){
        tab$G1P1[toKeep] <- df$G1P1
        tab$G0P0[toKeep] <- df$G0P0
        tab$G1P0[toKeep] <- df$G1P0
        tab$G0P1[toKeep] <- df$G0P1
      }
    }

    ## Subsequent:
    if(length(x$subsequent$sig.snps) > 1){
      df <- x$subsequent$sig.snps
      toKeep <- match(df$SNP.locus, locus)
      tab$p.value.3[toKeep] <- df$p.value
      tab$score.3[toKeep] <- df$score
      if("G1P1" %in% names(df)){
        tab$G1P1[toKeep] <- df$G1P1
        tab$G0P0[toKeep] <- df$G0P0
        tab$G1P0[toKeep] <- df$G1P0
        tab$G0P1[toKeep] <- df$G0P1
      }
    }
  }

  ## Add .CSV to filename:
  suff <- tolower(keepLastN(filename, 4))
  names(suff) <- NULL
  if(!identical(suff, ".csv")){
    filename <- paste(filename, ".csv", sep="")
  }

  ## Write table to CSV file:
  write.table(tab, file=filename, row.names = FALSE)

} # end write.treeWAS




#############
## treeWAS ##
#############


###################################################################################################################################

###################
## DOCUMENTATION ##
###################

#' Phylogenetic tree-based GWAS for microbes.
#'
#' This function implements a phylogenetic approach to genome-wide association studies (GWAS)
#' designed for use in bacteria and viruses.
#' The \code{treeWAS} approach measures the statistical association
#' between a phenotype of interest and the genotype at all loci,
#' seeking to identify significant associations, while accounting for the
#' confounding effects of clonal population structure and homologous recombination.
#'
#'
#' @param snps A matrix containing binary genetic data, with individuals in the rows
#'                and genetic loci in the columns and both rows and columns labelled.
#' @param phen A vector containing the phenotypic state of each individual, whose length is equal to the number of rows in
#'              \code{snps} and which is named with the same set of labels. The phenotype can be either
#'              binary (character or numeric), categorical (character or numeric), discrete or continuous (numeric).
#'              You must specify which type it is via the \code{phen.type} argument
#'              (otherwise treeWAS will assume binary or continuous).
#' @param tree A \code{phylo} object containing the phylogenetic tree; or, a character string,
#'                one of \code{"NJ"}, \code{"BIONJ"} (the default), or \code{"parsimony"};
#'                or, if NAs are present in the distance matrix, one of: \code{"NJ*"} or \code{"BIONJ*"},
#'                specifying the method of phylogenetic reconstruction.
#' @param phen.type An optional character string specifying whether the phenotypic variable should be treated
#'                  as either \code{"categorical"}, \code{"discrete"} or \code{"continuous"}.
#'                  If \code{phen.type} is \code{NULL} (the default), ancestral state reconstructions performed via ML
#'                  will treat any binary phenotype as discrete and any non-binary phenotype as continuous.
#'                  If \code{phen.type} is \code{"categorical"}, ML reconstructions and association tests will treat
#'                  values as nominal (not ordered) levels and not as meaningful numbers.
#'                  Categorical phenotypes must have >= 3 unique values (<= 5 recommended).
#'                  If \code{phen.type} is \code{"continuous"}, ML reconstructions will treat
#'                  values as meaningful numbers and may infer intermediate values.
#' @param n.subs A numeric vector containing the homoplasy distribution (if known, see details), or \code{NULL} (the default).
#' @param n.snps.sim An integer specifying the number of loci to be simulated for estimating the null distribution
#'                      (by default \code{10*ncol(snps)}). If memory errors arise during the analyis of a large dataset,
#'                      it may be necesary to reduce \code{n.snps.sim} from a multiple of 10 to, for example, 5x the number of loci.
#' @param chunk.size An integer indicating the number of \code{snps} loci to be analysed at one time.
#'                   This may be needed for large datasets or machines with insufficient memory.
#'                   Smaller values of \code{chunk.size} will increase the computational time required
#'                   (e.g., if \code{chunk.size = ncol(snps)/2}, treeWAS will take twice as long to complete).
#' @param mem.lim A logical or numeric value to set a memory limit for large datasets.
#'                If \code{FALSE} (the default), no limit is estimated and \code{chunk.size} is not changed.
#'                If \code{TRUE}, available memory is estimated with \code{memfree()} and \code{chunk.size} is reduced.
#'                A single numeric value can be used to set a memory limit (in GB) and \code{chunk.size} will be reduced accordingly.
#' @param test A character string or vector containing one or more of the following available tests of association:
#'              \code{"terminal"}, \code{"simultaneous"}, \code{"subsequent"}, \code{"cor"}, \code{"fisher"}.
#'              By default, the first three tests are run (see details).
#' @param correct.prop A logical indicating whether the \code{"terminal"} and \code{"subsequent"} tests will be corrected for
#'                     phenotypic class imbalance. Recommended if the proportion of individuals varies significantly across
#'                     the levels of the phenotype (if binary) or if the phenotype is skewed (if continuous).
#'                     If \code{correct.prop} is \code{FALSE} (the default), the original version of each test is run.
#'                     If \code{TRUE}, an alternate association metric based on the phi correlation coefficient
#'                     is calculated across the terminal and all (internal and terminal) nodes, respectively.
#' @param snps.reconstruction Either a character string specifying \code{"parsimony"} (the default) or \code{"ML"} (maximum likelihood)
#'                              for the ancestral state reconstruction of the genetic dataset,
#'                              or a matrix containing this reconstruction if it has been performed elsewhere
#'                              \emph{and} you provide the tree.
#' @param snps.sim.reconstruction A character string specifying \code{"parsimony"} (the default) or \code{"ML"} (maximum likelihood)
#'                                  for the ancestral state reconstruction of the simulated null genetic dataset.
#' @param phen.reconstruction Either a character string specifying \code{"parsimony"} (the default) or \code{"ML"} (maximum likelihood)
#'                              for the ancestral state reconstruction of the phenotypic variable,
#'                              or a vector containing this reconstruction if it has been performed elsewhere.
#' @param na.rm A logical indicating whether columns in \code{snps} containing more than 75\% \code{NA}s
#'                should be removed at the outset (TRUE, the default) or not (FALSE).
#' @param p.value A number specifying the base p-value to be set as the threshold of significance (by default, \code{0.01}).
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
#' @param plot.null.dist.pairs Either a logical indicating, for categorical phenotypes only, whether to plot
#'                             additional null distributions of association score statistics for each
#'                             pairwise comparison of phenotype levels (\code{TRUE}) or not (\code{FALSE}),
#'                             or the character string \code{"grid"} (the default) which will
#'                             print all of these plots in one grid (N pairs x N tests).
#' @param snps.assoc An optional character string or vector specifying known associated loci to be demarked in
#'                      results plots (e.g., from previous studies or if data is simulated); else \code{NULL}.
#' @param filename.plot An optional character string denoting the file location for
#'                        saving any plots produced (eg. "C:/Home/treeWAS_plots.pdf"); else \code{NULL}.
#' @param seed An optional integer to control the pseudo-randomisation process and allow
#'                for identical repeat runs of the function; else \code{NULL}.
#'
###########################################################
#'
#' @details
###########################################################
#' \strong{The treeWAS Approach}
#'
#' For a description of the approach adopted in our method,
#' please see either the \code{treeWAS} vignette
#' or the \code{treeWAS} \href{https://github.com/caitiecollins/treeWAS/wiki}{GitHub Wiki}.
#'
#' A more detailed description of our method can also be found in our paper, available in
#' \href{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005958}{PLOS Computational Biology}.
#'
###########################################################
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
#' (i.e., Denoting the state of the second allele as the inverse of the previous colummn), these should be removed.
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
#' \item{\code{terminal}}{The \code{terminal} test solves the following equation, for each genetic locus,
#'                        at the terminal nodes of the tree only:
#'                        \deqn{Terminal = | (1/Nterm)*(Pt*Gt - (1 - Pt)*Gt - Pt*(1 - Gt) + (1 - Pt)*(1 - Gt)) |}
#'                        The \code{terminal} test is a sample-wide test of association that
#'                        seeks to identify broad patterns of correlation between genetic loci and the phenotype,
#'                        without relying on inferences drawn from reconstructions of the ancestral states.}
#'
#' \item{\code{simultaneous}}{The \code{simultaneous} test solves the following equation, for each genetic locus,
#'                            across each branch in the tree:
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
#' \item{\code{subsequent}}{The \code{subsequent} test solves the following equation, for each genetic locus,
#'                            across each branch in the tree:
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
#'           The minimum p-value. P-values listed as zero can only truly be defined as being below this value.}
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
#'                 seed = 1)
#'
#' ## examine output:
#' print(out)
#'
#' }
#'
#'
###########################################################
#' @references Collins, C. and Didelot, X. "A phylogenetic method to perform genome-wide association studies in microbes
#' that accounts for population structure and recombination."
#' \href{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005958}{\emph{PLoS Comput. Biol.}},
#' vol. 14, p. e1005958, Feb. 2018.
#'
###########################################################
#'
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @import adegenet
#' @rawNamespace import(ape, except = zoom)
#' @importFrom Hmisc all.is.numeric
#' @importFrom phangorn midpoint
#' @importFrom scales rescale
#' @importFrom pryr object_size
#' @importFrom grDevices col2rgb dev.off heat.colors pdf rgb
#' @importFrom graphics arrows axis barplot box hist image lines par plot.new points rect text title
#' @importFrom stats anova as.formula chisq.test cor density dist ecdf fisher.test ftable glm lm mantelhaen.test
#' p.adjust quantile residuals rexp rnorm rpois
#' @importFrom utils str write.table
#'
#' @export

###################################################################################################################################
# @useDynLib phangorn, .registration = TRUE

# #' @import ape, except=zoom

# importFrom("grDevices", "col2rgb", "dev.off", "heat.colors", "pdf",
#            "rgb")
# importFrom("graphics", "arrows", "axis", "barplot", "box", "hist",
#            "image", "lines", "par", "plot.new", "points", "rect",
#            "text", "title")
# importFrom("stats", "anova", "as.formula", "cor", "density", "dist",
#            "ecdf", "fisher.test", "ftable", "glm", "lm",
#            "mantelhaen.test", "p.adjust", "quantile", "residuals",
#            "rexp", "rnorm", "rpois")
# importFrom("utils", "str", "write.table")

treeWAS <- function(snps,
                    phen,
                    tree = c("BIONJ", "NJ", "parsimony", "BIONJ*", "NJ*"),
                    phen.type = NULL,
                    n.subs = NULL,
                    n.snps.sim = ncol(snps)*10,
                    chunk.size = ncol(snps),
                    mem.lim = FALSE,
                    test = c("terminal", "simultaneous", "subsequent"),
                    correct.prop = FALSE,
                    snps.reconstruction = "parsimony",
                    snps.sim.reconstruction = "parsimony",
                    phen.reconstruction = "parsimony",
                    na.rm = TRUE,
                    p.value = 0.01,
                    p.value.correct = c("bonf", "fdr", FALSE),
                    p.value.by = c("count", "density"),
                    dist.dna.model = "JC69",
                    plot.tree = TRUE,
                    plot.manhattan = TRUE,
                    plot.null.dist = TRUE,
                    plot.dist = FALSE,
                    plot.null.dist.pairs = "grid",
                    snps.assoc = NULL, # (for manhattan plot)
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
    if("try-error" %in% class(try(match.arg(arg = tree,
                             choices = c("bionj", "nj", "parsimony", "nj*", "bionj*"),
                             several.ok = FALSE), silent=TRUE))){
      tree <- "bionj"
      cat("If tree is not a phylo object, please specify one of the following reconstruction methods:
          'NJ', 'BIONJ', 'parsimony', 'NJ*', 'BIONJ*'. Choosing 'BIONJ' by default.\n")
    }else{
      tree <- match.arg(arg = tree,
                        choices =  c("bionj", "nj", "parsimony", "nj*", "bionj*"),
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

  ##########
  ## CALL ##
  ##########
  ## Get arguments inputed to treeWAS for this run
  ## TO DO: update arguments to get actual values used (report both?).
  args <- NA

  ## store input args to return call at end of fn:
  # args <- match.call()
  args <- mget(names(formals()), sys.frame(sys.nframe()))


  ########################
  ## HANDLE SNPS & PHEN ##
  ########################
  if(!is.matrix(snps)) snps <- as.matrix(snps)
  # x <- snps
  n.snps <- ncol(snps)

  ## convert phenotype to factor
  phen.input <- phen
  phen <- as.factor(phen)
  # y <- phen

  ## set n.ind:
  n.ind <- length(phen)
  inds <- c(1:n.ind)

  #################
  ## HANDLE TREE ##
  #################

  ## RECONSTRUCTED TREE ##
  if("character" %in% class(tree)){

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
    if(length(snps.reconstruction) != 1){
      snps.reconstruction <- "parsimony"
      warning("Cannot use input snps.reconstruction if building tree internally.
              Performing a new snps reconstruction via parsimony.")
    }
    ## And phen.rec...
    if(length(phen.reconstruction) != 1){
      phen.reconstruction <- "parsimony"
      warning("Cannot use input phen.reconstruction if building tree internally.
              Performing a new phen reconstruction via parsimony.")
    }

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
    if(!"phylo" %in% class(tree)) tree <- as.phylo(tree)

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


  ####################################################################
  ################
  ## CHECK PHEN ##
  ################
  ####################################################################
  ## CHECK IF ANY PHEN MISSING:
  if(any(is.na(phen))){
    toRemove <- names(which(is.na(phen)))
    warning(c("The phenotypic variable for individual(s) ", list(toRemove), " is missing.
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
  ## CHECK IF/not BINARY:
  if(length(unique(as.vector(snps[!is.na(snps)]))) != 2){
    stop("snps must be a binary matrix")
  }else{
    ## Convert to numeric, binary (then logical):
    if(!is.numeric(snps) && !is.logical(snps)){
      if(!is.numeric(snps)){
        na.before <- length(which(is.na(snps)))
        if(!all.is.numeric(snps)){
          r.noms <- rownames(snps)
          c.noms <- colnames(snps)
          snps <- matrix(as.numeric(as.factor(snps))-1, nrow=nrow(snps), ncol=ncol(snps))
          snps <- matrix(as.logical(snps), nrow=nrow(snps), ncol=ncol(snps)) ## logical
          rownames(snps) <- r.noms
          colnames(snps) <- c.noms
        }else{
          r.noms <- rownames(snps)
          c.noms <- colnames(snps)
          snps <- matrix(as.numeric(as.character(snps)), nrow=nrow(snps), ncol=ncol(snps))
          snps <- matrix(as.logical(snps), nrow=nrow(snps), ncol=ncol(snps)) ## logical
          rownames(snps) <- r.noms
          colnames(snps) <- c.noms
        }
        na.after <- length(which(is.na(snps)))
        if(na.after > na.before){
          stop("NAs created in converting snps to numeric.")
        }
      }else{
        r.noms <- rownames(snps)
        c.noms <- colnames(snps)
        snps <- matrix(as.logical(snps), nrow=nrow(snps), ncol=ncol(snps)) ## logical
        rownames(snps) <- r.noms
        colnames(snps) <- c.noms
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

  # ## Remove any NON-POLYMORPHIC rows? (i.e., entirely 1 or 0)
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
  ## CHECK FOR NON-BINARY SNPS (?): ## DO THIS OUTSIDE OF THE treeWAS PIPELINE...
  ## ... TOO MANY OPPORTUNITIES FOR ERRORS IF DONE AUTOMATICALLY.
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
    snps <- snps[ord, ]
    ## check:
    if(!identical(as.character(rownames(snps)), as.character(tree$tip.label))){
      stop("Unable to rearrange snps such that rownames(snps) match content and order of tree$tip.label.
            Please check that these match.")
    }
  }

  ## (+ snps.reconstruction)
  if(is.matrix(snps.reconstruction)){
    ## REORDER SNPS/PHEN to match TREE LABELS:
    if(!is.null(rownames(snps.reconstruction))){
      if(!is.null(tree$node.label) && all(rownames(snps.reconstruction) %in% c(tree$tip.label, tree$node.label))){
        ord <- match(c(tree$tip.label, tree$node.label), rownames(snps.reconstruction))
        snps.reconstruction <- snps.reconstruction[ord, ]
      }
      if(!identical(rownames(snps.reconstruction), c(tree$tip.label, tree$node.label))){
        warning("Unable to rearrange snps.reconstruction such that rownames(snps.reconstruction)
                   match c(tree$tip.label, tree$node.label). Performing a new parsimonious reconstruction instead.")
        snps.reconstruction <- "parsimony"
      }
    }
  }


  ## REORDER PHEN TO MATCH TREE$TIP.LABEL
  if(!identical(as.character(names(phen)), as.character(tree$tip.label))){
    ord <- match(tree$tip.label, names(phen))
    phen <- phen[ord]
    ## check:
    if(!identical(as.character(names(phen)), as.character(tree$tip.label))){
      stop("Unable to rearrange phen such that names(phen) match content and order of tree$tip.label.
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
    ## CATEGORICAL, DISCRETE or CONTINUOUS: ##
    ## Convert phen to numeric:
    if(!is.numeric(phen)){
        if(all.is.numeric(phen)){
          phen <- as.numeric(as.character(phen))
        }else{
          phen.type <- "categorical"
          warning("phen has more than 2 levels but is not numeric.
             Setting phen.type to 'categorical'.")
        }
      phen <- as.numeric(as.factor(phen))
    }

    ## Set phen.rec.method: ##
    phen.rec.method <- "continuous"
    ## Get proportion unique:
    prop.u <- length(unique(phen))/length(phen)
    if(!is.null(phen.type)){
      if(phen.type == "discrete"){
        phen.rec.method <- "discrete"
        if(prop.u > 0.5){
          warning("Performing discrete reconstruction, although phen is ", round(prop.u, 2)*100, "% unique.
              Are you sure phen.type is 'discrete' not 'continuous'?\n", sep="")
        }
        if(n.levs < 6){
          warning("Your phen has only ", n.levs, " unique values.
          Are you sure phen.type is 'discrete' not 'categorical'?
          (If values should be treated as levels and not as meaningful numbers, set phen.type to 'categorical').\n", sep="")
        }
      }
      if(phen.type == "categorical"){
        phen.rec.method <- "discrete"
        if(n.levs > 5){
          warning("Association tests with categorical traits work best with fewer levels: your phen has ", n.levs, " levels.
          Are you sure phen.type is 'categorical' not 'discrete'?
          If categorical, consider eliminating some levels by collapsing similar ones or
          removing levels with few individuals.\n", sep="")
        }
      }
    } # end phen.type (categorical/discrete)
    if(is.null(phen.type) || if(!is.null(phen.type)){phen.type == "continuous"}){
      if(n.levs < 6){
        warning("Your phen is being treated as continuous, but it has only ", n.levs, " unique values.
          Set phen.type to 'discrete' or 'categorical' if these better describe your data.\n", sep="")
      }
    }
  } # end phen.type (categorical/discrete/continuous)
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
        phen.reconstruction <- "ml"
        if(phen.rec.method == "discrete" && n.levs == 2){phen.reconstruction <- "parsimony"}
        if(!is.null(phen.type)){if(phen.type == "categorical"){phen.reconstruction <- "parsimony"}}
      }else{
        ## If phen.rec provided and correct length:
        phen.reconstruction.ori <- phen.reconstruction
        levs <- unique(as.vector(unlist(phen.reconstruction)))
        n.levs <- length(levs[!is.na(levs)])
        ## Convert phen.reconstruction to numeric:
        if(!is.numeric(phen.reconstruction)){
          ## convert to numeric if possible:
          if(all.is.numeric(phen.reconstruction)){
            phen.reconstruction <- as.numeric(as.character(phen.reconstruction))
          }else{
            ## convert from factor:
            phen.reconstruction <- as.numeric(as.factor(phen.reconstruction))
          }
        }
        if(length(phen.reconstruction) > 1) names(phen.reconstruction) <- names(phen.reconstruction.ori)

        ##########
        ## REORDER PHEN.REC TO MATCH TREE$TIP.LABEL, NODE.LABEL:
        # if(!identical(as.character(names(phen.reconstruction)[1:length(phen)]), as.character(tree$tip.label))){
        #   ord <- match(tree$tip.label, names(phen.reconstruction)[1:length(phen)])
        #   phen.reconstruction <- phen.reconstruction[c(ord, (length(ord)+1):length(phen.reconstruction))]
        # }
        if(!identical(as.character(names(phen.reconstruction)), as.character(c(tree$tip.label, tree$node.label)))){
          if(!is.null(tree$node.label) && all(names(phen.reconstruction) %in% c(tree$tip.label, tree$node.label))){
            ord <- match(c(tree$tip.label, tree$node.label), names(phen.reconstruction))
            phen.reconstruction <- phen.reconstruction[ord]
          }
          ## check again:
          if(!identical(as.character(names(phen.reconstruction)), as.character(c(tree$tip.label, tree$node.label)))){
            warning("Unable to rearrange phen.reconstruction such that names(phen.reconstruction)
               match c(tree$tip.label, tree$node.label). Performing a new reconstruction instead.")
            phen.reconstruction <- "ml"
            if(phen.rec.method == "discrete" && n.levs == 2){phen.reconstruction <- "parsimony"}
            if(!is.null(phen.type)){if(phen.type == "categorical"){phen.reconstruction <- "parsimony"}}
          }
        }
        ##########
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

    if(!is.null(phen.rec)){

      ## PLOT TERMINAL PHEN + RECONSTRUCTED PHEN:
      plot_phen(phen.nodes = phen.rec, tree = tree,
                main.title = "Phylogenetic tree", align.tip.label = T, RTL = F)

    }else{

      ## PLOT TERMINAL PHEN ONLY:

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
    }

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
    n.subs <- get.fitch.n.mts(x=snps, tree=tree)
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


  #########################
  ## Handle memory limit ##
  #########################
  ml.est <- FALSE
  if(!is.null(mem.lim)){
    ## If mem.lim = FALSE ##
    ## --> work with chunk.size set by user.
    ## * May run out of memory
    ## * But, allows user to set a chunk.size that works for them,
    ##   eg. if mem.lim fails/makes poor choice of chunk.size given their PC/data.
    if(mem.lim == FALSE){
      mem.lim <- NULL
    }else{
      ## If mem.lim = TRUE ##
      ## --> automatically determine optimal max chunk.size, given available memory,
      ## regardless of input chunk.size value.
      if(mem.lim == TRUE){
        mem.lim <- NULL
        ## run memfree, or warn if it fails:
        if("try-error" %in% class(try(memfree(), silent=TRUE))){
          warning("Unable to determine amount of available memory.
                  Set mem.lim or chunk.size by hand.")
        }else{
          ## set mem.lim w memfree:
          mem.lim <- memfree()
          ## round down to be conservative and leave a little space:
          if(mem.lim > 1) mem.lim <- floor(mem.lim)
          ml.est <- TRUE
        }
      }else{
        ## If mem.lim = numeric ##
        ## --> Skip memfree estimate, use input mem.lim instead,
        ## but estimate n.chunks, chunk.size, regardless of input chunk.size value.
        ## * May run out of memory
        ## * But, allows user to set mem.lim in case memfree() fails for their OS.
        if(!is.numeric(mem.lim)){
          warning("mem.lim was not a number or a logical. Ignoring mem.lim and using default chunk.size.")
          mem.lim <- NULL
        }else{
          if(length(mem.lim) > 1){
            warning("mem.lim must be of length 1 if numeric. Ignoring mem.lim and using default chunk.size.")
            mem.lim <- NULL
          }
        }
      }
    } # end mem.lim check
    ## If mem.lim has been estimated or input...
    if(!is.null(mem.lim)){
      #######################################
      ## Update chunk.size, given mem.lim: ##
      #######################################

      ### TO DO: UPDATE MEMREQ ESTIMATE ON THE BASIS OF UNIQUE SNPS COLUMNS INSTEAD OF TOTAL...   %%   <---  (!)

      ## Check memory occuppied by snps:
      # require(pryr)
      memreq <- as.numeric(object_size(snps))/1000000000 # bytes --> GB

      ## Get n.snps.sim factor:
      fac <- n.snps.sim/ncol(snps)
      ## Update memreq:
      memreq <- memreq*fac

      ## Multiply by second factor to approximate total mem req'd by treeWAS:
      fac2 <- 70 ## (!!) NOTE---THIS IS AN ESTIMATE---IT MAY NEED TO BE IMPROVED... (!!)
      ## Update memreq:
      memreq <- memreq*fac2

      ## Get n.chunks (round up):
      nc <- ceiling(memreq/mem.lim)

      ## Get chunk.size (round up to nearest 1):
      if(nc > 1){
        chunk.size <- ceiling(ncol(snps)/nc)
      }else{
        chunk.size <- ncol(snps)
      }

      ## Print notice:
      if(ml.est == TRUE){
        cat(c("Updated chunk.size:", chunk.size,
              "\nNumber of chunks:", nc,
              "\nEstimated memory limit:", mem.lim, "GB",
              "\nEstimated memory required:", round(memreq, 2), "GB\n"))
      }else{
        cat(c("Updated chunk.size:", chunk.size,
              "\nNumber of chunks:", nc,
              "\nInput memory limit:", mem.lim, "GB",
              "\nEstimated memory required:", round(memreq, 2), "GB\n"))
      }
    }
  } # end mem.lim



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
      ## (What if chunk size is really at max memory??)
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

  ## MEMORY ISSUES:
  ## (1) CHECK MEMORY LIMITS AUTOMATICALLY & CHANGE CHUNK.SIZE DEFAULT?
  ## (2) ENABLE PARALLEL RUNS FOR CHUNK.SIZE FOR LOOP?
  ##### (+) Automatically check n.cores available, and if chunk.size == ncol(snps),
  ##### reduce chunk.size to ncol(snps)/n.cores...
  # require(pryr)
  # ## Check memory used in current R session (only):
  # mem_used()
  # ## Check memory occuppied by a given object:
  # object_size(snps)
  # memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",  intern=TRUE))
  # memfree/1000000 # convert from KB to GB (??)
  ## CHECK UNITS OF memfree
  ## + CHECK IF SYSTEM (and/or) COMMAND PATH WORK ON WINDOWS, MAC, ETC...

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


    if(length(CHUNKS) > 1){
      ## Print update notice:
      cat("########### Running chunk ", i, " of ", length(CHUNKS), "###########\n")
    }

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
          cat("Note: Updating n.snps.sim to match the number of loci after data cleaning.
              Input:", n.snps.sim.ori, " -->
              Updated:", n.snps.sim, "\n")

        }else{
          Nx <- c(2:200)
          if(any(n.snps*Nx == n.snps.sim)){
            n.snps.sim.ori <- n.snps.sim
            Nx <- Nx[which(Nx == (n.snps.sim/n.snps))]
            n.snps.sim <- ncol(snps)*Nx

            ## Print update notice:
            cat("Note: Updating n.snps.sim to match ", Nx, "x the number of loci after data cleaning.
                Input: ", n.snps.sim.ori, " -->
                Updated: ", n.snps.sim, "\n", sep="")

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

    snps.unique <- snps.index <- snps.sim.unique <- snps.sim.index <- snps.rec.index <- NULL

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
    if("list" %in% class(snps.sim)){
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


    ##########################
    ## Add NAs to snps.sim? ##  ### ###  ### ###  ### ###  ### ###  ### ###  ###
    ##########################
    ## Try distributing NAs proportionately (as with n.subs),
    ## but randomly* across each site
    ## (*could do 3-state (c(0,1,NA) reconstruction, in case neighbouring inds have NAs at a given site,
    ## which could affect scores; but, would have to modify snps.sim not to simulate NA ancestors)
    ## Get number of NAs per snps locus column:
    NA.tab <- sapply(c(1:ncol(snps)), function(e) length(which(is.na(snps[,e]))))
    NA.tab <- rep(NA.tab, round(ncol(snps.sim)/ncol(snps)))
    ss.nr <- nrow(snps.sim)
    ## Sample indices for each column
    col_inds <- lapply(seq_along(NA.tab), function(x) sample(ss.nr, NA.tab[x], replace = FALSE))
    snps.sim[unlist(col_inds)] <- NA

    ### ###  ### ###  ### ###  ### ###  ### ###  ### ###### ###  ### ###  ### ##

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

        ## CHECK row N:
        if(nrow(snps.reconstruction) != max(tree$edge[,2])){
          warning("The number of rows in the provided snps.reconstruction is not equal to the
                  total number of nodes in the tree. Performing a new parsimonious reconstruction instead.\n")
          snps.reconstruction <- snps.sim.reconstruction <- "parsimony"
        }

        ## CHECK column N:
        if(ncol(snps.reconstruction) != ncol(snps)){
          ## If user-provided:
          snps.rec <- snps.reconstruction
          ## Get UNIQUE snps.reconstruction
          snps.rec.complete <- snps.rec
          temp <- get.unique.matrix(snps.rec, MARGIN=2)
          snps.rec <- temp$unique.data
          snps.rec.index <- temp$index
          snps.reconstruction <- snps.rec

          ## CHECK index:
          if(!identical(snps.rec.index, snps.index)){
            warning("Careful-- snps and snps.rec do not reduce to the same set of unique columns.\n")
          }

        }
      }

      ## If not user-provided or checks failed, reconstruct ancestral states:
      if(!is.matrix(snps.reconstruction)){
        snps.rec <- asr(var = snps, tree = tree, type = snps.reconstruction, unique.cols = TRUE)
      }

      if(is.null(snps.rec.index)){
        snps.rec.index <- snps.index
      }

      #################################
      ## Reconstruct SIMULATED SNPs: ##
      #################################
      snps.sim.rec <- asr(var = snps.sim, tree = tree, type = snps.sim.reconstruction, unique.cols = TRUE)

     } # end reconstruction for tests 2 & 3


    ## !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ### !!! ###

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

    ## Add categorical flag:
    ctg <- FALSE
    if(!is.null(phen.type)){
      if(phen.type == "categorical"){
        ctg <- TRUE
      }
    }

    ## Run get.assoc.scores for each assoc.test
    # system.time(
    for(t in 1:length(TEST)){
      assoc.scores <- get.assoc.scores(snps = snps,
                                       snps.sim = snps.sim,
                                       phen = phen,
                                       tree = tree,
                                       test = TEST[[t]],
                                       correct.prop = correct.prop,
                                       categorical = ctg,
                                       snps.reconstruction = snps.rec,
                                       snps.sim.reconstruction = snps.sim.rec,
                                       phen.reconstruction = phen.rec,
                                       unique.cols = TRUE)
      ## EXPAND CORR.DAT & CORR.SIM:
      cd <- assoc.scores$corr.dat
      cs <- assoc.scores$corr.sim

      if(TEST[[t]] %in% c("simultaneous", "subsequent")){
        if(ncol(snps.rec) != ncol(snps)){
          cd <- cd[snps.rec.index]
        }else{
          cd <- cd[snps.index]
        }
      }else{
        cd <- cd[snps.index]
      }
      cs <- cs[snps.sim.index]
      names(cs) <- colnames(snps.sim.complete)

      names(cd) <- colnames(snps.complete)


      ## STORE DATA FOR EACH TEST:
      CORR.DAT[[i]][[t]] <- cd
      CORR.SIM[[i]][[t]] <- cs
      rm(assoc.scores)
    } # end for (t) loop
    # )

    ## EXPAND DATA BEFORE STORING FOR THIS CHUNK:
    snps.sim <- snps.sim.complete
    if(!is.null(snps.rec)){
      snps.rec <- snps.rec[,snps.rec.index]
      colnames(snps.rec) <- colnames(snps.complete)
    }
    if(!is.null(snps.sim.rec)){
      snps.sim.rec <- snps.sim.rec[,snps.sim.index]
      colnames(snps.sim.rec) <- colnames(snps.sim.complete)
    }

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
    ## add names to p.vals:
    names(p.vals) <- names(corr.dat)
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

      if(TEST[[i]] == "fisher"){
        ylab <- NULL # --> (uncorrected or) -log10 p-value
        log10 <- TRUE
        ttl <- paste("\n \n(", TEST[[i]], "test)")
      }else{
        ylab <- paste(TEST[[i]], "score", sep=" ")
        log10 <- FALSE
        ttl <- paste("\n \n(", TEST[[i]], "score)")
      }

      ## plot:
      manhattan.plot(p.vals = abs(corr.dat),
                     col = "funky",
                     transp = 0.25,
                     sig.thresh = sig.thresh,
                     thresh.col="red",
                     snps.assoc = NULL,
                     snps.assoc.col = "red",
                     jitter.amount = 0.00001,
                     min.p = NULL,
                     log10=log10,
                     ylab=ylab)
      ## Add subtitle:
      title(ttl, cex.main=0.9)
    } # end plot manhattan


    ####################################
    ## 4) (B,C) Plot the distribution ##
    ####################################

    if(plot.null.dist == TRUE){
      if(TEST[[i]] == "fisher"){
        ## Generate one histogram per test:
        plot_sig_snps(corr.dat = -log10(abs(corr.dat)),
                      corr.sim = -log10(abs(corr.sim)),
                      corr.sim.subset = NULL,
                      sig.corrs = -log10(abs(corr.dat[sig.snps])),
                      sig.snps = sig.snps.names,
                      sig.thresh = -log10(sig.thresh),
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
      }else{
        ## Generate one histogram per test:
        plot_sig_snps(corr.dat = abs(corr.dat),
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
      }
    } # end plot.null.dist

    if(plot.dist == TRUE){
      if(TEST[[i]] == "fisher"){
        ## Generate one histogram per test:
        plot_sig_snps(corr.dat = -log10(abs(corr.dat)),
                      corr.sim = -log10(abs(corr.sim)),
                      corr.sim.subset = NULL,
                      sig.corrs = -log10(abs(corr.dat[sig.snps])),
                      sig.snps = sig.snps.names,
                      sig.thresh = -log10(sig.thresh),
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
      }else{
        ## Generate one histogram per test:
        plot_sig_snps(corr.dat = abs(corr.dat),
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
      }
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
        ## BINARY
        ## store ind names:
        noms <- names(phen)
        ## If binary, convert phen to 0/1:
        phen <- as.numeric(as.factor(phen))
        phen <- rescale(phen, to=c(0,1)) # require(scales)
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
      }else{

        if(is.null(phen.type) || if(!is.null(phen.type)){phen.type != "categorical"}){
          ## CONTINUOUS
          ## Return df without G1P1 columns
          df <- data.frame(sig.snps,
                           sig.p.vals,
                           sig.corrs)
          names(df) <- c("SNP.locus",
                         "p.value",
                         "score")
        }else{
          ## CATEGORICAL
          if(!is.null(phen.type)){
            if(phen.type == "categorical"){
              phen.output <- phen
              phen <- phen.input
              if(length(toKeep) == 1){
                tb <- t(table(phen, snps.toKeep))
                cn <- paste(c("G0P","G1P"), rep(colnames(tb), each=2), sep="")
                tb <- t(as.matrix(as.vector(unlist(tb))))
                colnames(tb) <- cn
              }else{
                tb <- t(table(phen, snps.toKeep[,1]))
                cn <- paste(c("G0P","G1P"), rep(colnames(tb), each=2), sep="")
                tb <- t(sapply(c(1:ncol(snps.toKeep)),
                               function(e) as.vector(unlist(t(table(phen, snps.toKeep[,e]))))))
                colnames(tb) <- cn
              }
              ## return df with k*2 extra columns
              df <- data.frame(sig.snps,
                               sig.p.vals,
                               sig.corrs)
              names(df) <- c("SNP.locus",
                             "p.value",
                             "score")
              df <- cbind(df, tb)
            }
          }
        }
        ## NOTE: Could return df with rownames=# rather than rownames=sig.snps.names.. ####   (??)   ####
      }
    }else{
      df <- "No significant SNPs found."
    }

    ## 0 p.vals
    #   min.p <- paste("p-values listed as 0 are <",
    #                  1/length(corr.sim), sep=" ")
    min.p <- 1/length(corr.sim)
    names(min.p) <- c("p-values listed as 0 are less than:")



    ## Get output:
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



  ###########################
  ## Get COMBINED results: ##
  ###########################
  ## get unique(sig.snps) for terminal, simultaneous, subsequent results combined.
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


  ##############################################################################
  #######################################
  ## PAIRWISE categorical assoc tests: ##
  #######################################
  ## Post-hoc pairwise assoc tests* btw each sig. SNP & each comparison of 2 levels of the phenotype.
  ## *Score 1 (phi), Score 2 (by phen pair), Score 3 (phi), chi-squared p-values.
  if(!is.null(phen.type)){
    if(phen.type == "categorical"){
      snps.sig <- sort(as.numeric(RES[[1]]$treeWAS.combined))

      phen.output <- phen
      phen <- phen.input
      ## ensure phen.rec levels match input levels:
      if(!all(unique(phen.rec[!is.na(phen.rec)]) %in% unique(phen[!is.na(phen)]))){
        pr <- phen.rec
        pn <- as.numeric(as.factor(phen))
        prn <- as.numeric(as.factor(phen.rec))
        levs <- unique(pn[!is.na(pn)])
        if(all(unique(prn[!is.na(prn)]) %in% levs)){
          for(l in 1:length(levs)){
            lev.c <- as.character(phen[which(pn == levs[l])[1]])
            prn[which(prn == levs[l])] <- lev.c
          }
          pr <- prn
          names(pr) <- names(phen.rec)
          phen.rec <- pr
        }
        #
      }

      ###################################
      ## PAIRWISE SCORES on SIG.SNPS:  ##
      ###################################
      if(length(snps.sig) > 0){

        print(paste("Running pairwise tests @", Sys.time()))

        ## Get phen pairs:
        pairs <- t(combn(unique(phen.rec[!is.na(phen.rec)]), m=2))

        ## Print plots in grid?
        if(!is.null(plot.null.dist.pairs)){
          if(plot.null.dist.pairs == "grid"){
            mfrow.ori <- par()$mfrow
            par(mfrow=c(nrow(pairs), length(TEST)))
          }
        }

        SIG.LIST <- res.pairs <- PW.DF <- list()
        #############################################
        ## FOR (p) LOOP:
        for(p in 1:nrow(pairs)){

          ## Get phen diffs: ##
          ph <- phen
          ph[which(!ph %in% pairs[p,])] <- NA
          ph <- as.numeric(as.factor(as.character(ph)))-1

          pr <- phen.rec
          pr[which(!pr %in% pairs[p,])] <- NA
          pr <- as.numeric(as.factor(as.character(pr)))-1

          r.toKeep <- which(!is.na(ph))
          r.toKeep.rec <- which(!is.na(pr))

          CORR.DAT <- CORR.SIM <- sig.list <- SIG.LIST[[p]] <- PW.DF[[p]] <- list()
          ## FOR (t) LOOP:
          for(t in 1:length(TEST)){

            if(TEST[[t]] == "simultaneous"){
              r.toKeep <- c(1:length(ph))
              r.toKeep.rec <- c(1:length(pr))
              }

            assoc.scores <- get.assoc.scores(snps = snps[r.toKeep, snps.sig],
                                             snps.sim = snps.sim[r.toKeep, ],
                                             phen = ph[r.toKeep],
                                             tree = tree,
                                             test = TEST[[t]],
                                             correct.prop = TRUE,
                                             categorical = FALSE,
                                             snps.reconstruction = snps.rec[r.toKeep.rec, snps.sig],
                                             snps.sim.reconstruction = snps.sim.rec[r.toKeep.rec, ],
                                             phen.reconstruction = pr[r.toKeep.rec],
                                             unique.cols = TRUE)

            CORR.DAT[[t]] <- assoc.scores$corr.dat
            CORR.DAT[[t]][which(is.na(CORR.DAT[[t]]))] <- 0
            CORR.SIM[[t]] <- assoc.scores$corr.sim
            CORR.SIM[[t]][which(is.na(CORR.SIM[[t]]))] <- 0

            sig.list[[t]] <- get.sig.snps(corr.dat = CORR.DAT[[t]],
                                          corr.sim = CORR.SIM[[t]],
                                          snps.names = colnames(snps)[snps.sig],
                                          test = TEST[[t]],
                                          p.value = p.value,
                                          p.value.correct = p.value.correct,
                                          p.value.by = p.value.by)

            SIG.LIST[[p]][[t]] <- sig.list[[t]]

            pw.df <- data.frame(sig.list[[t]]$corr.dat, sig.list[[t]]$p.vals)
            names(pw.df) <- c(TEST[[t]], paste("p", TEST[[t]], sep="."))
            PW.DF[[p]][[t]] <- pw.df

            corr.dat <- sig.list[[t]]$corr.dat
            corr.sim <- sig.list[[t]]$corr.sim
            sig.thresh <- sig.list[[t]]$sig.thresh
            sig.snps <- sig.list[[t]]$sig.snps
            sig.snps.names <- sig.list[[t]]$sig.snps.names

            ## (*) ADD ARG:
            if(plot.null.dist.pairs == TRUE || plot.null.dist.pairs == "grid"){

              ## Generate one histogram per test:
              plot_sig_snps(corr.dat = abs(corr.dat),
                            corr.sim = abs(corr.sim),
                            corr.sim.subset = NULL,
                            sig.corrs = abs(corr.dat[sig.snps]),
                            sig.snps = sig.snps.names,
                            sig.thresh = sig.thresh,
                            test = TEST[[t]],
                            sig.snps.col = "black",
                            hist.col = rgb(0.5,0,1,0.5), # purple
                            hist.subset.col = rgb(1,0,0,0.5),
                            thresh.col = "red",
                            snps.assoc = NULL,
                            snps.assoc.col = "blue",
                            bg = "lightgrey",
                            grid = TRUE,
                            freq = FALSE,
                            plot.null.dist = TRUE,
                            plot.dist = FALSE,
                            main.title = paste("Null distribution \n(", TEST[[t]], ": ",
                                               paste0(pairs[p,], collapse="_"), ")"
                                               , sep=""))
            } # end plot_sig_snps

          } # end for (t) loop
          #############################################
          names(SIG.LIST[[p]]) <- test

          ## get counts of n.pw.tests sig per snps.sig for phen.pair p:
          pw.sig.snps <- as.vector(unlist(sapply(c(1:length(TEST)), function(e) sig.list[[e]]$sig.snps.names)))
          if(length(pw.sig.snps) > 0){
            pw.sig <- as.vector(table(factor(pw.sig.snps, levels=as.character(snps.sig))))
          }else{
            pw.sig <- rep(0, length(snps.sig))
          }

          ## Get score & p-values for t tests for this phen.pair:
          pw.df <- do.call(cbind, PW.DF[[p]])

          ## Get G0P1 table for this phen.pair:
          snps.toKeep <- snps[r.toKeep, snps.sig]
          if(length(snps.sig) == 1){
            tb <- t(table(phen[r.toKeep], snps.toKeep))
            cn <- paste(c("G0P","G1P"), rep(colnames(tb), each=2), sep="")
            tb <- t(as.matrix(as.vector(unlist(tb))))
            colnames(tb) <- cn
          }else{
            tb <- t(table(phen[r.toKeep], snps.toKeep[,1]))
            cn <- paste(c("G0P","G1P"), rep(colnames(tb), each=2), sep="")
            tb <- t(sapply(c(1:ncol(snps.toKeep)),
                           function(e) as.vector(unlist(t(table(phen[r.toKeep], snps.toKeep[,e]))))))
            colnames(tb) <- cn
          }


          res.pairs[[p]] <- cbind(data.frame("SNP.locus" = snps.sig,
                                       "pw.sig" = pw.sig),
                                       pw.df, tb)

        } # end for (p) loop
        #############################################

        ## Name pw list by phen pairs:
        names(SIG.LIST) <- sapply(c(1:nrow(pairs)), function(p) paste0(pairs[p,], collapse="_"))
        names(res.pairs) <- sapply(c(1:nrow(pairs)), function(p) paste0(pairs[p,], collapse="_"))

        ## Print plots in grid?
        if(!is.null(plot.null.dist.pairs)){
          if(plot.null.dist.pairs == "grid"){
            par(mfrow=mfrow.ori)
          }
        }

        ## ADD NAMES to res.pairs, sig.list by phen : pair, test...

        ## Enable manhattan plots etc for pairwise tests...
        ## Build summary plot (tree+grid) fn?



        ## Return res.pairs sig tables & SIG.LIST tests data with output:
        PW <- list("pairwise.combined" = res.pairs, "pairwise.tests" = SIG.LIST)
      }else{
        ## No sig. snps:
        PW <- list("pairwise.combined" = NULL, "pairwise.tests" = NULL)
      }

      ## Add as last item of RES list.
      RES[[(length(RES)+1)]] <- PW

      phen <- phen.output
    }
  } # end categorical sig.snps
  ##############################################################################

  #######################
  ## end filename.plot ##
  #######################
  if(!is.null(filename.plot)){
    dev.off()
  }

  phen <- phen.curr

  ## Also return data (in the form we were working with):
  dat <- list("snps" = snps,
              "snps.reconstruction" = snps.rec,
              "snps.sim" = snps.sim,
              "snps.sim.reconstruction" = snps.sim.rec,
              "phen" = phen,
              "phen.reconstruction" = phen.rec,
              "tree" = tree,
              "n.subs" = n.subs,
              "call" = args)

  ## Assign data to last element of RES:
  RES[[(length(RES)+1)]] <- dat

  ## name elements of RES:
  if(is.null(phen.type) || if(!is.null(phen.type)){phen.type != "categorical"}){
    names(RES) <- c("treeWAS.combined", test, "dat")
  }else{
    names(RES) <- c("treeWAS.combined", test, "pairwise", "dat")
  }

  results <- RES

  class(results) <- "treeWAS"

  return(results)

} # end treeWAS





##########################################################################
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
