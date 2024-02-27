###################
## treeWAS DATA: ##
###################

################################################################################################################################################
################################################################################################################################################
########################################################################
###################
## DOCUMENTATION ##
###################
#' Nsubs per site with no recombination (R = 0; M = 0.01).
#'
#' This vector contains a homoplasy distribution,
#' representing the relative number of substitutions per site
#' that occurred along a phylogenetic tree when evolution was
#' simulated with a mutation rate of M = 0.01 and a recombination rate of R = 0
#' (\emph{r/m} = 0).
#'
#' A per-site mutation rate of M = 0.01
#' indicates that each site, on average, undergoes 0.01
#' substitutions due to mutation along the phylogenetic tree.
#' A per-site recombination rate of R = 0
#' indicates that no recombination occurred.
#'
#' Each element of the vector indicates the number of genetic loci
#' that have undergone the number of substitutions indicated by the name of that element (Nsub = i).
#'
#' If visualised as a bar plot (with \code{barplot(dist_0)}),
#' one would see that the Nsub distribution is arranged as if it were the counts of a histogram
#' with index names along the x-axis, corresponding to Nsub (the number of substitutions per site),
#' and cell counts along the y-axis, showing Nloci (the number of genetic sites undergoing Nsub=i substitutions along the tree).
#'
#' @docType data
#' @name dist_0
#' @usage data(dist_0)
#' @format A named vector of length 4.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL

########################################################################
###################
## DOCUMENTATION ##
###################
#' Nsubs per site with limited recombination (R = 0.01; M = 0.01).
#'
#' This vector contains a homoplasy distribution,
#' representing the relative number of substitutions per site
#' that occurred along a phylogenetic tree when evolution was
#' simulated with a mutation rate of M = 0.01 and a recombination rate of R = 0.01
#' (\emph{r/m} = 1).
#'
#' A per-site mutation rate of M = 0.01
#' indicates that each site, on average, undergoes 0.01
#' substitutions due to mutation along the phylogenetic tree.
#' A per-site recombination rate of R = 0.01
#' indicates that each site, on average, undergoes 0.01
#' substitutions per site due to within-species recombination.
#'
#' Each element of the vector indicates the number of genetic loci
#' that have undergone the number of substitutions indicated by the name of that element (Nsub = i).
#'
#' If visualised as a bar plot (with \code{barplot(dist_0.01)}),
#' one would see that the Nsub distribution is arranged as if it were the counts of a histogram
#' with index names along the x-axis, corresponding to Nsub (the number of substitutions per site),
#' and cell counts along the y-axis, showing Nloci (the number of genetic sites undergoing Nsub=i substitutions along the tree).
#'
#' @docType data
#' @name dist_0.01
#' @usage data(dist_0.01)
#' @format A named vector of length 15.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL

########################################################################
###################
## DOCUMENTATION ##
###################
#' Nsubs per site with recombination (R = 0.05; M = 0.01).
#'
#' This vector contains a homoplasy distribution,
#' representing the relative number of substitutions per site
#' that occurred along a phylogenetic tree when evolution was
#' simulated with a mutation rate of M = 0.01 and a recombination rate of R = 0.05
#' (\emph{r/m} = 5).
#'
#' A per-site mutation rate of M = 0.01
#' indicates that each site, on average, undergoes 0.01
#' substitutions due to mutation along the phylogenetic tree.
#' A per-site recombination rate of R = 0.05
#' indicates that each site, on average, undergoes 0.05
#' substitutions per site due to within-species recombination.
#'
#' Each element of the vector indicates the number of genetic loci
#' that have undergone the number of substitutions indicated by the name of that element (Nsub = i).
#'
#' If visualised as a bar plot (with \code{barplot(dist_0.05)}),
#' one would see that the Nsub distribution is arranged as if it were the counts of a histogram
#' with index names along the x-axis, corresponding to Nsub (the number of substitutions per site),
#' and cell counts along the y-axis, showing Nloci (the number of genetic sites undergoing Nsub=i substitutions along the tree).
#'
#' @docType data
#' @name dist_0.05
#' @usage data(dist_0.05)
#' @format A named vector of length 26.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL


########################################################################
###################
## DOCUMENTATION ##
###################
#' Nsubs per site with recombination (R = 0.1; M = 0.01).
#'
#' This vector contains a homoplasy distribution,
#' representing the relative number of substitutions per site
#' that occurred along a phylogenetic tree when evolution was
#' simulated with a mutation rate of M = 0.01 and a recombination rate of R = 0.1
#' (\emph{r/m} = 10).
#'
#' A per-site mutation rate of M = 0.01
#' indicates that each site, on average, undergoes 0.01
#' substitutions due to mutation along the phylogenetic tree.
#' A per-site recombination rate of R = 0.1
#' indicates that each site, on average, undergoes 0.1
#' substitutions per site due to within-species recombination.
#'
#' Each element of the vector indicates the number of genetic loci
#' that have undergone the number of substitutions indicated by the name of that element (Nsub = i).
#'
#' If visualised as a bar plot (with \code{barplot(dist_0.1)}),
#' one would see that the Nsub distribution is arranged as if it were the counts of a histogram
#' with index names along the x-axis, corresponding to Nsub (the number of substitutions per site),
#' and cell counts along the y-axis, showing Nloci (the number of genetic sites undergoing Nsub=i substitutions along the tree).
#'
#' @docType data
#' @name dist_0.1
#' @usage data(dist_0.1)
#' @format A named vector of length 31.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL

########################################################################
###################
## DOCUMENTATION ##
###################
#' Nsubs per site with considerable recombination (R = 0.2; M = 0.01).
#'
#' This vector contains a homoplasy distribution,
#' representing the relative number of substitutions per site
#' that occurred along a phylogenetic tree when evolution was
#' simulated with a mutation rate of M = 0.01 and a recombination rate of R = 0.2
#' (\emph{r/m} = 20).
#'
#' A per-site mutation rate of M = 0.01
#' indicates that each site, on average, undergoes 0.01
#' substitutions due to mutation along the phylogenetic tree.
#' A per-site recombination rate of R = 0.2
#' indicates that each site, on average, undergoes 0.2
#' substitutions per site due to within-species recombination.
#'
#' Each element of the vector indicates the number of genetic loci
#' that have undergone the number of substitutions indicated by the name of that element (Nsub = i).
#'
#' If visualised as a bar plot (with \code{barplot(dist_0.2)}),
#' one would see that the Nsub distribution is arranged as if it were the counts of a histogram
#' with index names along the x-axis, corresponding to Nsub (the number of substitutions per site),
#' and cell counts along the y-axis, showing Nloci (the number of genetic sites undergoing Nsub=i substitutions along the tree).
#'
#' @docType data
#' @name dist_0.2
#' @usage data(dist_0.2)
#' @format A named vector of length 30.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL

########################################################################
###################
## DOCUMENTATION ##
###################
#' Nsubs per site with considerable recombination (R = 0.25; M = 0.01).
#'
#' This vector contains a homoplasy distribution,
#' representing the relative number of substitutions per site
#' that occurred along a phylogenetic tree when evolution was
#' simulated with a mutation rate of M = 0.01 and a recombination rate of R = 0.25
#' (\emph{r/m} = 25).
#'
#' A per-site mutation rate of M = 0.01
#' indicates that each site, on average, undergoes 0.01
#' substitutions due to mutation along the phylogenetic tree.
#' A per-site recombination rate of R = 0.25
#' indicates that each site, on average, undergoes 0.25
#' substitutions per site due to within-species recombination.
#'
#' Each element of the vector indicates the number of genetic loci
#' that have undergone the number of substitutions indicated by the name of that element (Nsub = i).
#'
#' If visualised as a bar plot (with \code{barplot(dist_0.25)}),
#' one would see that the Nsub distribution is arranged as if it were the counts of a histogram
#' with index names along the x-axis, corresponding to Nsub (the number of substitutions per site),
#' and cell counts along the y-axis, showing Nloci (the number of genetic sites undergoing Nsub=i substitutions along the tree).
#'
#' @docType data
#' @name dist_0.25
#' @usage data(dist_0.25)
#' @format A named vector of length 34.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL

################################################################################################################################################
################################################################################################################################################
########################################################################
###################
## DOCUMENTATION ##
###################
#' A binary phenotype.
#'
#' This vector specifies the phenotype of each individual.
#' In this case, the phenotype is a binary variable.
#' Because the phenotypic vector is encoded as a factor
#' with two possible phenotypic states, "A" and "B",
#' which may be represented by the numeric values 1 and 2 (as in \code{str(phen)}).
#'
#'
#' Each individual in the sample is represented by a unique identifier (name)
#' which corresponds to the name of one element of the phenotypic vector.
#' Each element of the phenotypic vector gives the phenotypic value of the named individual.
#'
#'
#' @docType data
#' @name phen
#' @usage data(phen)
#' @format A named vector of length 100.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL

########################################################################
###################
## DOCUMENTATION ##
###################
#' The ancestral state reconstruction of a binary phenotype.
#'
#' This vector contains the terminal and ancestral states of a binary phenotypic variable (\code{data(phen)}).
#' The observed phenotypic states of sampled individuals
#' (i.e., those represented at the terminal nodes of a phylogenetic tree)
#' are presented first, in elements 1:N (here 1:100).
#' The unobserved ancestral states of the phenotype at internal nodes have been
#' inferred via ancestral state reconstruction, using \code{asr(phen, tree)}.
#'
#' Like the original phenotypic vector (\code{data(phen)}),
#' \code{phen.reconstruction} is a binary variable that is encoded as a factor
#' with two possible phenotypic states, "A" and "B",
#' which may be represented by the numeric values 1 and 2 (as in \code{str(phen.reconstruction)}).
#'
#' Each individual in the sample is represented by a unique identifier (name)
#' which corresponds to the name of one element of the phenotypic vector.
#' (Internal node names have been generated during ancestral state reconstruction.)
#' Each element of the phenotypic vector gives the phenotypic value of the named individual.
#'
#'
#' @docType data
#' @name phen.reconstruction
#' @usage data(phen.reconstruction)
#' @format A named vector of length 199.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL

########################################################################
###################
## DOCUMENTATION ##
###################
#' A continuous phenotype.
#'
#' This vector specifies the phenotype of each individual.
#' In this case, the phenotype is a continuous numeric value.
#'
#' Each individual in the sample is represented by a unique identifier (name)
#' which corresponds to the name of one element of the phenotypic vector.
#' Each element of the phenotypic vector gives the phenotypic value of the named individual.
#'
#' Note that, due to some skew in the distribution of this continuous variable,
#' it may be useful to transform the phenotype by rank prior to analysis by treeWAS,
#' as in \code{data(phen.cont.rank)} (see the treeWAS vignette).
#' % (see \code{vignette("treeWAS")}).
#'
#'
#' @docType data
#' @name phen.cont
#' @usage data(phen.cont)
#' @format A named vector of length 533.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL

########################################################################
###################
## DOCUMENTATION ##
###################
#' A rank-transformed continuous phenotype.
#'
#' This vector specifies the phenotype of each individual.
#' In this case, the phenotype is a rank, that has been derived by
#' rank-ordering the elements of the original continuous phenotype (\code{data(phen.cont)})
#' from lowest to highest.
#' Transforming by rank prior to analysis by treeWAS can be useful
#' for continuous phenotypic variables that are highly skewed or contain significant outliers
#' (see the treeWAS vignette).
#' % (see \code{vignette("treeWAS")}).
#'
#' Each individual in the sample is represented by a unique identifier (name)
#' which corresponds to the name of one element of the phenotypic vector.
#' Each element of the phenotypic vector gives the phenotypic value of the named individual.
#'
#'
#' @docType data
#' @name phen.cont.rank
#' @usage data(phen.cont.rank)
#' @format A named vector of length 533.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL


########################################################################
###################
## DOCUMENTATION ##
###################
#' Phenotypic tree-colouring schemes.
#'
#' A list containing the colour values that \code{plot_phen} generates to represent
#' the states and substitutions of the phenotypic variable (\code{data(phen)})
#' along the phylogenetic tree (\code{data(tree)}), with \code{plot_phen(tree, phen.nodes=phen)}.
#' You are unlikely to have to interact with this list,
#' as the colours are automatically plotted by the \code{plot_phen} function.
#'
#' The five elements of this list give the colour schemes used to indicate the phenotypic state at:
#' edge.labels, edges, all.nodes, internal.nodes, and tip.labels.
#'
#'
#' @docType data
#' @name phen.plot.col
#' @usage data(phen.plot.col)
#' @format A list of length 5.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL



########################################################################
###################
## DOCUMENTATION ##
###################
#' A genetic data matrix.
#'
#' This binary matrix contains the allelic states of genetic variables,
#' typically single-nucleotide polymorphisms (SNPs) (or the presence/absence states of accessory genes),
#' showing individuals in the rows and genetic loci in the columns.
#'
#' Each individual in the sample is represented by a unique identifier (name)
#' which corresponds to the name of one row of the snps matrix.
#' Each genetic locus is also required to have a unique name.
#'
#' In this \code{snps} matrix, redundant columns are present for biallelic loci,
#' denoting the state of the second allele as the inverse of the previous column
#' (e.g., compare locus 1.g and locus 1.a).
#' These biallelic sites can be condensed into a more efficient binary form
#' by using \code{get.binary.snps(snps)}
#' (see the treeWAS vignette).
#' % (see vignette("treeWAS")).
#'
#'
#' @docType data
#' @name snps
#' @usage data(snps)
#' @format A binary matrix with 100 rows and 20,003 columns.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL

########################################################################
###################
## DOCUMENTATION ##
###################
#' The ancestral state reconstruction of a genetic data matrix.
#'
#' This binary matrix contains the terminal and ancestral allelic states of a set of genetic variables
#' (for the original genetic data matrix, see: \code{data(snps)}),
#' showing individuals in the rows and genetic loci in the columns.
#' The observed genotypic states of sampled individuals
#' (i.e., those represented at the terminal nodes of a phylogenetic tree)
#' are presented first, in elements 1:N (here 1:100).
#' These rows of the matrix are identical to the input \code{snps} matrix (see: \code{data(snps)}).
#' The unobserved ancestral states of the genotype at internal nodes have been
#' inferred via ancestral state reconstruction, using \code{asr(snps, tree)}.
#'
#' Each individual in the sample is represented by a unique identifier (name)
#' which corresponds to the name of one row of the snps matrix.
#' (Internal node names have been generated during ancestral state reconstruction.)
#' Each genetic locus is also required to have a unique name.
#'
#' In this \code{snps.reconstruction} matrix, redundant columns are present for biallelic loci,
#' denoting the state of the second allele as the inverse of the previous column
#' (e.g., compare locus 1.g and locus 1.a).
#' These biallelic sites can be condensed into a more efficient binary form
#' by using \code{get.binary.snps(snps)}
#' (see the treeWAS vignette).
#' % (see vignette("treeWAS")).
#'
#'
#' @docType data
#' @name snps.reconstruction
#' @usage data(snps.reconstruction)
#' @format A binary matrix with 199 rows and 20,003 columns.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL

########################################################################
###################
## DOCUMENTATION ##
###################
#' The phenotypically-associated sites in the \code{snps} matrix.
#'
#' This vector specifies the identities (names) of the loci in the genetic data matrix (see: \code{data(snps)})
#' that have been simulated along the phylogenetic tree (see: \code{data(tree)})
#' to be in statistical association with the phenotype (see: \code{data(phen)}).
#' Comparing this vector of snps column names to the set of snps loci identified by treeWAS
#' allows us to evaluate the performance of the treeWAS GWAS method.
#' After applying treeWAS to the components of this dataset,
#' using: \code{treeWAS(snps, phen, tree)},
#' we can assess the ability of treeWAS to recover these "known" associated sites
#' via any of its three association scores
#' (see the treeWAS vignette).
#' % (see vignette("treeWAS")).
#'
#'
#' @docType data
#' @name snps.assoc
#' @usage data(snps.assoc)
#' @format A vector of length 10.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL

########################################################################
###################
## DOCUMENTATION ##
###################
#' A phylogenetic tree.
#'
#' This phylogenetic tree is a phylo object (see \code{vignette("Trees", package="phangorn")})
#' connecting the individuals represented in
#' the rows of the example genetic data matrix (see: \code{data(snps)})
#' and the elements of the example phenotypic vector (see: \code{data(phen)}).
#'
#' In this case, the tree was generated via simulation and used to simulate the genotypic and phenotypic data.
#' In a typical empirical analysis, however, a phylogenetic tree would represent the inferred
#' ancestral relationships between individuals, and it would be estimated from the available genetic data.
#' For example, such a phylogeny could be reconstructed using \code{tree.reconstruct(snps, method="NJ")},
#' or automatically generated within treeWAS, according to the \code{tree} argument, as in:
#' \code{treeWAS(snps, phen, tree="NJ")}
#' (see the treeWAS vignette).
#' % (see vignette("treeWAS")).
#'
#'
#' @docType data
#' @name tree
#' @usage data(tree)
#' @format A phylo object with 100 terminal nodes and 99 internal nodes.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL

########################################################################
###################
## DOCUMENTATION ##
###################
#' Example output of treeWAS.
#'
#' This "treeWAS" class object is a list containing the output of a treeWAS analysis.
#' This GWAS analysis was performed to identify associations between
#' loci in the example genetic data matrix (see: \code{data(snps)})
#' and phenotypic states in the example phenotypic vector (see: \code{data(phen)}),
#' along the phylogenetic tree (see: \code{data(tree)}).
#'
#'
#' This \code{treeWAS} output was returned by the function:
#' \code{treeWAS(snps, phen, tree)}.
#' treeWAS output contains elements describing
#' the significant associations identified by each of the
#' three association scores applied to all genetic loci.
#' Additional elements of the output return all data that was
#' used in the GWAS analysis, including both data input to treeWAS
#' and all relevant data generated by treeWAS.
#'
#' For a detailed description of the elements of this output,
#' please scroll down to the "Value" section of the \code{treeWAS} function documentation,
#' which can be accessed with: \code{?treeWAS}.
#' More information can also be found in the treeWAS vignette.
#' % (see vignette("treeWAS")).
#'
#'
#' @docType data
#' @name treeWAS.example.out
#' @usage data(treeWAS.example.out)
#' @format A treeWAS class object, comprising a list of length 5.
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#' @keywords data datasets
NULL
