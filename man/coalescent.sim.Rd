% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coalescent.sim.R
\name{coalescent.sim}
\alias{coalescent.sim}
\title{Simulate a tree, phenotype, and genetic data.}
\usage{
coalescent.sim(
  n.ind = 100,
  n.snps = 10000,
  n.subs = 1,
  n.snps.assoc = 0,
  assoc.prob = 100,
  n.phen.subs = 15,
  phen = NULL,
  plot = TRUE,
  heatmap = FALSE,
  reconstruct = FALSE,
  dist.dna.model = "JC69",
  grp.min = 0.25,
  row.names = TRUE,
  set = 1,
  tree = NULL,
  coaltree = TRUE,
  s = 20,
  af = 10,
  filename.plot = NULL,
  seed = NULL
)
}
\arguments{
\item{n.ind}{An integer specifying the number of individual genomes to simulate
(ie. the number of terminal nodes in the tree).}

\item{n.snps}{An integer specifying the number of genetic loci to simulate.}

\item{n.subs}{Either an integer or a vector (containing a distribution) that is
used to determine the number of substitutions
to occur on the phylogenetic tree for each genetic locus (see details).}

\item{n.snps.assoc}{An optional integer specifying the number of genetic loci}

\item{assoc.prob}{An optional integer (> 0, <= 100) specifying the strength of the
association between the n.snps.assoc loci and the phenotype (see details).}

\item{n.phen.subs}{An integer specifying the expected number of phenotypic
substitutions to occur on the phylogenetic tree (through the same process as
the n.subs parameter when n.subs is an integer (see details)).}

\item{phen}{An optional vector containing a phenotype for each of the
n.ind individuals if no phenotypic simulation is desired.}

\item{plot}{A logical indicating whether to generate a plot of the phylogenetic tree (\code{TRUE}) or not (\code{FALSE}, the default).}

\item{heatmap}{A logical indicating whether to produce a heatmap of the genetic distance
between the simulated genomes of the n.ind individuals.}

\item{reconstruct}{Either a logical indicating whether to attempt to reconstruct
a phylogenetic tree using the simulated genetic data, or one of c("UPGMA", "nj", "ml")
to specify that tree reconstruction is desired by one of these three methods
(Unweighted Pair Group Method with Arithmetic Mean, Neighbour-Joining, Maximum-Likelihood).}

\item{dist.dna.model}{A character string specifying the type of model to use in reconstructing the phylogenetic tree for
calculating the genetic distance between individual genomes, only used if \code{tree} is
a character string (see ?dist.dna).}

\item{grp.min}{An optional number between 0.1 and 0.9 to control the proportional size of the smaller phenotypic group.}

\item{row.names}{An optional vector containing row names for the individuals to be simulated.}

\item{set}{An integer (1, 2, or 3) required to select the method of generating associated loci if \code{n.snps.assoc} is not zero.}

\item{coaltree}{A logical indicating whether to generate a coalescent tree (\code{TRUE}, the default),
or an rtree-type tree (\code{FALSE}, see ?rtree).}

\item{s}{If \code{set} is 3, the \code{s} parameter controls a baseline number of substitutions to be
experienced by the phenotype and associated loci: by default, 20.}

\item{af}{If \code{set} is 3, the \code{af} parameter provides an association factor,
controlling the preference for association over non-association at associated loci:  by default, 10 (for a 10x preference).}

\item{filename.plot}{An optional character string denoting the file location for saving any plots produced; else \code{NULL}.}

\item{seed}{An optional integer to control the pseudo-randomisation process and allow for identical repeat runs of the function;
else \code{NULL}.}
}
\description{
This funtion allows the user to simulate a phylogenetic tree, as well as
phenotypic and genetic data, including associated and unassociated loci.
}
\details{
\strong{Homoplasy Distribution}

The homoplasy distribution contains the number of substitutions per site.

If the value of the \code{n.subs} parameter is set to an integer, this integer is
used as the parameter of a Poisson distribution from which the number of substitutions to
occur on the phylogenetic tree is drawn for each of the \code{n.snps} simulated genetic loci.

The \code{n.subs} argument can also be used to provide a distribution
to define the number of substitutions per site.

It must be in the form of a \emph{named} vector (or table), or a vector in which the \emph{i}'th element
contains the number of \emph{loci} that have been estimated to undergo \emph{i} substitutions on the tree.
The vector must be of length \emph{max n.subs}, and "empty" indices must contain zeros.
For example: the vector \code{n.subs = c(1833, 642, 17, 6, 1, 0, 0, 1)},
could be used to define the homoplasy distribution for a dataset with 2500 loci,
where the maximum number of substitutions to be undergone on the tree by any locus is 8,
and no loci undergo either 6 or 7 substitutions.


\strong{Association Probability}

The \code{assoc.prob} parameter is only functional when \code{set} is set to 1.
If so, \code{assoc.prob} controls the strength of association through a process analagous to dilution.
All \code{n.snps.assoc} loci are initially simulated to undergo a substitution
every time the phenotype undergoes a substitution (ie. perfect association).
The assoc.prob parameter then acts like a dilution factor, removing \code{(100 - assoc.prob)\%}
of the substitutions that occurred during simulation under perfect association.
}
\examples{
\dontrun{
## load example homoplasy distribution
data(dist_0)
str(dist_0)

## simulate a matrix with 10 associated loci:
dat <- coalescent.sim(n.ind = 100,
                        n.snps = 1000,
                        n.subs = dist_0,
                        n.snps.assoc = 10,
                        assoc.prob = 90,
                        n.phen.subs = 15,
                        phen = NULL,
                        plot = TRUE,
                        heatmap = FALSE,
                        reconstruct = FALSE,
                        dist.dna.model = "JC69",
                        grp.min = 0.25,
                        row.names = NULL,
                        coaltree = TRUE,
                        s = NULL,
                        af = NULL,
                        filename = NULL,
                        set = 1,
                        seed = 1)

## examine output:
str(dat)

## isolate elements of output:
snps <- dat$snps
phen <- dat$phen
snps.assoc <- dat$snps.assoc
tree <- dat$tree
}
}
\author{
Caitlin Collins \email{caitiecollins@gmail.com}
}
