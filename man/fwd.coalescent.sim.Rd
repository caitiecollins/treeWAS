% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fwd.coalescent.sim.R
\name{fwd.coalescent.sim}
\alias{fwd.coalescent.sim}
\title{Short one-phrase description.}
\usage{
fwd.coalescent.sim(
  n.ind = 100,
  n.snps = 10000,
  n.subs = 1,
  n.snps.assoc = 10,
  n.subs.assoc = 15,
  p = 1,
  heatmap = FALSE,
  reconstruct = FALSE,
  dist.dna.model = "JC69",
  seed = 1
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

\item{heatmap}{A logical indicating whether to produce a heatmap of the genetic distance
between the simulated genomes of the n.ind individuals.}

\item{reconstruct}{Either a logical indicating whether to attempt to reconstruct
a phylogenetic tree using the simulated genetic data, or one of c("UPGMA", "nj", "ml")
to specify that tree reconstruction is desired by one of these three methods
(Unweighted Pair Group Method with Arithmetic Mean, Neighbour-Joining, Maximum-Likelihood).}

\item{seed}{An optional integer controlling the pseudo-random process of simulation. Two
instances of coalescent.sim with the same seed and arguments will produce identical output.}

\item{assoc.prob}{An optional integer (> 0, <= 100) specifying the strength of the
association between the n.snps.assoc loci and the phenotype (see details).}

\item{n.phen.subs}{An integer specifying the expected number of phenotypic
substitutions to occur on the phylogenetic tree (through the same process as
the n.subs parameter when n.subs is an integer (see details)).}

\item{phen}{An optional vector containing a phenotype for each of the
n.ind individuals if no phenotypic simulation is desired.}
}
\description{
Longer proper discription of function...
}
\details{
#### n.subs ####
If the value of the n.subs parameter is set to an integer, this integer is
used as the parameter of a Poisson distribution from which the number of substitutions to
occur on the phylogenetic tree is drawn for each of the n.snps simulated genetic loci.
If n.subs is a vector containing a distribution, this is used directly (in proportion to n.snps)
to define the number of substitutions per site. For example, if n.subs=c(3000, 900, 70, 20, 0, 10)
and n.snps=8000, then 6000 simulated sites will undergo exactly
one substitution somewhere on the phylogenetic tree, 1800 will undergo two,
140 three, 40 four, 0 five, and 20 six.
#### assoc.prob ####
The assoc.prob parameter controls the strength of association through a process analagous to dilution.
All n.snps.assoc loci are initially simulated to undergo a substitution
every time the phenotype undergoes a substitution (ie. perfect association).
The assoc.prob parameter then acts like a dilution factor, removing (100 - assoc.prob)%
of the substitutions that occurred during simulation under perfect association.
}
\author{
Caitlin Collins \email{caitiecollins@gmail.com}
}
