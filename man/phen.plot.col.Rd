% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{phen.plot.col}
\alias{phen.plot.col}
\title{Phenotypic tree-colouring schemes.}
\format{
A list of length 5.
}
\usage{
data(phen.plot.col)
}
\description{
A list containing the colour values that \code{plot_phen} generates to represent
the states and substitutions of the phenotypic variable (\code{data(phen)})
along the phylogenetic tree (\code{data(tree)}), with \code{plot_phen(tree, phen.nodes=phen)}.
You are unlikely to have to interact with this list,
as the colours are automatically plotted by the \code{plot_phen} function.
}
\details{
The five elements of this list give the colour schemes used to indicate the phenotypic state at:
edge.labels, edges, all.nodes, internal.nodes, and tip.labels.
}
\author{
Caitlin Collins \email{caitiecollins@gmail.com}
}
\keyword{data}
\keyword{datasets}
