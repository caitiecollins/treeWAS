% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeWAS.R
\name{print.treeWAS}
\alias{print.treeWAS}
\title{Print \code{treeWAS} output.}
\usage{
\method{print}{treeWAS}(x, sort.by.p = FALSE, digits = 3)
}
\arguments{
\item{x}{The output returned by \code{treeWAS}.}

\item{sort.by.p}{A logical indicating whether to sort the results by decreasing p-value (\code{TRUE})
or by locus (\code{FALSE}, the default).}
}
\description{
Print the results of \code{treeWAS}, excluding longer data elements within the output.
}
\author{
Caitlin Collins \email{caitiecollins@gmail.com}
}
