% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cauchycombine.R
\name{CombinePValues}
\alias{CombinePValues}
\title{Cauchy combination}
\usage{
CombinePValues(pvalues, weights = NULL)
}
\arguments{
\item{pvalues}{A p-value matrix}

\item{weights}{A weight matrix, which should have the same shape with the \code{pvalues}. The row sum should be 1.}
}
\value{
Combined p-values.
}
\description{
Using Cauchy method to cambine different p-values.
}
