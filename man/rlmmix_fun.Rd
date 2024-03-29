% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rlmmix_fun.R
\name{rlmmix_fun}
\alias{rlmmix_fun}
\title{Robust DAA of Microbiome Sequencing Data}
\usage{
rlmmix_fun(
  otu_tab,
  meta,
  formula,
  adaptive = TRUE,
  imputation = FALSE,
  pseudo_cnt = 0.5,
  corr_cut = 0.1,
  prev_cut = 0,
  lib_cut = 1,
  hyper_para_vec = NULL,
  adj_method = "BH",
  alpha = 0.05
)
}
\arguments{
\item{otu_tab}{Data frame or matrix representing observed OTU table. Row: taxa; column: samples.
NAs are not expected in OTU tables so are not allowed in function \code{rlm_fun}.}

\item{meta}{Data frame of covariates. The rows of \code{meta} correspond to the columns of \code{otu_tab}.
NAs are allowed. If there are NAs, the corresponding samples will be removed in the analysis.}

\item{formula}{Character. For example: \code{formula = '~x1 + (1|id)'}.}

\item{adaptive}{TRUE or FALSE. The default is TRUE. If TRUE, the parameter \code{imputation} will be treated as FALSE no matter
what it is actually set to be. Then the significant correlations between the sequencing depth and explanatory variables will be tested
via the linear regression between the log of the sequencing depths and \code{formula}. If any p-value is smaller than or equal to
\code{corr_cut}, the imputation approach will be used; otherwise, the pseudo-count approach will be used.}

\item{imputation}{TRUE or FALSE. The default is FALSE. If TRUE, then we use the imputation approach, i.e., zeros in \code{otu_tab} will be
imputed using the formula in the referenced paper.}

\item{pseudo_cnt}{A positive real value. The default is 0.5. If \code{adaptive} and \code{imputation} are both FALSE,
then we use the pseudo-count approach, i.e., we add \code{pseudo_cnt} to each value in \code{otu_tab}.}

\item{corr_cut}{A real value between 0 and 1; significance level of correlations between the sequencing depth and
explanatory variables. The default is 0.1.}

\item{prev_cut}{A real value between 0 and 1; taxa with prevalence (percentage of nonzeros)
less than \code{prev_cut} are excluded. The default is 0 (no taxa will be excluded).}

\item{lib_cut}{A non-negative real value; samples with less than \code{lib_cut} read counts are excluded.
The default is 1 (no samples will be excluded).}

\item{hyper_para_vec}{A numeric vector for hyperparameters used in regression loss function.
The default is c(1.345, 2.28).}

\item{adj_method}{A Character; p-value adjusting approach. See R function \code{p.adjust}. The default is 'BH'.}

\item{alpha}{A real value between 0 and 1; significance level of differential abundance. The default is 0.05.}
}
\value{
A list with the elements
\item{p_value_mat}{A matrix of p-values. Each column corresponding to a hyperparameter; each row corresponding to a taxa.}
\item{combine_pvalue}{A numeric vector. Combined p-value gotten by Cauchy combination method.}
\item{index_select}{A index vector. The index of significant taxa.}
}
\description{
Robust differential abundance analysis of Microbiome sequencing data based on M-estimation method with random effect.
}
