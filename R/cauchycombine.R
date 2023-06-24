#' Cauchy combination
#'
#' Using Cauchy method to cambine different p-values.
#'
#' @import stats
#'
#' @param pvalues A p-value matrix
#' @param weights A weight matrix, which should have the same shape with the \code{pvalues}. The row sum should be 1.
#'
#' @return Combined p-values.
#'
CombinePValues <- function(pvalues, weights = NULL) {
  if (!is.matrix(pvalues)) {
    pvalues <- as.matrix(pvalues)
  }
  ## to avoid extremely values
  pvalues[which(pvalues == 0)] <- 5.55e-17
  pvalues[which((1 - pvalues) < 1e-3)] <- 0.99

  num_pval <- ncol(pvalues)
  num_gene <- nrow(pvalues)
  if (is.null(weights)) {
    weights <- matrix(rep(1.0 / num_pval, num_pval * num_gene), ncol = num_pval)
  } # end fi
  if ((nrow(weights) != num_gene) || (ncol(weights) != num_pval)) {
    stop("the dimensions of weights does not match that of combined pvalues")
  } # end fi

  Cstat <- tan((0.5 - pvalues) * pi)

  wCstat <- weights * Cstat
  Cbar <- apply(wCstat, 1, sum)
  # combined_pval <- 1.0/2.0 - atan(Cbar)/pi
  combined_pval <- 1.0 - pcauchy(Cbar)
  combined_pval[which(combined_pval <= 0)] <- 5.55e-17
  return(combined_pval)
}
