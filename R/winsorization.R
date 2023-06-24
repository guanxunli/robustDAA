#' Winsorization
#'
#' @param out_tab Data frame or matrix representing observed OTU table. Row: taxa; column: samples.
#' @param quan A real value between 0 and 1; winsorization cutoff (quantile) for \code{otu.tab}, e.g., 0.97.
#'
#' @return OTU table after winsorization.
#'
#' @export
#'
winsor.fun <- function(out_tab, quan) {
  N <- colSums(out_tab)
  P <- t(t(out_tab) / N)
  cut <- apply(P, 1, quantile, quan)
  Cut <- matrix(rep(cut, ncol(out_tab)), nrow(out_tab))
  ind <- P > Cut
  P[ind] <- Cut[ind]
  out_tab <- round(t(t(P) * N))
  return(out_tab)
}
