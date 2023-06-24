#' Robust DAA of Microbiome Sequencing Data
#'
#' Robust differential abundance analysis of Microbiome sequencing data based on quantile regression.
#'
#' @import stats
#' @importFrom quantreg rq
#' @importFrom modeest mlv
#'
#' @param otu_tab Data frame or matrix representing observed OTU table. Row: taxa; column: samples.
#' NAs are not expected in OTU tables so are not allowed in function \code{rlm_fun}.
#' @param meta Data frame of covariates. The rows of \code{meta} correspond to the columns of \code{otu_tab}.
#' NAs are allowed. If there are NAs, the corresponding samples will be removed in the analysis.
#' @param formula Character. For example: \code{formula = '~x1*x2+x3'}. At least one fixed effect is required.
#' @param adaptive TRUE or FALSE. The default is TRUE. If TRUE, the parameter \code{imputation} will be treated as FALSE no matter
#' what it is actually set to be. Then the significant correlations between the sequencing depth and explanatory variables will be tested
#' via the linear regression between the log of the sequencing depths and \code{formula}. If any p-value is smaller than or equal to
#' \code{corr_cut}, the imputation approach will be used; otherwise, the pseudo-count approach will be used.
#' @param imputation TRUE or FALSE. The default is FALSE. If TRUE, then we use the imputation approach, i.e., zeros in \code{otu_tab} will be
#' imputed using the formula in the referenced paper.
#' @param pseudo_cnt A positive real value. The default is 0.5. If \code{adaptive} and \code{imputation} are both FALSE,
#' then we use the pseudo-count approach, i.e., we add \code{pseudo_cnt} to each value in \code{otu_tab}.
#' @param corr_cut A real value between 0 and 1; significance level of correlations between the sequencing depth and
#' explanatory variables. The default is 0.1.
#' @param prev_cut A real value between 0 and 1; taxa with prevalence (percentage of nonzeros)
#' less than \code{prev_cut} are excluded. The default is 0 (no taxa will be excluded).
#' @param lib_cut A non-negative real value; samples with less than \code{lib_cut} read counts are excluded.
#' The default is 1 (no samples will be excluded).
#' @param tau_vec A numeric vector for quantile to use in quantile regression. The default is seq(0.25, 0.75, by = 0.05).
#' @param adj_method A Character; p-value adjusting approach. See R function \code{p.adjust}. The default is 'BH'.
#' @param alpha A real value between 0 and 1; significance level of differential abundance. The default is 0.05.
#'
#' @return A list with the elements
#' \item{p_value_mat}{A matrix of p-values. Each column corresponding to a hyperparameter; each row corresponding to a taxa.}
#' \item{combine_pvalue}{A numeric vector. Combined p-value gotten by Cauchy combination method.}
#' \item{index_select}{A index vector. The index of significant taxa.}
#'
#' @export
#'
qr_fun <- function(otu_tab, meta, formula, adaptive = TRUE, imputation = FALSE,
                   pseudo_cnt = 0.5, corr_cut = 0.1, prev_cut = 0, lib_cut = 1,
                   tau_vec = NULL, adj_method = "BH", alpha = 0.05) {
  #### check input
  if (any(is.na(otu_tab))) {
    stop("The OTU table contains NAs! Please remove!\n")
  }
  allvars <- colnames(meta)
  if (is.null(tau_vec)) {
    tau_vec <- seq(0.25, 0.75, by = 0.05)
  }
  n_tau <- length(tau_vec)

  ## preprocessing
  keep.sam <- which(colSums(otu_tab) >= lib_cut & rowSums(is.na(meta)) == 0)
  Y <- otu_tab[, keep.sam]
  Z <- as.data.frame(meta[keep.sam, ])
  colnames(Z) <- allvars

  n <- ncol(Y)
  keep.tax <- which(rowSums(Y > 0) / n >= prev_cut)
  Y <- Y[keep.tax, ]
  m <- nrow(Y)

  ## some samples may have zero total counts after screening taxa
  if (any(colSums(Y) == 0)) {
    ind <- which(colSums(Y) > 0)
    Y <- Y[, ind]
    Z <- as.data.frame(Z[ind, ])
    names(Z) <- allvars
    keep.sam <- keep.sam[ind]
    n <- ncol(Y)
  }

  ## scaling numerical variables
  ind <- sapply(seq_len(ncol(Z)), function(i) is.numeric(Z[, i]))
  Z[, ind] <- scale(Z[, ind])
  p <- ncol(Z) + 1

  ## handling zeros
  if (any(Y == 0)) {
    N <- colSums(Y)
    if (adaptive) {
      logN <- log(N)
      tmp <- lm(as.formula(paste0("logN", formula)), Z)
      corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
      if (any(corr.pval <= corr_cut)) {
        # cat("Imputation approach is used.\n")
        imputation <- TRUE
      } else {
        # cat("Pseudo-count approach is used.\n")
        imputation <- FALSE
      }
    }
    if (imputation) {
      N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
      N.mat[Y > 0] <- 0
      tmp <- N[max.col(N.mat)]
      Y <- Y + N.mat / tmp
    } else {
      Y <- Y + pseudo_cnt
    }
  }

  ## CLR transformation
  logY <- log2(Y)
  W <- t(logY) - colMeans(logY)

  ## quantile regression
  options(warn = -1)
  formula_use <- as.formula(paste0("w_tmp", formula))
  pvalue_mat <- matrix(NA, nrow = m, ncol = n_tau)
  for (iter_tau in seq_len(n_tau)) {
    tau <- tau_vec[iter_tau]
    alpha_vec <- numeric(m)
    sd_alpha <- numeric(m)
    for (iter_m in seq_len(m)) {
      w_tmp <- W[, iter_m]
      fit <- rq(formula_use, tau = tau, data = Z)
      alpha_vec[iter_m] <- coef(fit)["u1"]
      summary_fit <- summary(fit, se = "boot", R = 500)
      sd_alpha[iter_m] <- coef(summary_fit)["u1", 2]
    }
    #### mode correction
    # bias <- median(alpha_vec)
    bias <- mlv(sqrt(n) * alpha_vec,
                         method = "meanshift", kernel = "gaussian"
    ) / sqrt(n)
    alpha_correct <- alpha_vec - bias

    #### test
    ## t-test
    stat <- alpha_correct / sd_alpha
    p_value <- 2 * pt(-abs(stat), df = n - p)
    pvalue_mat[, iter_tau] <- p_value
  }
  options(warn = 0)
  weights <- matrix(rep(1.0 / n_tau, n_tau * m), ncol = n_tau)
  combine_pvalue <- CombinePValues(pvalue_mat, weights)
  p_adj <- p.adjust(combine_pvalue, "BH")
  index_select <- which(p_adj < alpha)
  #### return results
  return(list(
    p_value_mat = pvalue_mat, combine_pvalue = combine_pvalue,
    index_select = index_select
  ))
}
