#' Robust DAA of Microbiome Sequencing Data
#'
#' Robust differential abundance analysis of Microbiome sequencing data based on M-estimation method.
#'
#' @import stats
#' @importFrom MASS rlm
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
#' @param reg_method Character. The regression method. It can be "psi.huber" for the Huber method and "psi.bisquare" for the Bisquare method.
#' The default is "psi.huber".
#' @param hyper_para_vec A numeric vector for hyperparameters used in regression loss function.
#' The default is exp(seq(log(1.345), log(5), length = 10)) for Huber method.
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
rlm_fun <- function(otu_tab, meta, formula, adaptive = TRUE, imputation = FALSE,
                    pseudo_cnt = 0.5, corr_cut = 0.1, prev_cut = 0, lib_cut = 1,
                    reg_method = "psi.huber", hyper_para_vec = NULL,
                    adj_method = "BH", alpha = 0.05) {
  #### check input
  if (any(is.na(otu_tab))) {
    stop("The OTU table contains NAs! Please remove!\n")
  }
  allvars <- colnames(meta)

  ## pre processing
  keep.sam <- which(colSums(otu_tab) >= lib_cut & rowSums(is.na(meta)) == 0)
  otu_tab <- otu_tab[, keep.sam]
  meta <- as.data.frame(meta[keep.sam, ])
  colnames(meta) <- allvars

  n <- ncol(otu_tab)
  keep.tax <- which(rowSums(otu_tab > 0) / n >= prev_cut)
  otu_tab <- otu_tab[keep.tax, ]
  m <- nrow(otu_tab)

  ## some samples may have zero total counts after screening taxa
  if (any(colSums(otu_tab) == 0)) {
    ind <- which(colSums(otu_tab) > 0)
    otu_tab <- otu_tab[, ind]
    meta <- as.data.frame(meta[ind, ])
    names(meta) <- allvars
    keep.sam <- keep.sam[ind]
    n <- ncol(otu_tab)
  }

  ## scaling numerical variables
  ind <- sapply(seq_len(ncol(meta)), function(i) is.numeric(meta[, i]))
  meta[, ind] <- scale(meta[, ind])
  p <- ncol(meta) + 1

  ## handling zeros
  if (any(otu_tab == 0)) {
    N <- colSums(otu_tab)
    if (adaptive) {
      logN <- log(N)
      tmp <- lm(as.formula(paste0("logN", formula)), meta)
      corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
      if (any(corr.pval <= corr_cut)) {
        imputation <- TRUE
      } else {
        imputation <- FALSE
      }
    }
    if (imputation) {
      N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
      N.mat[otu_tab > 0] <- 0
      tmp <- N[max.col(N.mat)]
      otu_tab <- otu_tab + N.mat / tmp
    } else {
      otu_tab <- otu_tab + pseudo_cnt
    }
  }

  ## CLR transformation
  logY <- log2(otu_tab)
  W <- t(logY) - colMeans(logY)

  ## robust linear regression
  options(warn = -1)
  formula_use <- as.formula(paste0("w_tmp", formula))
  ## Huber method
  if (reg_method == "psi.huber") {
    if (is.null(hyper_para_vec)) {
      hyper_para_vec <- exp(seq(log(1.345), log(5), length = 10))
    }
    n_hyper <- length(hyper_para_vec)
    pvalue_mat <- matrix(NA, nrow = m, ncol = n_hyper)
    ## different hyperparameter
    for (iter_hyper in seq_len(n_hyper)) {
      hyper_para <- hyper_para_vec[iter_hyper]
      alpha_vec <- numeric(m)
      sd_alpha <- numeric(m)
      for (iter_m in seq_len(m)) {
        w_tmp <- W[, iter_m]
        fit <- rlm(formula_use, data = meta, psi = "psi.huber", maxit = 50, k = hyper_para)
        alpha_vec[iter_m] <- coef(fit)["u1"]
        summary_fit <- summary(fit)
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
      pvalue_mat[, iter_hyper] <- p_value
    }
  } else if (reg_method == "psi.bisquare") {
    ## bisquare method
    if (is.null(hyper_para_vec)) {
      hyper_para_vec <- exp(seq(log(4.685), log(20), length = 10))
    }
    n_hyper <- length(hyper_para_vec)
    pvalue_mat <- matrix(NA, nrow = m, ncol = n_hyper)
    ## different hyperparameter
    for (iter_hyper in seq_len(n_hyper)) {
      hyper_para <- hyper_para_vec[iter_hyper]
      alpha_vec <- numeric(m)
      sd_alpha <- numeric(m)
      for (iter_m in seq_len(m)) {
        w_tmp <- W[, iter_m]
        fit <- rlm(formula_use, data = meta, psi = "psi.bisquare", maxit = 50, c = hyper_para)
        alpha_vec[iter_m] <- coef(fit)["u1"]
        summary_fit <- summary(fit)
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
      pvalue_mat[, iter_hyper] <- p_value
    }
  }
  options(warn = 0)
  weights <- matrix(rep(1.0 / n_hyper, n_hyper * m), ncol = n_hyper)
  combine_pvalue <- CombinePValues(pvalue_mat, weights)
  p_adj <- p.adjust(combine_pvalue, "BH")
  index_select <- which(p_adj < alpha)
  #### return results
  return(list(
    p_value_mat = pvalue_mat, combine_pvalue = combine_pvalue,
    index_select = index_select
  ))
}
