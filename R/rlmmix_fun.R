#' Robust DAA of Microbiome Sequencing Data
#'
#' Robust differential abundance analysis of Microbiome sequencing data based on M-estimation method with random effect.
#'
#' @import stats
#' @import lme4
#' @import robustlmm
#' @importFrom pbkrtest get_Lb_ddf
#' @importFrom modeest mlv
#'
#' @param otu_tab Data frame or matrix representing observed OTU table. Row: taxa; column: samples.
#' NAs are not expected in OTU tables so are not allowed in function \code{rlm_fun}.
#' @param meta Data frame of covariates. The rows of \code{meta} correspond to the columns of \code{otu_tab}.
#' NAs are allowed. If there are NAs, the corresponding samples will be removed in the analysis.
#' @param formula Character. For example: \code{formula = '~x1 + (1|id)'}.
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
#' @param hyper_para_vec A numeric vector for hyperparameters used in regression loss function.
#' The default is c(1.345, 2.28).
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
rlmmix_fun <- function(otu_tab, meta, formula, adaptive = TRUE, imputation = FALSE,
                    pseudo_cnt = 0.5, corr_cut = 0.1, prev_cut = 0, lib_cut = 1,
                    hyper_para_vec = NULL, adj_method = "BH", alpha = 0.05) {
  #### check input
  if (any(is.na(otu_tab))) {
    stop("The OTU table contains NAs! Please remove!\n")
  }
  allvars <- colnames(meta)

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

  ## check random effect
  if (grepl("\\(", formula)) {
    random_effect <- TRUE
  } else {
    random_effect <- FALSE
  }

  ## handling zeros
  if (any(Y == 0)) {
    N <- colSums(Y)
    if (adaptive) {
      logN <- log(N)
      if (random_effect) {
        tmp <- lmerTest::lmer(as.formula(paste0("logN", formula)), Z)
      } else {
        tmp <- lm(as.formula(paste0("logN", formula)), Z)
      }
      corr_pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
      if (any(corr_pval <= corr_cut)) {
        # cat("Imputation approach is used.\n")
        imputation <- TRUE
      } else {
        # cat("Pseudo-count approach is used.\n")
        imputation <- FALSE
      }
    }
    if (imputation) {
      Nmat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
      Nmat[Y > 0] <- 0
      tmp <- N[max.col(Nmat)]
      Y <- Y + Nmat / tmp
    } else {
      Y <- Y + pseudo_cnt
    }
  }

  ## CLR transformation
  logY <- log2(Y)
  W <- t(logY) - colMeans(logY)

  ## robust linear regression
  options(warn = -1)
  formula_use <- as.formula(paste0("w_tmp", formula))
  if (is.null(hyper_para_vec)) {
    hyper_para_vec <- c(1.345, 2.28)
  }
  n_hyper <- length(hyper_para_vec)
  pvalue_mat <- matrix(NA, nrow = m, ncol = n_hyper)
  alpha_mat <- matrix(NA, nrow = m, ncol = n_hyper)
  sd_mat <- matrix(NA, nrow = m, ncol = n_hyper)
  df_mat <- matrix(NA, nrow = m, ncol = n_hyper)
  for (iter_m in seq_len(m)) {
    ## fit robust rlm
    w_tmp <- W[, iter_m]
    rfit <- suppressMessages(robustlmm::rlmer(formula_use, data = Z))
    ## parameter
    summary_rfit <- summary(rfit)
    alpha_mat[iter_m, 1] <- coef(summary_rfit)[2, 1]
    sd_mat[iter_m, 1] <- coef(summary_rfit)[2, 2]
    ## degree of freedom
    rfit_df <- rfit
    class(rfit_df) <- "lmerMod" # to use get_Lb_ddf
    df_KR <- get_Lb_ddf(rfit_df, lme4::fixef(rfit_df))
    df_mat[iter_m, 1] <- df_KR
    ## different hyperparameter
    for (iter_hyper in c(2:n_hyper)) {
      hyper_para <- hyper_para_vec[iter_hyper]
      rfit_hyper <- update(rfit, rho.sigma.e = psi2propII(smoothPsi, k = hyper_para),
                           rho.sigma.b = psi2propII(smoothPsi, k = hyper_para))
      ## parameter
      summary_rfit <- summary(rfit_hyper)
      alpha_mat[iter_m, iter_hyper] <- coef(summary_rfit)[2, 1]
      sd_mat[iter_m, iter_hyper] <- coef(summary_rfit)[2, 2]
      ## degree of freedom
      rfit_df <- rfit_hyper
      class(rfit_df) <- "lmerMod" # to use get_Lb_ddf
      df_KR <- get_Lb_ddf(rfit_df, lme4::fixef(rfit_df))
      df_mat[iter_m, iter_hyper] <- df_KR
    }
  }
  options(warn = 0)

  for (iter_hyper in seq_len(n_hyper)) {
    alpha_vec <- alpha_mat[, iter_hyper]
    sd_alpha <- sd_mat[, iter_hyper]
    df_vec <- df_mat[, iter_hyper]
    #### mode correction
    bias <- mlv(sqrt(n) * alpha_vec,
                         method = "meanshift", kernel = "gaussian"
    ) / sqrt(n)
    alpha_correct <- alpha_vec - bias
    #### test
    ## t-test
    stat <- alpha_correct / sd_alpha
    p_value <- 2 * pt(-abs(stat), df = df_vec)
    pvalue_mat[, iter_hyper] <- p_value
  }
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
