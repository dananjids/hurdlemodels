#' A function to calculate important sampling cross validation log rpp
#'
#' @param log_cdf log cdf.
#' @param log_pmf log pmf.
#' @export
#'
iscv_logrpp <- function(log_cdf, log_pmf){
  mc_used <- dim(log_pmf)[1]
  n <- dim(log_pmf)[2]
  log_u <- log(matrix(runif(n), mc_used, n, byrow = TRUE))

  Log_Add_Exps <- function (log.cdf, log.pmf)
  {
    lM <- pmax(log.cdf,log.pmf)
    lM + log(exp (log_cdf-lM) + exp(log.pmf-lM))
  }

  log_pv<-Log_Add_Exps(log.cdf=log_cdf, log.pmf=log_pmf+log_u)

  logrpp_iscv<-apply(log_pv - log_pmf, 2, log_sum_exp) -
    apply(-log_pmf, 2, log_sum_exp)

  id.large <- which(exp(logrpp_iscv) == 1.000000e+00)
  id.small <- which(exp(logrpp_iscv) == 0.000000e+00)
  logrpp_iscv[id.large] <- log(9e-5)
  logrpp_iscv[id.small] <- log(1e-5)

  logrpp_iscv
}
