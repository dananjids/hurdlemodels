#' A function to calculate posterior log rpp
#'
#' @param log_cdf log cdf.
#' @param log_pmf log pmf.
#' @export
#'
post_logrpp <- function(log_cdf, log_pmf){
  mc_used <- dim(log_pmf)[1]
  n <- dim(log_pmf)[2]
  log_u <- log(matrix(runif(n), mc_used, n, byrow = TRUE))
  Log_Add_Exps <- function (log.cdf, log.pmf)
  {
    lM <- pmax(log.cdf,log.pmf)
    lM + log(exp (log_cdf-lM) + exp(log.pmf-lM))
  }
  log_pv<-Log_Add_Exps(log.cdf=log_cdf, log.pmf=log_pmf+log_u)

  log_mean_exp <- function (lx)
  {
    log_sum_exp (lx) - log(length (lx))
  }
  logrpp_post<- apply(log_pv, 2, log_mean_exp)

  logrpp_post
}
