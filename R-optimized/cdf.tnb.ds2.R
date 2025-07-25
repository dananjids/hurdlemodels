#' A function to calculate cdf of Truncated Negative Binomial
#'
#' @param y y values.
#' @param mu mu parameter of TNB distribution.
#' @param size size parameter of TNB distribution.
#' @export

cdf.tnb.ds2 <- function(y, mu, size, lower.tail = FALSE, log.p = FALSE) {

  n.y <- max(length(y), length(mu),length(size))
  y_mat <- matrix(rep(y, length=n.y), ncol = length(y), byrow = T)

  log_upper_prob <- matrix(rep(0, length=n.y), ncol = length(y))

  ## computing upper tail probabilities for positive y only (y>0)
  log_upper_prob <- pnbinom(q=y_mat, mu=mu,size=size, lower.tail = FALSE, log.p = TRUE) -
    pnbinom(q=0,mu=mu,size=size, lower.tail = FALSE, log.p = TRUE)

  upper_prob <- exp(log_upper_prob)
  log_lower_prob <- log_diff_exp(0, log_upper_prob)
  lower_prob <- exp (log_lower_prob)

  if (lower.tail==FALSE)
  {
    if(log.p) log_upper_prob
    else upper_prob
  } else
  {
    if(log.p) log_lower_prob
    else lower_prob
  }
}
