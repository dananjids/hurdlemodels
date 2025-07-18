#' A function to calculate pdf of Truncated Negative Binomjial
#'
#' @param y y values.
#' @param mu mu parameter of TNB distribution.
#' @param size size parameter of TNB distribution.
#' @export

pdf.tnb.ds2 <- function(y, mu, size, log.p = FALSE)
{

  n.y <- max(length(y), length(mu),length(size))
  y_mat <- matrix(rep(y, length=n.y), ncol = length(y), byrow = T)

  log_prob <- matrix(rep(-Inf, length=n.y), ncol = length(y))
  #if(count_only) i.py <- which(y>0) else i.py <- 1:n
  #i.py <- which (y>0) ## for y > 0

  ## computing log probabilities for positive y only (y>0)
  log_prob <- dnbinom(x=y_mat,mu=mu,size=size, log = TRUE) -
    pnbinom(q=0,mu=mu,size=size, lower.tail = FALSE, log.p = TRUE)

  prob <- exp(log_prob)

  if(log.p) log_prob
  else prob

}
