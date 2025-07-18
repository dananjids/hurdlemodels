#' A function to calculate predictive log truncated negative binomial pmf and cdf of 'brm' fit.
#'
#' @param fit A `brm` fit.
#' @export log.pred.dist.TNB

log.pred.dist.TNB2 <- function(fit){

  # summary.fit <- summary(fit)
  # chains <- summary.fit$chains
  # iter <- summary.fit$iter
  # warmup <- summary.fit$warmup
  chains <- fit$fit@sim$chains
  iter <- fit$fit@sim$iter
  warmup <- fit$fit@sim$warmup
  mc_used <- chains*(iter - warmup)

  data <- fit$data
  response <- fit$formula$resp
  model.data <- model.frame(fit$formula, data=data)
  model.var <- names(model.data)
  if(!response %in% model.var) response <- strsplit(as.character(fit$formula$formula), "~")[[2]]
  sim.y <- model.data[, response]
  n <- length(sim.y)
  count.id <- which(sim.y > 0)
  zero_id <- which(sim.y == 0)

  mu <- posterior.pred.ds(fit, dpar = "mu")
  shape <- posterior.pred.ds(fit, dpar = "shape")

  lpmf_hat <- matrix(NA, mc_used, n)
  lcdf_hat <- matrix(NA, mc_used, n)

  lpmf_hat <- pdf.tnb.ds2(sim.y, mu = mu, size = shape, log = TRUE)
  lcdf_hat <- cdf.tnb.ds2(sim.y, mu = mu, size = shape, lower.tail = FALSE, log.p = TRUE)

  # pred_dist <- list(lpmf_hat = lpmf_hat, lcdf_hat = lcdf_hat, zero_id = zero_id)
  return(list(lpmf_hat = lpmf_hat, lcdf_hat = lcdf_hat, zero_id = zero_id))
}
