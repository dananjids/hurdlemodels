#' A function to calculate posterior predictive log truncated negative binomial pmf and cdf of 'brm' fit.
#'
#' @param fit A `brm` fit.
#' @export log.pred.dist.TNB

log.pred.dist.TNB <- function(fit){

  summary.fit <- summary(fit)
  chains <- summary.fit$chains
  iter <- summary.fit$iter
  warmup <- summary.fit$warmup
  mc_used <- chains*(iter - warmup)

  data <- fit$data
  response <- fit$formula$resp
  model.data <- model.frame(fit$formula, data=data)
  model.var <- names(model.data)
  if(!response %in% model.var) response <- strsplit(as.character(fit$formula$formula), "~")[[2]]
  sim.y <- as.matrix(model.data[, response])
  n <- length(sim.y)
  count.id <- which(sim.y > 0)
  zero_id <- which(sim.y == 0)

  mu <- posterior.pred.ds(fit, dpar = "mu")
  shape <- posterior.pred.ds(fit, dpar = "shape")

  lpmf_hat <- matrix(NA, mc_used, n)
  lcdf_hat <- matrix(NA, mc_used, n)

  for (i in count.id){
    lpmf_hat[,i] <- pdf.tnb.ds(sim.y[i], mu = mu[,i], size = shape[,i], log = TRUE)
    lcdf_hat[,i] <- cdf.tnb.ds(sim.y[i], mu = mu[,i], size = shape[,i], lower.tail = FALSE, log.p = TRUE)
  }

  pred_dist <- list(lpmf_hat = lpmf_hat, lcdf_hat = lcdf_hat, zero_id = zero_id)
  return(pred_dist)
}
