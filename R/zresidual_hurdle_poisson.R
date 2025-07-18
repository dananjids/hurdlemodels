#' A function to calculate z-residuals of a 'brm' fit.
#'
#' @param fit A `brm` fit.
#' @param type The component which the z-residuals should be calculated.
#' @export
#'
zresidual_hurdle_poisson <- function(fit,  type , method = "iscv", nrep = 1){

  data <- fit$data
  response <- fit$formula$resp
  model.data <- model.frame(fit$formula, data=data)
  model.var <- names(model.data)
  if(!response %in% model.var) response <- strsplit(as.character(fit$formula$formula), "~")[[2]]
  sim.y <- as.matrix(model.data[, response])
  n <- length(sim.y)
  zero_id <- which(sim.y == 0)

  # Argument should be one of the element in type_list
  type_list <- c("zero", "TP", "HP")
  names(type_list) <- c("zero", "count", "hurdle")

  ldist <- get(paste0("log.pred.dist.", type_list[type]))(fit)
  lpmf <- ldist$lpmf_hat
  lcdf <- ldist$lcdf_hat

  # Argument should be one of the element in rpp_list
  rpp_list <- c(iscv = "iscv_logrpp", post = "post_logrpp")
  names(rpp_list) <- c("iscv", "post")

  z_res<- matrix(NA, ncol = nrep, nrow = n)

  for (i in 1:nrep) {
    rpp <- get(rpp_list[[method]])(lcdf, lpmf)
    z_res[, i] <- -qnorm(rpp, log.p = T)
  }

  colnames(z_res) <- paste0("zresidual", 1:nrep)

  fitted.type <- c(zero="zero", count="mu")
  if(type=="hurdle"){
    fitted.value <- colMeans(posterior_epred(fit))
  } else if(type=="zero"){
    fitted.value <- 1-colMeans(posterior.pred.ds(fit, dpar = fitted.type[type]))
  } else {
    fitted.value <- colMeans(posterior.pred.ds(fit, dpar = fitted.type[type]))
  }

  attributes(z_res) <- c(attributes(z_res),list(
    type = type,
    zero_id = zero_id,
    log_pmf = lpmf,
    log_cdf = lcdf,
    covariates = subset(fit$data, select = -get(response)),
    fitted.value = fitted.value
  ))

  class(z_res) <- c("zresid", class(z_res))

  return(z_res)
}
