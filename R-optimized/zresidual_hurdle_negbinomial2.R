#' A function to calculate z-residuals of a 'brm' fit.
#' This function is to be used when the user needs to calculate the z-residuals of TNB/HNB
#'
#' @param fit A `brm` fit.
#' @param type The component which the z-residuals should be calculated.
#' @export
#'
zresidual_hurdle_negbinomial2 <- function(fit,  type , method = "iscv", nrep = 1){

  data <- fit$data
  response <- fit$formula$resp
  model.data <- model.frame(fit$formula, data=data)
  model.var <- names(model.data)
  if(!response %in% model.var) response <- strsplit(as.character(fit$formula$formula), "~")[[2]]
  sim.y <- as.matrix(model.data[, response])
  n <- length(sim.y)
  #if(type == "count") id <- which(sim.y > 0) else id <- 1:n
  zero_id <- which(sim.y == 0)

  # Argument should be one of the element in type_list
  type_list <- c(zero = "zero2", count = "TNB2", hurdle = "HNB2")
  rpp_list <- c(iscv = "iscv_logrpp", post = "post_logrpp")

  ldist <- get(paste0("log.pred.dist.", type_list[type]))(fit)
  lpmf <- ldist$lpmf_hat
  lcdf <- ldist$lcdf_hat

  #if(count_only) z_res<- matrix(NA, ncol = nrep, nrow = dim(lpmf)[2])
  z_res<- matrix(NA, ncol = nrep, nrow = n)
  rpp_fun <- get(rpp_list[[method]])
  z_res <- matrix(-qnorm(replicate(nrep, rpp_fun(lcdf, lpmf)), log.p = T),
                         nrow = n, ncol = nrep)
  # for (i in 1:nrep) {
  #   rpp <- get(rpp_list[[method]])(lcdf, lpmf)
  #   z_res[, i] <- -qnorm(rpp, log.p = T)
  # }

  #if(type == "count") z_res <- z_res[-zero_id,]

  colnames(z_res) <- paste0("zresidual", 1:nrep)

  attributes(z_res) <- c(attributes(z_res),list(
    type = type,
    #count_only = count_only,
    zero_id = zero_id,
    log_pmf = lpmf,
    log_cdf = lcdf,
    covariates = fit$data[names(fit$data) != response],
    fitted.value = predict(fit, type = "conditional")[,1]
  ))

  class(z_res) <- c("zresid", class(z_res))

  return(z_res)
}
