#' A function to calculate Bartlett of Zresidual
#'
#' @param Zresidual A z-residual.
#' @param X Linear predictor or covariate
#' @param k.bl Number of bins if applicable
#' @export bartlett.test.zresid

bartlett.test.zresid <- function (Zresidual, X = c("fitted.value", "covariate"), k.bl=10)
{
  if (missing(X))
    X = "fitted.value"

  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10

  if (X == "fitted.value") {
    fitted.value <- attr(Zresidual, "fitted.value")
    bl.pv<-rep(0,ncol(Zresidual))
    for(j in 1:ncol(Zresidual)){
      #fitted.value <- fitted.value[which(!is.na(Zresidual[,j]))]
      bl.pv[j]<-test.var.bartl(Zresidual[,j], fitted.value, k.bl)
    }
    bl.pv
  }
  if (X != "fitted.value") {
    fitted.value <- attr(Zresidual, "covariates")
    if(X == "covariate"){
      i<-1
      cat("To plot against other covariates, set X to be the covariate name. Please copy one of the covariate name:",
          variable.names(fitted.value))

    } else if(X %in% variable.names(fitted.value)){
      cov.name<-variable.names(fitted.value)
      i<- which(cov.name==X)
    } else{stop(paste0("X must be the one of covariate name: ", variable.names(fitted.value),". ")) }

    bl.pv<-rep(0,ncol(Zresidual))
    for(j in 1:ncol(Zresidual)){
      id.na <- which(is.na(Zresidual[,j]))
      count.id <- which(!is.na(Zresidual[,j]))
      new.Zresidual <- Zresidual[count.id, j]
      bl.pv[j]<-test.var.bartl(new.Zresidual, fitted.value[,i][count.id], k.bl)
    }
    bl.pv
  }
  bl.pv
}
