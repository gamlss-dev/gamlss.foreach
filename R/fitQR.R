################################################################################
################################################################################
################################################################################
################################################################################
fitQR <- function(formula, data, quantile=c(0.5), ...)
{
################################################################################  
#  require(gamlss.foreach)
  m <- NULL
################################################################################  
    fnu <- function(quantile) 
  {
    if (quantile<0 || quantile>1) stop("quantile should be between 0 and 1")
    ((1-quantile)/quantile)^.5 
  }
################################################################################ 
  if (any(quantile<0) || any(quantile>1)) stop("quantile should be between 0 and 1")
models <- list()
     R <- length(quantile)
models <- foreach(m = 1:R)%dopar%{
  m1 <- gamlss(formula=formula, data=data, family=gamlss.dist::SEP3, 
               sigma.fix=TRUE, sigma.start=1,
               nu.fix=TRUE, nu.start=fnu(quantile[m]),  
               tau.fix=TRUE, tau.start=1, trace=FALSE,
               ...)
  m1
}
models  
}
######################################################################################
######################################################################################
