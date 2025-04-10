################################################################################
################################################################################
################################################################################
################################################################################
# Comments
# i) the function fitQR() is the latest and its in gamlss.foreach
# ii) the function qreg() is the one I used in GAMLSS_quantile_regression.Rnw
# iii) fitQR() is general in the sense it can take any formula including smoothing 
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
################################################################################
################################################################################
################################################################################
################################################################################
qreg <- function(formula, data, alpha=0.5, weights,  ...)
{
  fnu <- function(alpha) 
  {
    if (alpha<0 || alpha>1) stop("alpha should be between 0 and 1")
    ((1-alpha)/alpha)^.5 
  }
 # if (missing(weights)) weights <- rep(1,dim(data)[1])
  if (alpha<0 || alpha>1) stop("alpha should be between 0 and 1")
  m1 <- gamlss(formula=formula, data=data,
               family=gamlss.dist::SEP3, 
               sigma.fix=TRUE, sigma.start=1,
               nu.fix=TRUE, nu.start=fnu(alpha),  
               tau.fix=TRUE, tau.start=1, c.crit=0.01,...)
  m1
}
################################################################################
################################################################################
################################################################################
################################################################################






