###########################################################################
###########################################################################
###########################################################################
# A set of function for doing bayesian and not parametric boosting 
require(gamlss.foreach)
###########################################################################
###########################################################################
###########################################################################
# functions
# i)    coefAll1()         which probably should go in the main gamlss
# ii)   BayesianBoost()    at the moment only coef and parameters are saved 
#                          not gammas or the fitted smoothers  
#                          TO DO 
#        tidy up what to save maybe I should have 
#         options save = c("coef", "parameters", "lambda", "gammas", "smoothers")
# iii)  summary.Bayesian.boot()   summary for the coefficients
#                         
#       there is a summary for parameters 27-10-2022
#       plot.Bayesian.boot()      plotting the distribution of the coefficients  
#                          TO DO
#                     maybe should be on option pages=1 see gamlss.ggplots                           
# iii)  print.Bayesian.boot()     information about boost
# 
# vi)   NonParametricBoot()      This only saves coef 
#                         TO DO 
#                         saved at least parameters  
#       summary.NonParametricBoot()
#       plot.NonParametricBoot()
# v)    print.NonParametricBoot()
# vi)   ratioBoot()
# vii)  momentBaysianBoot
# viii) moment Boot
###########################################################################
# TODO
# # save original coefficients OK DONE
# check Kriton parallelism  
# add summary on the Bayesian.boot and classical.boot classes
###########################################################################
###########################################################################
###########################################################################
# FUNCTION coefAll()
# Function 1
########################################################################### 
# this probably should go in the main gamlss 
# It went but it does not get lambdas yet 
# This is with lambda
########################################################################### 
# TO DO exceptions with different smoothers 
# check 
# pbm 
# pbc 
# pbz
# cs
# scs 
# ri
# pcat
# gmrf
# pvc
# ga
# ba
# nn
# tr
# lo
# ma
# pc
# pcr
# gnet
####################################################################################
####################################################################################
####################################################################################  
coefAll1 <- function(object, lambdas= TRUE, deviance=FALSE, ...)
{
  out <- list()
if ("mu" %in% object$par) #
 {
    out$mu <- coef(object, "mu")  
    if (!is.null(object$mu.s)&&lambdas)
    {
      M <- dim(object$mu.s)[2]
      outLambda <- rep(0,M)
      for (i in 1:M)
      {
        outLambda[i] <-  getSmo(object, which=i)$lambda
      }
      out$mu <- c(out$mu, lambda=outLambda) 
    }
 }
if ("sigma" %in% object$par)#
 {
  out$sigma <- coef(object, "sigma")
  if (!is.null(object$sigma.s)&&lambdas)
    {
      M <- dim(object$sigma.s)[2]
      outLambda <- rep(0,M)
      for (i in 1:M)
      {
        outLambda[i] <-  getSmo(object, parameter="sigma", which=i)$lambda
      }
      out$sigma <- c(out$sigma, lambda=outLambda) 
    }
 }  
if ("nu" %in% object$par) #
 {
    out$nu <- coef(object, "nu")  
  if (!is.null(object$nu.s)&&lambdas)
    {
      M <- dim(object$nu.s)[2]
      outLambda <- rep(0,M)
      for (i in 1:M)
      {
        outLambda[i] <-  getSmo(object, parameter="nu", which=i)$lambda
      }
      out$nu <- c(out$nu, lambda=outLambda) 
    }
 }  
if ("tau" %in% object$par) #
  {
    out$tau <- coef(object, "tau")  
  if (!is.null(object$tau.s)&&lambdas)
    {
      M <- dim(object$tau.s)[2]
      outLambda <- rep(0,M)
      for (i in 1:M)
      {
        outLambda[i] <-  getSmo(object, parameter="tau", which=i)$lambda
      }
      out$tau <- c(out$tau, lambda=outLambda) 
    }
  } 
  if ("tau" %in% object$par)
    out$tau <-  coef(object, "tau")
  if (deviance) out$deviance <- deviance(object)
  return(out)
}
#########################################################################################
#########################################################################################
#########################################################################################
# this should work with any parametric GAMLSS models
# FUNCTION BayesianBoot Rubin (1981) 
###########################################################################
###########################################################################
###########################################################################
BayesianBoot <- function(obj, data=NULL, B=100, newdata=NULL)
{
# local function ----------------------------------------------------------
# Rubin 1981 generation of weights 
rFun  <- function(n)
{
     U <- runif(n-1)
    oU <- U[order(U)]
    oU <- c(0,oU,1)
    g  <- diff(oU)
    g
}
#-------------------------------------------------------------------------
#pb <- tkProgressBar(max=B)
#progress <- function(i) setTkProgressBar(pb, i)
#progress <- function(n) cat(sprintf("task %d is complete\n", n))
#------------------------------------------------------------------------
#opts <- list(progress=progress)
#--------------------------------- --------------------------------------- 
if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
     data <- if ("data"%in%names(obj$call)) eval(obj$call$data)
               else if (!is.null(data)) data
               else stop("data are not defined")  
if (!is.data.frame(data)) stop("data is not a data.frame")
     N <- nN <- dim(data)[1]     # length of data
if (!is.null(newdata))
{
  if(!is.data.frame(newdata))  stop("newdata should be a data.frame") 
    nN <- dim(newdata)[1]
} else newdata <- data  
## those are the original parameters and saved for the output
     fitted_par <- predictAll(obj, newdata=newdata, output="matrix")
## Not that if no newdata you want to predict in the original data to get 
## estimates CI for the fitted values 
        thecall <- as.call(obj$call) # the call
thecall$control <- gamlss.control(trace=FALSE) # suppressing tracing
## getting the length of the coefficients 
      orig.coef <- unlist(coefAll1(obj, deviance=TRUE))
              Q <- length(orig.coef)
## FOREACH HERE              
      res <- foreach(i = 1:B, .packages="gamlss", .errorhandling = "remove",
                     .inorder = FALSE)%dopar%{
             bootstrap_weights <- rFun(N)*N #
               thecall$weights <- bootstrap_weights
                       model_1 <- eval(thecall) # how to catch?
      list(coef=coefAll1(model_1, deviance = TRUE), par=predictAll(model_1, newdata=newdata, output="matrix"))# output
                     }
##   FOREACH ENDS    
           newB <- length(res)
if (newB!=B) warning(cat(B-newB, "simulations failed only", newB, "were finished \n"))  
           COEF <- matrix(0, ncol=Q, nrow=newB)
if ("mu" %in% obj$par)       MU <- matrix(0, ncol=newB, nrow=nN)
if ("sigma" %in% obj$par) SIGMA <- matrix(0, ncol=newB, nrow=nN)
if ("nu" %in% obj$par)       NU <- matrix(0, ncol=newB, nrow=nN)
if ("tau" %in% obj$par)     TAU <- matrix(0, ncol=newB, nrow=nN)
## maybe I should do this with foreach      
for (i in 1:newB)
    {
      COEF[i,] <- unlist(res[[i]][["coef"]])
      if ("mu" %in% obj$par)        MU[,i] <- unlist(res[[i]][["par"]])[,"mu"]
      if ("sigma" %in% obj$par)  SIGMA[,i] <- unlist(res[[i]][["par"]])[,"sigma"]
      if ("nu" %in% obj$par)        NU[,i] <- unlist(res[[i]][["par"]])[,"nu"]
      if ("tau" %in% obj$par)      TAU[,i] <- unlist(res[[i]][["par"]])[,"tau"]
    }
parameters <- switch(length(obj$par),
                     list(mu=MU),
                     list(mu=MU, sigma=SIGMA),
                     list(mu=MU, sigma=SIGMA, nu=NU),
                     list(mu=MU, sigma=SIGMA, nu=NU, tau=TAU))
colnames(COEF) <-   names(unlist(res[[1]][["coef"]]))
# is any column of COEF has all valus NA the omit it 
if (any(is.na(COEF)))
{
  ind <- rep(0,Q)
for (i in 1:Q)  ind[i] <- ifelse(all(is.na(COEF[,i])),-i, ind[1])
 COEF <- COEF[,ind]
}
#################################################################
#                           OUTPUT                              #
#################################################################  
    out <- list(boot = COEF, 
                   B = B,
               trueB = newB,
               param = parameters,
           orig.coef = orig.coef,
          orig.param = fitted_par, 
           orig.call = obj$call) 
class(out) <- list("Bayesian.boot")
out 
}
#########################################################################
#########################################################################
#########################################################################
# summary 
summary.Bayesian.boot <- function(object, 
                    what=c("coef", "mu", "sigma", "nu", "tau"),
                    show = TRUE,
                    ...)
{
  what <- match.arg(what)
if (what=="coef")  
 {
if (show)   cat("\n Coefficients from Baysian.boot() taken from:", object$trueB, "simulations \n")  
   tcoef <-  t(apply(object$boot,2,"quantile", probs=c(0.025, 0.50, 0.975), na.rm=TRUE))
   mcoef <-  apply(object$boot,2,"mean", na.rm=TRUE )
     out <- cbind(mean=mcoef, tcoef, mle=object$orig.coef[names(object$orig.coef)%in%names(mcoef)]) 
     if (show) return(printCoefmat(out,...)) else return(invisible(out))
 }
if (what=="mu")  
 {
   if ("mu" %in% colnames(object$orig.param))
   {
     if (show)      cat("\n parameter mu from Baysian.boot() taken from:", object$trueB, "simulations \n")  
     tmu <-  t(apply(object$param[["mu"]],1,"quantile", probs=c(0.025, 0.50, 0.975), na.rm=TRUE))
     mmu <-  apply(object$param[["mu"]],1,"mean", na.rm=TRUE )
     out <- cbind(mean=mmu, tmu, mu=object$orig.param[,"mu"]) 
     #printCoefmat(out,...)
if (show) return(printCoefmat(out,...)) else return(invisible(out))
   } else stop("mu is not a parameter in the model")
 }  
if (what=="sigma")  
 {
 if ("sigma" %in% colnames(object$orig.param))
  {
   if (show)   cat("\n parameter sigma from Baysian.boot() taken from:", object$trueB, "simulations \n")  
    tsigma <-  t(apply(object$param[["sigma"]],1,"quantile", probs=c(0.025, 0.50, 0.975), na.rm=TRUE))
    msigma <-  apply(object$param[["sigma"]],1,"mean", na.rm=TRUE )
       out <- cbind(mean=msigma, tsigma, sigma=object$orig.param[,"sigma"]) 
       if (show) return(printCoefmat(out,...)) else return(invisible(out))
  } else stop("sigma is not a parameter in the model")
 }    
if (what=="nu")  
 {
 if ("nu" %in% colnames(object$orig.param))
  {
   if (show)    cat("\n parameter nu from Baysian.boot() taken from:", object$trueB, "simulations \n")  
    tnu <-  t(apply(object$param[["nu"]],1,"quantile", probs=c(0.025, 0.50, 0.975), na.rm=TRUE))
    mnu <-  apply(object$param[["nu"]],1,"mean", na.rm=TRUE )
    out <- cbind(mean=mnu, tnu, nu=object$orig.param[,"nu"]) 
    if (show) return(printCoefmat(out,...)) else return(invisible(out))
   } else stop("nu is not a parameter in the model")
 }
if (what=="tau")  
 {
  if ("nu" %in% colnames(object$orig.param))
  {
    if (show)  cat("\n parameter tau from Baysian.boot() taken from:", object$trueB, "simulations \n")  
    ttau <-  t(apply(object$param[["tau"]],1,"quantile", probs=c(0.025, 0.50, 0.975), na.rm=TRUE))
    mtau <-  apply(object$param[["tau"]],1,"mean", na.rm=TRUE )
     out <- cbind(mean=mtau, ttau, tau=object$orig.param[,"tau"]) 
     if (show) return(printCoefmat(out,...)) else return(invisible(out))
   } else stop("tau is not a parameter in the model") 
 }  
}
##########################################################################
##########################################################################
##########################################################################
# plot (only we plot coefficients)
plot.Bayesian.boot <- function(x, coef=1, add=FALSE,  ...)
{
  if (add) {lines(density(x$boot[,coef]), ... )}
  else {
    plot(density(x$boot[,coef]), main=colnames(x$boot)[coef], ... )
    rug(x$boot[,coef])
  }
  abline(v=x$orig.coef[coef], col="gray")  
  abline(v=quantile(x$boot[,coef], c(0.025, .5, .975)), col=gray(.9))
}
###########################################################################
###########################################################################
###########################################################################
#  iii) print.Bayesian.print
print.Bayesian.boot <- function(x,...)
{
  cat("\n Results from Baysian.boot() taken from:", x$trueB, " simulations \n")
  cat("\n The simulated coefficient values are in: object$boot", "\n") 
  cat("\n The simulated distribution parametets values are in: object$param", "\n") 
  cat("\n Simulation was done from the original call: \n")
  ccc<-x$orig.call; print(ccc)
}  
######################################################################################
######################################################################################
######################################################################################