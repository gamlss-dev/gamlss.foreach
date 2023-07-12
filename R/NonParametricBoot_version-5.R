###########################################################################
###########################################################################
###########################################################################
# FUNCTION classical Boot
# mikis revision 27-10-2022
########################################################################### 
NonParametricBoot <- function(obj, data=NULL, B=100, newdata=NULL)
{
if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
  datA <- if ("data"%in%names(obj$call)) eval(obj$call$data)
          else if (!is.null(data)) data
          else stop("data are not defined")  
if (!is.data.frame(datA)) stop("data is not a data.frame")
     N <- nN <- dim(datA)[1]     # length of data
     R <-  dim(datA)[2] 
if (!is.null(newdata))
  {
    if(!is.data.frame(newdata))  stop("newdata should be a data.frame") 
    nN <- dim(newdata)[1]
  } else newdata <- datA  
## those are the original parameters and saved for the output
     fitted_par <- predictAll(obj, newdata=newdata, output="matrix")
## Not that if no newdata you want to predict in the original data to get 
## estimates CI for the fitted values 
        thecall <- as.call(obj$call) # the call
#thecall$control <- gamlss.control(trace=FALSE) # supressing tracing
## getting the length of the coefficients 
      orig.coef <- unlist(coefAll1(obj, deviance=TRUE))
     names.coef <-  names(orig.coef)
              Q <- length(orig.coef) 
## FOREACH HERE   
  res <- foreach(i = 1:B, .packages="gamlss", .errorhandling = "remove",
                    .export = c("coefAll1"), .inorder = FALSE)%dopar%
    {
               ii <- sample(N,N,replace=T) 
        boot_data <- as.data.frame(datA[ii,])
names( boot_data) <- names(datA)
             mooo <- update(obj, data= boot_data, trace=F ) 
  list(coef=coefAll1(mooo, deviance = TRUE), par=predictAll(mooo, data=boot_data, newdata=newdata, output="matrix")) 
    }
##   FOREACH ENDS  
            newB <- length(res)
if (newB!=B) warning(cat(B-newB, "simulation failed only", newB, "were finished \n"))        
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
# #################################################################
# #                           OUTPUT                              #
# #################################################################  
out <- list(boot = COEF, 
               B = B,
           trueB = newB,
           param = parameters,
       orig.coef = orig.coef,   
      orig.param = fitted_par, 
       orig.call = obj$call) 
class(out) <- list("NonParametric.Boot")
out 
}
###########################################################################
###########################################################################
########################################################################### 
summary.NonParametric.Boot <- function(object, 
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
if (show) return(printCoefmat(out,...)) else return(invisible(out))
    } else stop("mu is not a parameter in the model")
  }  
  if (what=="sigma")  
  {
    if ("sigma" %in% colnames(object$orig.param))
    {
      if (show) cat("\n parameter sigma from Baysian.boot() taken from:", object$trueB, "simulations \n")  
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
      if (show) cat("\n parameter nu from Baysian.boot() taken from:", object$trueB, "simulations \n")  
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
      if (show)   cat("\n parameter tau from Baysian.boot() taken from:", object$trueB, "simulations \n")  
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
########################################################################################
########################################################################################
########################################################################################
# # plot 
plot.NonParametric.Boot <- function(x, par=1, add=FALSE, ...)
{
  if (add) {lines(density(x$boot[,par]), ... )}
  else {
    plot(density(x$boot[,par]), main=colnames(x$boot)[par], ... )
    rug(x$boot[,par])
  }
  
  abline(v=x$orig.coef[par], col="gray")  
  abline(v=quantile(x$boot[,par], c(0.025, .5, .975)), col=gray(.9))
}
########################################################################## 
########################################################################## 
##########################################################################
print.NonParametric.Boot <- function(x,...)
{
  cat("\n Results from classical.boot() taken from:", x$trueB, " simulations \n")
  cat("\n The simulated coefficient values are in: object$boot", "\n")
  cat("\n The simulated distribution parametets values are in: object$param", "\n") 
  cat("\n Simulation was done from the original call: \n")
  ccc<-x$orig.call; print(ccc)
}  

##########################################################################
##########################################################################
##########################################################################
