################################################################################
################################################################################
################################################################################
################################################################################
# function 1
#require(gamlss)
gamlss_fits <- function (formula = NULL, 
                   sigma.formula = ~1, 
                      nu.formula = ~1, 
                     tau.formula = ~1, 
                            data = NULL,
                         newdata = NULL, 
                          family = NO,
                         weights,
                         control = gamlss.control(trace = FALSE), 
                             ...) 
{
################################################################################  
rqres <- function(pfun = "pNO", type = c("Continuous", "Discrete", 
         "Mixed"), censored = NULL, ymin = NULL, mass.p = NULL, 
          prob.mp = NULL, y = y, ...) {
  }
  body(rqres) <- eval(quote(body(rqres)), envir = getNamespace("gamlss"))
################################################################################  
  sys_call <- sys.call()
if (is.null(data)) 
    stop("data should be set here")
if (is.null(newdata)) 
    stop("newdata should be set")
    # dataTraining <- data
    #    dataValid <- newdata
            Data <- concat(data, newdata, names=c("train", "test"))
           fname <- gamlss.dist::as.gamlss.family(family)
            dfun <- paste("d", fname$family[[1]], sep = "")
            pfun <- paste("p", fname$family[[1]], sep = "")
            lpar <- length(fname$parameters)
           dtype <- fname$type
if (is.null(formula)) 
    stop("no formula is set in gamlssVGD")
if (is.null(data)) 
    stop("the data argument is needed in gamlssVGD")
    #        browser()
    # m1 <-  eval(expression(gamlss::gamlss(formula = formula, sigma.formula = sigma.formula, 
    #               nu.formula = nu.formula, tau.formula = tau.formula, 
    #               weights=ftest,
    #               data = Data, family = family, control = control)), envir=.GlobalEnv)
    #  nfitted <- predictAll(m1, newdata = newdata, data = data, output="matrix")
    
   #  dim1newdata <- dim(dataValid)[1]
#   }
#   else {
    m1 <- gamlss(formula = formula, sigma.formula = sigma.formula,
                 nu.formula = nu.formula, tau.formula = tau.formula,
                 weights=weights,
                 data = data, family = family, control = control,
                 ...)
    nfitted <- predictAll(m1, newdata = newdata, data = data)
    dim1newdata <- dim(newdata)[1]
#   }
#   if (fname$family[1] %in% .gamlss.bi.list) {
#     if (NCOL(nfitted$y) == 1) {
#       y1 <- nfitted$y
#       bd <- nfitted$bd
#     }
#     else {
#       bd <- nfitted$y[, 1] + nfitted$y[, 2]
#       y1 <- nfitted$y[, 1]
#     }
#   }
#   else {
#     y1 <- nfitted$y
#   }
#   if (lpar == 1) {
#     if (fname$family[[1]] %in% .gamlss.bi.list) {
#       devi <- call(dfun, x = y1, mu = nfitted$mu, bd = bd, 
#                    log = TRUE)
#       ures <- call("rqres", pfun = pfun, type = dtype, 
#                    ymin = fname$rqres[[1]][["ymin"]], y = y1, bd = bd, 
#                    mu = nfitted$mu)
#     }
#     else {
#       devi <- call(dfun, x = y1, mu = nfitted$mu, log = TRUE)
#       ures <- call("rqres", pfun = pfun, type = dtype, 
#                    ymin = fname$rqres[[1]][["ymin"]], y = y1, mu = nfitted$mu)
#     }
#   }
#   else if (lpar == 2) {
#     if (fname$family[[1]] %in% .gamlss.bi.list) {
#       devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma, 
#                    bd = bd, log = TRUE)
#       ures <- call("rqres", pfun = pfun, type = dtype, 
#                    ymin = fname$rqres[[1]][["ymin"]], y = y1, mu = nfitted$mu, 
#                    sigma = nfitted$sigma, bd = bd)
#     }
#     else {
#       devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma, 
#                    log = TRUE)
#       ures <- call("rqres", pfun = pfun, type = dtype, 
#                    ymin = fname$rqres[[1]][["ymin"]], y = y1, mu = nfitted$mu, 
#                    sigma = nfitted$sigma)
#     }
#   }
#   else if (lpar == 3) {
#     if (fname$family[[1]] %in% .gamlss.bi.list) {
#       devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma, 
#                    nu = nfitted$nu, bd = bd, log = TRUE)
#       ures <- call("rqres", pfun = pfun, type = dtype, 
#                    ymin = fname$rqres[[1]][["ymin"]], y = y1, mu = nfitted$mu, 
#                    sigma = nfitted$sigma, nu = nfitted$nu, bd = bd)
#     }
#     else {
#       devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma, 
#                    nu = nfitted$nu, log = TRUE)
#       ures <- call("rqres", pfun = pfun, type = dtype, 
#                    ymin = fname$rqres[[1]][["ymin"]], y = y1, mu = nfitted$mu, 
#                    sigma = nfitted$sigma, nu = nfitted$nu)
#     }
#   }
#   else {
#     if (fname$family[[1]] %in% .gamlss.bi.list) {
#       devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma, 
#                    nu = nfitted$nu, tau = nfitted$tau, bd = bd, 
#                    log = TRUE)
#       ures <- call("rqres", pfun = pfun, type = dtype, 
#                    ymin = fname$rqres[[1]][["ymin"]], y = y1, mu = nfitted$mu, 
#                    sigma = nfitted$sigma, nu = nfitted$nu, tau = nfitted$tau, 
#                    bd = bd)
#     }
#     else {
#       devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma, 
#                    nu = nfitted$nu, tau = nfitted$tau, log = TRUE)
#       ures <- call("rqres", pfun = pfun, type = dtype, 
#                    ymin = fname$rqres[[1]][["ymin"]], y = y1, mu = nfitted$mu, 
#                    sigma = nfitted$sigma, nu = nfitted$nu, tau = nfitted$tau)
#     }
#   }
#   Vresid <- eval(ures)
#   dev.incr <- -2 * eval(devi)
#   dev <- sum(dev.incr)
#   m1$VGD <- dev
#   m1$IncrVGD <- dev.incr
#   m1$predictError <- dev/dim1newdata
#   m1$residVal <- Vresid
#   m1$call$family <- sys_call$family
#   class(m1) <- c("gamlssVGD", "gamlss", "gam", "glm", "lm")
#   m1
 } 

################################################################################
################################################################################
################################################################################
################################################################################
# function 2
# gamlss_Kfold <- function(formula = NULL, 
#                   sigma.formula = ~1, 
#                      nu.formula = ~1, 
#                     tau.formula = ~1, 
#                            data = NULL, 
#                          family = NO, 
#                         control = gamlss.control(trace = FALSE), 
#                               K = 10, 
#                        set.seed = 123, ...
#     ) 
# {
#       sys_call <- sys.call()
# if (is.null(data)) 
#       stop("data should be set here")
#              N <- dim(data)[1]
#       set.seed <- set.seed
#           rand <- data_Kfold_weights(data, K=K)
#         K.fold <- K
#     devianceCV <- rep(0, K.fold)
#        residCV <- rep(0, N)
#      #    i <- sort(unique(rand))
#            fn <- function(i, ...) {
#       cat("fold ", i, "\n", sep = "")
#       learn <- gamlssVGD(formula = formula, sigma.formula = sigma.formula, 
#                          nu.formula = nu.formula, tau.formula = tau.formula, 
#                          family = family, control = control, data = data[rand != 
#                                                                            i, ], newdata = data[rand == i, ], ...)
#       list(CV = learn$VGD, resid = learn$residVal)
#     }
#     CV <- if (ncpus > 1L && (have_mc || have_snow)) {
#       if (have_mc) {
#         parallel::mclapply(i, fn, mc.cores = ncpus)
#       }
#       else if (have_snow) {
#         list(...)
#         if (is.null(cl)) {
#           cl <- parallel::makeForkCluster(ncpus)
#           if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
#             parallel::clusterSetRNGStream(cl)
#           res <- parallel::parLapply(cl, i, fn)
#           parallel::stopCluster(cl)
#           res
#         }
#         else t(parallel::parLapply(cl, i, fn))
#       }
#     }
#     else lapply(i, fn)
#     cv <- rep(0, K.fold)
#     for (i in 1:K.fold) {
#       residCV[rand == i] <- CV[[i]]$resid
#       cv[i] <- CV[[i]]$CV
#     }
#     out <- list(CV = sum(cv), allCV = cv, residCV = residCV)
#     class(out) <- "gamlssCV"
#     out
#   }  
#   if (!is(obj,"gamlss")) stop("the obj should be a gamlss object")
#   if (!grepl("$", deparse(substitute(as.time)), fixed=TRUE))
#   {
#     as.time <- eval(substitute(as.time), envir=as.environment(data))
#   }
#      pdf <- obj$family[1]
#   binom  <- pdf%in%gamlss::.gamlss.bi.list # whether binomial
#       fn <- eval(parse(text=pdf))()$G.dev.inc # get pdf
#     lpar <- eval(parse(text=pdf))()$nopar
#   if (binom) {bd <- obj$bd ; Y <- obj$y}
#        N <- nrow(data)
#     maxK <- N-window
#        i <- 1 :  maxK
#     yhat <- foreach(i = 1 :  maxK, .packages="gamlss", .combine = rbind) %dopar% {
#      um1 <- try(update(obj, data = data[i:(window+(i-1)), ], 
#                     mu.start = obj$mu.fv[i:(window+(i-1))],
#                  sigma.start = obj$sigma.fv[i:(window+(i-1))],
#                     nu.start = obj$nu.fv[i:(window+(i-1))],
#                    tau.start = obj$tau.fv[i:(window+(i-1))],
#                    trace=FALSE))
#      if (any(class(um1)%in%"try-error")) stop("the algorithm failed in iteration", i, "\n")
#      suppressMessages(pp<- predictAll(um1, data = data[i:(window+(i-1)), ],  newdata = data[window+i,], output="matrix"))
#    # if) 
#     if (binom)
#     {
#       DevIncr <- switch(lpar, 
#                         fn( Y[i], mu = pp[,"mu"], bd=bd[i]),   # 1
#                         fn( Y[i], mu = pp[,"mu"],              # 2
#                                sigma = pp[,"sigma"], bd=bd[i]),                        
#                         fn( Y[i], mu = pp[,"mu"],              # 3
#                                sigma = pp[,"sigma"],
#                                   nu = pp[,"nu"], bd=bd[i]),
#                         fn( Y[i], mu = pp[,"mu"],              # 4
#                                sigma = pp[,"sigma"],  
#                                   nu = pp[,"nu"],
#                                  tau = pp[,"tau"],bd=bd[i]))
#     } else
#     {
#       DevIncr <- switch(lpar, 
#                         fn( obj$y[i], mu =  pp[,"mu"]),       # 1
#                         fn( obj$y[i], mu =  pp[,"mu"],        # 2
#                             sigma =  pp[,"sigma"]),  # or fitted(obj, "sigma")                    
#                         fn( obj$y[i], mu =  pp[,"mu"],        # 3
#                             sigma =  pp[,"sigma"],
#                             nu = pp[,"nu"]),
#                         fn( obj$y[i], mu =  pp[,"mu"],        # 4
#                                    sigma =  pp[,"sigma"],  
#                                      nu = pp[,"nu"],
#                                     tau = pp[,"tau"]))
#     }  
#     yhat <- cbind(pp, GDincr=DevIncr) 
#     }
#           cn <- colnames(yhat)
#            M <- matrix(unlist(yhat), ncol=length(cn))
# colnames(M)  <- cn
# rownames(M)  <- if (!is.null(as.time)) as.time[(window+1):N]  else rownames(data[(window+1):N,])
#   return(M)
# }

################################################################################
################################################################################
################################################################################
################################################################################
# gamlssCV <- function (formula = NULL, sigma.formula = ~1, nu.formula = ~1, 
#           tau.formula = ~1, data = NULL, family = NO, control = gamlss.control(trace = FALSE), 
#           K.fold = 10, set.seed = 123, rand = NULL, parallel = c("no", 
#                                                                  "multicore", "snow"), ncpus = 1L, cl = NULL, ...) 
# {
#   if (missing(parallel)) 
#     parallel <- "no"
#   parallel <- match.arg(parallel)
#   have_mc <- have_snow <- FALSE
#   if (parallel != "no" && ncpus > 1L) {
#     if (parallel == "multicore") 
#       have_mc <- .Platform$OS.type != "windows"
#     else if (parallel == "snow") 
#       have_snow <- TRUE
#     if (!have_mc && !have_snow) 
#       ncpus <- 1L
#     loadNamespace("parallel")
#   }
#   sys_call <- sys.call()
#   if (is.null(data)) 
#     stop("data should be set here")
#   N <- dim(data)[1]
#   set.seed <- set.seed
#   rand <- if (is.null(rand)) 
#     sample(K.fold, N, replace = TRUE)
#   else rand
#   if (length(rand) != N) 
#     stop("the length of the rand should be equal to data")
#   K.fold <- length(unique(rand))
#   CV <- rep(0, K.fold)
#   residCV <- rep(0, N)
#   i <- sort(unique(rand))
#   fn <- function(i, ...) {
#     cat("fold ", i, "\n", sep = "")
#     learn <- gamlssVGD(formula = formula, sigma.formula = sigma.formula, 
#                        nu.formula = nu.formula, tau.formula = tau.formula, 
#                        family = family, control = control, data = data[rand != 
#                                                                          i, ], newdata = data[rand == i, ], ...)
#     list(CV = learn$VGD, resid = learn$residVal)
#   }
#   CV <- if (ncpus > 1L && (have_mc || have_snow)) {
#     if (have_mc) {
#       parallel::mclapply(i, fn, mc.cores = ncpus)
#     }
#     else if (have_snow) {
#       list(...)
#       if (is.null(cl)) {
#         cl <- parallel::makeForkCluster(ncpus)
#         if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
#           parallel::clusterSetRNGStream(cl)
#         res <- parallel::parLapply(cl, i, fn)
#         parallel::stopCluster(cl)
#         res
#       }
#       else t(parallel::parLapply(cl, i, fn))
#     }
#   }
#   else lapply(i, fn)
#   cv <- rep(0, K.fold)
#   for (i in 1:K.fold) {
#     residCV[rand == i] <- CV[[i]]$resid
#     cv[i] <- CV[[i]]$CV
#   }
#   out <- list(CV = sum(cv), allCV = cv, residCV = residCV)
#   class(out) <- "gamlssCV"
#   out
# }
################################################################################
################################################################################
################################################################################
################################################################################
concat <- function(..., names = NULL) 
  {
  tmp <- list(...)
  if (is.null(names)) 
    names <- names(tmp)
  if (is.null(names)) 
    names <- sapply(as.list(match.call()), deparse)[-1]
  if (any(sapply(tmp, is.matrix) | sapply(tmp, is.data.frame))) {
    len <- sapply(tmp, function(x) c(dim(x), 1)[1])
    len[is.null(len)] <- 1
    data <- rbind(...)
  }
  else 
    {
    len <- sapply(tmp, length)
    data <- unlist(tmp)
  }
 # names1 <- if(grep("test", names[1])==1) c("test", "train") else c("train","test")
   f <- factor(rep(names, len), levels = names)
  namelist <- model.matrix(~f-1)
  return(data.frame(data, namelist))
}
################################################################################
################################################################################
################################################################################
################################################################################
fittedAll <- function (obj, ...) 
{
  out <- list()
  if ("mu" %in% obj$par) 
    out$mu <- fitted(obj, "mu")
  if ("sigma" %in% obj$par) 
    out$sigma <- fitted(obj, "sigma")
  if ("nu" %in% obj$par) 
    out$nu <- fitted(obj, "nu")
  if ("tau" %in% obj$par) 
    out$tau <- fitted(obj, "tau")
if (matrix)   
  return(out)
}
################################################################################
################################################################################
################################################################################
################################################################################
# fittedAll_list <- function (obj, ...) 
# {
#   out <- list()
#   if ("mu" %in% obj$par) 
#     out$mu <- fitted(obj, "mu")
#   if ("sigma" %in% obj$par) 
#     out$sigma <- fitted(obj, "sigma")
#   if ("nu" %in% obj$par) 
#     out$nu <- fitted(obj, "nu")
#   if ("tau" %in% obj$par) 
#     out$tau <- fitted(obj, "tau")
#   if (matrix)   
#     return(out)
# }
################################################################################
################################################################################
################################################################################
################################################################################
# fittedAll_list <- function (obj, ...) 
# {
#   fitted <- predictAll(obj, )
#   
# }
################################################################################
################################################################################
################################################################################
################################################################################
# function to perform a K-fold cross validation to a `gamlss` object
# the output is `CV_gamlss`
cross_val <- function(model, data, K.fold=10, rand=NULL, print=TRUE)
{
if (!is.gamlss(model))  stop(paste("This is not an gamlss object", "\n", ""))
     datA <- if ("data"%in%names(model$call)) eval(model$call$data)
            else if (!is.null(data)) data
           else stop("data are not defined")  
if (!is.data.frame(datA)) stop("data is not a data.frame") 
         N <- nN <- dim(datA)[1]     # length of data
         R <-  dim(datA)[2]          # number of variables 
#fitted_par <- predictAll(obj, output="matrix")
## estimates CI for the fitted values 
  residuals <- rep(0, N)
   dev_Incr <- rep(0, N)
      Index <- 1:N
rand_sample <- if (is.null(rand)) sample(K.fold, N, replace=TRUE) else rand 
## FOREACH HERE   
 res <- foreach(i = 1:K.fold, .packages="gamlss", .errorhandling = "remove",
         .inorder = FALSE)%dopar%
   {
if (print) cat("the", i, "fold", "\n" )     
            ii <- Index[rand_sample==i]
       CV_data <- as.data.frame(datA[ii,])
names(CV_data) <- names(datA)
          mooo <- update(model, data= CV_data, trace=F ) 
list(resid=resid(mooo), dev.Incr = devianceIncr(mooo)) 
   } 
## FOREACH FINISHED
for (i in 1:K.fold)
{
residuals[rand_sample==i] <- res[[i]]$resid
 dev_Incr[rand_sample==i] <- res[[i]]$dev.Incr
}
################################################################################
#                            OUTPUT                                            #
################################################################################  
 out <- list(model=model, 
             K = K.fold,
             orig.resid = resid(model),   
             orig.incr = devianceIncr(model), 
              cv.resid = residuals,   
               cv.incr = dev_Incr, 
          cv.deviance  = sum(dev_Incr),
             orig.call = model$call) 
 class(out) <- list("CV_gamlss")
 invisible(out)
}
# the function sves 
################################################################################
################################################################################
################################################################################
################################################################################
# function to perform K-fold cross validation on a `gamlss2` object
# the output is a `CV_gamlss` object
cross_val2 <- function(model, data, K.fold=10, rand=NULL, print=TRUE)
{
if (!is(model,"gamlss2"))  stop(paste("This is not an gamlss2 object", "\n", ""))
  datA <- if ("data"%in%names(model$call)) eval(model$call$data)
          else if (!is.null(data)) data
          else stop("data are not defined")  
if (!is.data.frame(datA)) stop("data is not a data.frame") 
         N <- nN <- dim(datA)[1]     # length of data
 # R <-  dim(datA)[2]          # number of variables 
  #fitted_par <- predictAll(obj, output="matrix")
  ## estimates CI for the fitted values 
  residuals <- rep(0, N)
   dev_Incr <- rep(0, N)
      Index <- 1:N
rand_sample <- if (is.null(rand)) sample(K.fold, N, replace=TRUE) else rand 
## FOREACH HERE   
  res <- foreach(i = 1:K.fold, .packages="gamlss", .errorhandling = "remove",
                 .inorder = FALSE)%dopar%
    {
      if (print) cat("the", i, "fold", "\n" )     
      ii <- Index[rand_sample==i]
      CV_data <- as.data.frame(datA[ii,])
      names(CV_data) <- names(datA)
      mooo <- update(model, data= CV_data, trace=F ) 
      list(resid=resid(mooo), dev.Incr = deviance(mooo, sum=FALSE)) 
    } 
## FOREACH FINISHED
for (i in 1:K.fold)
  {
    residuals[rand_sample==i] <- res[[i]]$resid
    dev_Incr[rand_sample==i] <- res[[i]]$dev.Incr
  }
################################################################################
#                            OUTPUT                                            #
################################################################################  
  out <- list(model=model, 
              K = K.fold,
              orig.resid = resid(model),   
               orig.incr = deviance(model, sum=FALSE), 
                cv.resid = residuals,   
                 cv.incr = dev_Incr, 
            cv.deviance  = sum(dev_Incr),
               orig.call = model$call) 
  class(out) <- list("CV_gamlss")
  invisible(out)
}
################################################################################
################################################################################
################################################################################
################################################################################
# this is a general function to obtained the deviance increments for any 
# `gamlss2` object. It takes a `gamlss2` object and creates a vector of deviance 
# increments   
# This function is not needed any more because deviavce(., sum=FALSE) 
# is doing similar job
# deviance_Incr <- function(model, newdata=NULL)
# {
#   if (!is(model, "gamlss2")) stop("This works for gamlss2 objects only")
#   DINC <- if (is.null(newdata))
#   {
#     -2*model$family$pdf(model.response(model.frame(model)), fitted(model, type="parameter"), 
#                         log=TRUE)  
#   } else 
#   {
#    -2*model$family$pdf(y=unlist(newdata[response_name(model)]), par=predict(model, type="parameter", newdata= newdata), log=TRUE)
#   } 
#   DINC
# }
################################################################################
################################################################################
################################################################################
################################################################################
# this function takes a `CV_gamlss` or `CV_gamlss2` object and extract 
# CV-residuals or CV_deviance increments  
resid.CV_gamlss <- function(obj, type=c("resid", "incr"))
{
if (!(is(obj, "CV_gamlss"))) stop("not a CV_gamlss object")
  type <- match.arg(type)
res <-   if (type=="resid") obj$cv.resid else obj$cv.incr 
return(res)
}
################################################################################
################################################################################
################################################################################
################################################################################
# this function takes a gamlss or gamlss2 object and creates R reprications
mult_cross_val <- function(model, R=10, data, K.fold=10)
{
  gamlssCVcall <- match.call()
  orig_gamlss <- if (is(model,"gamlss")) TRUE 
                 else FALSE
        CVlist <- list()
           dev <- rep(0, length = R)
for (i in (1:R)) {
    m1 <- try(if (orig_gamlss) cross_val(model = model, K.fold=K.fold, data = data) 
              else            cross_val2(model = model, K.fold=K.fold, data = data)) 
    if (any(class(m1) %in% "try-error")) {
        cat("model=", i, "failed", "\n")
      dev[i] <- NA
      CVlist[[i]] <- NA
      next
    }
    CVlist[[i]] <- m1
    dev[i] <- m1$cv.deviance
    cat("model=", i, "\n")
  }
  class(CVlist) <- "mult_CV_gamlss" 
  return(CVlist)
}
################################################################################
################################################################################
################################################################################
################################################################################
# we need summary and plot for a   "mult_CV_gamlss" object 
from_list2mat <- function(object, type=c("resid", "incr"), ...)
{
  resid_CV_gamlss <- function(obj, type=c("resid", "incr"))
  {
    if (!(is(obj, "CV_gamlss"))) stop("not a CV_gamlss object")
    type <- match.arg(type)
    res <-   if (type=="resid") obj$cv.resid else obj$cv.incr 
    return(res)
  }  
  type <- match.arg(type)  
    pp <- sapply(as.list(object), resid_CV_gamlss, type)
pp
}
################################################################################
################################################################################
################################################################################
################################################################################

# matplot(summary(PP, type="incr"))
# matplot(summary(PP, type="resid"))
# pairs(summary(PP, type="resid"))
#install.packages("matrixStats")
#library(matrixStats)

#rowSums2(summary(PP, type="resid"))
# plot(rowMeans2(summary(PP, type="resid")))
# plot(rowMedians(summary(PP, type="resid")))
# rowSds(summary(PP, type="resid"))
# rowMaxs(summary(PP, type="resid"))
# rowMins(summary(PP, type="resid"))
# 
# 
# plotting i ggplots
#  mm = rowMeans2(summary(PP, type="resid"))
# sds = rowSds(summary(PP, type="resid"))
# library(ggplot2)
# 
# df <- data.frame(
#   index = 1:3082,
#   mean = mm,
#   sd = sds
# )
# 
# ggplot(df, aes(index, mean)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = mean - sd,
#                     ymax = mean + sd),
#                 width = 0.2) +
#   labs(x = "Group", y = "Mean")
################################################################################
################################################################################
################################################################################
################################################################################