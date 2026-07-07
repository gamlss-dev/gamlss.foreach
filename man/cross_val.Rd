\name{cross_val}
\alias{cross_val}
\alias{cross_val2}
\alias{mult_cross_val}
\alias{from_list2mat}


\title{
Cross Validation for GAMLSS Models
}
\description{
There are several functions to help with K-th fold Cross Validation 

1) \code{cross_val()}: the function takes a \code{gamlss} object and produce an  \code{CV_gamlss} object with components:

      \itemize{

  \item \code{model}: the original fitted  \code{gamlss} model

  \item \code{orig.resid}: the residuals of the original model 

  \item \code{orig.incr}:  the deviance inrements  of the original model 

 \item \code{cv.resid}: the residuals of the cross validation model

  \item \code{cv.incr}:  the deviance inrements  of the cross validation model

  \item \code{cv.deviance}:  the deviance of the cross validation model i.e. the sum of
         \code{cv.incr}
}         
2) \code{cross_val2()}: the function takes a \code{gamlss2} object and produce an  \code{CV_gamlss2} object with identical components  as above.


3) The function \code{mult_cross_val()} cross validates a \code{gamlss} or \code{gamlss2} fitted model several times to produce a list of mupliple cross valadations. 

4) The function \code{from_list2mat()} takes a \code{mult_CV_gamlss} object which is a lsist of several cross validated \code{gamlss} objects and creates a matrix of all of the residuals from each object.  The matrix can be used then for summary statsirics or plots.
}
\usage{
cross_val(model, data, K.fold = 10, rand = NULL, print = TRUE)

cross_val2(model, data, K.fold = 10, rand = NULL, print = TRUE)

mult_cross_val(model, R = 10, data, K.fold = 10)

from_list2mat(object, type = c("resid", "incr"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A \code{gamlss} or \code{gamlss2} model}
  \item{data}{The data usen in the analysis of the \code{model} }
  \item{K.fold}{How many splits of the data we should do}
  \item{rand}{an alternative way to do the splits of the data (for reprodusability)}
  \item{print}{whether to prin that finish the the ith section of the data}
  \item{R}{how many times to repeat the cross validation in \code{mult_cross_val}}
  \item{object}{An object created from \code{mult_cross_val}}
  \item{type}{whether the resulting matrix should residuals or deviance increments}
  \item{...}{futher arguments}
}
\details{
The functions \code{cross_val()} and \code{cross_val2()} produce a  \code{CV_gamlss} object 

The function \code{mult_cross_val()} produces multiple cross validations and creates the object
 \code{mult_CV_gamlss}.
 
The function  \code{mult_cross_val()} produces multiple cross validations and creates the object
 \code{mult_CV_gamlss}.
}
\value{
The functions \code{cross_val()} and \code{cross_val2()} produce a  \code{CV_gamlss} object with a \code{resid()} as a generic function. 


 The function \code{mult_cross_val()} produces multiple cross validations and creates the object 
 \code{mult_CV_gamlss}which has as generic the function \code{summary()}. The information in produce with \code{summary()} can be use in plotting or summarize the performance of the cross validation.
}
\references{

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape, (with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/v23/i07/}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  

(see also \url{https://www.gamlss.com/}).
}

\author{Mikis Stasinopoulos \email{d.stasinopoulos@gre.ac.uk} and  Bob Rigby}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\seealso{
\code{\link[gamlss]{gamlssCV}}, \code{\link[gamlss]{CV}}}


\examples{

\donttest{
library(gamlss)
library(doParallel)
m1 <- gamlss(rent~area+I(yearc^2)+location, data=rent99, family=BCT)

registerDoParallel(cores = 10)
mCV=cross_val(m1)
stopImplicitCluster()
plot(resid(mCV))
registerDoParallel(cores = 10)
MMCV<- mult_cross_val(m1) 
stopImplicitCluster()

matplot(from_list2mat(MMCV, type="incr"))
matplot(from_list2mat(MMCV, type="resid"))
pairs(from_list2mat(MMCV, type="resid"))
library(matrixStats)
  mm = rowMeans2(summary(PP, type="resid"))
 sds = rowSds(summary(PP, type="resid"))
 library(ggplot2)
# 
 df <- data.frame(
   index = 1:3082,
   mean = mm,
   sd = sds
 )
ggplot(df, aes(index, mean)) +
   geom_point(size = 3) +
   geom_errorbar(aes(ymin = mean - sd,
                     ymax = mean + sd),
                 width = 0.2) +
   labs(x = "Group", y = "Mean")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
 \keyword{regression}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
