\name{fitQR}
\alias{fitQR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting Quantile Regression Using the SEP3 Distribution}
\description{
This is a function to demonstrate that quantile regression can be fitted using the \code{SEP3} distribution in GAMLSS. The function is done for demonstration purpose and it is not trying to compete with a proper implementation of quantile regression like
Fasiolo (2021). 
}
\usage{
fitQR(formula, data, quantile = c(0.5), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula for the quantile regression. It fit the \code{mu} for the \code{SEP3} distribution.} 
  \item{data}{the data frame}
  \item{quantile}{the required quantiles}
  \item{\dots}{arguments to pass to the \code{gamlss()}{} function}
}
\details{
The function is fitting multiple GAMLSS fits using the \code{SEP3} distribution 
with  \code{sigma} parameter fixed at 1, \code{tau} parameter fixed at 1 
and the \code{nu} parameters  set to \code{((1-quantile)/quantile)^.5 }.  
}
\value{
Multiple fitted GAMLSS objects.
}
\references{

Fasiolo M, Wood SN, Zaffran M, Nedellec R, Goude Y (2021). "qgam: Bayesian
  Nonparametric Quantile Regression Modeling in R." \emph{Journal of Statistical
  Software}, 100(9), 1-31. doi:10.18637/jss.v100.i09
  <https://doi.org/10.18637/jss.v100.i09>.
  
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape, (with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby, R. A., Stasinopoulos, D. M.,  Heller, G. Z.,  and De Bastiani, F. (2019)
	\emph{Distributions for modeling location, scale, and shape: Using GAMLSS in R}, Chapman and Hall/CRC, \doi{10.1201/9780429298547}. An older version can be found in \url{https://www.gamlss.com/}.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \doi{10.18637/jss.v023.i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  \doi{10.1201/b21973}

Stasinopoulos, M. D., Rigby, R. A., and De Bastiani F., (2018) GAMLSS: a distributional regression approach, \emph{Statistical Modelling}, Vol. \bold{18}, 
pp, 248-273, SAGE Publications Sage India: New Delhi, India. \doi{10.1177/1471082X18759144}




(see also \url{https://www.gamlss.com/}).  
}
\author{Mikis Stasinopoulos}


\seealso{
\code{\link{SEP3}}
}
\examples{
dbmi10 <- dbbmi[db$age>10,]
dim(dbmi10)
registerDoParallel(cores = 2)
m01 <-fitQR(bmi~age, data=dbmi10, quantile = c(0.01, 0.1,.25,.5,.75,0.90,0.99))
stopImplicitCluster()
require(ggplot2)
ggplot(data=dbmi10, aes(x=age, bmi))+
  geom_point(col="steelblue4")+
  geom_line(aes(x=age, y=fitted(m01[[1]])), lwd=1, col="yellow")+
  geom_line(aes(x=age, y=fitted(m01[[2]])), lwd=1, col="green")+
  geom_line(aes(x=age, y=fitted(m01[[3]])), lwd=1, col="red")+
  geom_line(aes(x=age, y=fitted(m01[[4]])), lwd=1, col="black")+
  geom_line(aes(x=age, y=fitted(m01[[5]])), lwd=1, col="red")+
  geom_line(aes(x=age, y=fitted(m01[[6]])), lwd=1, col="green")+
  geom_line(aes(x=age, y=fitted(m01[[7]])), lwd=1, col="yellow")
}

\keyword{regression}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
