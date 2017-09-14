\name{gamlss.inf-package}
\alias{gamlss.inf-package}
\alias{gamlss.inf}
\docType{package}
\title{
Models for Mixed (Inflated and Adjusted) Response Variables.
}
\description{
This package allows mixed distribution fitting in GAMLSS. A mixed distribution is  one containing both continuous and discrete parts, see Chapter 5 of Stasinopoulos et al. (2017).   There are some mixed distribution in the GAMLSS implementation in \R{} like the \code{BEINF, BEINF0, BEINF1} for data defined on [0,1] or \code{ZAGA, ZAIG} for data defined on a positive real line but the choice is very limited. This package  enhance the availability of mixed distribution within the GAMLSS framework.

For historical reasons the authors use the terminology "Inflated" for models on [0,1],
"Adjusted" for models on [0,Inf]. We will follow the same terminology here.
So this package allows the fit of an inflated GAMLSS model when the response variable is defined in the intervals [0,1), (0,1] and [0,1] and  the fit of zero adjusted models when the response variable is defined in the positive real line, (but where there are zeros in the data). 

For models with inflated proportion response variables the package provides up to two extra parameters, a mass point at zero and a mass point at one. Adding an extra inflation point at zero (or at one), is equivalent to fit two separate GAMLSS models, a GAMLSS model with a continuous distribution defined at the interval (0,1), and a logit model for zero (or ones). 
When both zero and one are present, i.e. [0,1],  a multinomial model  is needed  to fit the non-(0,1) part. 

For the zero adjusted models with a response defined on the positive real line (but where zeros exist in the data), the actual fitting  can be achieved by fitting two separate GAMLSS models one with a distribution on the real positive line and one binomial model for the zeros and non-zeros.      

This package uses the two models fitting procedures but the resulting fitted object behaves like a typical GAMLSS object so a lot of standard GAMLSS diagnostics can be used with it.   It also allow transformed or truncated \code{gamlss.family} distributions to be used  for the continuous part of the distribution therefore giving a great flexibility of the type of mixed distribution to be used. 


The function \code{gamlssInf0to1()} and \code{gamlssAdj()} can fit the Inflated and Adjusted  models, respectively, and they are described in the two vignettes available with the package. 

}
\details{
\tabular{ll}{
Package: \tab gamlss.inf\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-09-2\cr
License: \tab GPL-2 | GPL-3\cr
}
}
\author{
Marco Enea <marco.enea@unipa.it>, Mikis Stasinopoulos <mikis.stasinopoulos@gamlss.org>, 
Bob Rigby and Abu Hossain \cr 

Maintainer: Marco Enea <marco.enea@unipa.it> 

}
\references{
Hossain, A., Stasinopoulos, M., Rigby, R. and Enea, M. (2015). Centile estimation for a proportion response variable. \emph{Statistics in Medicine},  doi: 10.1002/sim.6748. \cr

Ospina R. and Ferrari S. L. P. (2010) Inflated beta distributions, \emph{Statistical Papers}, 
    \bold{23},  111-126. 

Rigby R.A. and Stasinopoulos D.M. (2005). Generalized additive models for location, scale 
and shape,(with discussion), Appl. Statist., 54, part 3, pp 507-554. \cr

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R. Journal of Statistical Software, Vol. 23, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}. 

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC. 
\url{https://www.crcpress.com/Flexible-Regression-and-Smoothing-Using-GAMLSS-in-R/Stasinopoulos-Rigby-Heller-Voudouris-Bastiani/p/book/9781138197909}.

}

\seealso{
\code{\link[gamlss]{gamlss}}
}

\keyword{ package }

\examples{

# An artificial example using simulated data 0n [0,1)

# Firstly, we use function gen.Family() to create 
# the logit skew student t (logitSST) distribution which is defined in the (0,1) 
# interval. Then we use function gen.Inf0to1() to create the 0-inflated logitSST 
# distribution defined in [0,1).

gen.Family("SST", "logit")
gen.Inf0to1("logitSST","Zero") 

#now we can generate the data and run the model 

set.seed(10)
Y <- rlogitSSTInf0(500,mu=0.5,sigma=0.7,nu=0.5,tau=5,xi0=0.5,log=FALSE)
m1 <- gamlssInf0to1(y=Y,mu.formula=~1,
                    sigma.formula=~1,
                    nu.formula=~1,
                    tau.formula=~1,
                    xi0.formula=~1,
                    family=logitSST,
                    trace = TRUE)
summary(m1)

}
