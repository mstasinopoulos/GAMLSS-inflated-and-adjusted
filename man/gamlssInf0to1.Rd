\name{gamlssInf0to1}
\alias{gamlssInf0to1}

\title{
GAMLSS model for a proportion response variable with point(s) mass at 0 and or 1.
}
\description{
Function \code{gamlssInf0to1()} allows to fit inflated gamlss models when the response variable distribution is defined in the intervals [0,1), (0,1] and [0,1]. 
The \code{gamlssInf0to1} model for inflated proportion variables is a \code{gamlss} model provided of up to two extra parameters for the mass point(s). In the case of inflation point at zero (one), this is equivalent to fit two separate models, a gamlss model for the (0,1) part, and a logit model for zero (one) vs non-zero (non-one) part. When both zero and one are present, a multinomial model is involved to fit the non-(0,1) part.
}
\usage{
gamlssInf0to1(y = NULL, mu.formula = ~1, sigma.formula = ~1, 
             nu.formula = ~1,tau.formula = ~1, 
             xi0.formula = ~1,xi1.formula = ~1, data = NULL, 
             family = BE, weights = rep(1, length(Y_)), 
             trace = FALSE, ...)
}

\arguments{
  \item{y}{ the proportion response variable with inflation at zero and/or one}
  \item{mu.formula}{
  a model formula for \code{mu}
}
  \item{sigma.formula}{
  a model formula for  \code{sigma} 
}
  \item{nu.formula}{
  a model formula for \code{nu} 
}
  \item{tau.formula}{
a model formula for \code{tau} 
}
  \item{xi0.formula}{
a model formula for the probability at zero
}
  \item{xi1.formula}{
a model formula for the probability at one
}
  \item{data}{
a data frame containing the variables occurring in the formula.
}
  \item{family}{
any \code{gamlss} distribution family defined in (0,1)
}
  \item{weights}{a vector of weights as in \link{gamlss} 
}
  \item{trace}{
logical, if TRUE information on model estimation will be printed during the fitting
}
  \item{\dots}{
 for extra parameters
}
}
\details{
The default family is a Beta distribution (BE), but other (0,1) distributions can be used, e.g.  those generated from existing continuous gamlss family distributions by using \code{\link{gen.Family}} with link "logit". 
}
\value{
returns a \code{gamlssInf0to1} object which has its own methods  
 }
\references{
Hossain, A., Stasinopoulos, M., Rigby, R. and Enea, M. (2015). Centile estimation for a proportion response variable. \emph{Statistics in Medicine},  doi: 10.1002/sim.6748. \cr

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{
Mikis Stasinopoulos, Robert Rigby, Abu Hossain and Marco Enea
}

\seealso{
\code{\link{gamlss.family}}, \code{\link{BEINF}}, \code{\link{BE}}, \code{\link{BEo}}, \code{\link{BEZI}}, \code{\link{BEOI}}
}
\examples{

# 1. An artificial example using simulated data

# Firstly, we use function gen.Family() to create the logit skew 
# student t (logitSST) distribution defined in the (0,1) interval, 
# and function gen.Inf0to1() to create the 0-inflated logitSST 
# distribution defined in [0,1).

gen.Family("SST", "logit")
gen.Inf0to1("logitSST","Zero") 

#now we can generate the data and run the model 
set.seed(10)
Y <- rlogitSSTInf0(500,mu=0.5,sigma=0.7,nu=0.5,tau=5,xi0=0.5,log=FALSE)
dat <- data.frame(Y)
dat$x <- rnorm(500)
m1 <- gamlssInf0to1(y=Y,mu.formula=~x, sigma.formula=~x,
                    nu.formula=~x, tau.formula=~x,
                    xi0.formula=~x,data=dat, family=logitSST)
summary(m1)

# 2. Example of equivalent gamlss models for an inflated-at-1 Beta distribution 

Y <- rBEINF1(500,mu=0.5,sigma=0.7,nu=0.5)
m2 <- gamlss(Y~1,sigma.formula=~1,nu.formula=~1,family=BEINF1)
m3.1 <- gamlss(Y[Y<1]~1,sigma.formula=~1,family=BE)
m3.2 <- gamlss(I(Y==1)~1,family=BI)
m4 <- gamlssInf0to1(Y,mu.formula=~1,sigma.formula=~1,xi1=~1,family=BE)
stopifnot(all.equal(deviance(m2),(deviance(m3.1)+deviance(m3.2))), 
          all.equal(deviance(m2),deviance(m4)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}

