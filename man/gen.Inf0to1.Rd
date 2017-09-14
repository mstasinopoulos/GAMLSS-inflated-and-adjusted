\name{gen.Inf0to1}
\alias{gen.Inf0to1}
\alias{Inf0to1.d}
\alias{Inf0to1.p}
\alias{Inf0to1.q}
\alias{Inf0to1.r}
\alias{plotInf0to1}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to generate inflated 0-to-1 distributions from existing continuous gamlss.family distributions defined in (0,1).
}
\description{
There are six functions here. Only the function \code{gen.Inf0to1()} should be used. The remaing five functions will be automatically created once \code{gen.Inf0to1()} has been run.
}
\usage{
gen.Inf0to1(family = "BE", type.of.Inflation = c( "Zero&One", "Zero", "One"), 
            ...)
Inf0to1.d(family = "BE", type.of.Inflation = c( "Zero&One", "Zero", "One"), 
            ...) 
Inf0to1.p(family = "BE", type.of.Inflation = c( "Zero&One", "Zero", "One"), 
            ...) 
Inf0to1.q(family = "BE", type.of.Inflation = c( "Zero&One", "Zero", "One"), 
            ...) 
Inf0to1.r(family = "BE", type.of.Inflation = c( "Zero&One", "Zero", "One"), 
            ...) 
plotInf0to1(family = "BE", type.of.Inflation = c( "Zero&One", "Zero", "One"), 
            ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{family}{a continuous (0,1) distribution (extremes not included) \code{gamlss.family} distribution}
  \item{type.of.Inflation }{the type of inflation}
  \item{\dots}{for passing extra arguments}
}
\details{
Functions \code{Inf0to1.d}, \code{Inf0to1.p}, \code{Inf0to1.q} and \code{Inf0to1.r} allow to create the density function, distribution function, quantile function and random generation, respectively. Function \code{plotInf0to1} can be used to create the plot the distributions.

Alternatively, the function \code{gen.Inf0to1} creates the all the standard \code{d,p,q,r} functions plus the plotting function. 

For example, let us take the case of the logit SST distribution with inflation at 1. First generate the \code{"logitSST"} distribution by using \code{gen.Family("SST", "logit")}, and then, by use \code{gen.Inf0to1("logitSST", "One")}.
The functins  \code{dlogitSSTInf1}, \code{plogitSSTInf1}, \code{qlogitSSTInf1}, \code{rlogitSSTInf1} and \code{plotlogitSSTInf1} will be automatically generated.
Note that \code{gen.Inf0to1} never creates a fitting function of the type "logitSSTInf1", but the existing \code{logitSST} must be specified instead as an  argument \code{family} of function  \code{gamlssInf0to1()}. 
}
\value{
 The function \code{gen.Inf0to1} returns  the \code{d},  \code{p}, \code{q} and \code{r} functions plus the plotting function.  
}
\references{
 Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC. 
\url{https://www.crcpress.com/Flexible-Regression-and-Smoothing-Using-GAMLSS-in-R/Stasinopoulos-Rigby-Heller-Voudouris-Bastiani/p/book/9781138197909}.

}
\author{
Mikis Stasinopoulos \email{mikis.stasinopoulos@gamlss.org}, Bob Rigby and Marco Enea 
}


\examples{
# 1.
gen.Inf0to1("BE","Zero&One")
rBETAInf0to1 <- Inf0to1.r("BE","Zero&One") 
all.equal(rBETAInf0to1, rBEInf0to1) 
plotBEInf0to1()
plotBEInf0to1(mu=0.3,sigma=0.35,xi0=0.5,xi1=0.3)

# 2.
gen.Family("SST", "logit")
gen.Inf0to1("logitSST","One") 
set.seed(30)

args(rlogitSSTInf1)
y <- rlogitSSTInf1(1000,mu=0.2,sigma=0.5,nu=1,tau=5,xi1=0.2)
quantile(y,c(0.1,0.25,0.5,0.75,0.9))

args(qlogitSSTInf1)
qlogitSSTInf1(p=c(0.1,0.25,0.5,0.75,0.9),mu=0.2,sigma=0.5,nu=1,tau=5,xi1=0.2)
plotlogitSSTInf1(mu=0.2,sigma=0.5,nu=1,tau=5,xi1=0.2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{regression}
