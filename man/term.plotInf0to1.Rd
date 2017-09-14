\name{term.plotInf0to1}
\alias{term.plotInf0to1}

\title{Plot regression terms for a specified parameter of a fitted gamlssInf0to1 object }
\description{
 \code{term.plotInf0to1} produces term plots for a specified parameter from a    
  gamlssinf0to1 object. 
}
\usage{
term.plotInf0to1(object, parameter = c("mu", "sigma", "nu", "tau", "xi0", "xi1"),...) 
                                   
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a gamlssinf0to1 fitted model}
  \item{parameter}{ which distribution (or inflation) parameter is required, default \code{parameter="mu"}}
 \item{\dots}{ extra arguments, the same of \code{\link[gamlss]{term.plot}} (except 'what')  }
}
\details{ see function \code{\link[gamlss]{term.plot}}}
 
\value{
  A plot of fitted terms.
}

\references{
Hossain, A., Stasinopoulos, M., Rigby, R. and Enea, M. (2015). Centile estimation for a proportion response variable. \emph{Statistics in Medicine},  doi: 10.1002/sim.6748. \cr

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
 }

\author{Abu Hossain, Mikis Stasinopoulos \email{mikis.stasinopoulos@gamlss.org}, Bob Rigby and Marco Enea }

\seealso{  \code{\link{gamlssInf0to1}}}
\examples{

gen.Family("SST", "logit")
gen.Inf0to1("logitSST","Zero") 
set.seed(10)
Y <- rlogitSSTInf0(500,mu=0.5,sigma=0.7,nu=0.5,tau=5,xi0=0.5,log=FALSE)
dat <- data.frame(Y)
dat$x <- rnorm(500)
m1 <- gamlssInf0to1(y=Y,mu.formula=~x, sigma.formula=~x,
                    nu.formula=~x, tau.formula=~x,
                    xi0.formula=~x,data=dat, family=logitSST)

# term.plot for the mu parameter
term.plot(m1$dist,parameter="mu")
term.plotInf0to1(m1,parameter="mu",col.shaded = 2)

# term.plot for the binomial parameter
term.plot(m1$multinom,parameter="mu")
term.plotInf0to1(m1,parameter="xi0",col.shaded = 2)

}
\keyword{regression}
