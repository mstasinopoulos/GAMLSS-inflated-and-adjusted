\name{summary.gamlssinf0to1}
\alias{summary.gamlssinf0to1}
\alias{summary.gamlssZadj}

\title{Summarizes an inflated GAMLSS fitted model}
\description{These are specific methods for the generic function \code{summary} which summarize objects retuned by \code{gamlssinf0to1} or \code{gamlssZadj}. 
}
\usage{
\method{summary}{gamlssinf0to1}(object, type = c("vcov", "qr"),
       robust=FALSE, save = FALSE, hessian.fun = c("R", "PB"), 
       digits = max(3, getOption("digits") - 3),...)

\method{summary}{gamlssZadj}(object, type = c("vcov", "qr"),
       robust=FALSE, save = FALSE, hessian.fun = c("R", "PB"), 
       digits = max(3, getOption("digits") - 3),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a gamlssinf0to1 or gamlssZadj fitted model}
  \item{type}{the default value \code{vcov} uses the \code{vcov()} method for gamlss to get the 
  variance-covariance  matrix of the estimated beta coefficients, see details below. 
  The alternative \code{qr} is the original method used in gamlss to 
  estimated the standard errors but it is not reliable since it do not take into the account the inter-correlation between
  the distributional parameters \code{mu}, \code{sigma}, \code{nu} and \code{tau}, while the inflation  parameters \code{xi0} and \code{xi1} are uncorrelated anywway. 
 }
  \item{robust}{whether robust (sandwich) standard errors are required}
  \item{save}{whether to save the environment of the function so to have access to its values}
  \item{hessian.fun}{whether when calculate the Hessian should use the "R" function \code{optimHess()} or a function based on Pinheiro and Bates \code{nlme} package, "PB".}
 \item{digits}{the number of digits in the output}
\item{\dots}{for extra arguments}
}
\details{Using the  default value \code{type="vcov"}, the \code{vcov()} method is used to get  
  the variance covariance matrix (and consequently the standard errors) of the beta parameters. 
  The variance covariance matrix is  calculated using the inverse of the numerical second derivatives
  of the observed information matrix. This is a more reliable method since it take into the account the 
  inter-correlation between the all the parameters. The \code{type="qr"} assumes that the parameters are fixed 
  at the estimated values. Note that both methods are not appropriate and should be used with caution if smoothing 
  terms are used in the fitting. }
 
\value{
 Print summary of a gamlssinf0to1 or a gamlssZadj object
}

\references{
Houssain, A., Stasinopoulos, M., Rigby, R. and Enea, M. (2015). Centile estimation for a proportion response variable. Accepted for publication on \emph{Statistics in Medicine}.  \cr

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2003) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
 }

\author{Abu Hossain, Mikis Stasinopoulos \email{mikis.stasinopoulos@gamlss.org}, Bob Rigby and Marco Enea }

\seealso{  \code{\link{gamlssInf0to1}}, \code{\link{gamlssZadj}}}
\examples{

# 1. An artificial example using simulated data for the gamlssInf0to1 example
# Firstly, we use function gen.Family() to create the logit skew student t 
# (logitSST) distribution defined in the (0,1) interval, and function 
# gen.Inf0to1() to create the 0-inflated logitSST distribution defined 
# in [0,1).

gen.Family("SST", "logit")
gen.Inf0to1("logitSST","Zero") 

#now we can generate the data and run the model 
set.seed(100)
Y <- rlogitSSTInf0(500,mu=0.5,sigma=0.7,nu=0.5,tau=5,xi0=0.5,log=FALSE)
dat <- data.frame(Y)
dat$x <- rnorm(500)
m1 <- gamlssInf0to1(y=Y,mu.formula=~x, sigma.formula=~1,
                    nu.formula=~1, tau.formula=~1,
                    xi0.formula=~pb(x),data=dat, family=logitSST)
summary(m1)

# 2. The gamlssZadj example 
set.seed(3210)
x <- (runif(1000)*4)-2
data(sda)
fmu <- splinefun(sda$x, sda$mu)
fsigma <- splinefun(sda$x, sda$sigma)
fnu <- function(x)
  {f <- splinefun(sda$x, sda$nu)
f(x)/6
}
y0 <- rZAGA(1000, mu=fmu(x), sigma=fsigma(x), nu=fnu(x))
da <- data.frame(y0,x)
t0p <- gamlssZadj(y=y0, mu.fo=~pb(x), sigma.fo=~pb(x),data=da,
                  trace=TRUE, xi0.fo=~pb(x), family="GA")
summary(t0p)
}
\keyword{regression}
