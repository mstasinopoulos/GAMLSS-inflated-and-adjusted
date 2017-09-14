\name{term.plotZadj}
\alias{term.plotZadj}

\title{Plot regression terms for a specified parameter of a fitted gamlssZadj object }
\description{
 \code{term.plotZadj} produces term plots for a specified parameter from a    
  gamlssZadjobject. 
}
\usage{
term.plotZadj(object, parameter = c("mu", "sigma", "nu", "tau", "xi0"),...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a gamlssZadj fitted model}
  \item{parameter}{ which distribution (or inflation) parameter is required, default \code{parameter="mu"}}
 \item{\dots}{ extra arguments, the same of \link[gamlss]{term.plot} (except 'what')  }
}
\details{ see function \link[gamlss]{term.plot}}
 
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

\seealso{  \code{\link{gamlssZadj}}}
\examples{

set.seed(3210)
x <- (runif(1000)*4)-2
data(sda)
fmu <- splinefun(sda$x, sda$mu)
curve(fmu, -2,2)
fsigma <- splinefun(sda$x, sda$sigma)
curve(fsigma, -2,2)
fnu <- function(x)
  {f <- splinefun(sda$x, sda$nu)
f(x)/6
}
curve(fnu, -2,2)
set.seed(321)
y0 <- rZAGA(1000, mu=fmu(x), sigma=fsigma(x), nu=fnu(x))
da <- data.frame(y0,x)
g0p <- gamlss(y0~pb(x), sigma.fo=~pb(x), nu.fo=~pb(x), data=da, family=ZAGA)
t0p <- gamlssZadj(y=y0, mu.fo=~pb(x), sigma.fo=~pb(x),data=da,
                  trace=TRUE, xi0.fo=~pb(x), family="GA")

# term.plot for the mu parameter
term.plot(g0p);title("gamlss")
term.plot(t0p$dist,"mu");title("gamlssZadj")
term.plotZadj(t0p,"mu",col.shaded = 3);title("gamlssZadj")



# term.plot for the sigma parameter
term.plot(g0p, "sigma");title("gamlss")
term.plot(t0p$dist,"sigma");title("gamlssZadj")
term.plotZadj(t0p,"sigma",col.shaded = 3);title("gamlssZadj")

# term.plot for the binomial parameter
term.plot(g0p, "nu");title("gamlss")
term.plot(t0p$binom,"mu");title("gamlssZadj")
term.plotZadj(t0p,"xi0",col.shaded = 3);title("gamlssZadj")


}
\keyword{regression}
