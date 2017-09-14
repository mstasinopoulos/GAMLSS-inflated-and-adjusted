\name{gen.Zadj}
\alias{gen.Zadj}
\alias{Zadj.d}
\alias{Zadj.p}
\alias{Zadj.q}
\alias{Zadj.r}
\alias{plotZadj}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to generate zero adjusted distributions from existing continuous gamlss.family distributions defined on positve real line.}
\description{
There are six functions here. Only the function \code{gen.Zadj()} should be used. The remaing four functions will be automatically created once \code{gen.Zadj()} has been run.
}
\usage{
gen.Zadj(family = "GA", ...)
Zadj.d(family = "GA", ...)
Zadj.p(family = "GA", ...)
Zadj.q(family = "GA", ...)
Zadj.r(family = "GA", ...)
plotZadj(family = "GA", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{family}{a continuous positive rael line  distribution}
  \item{\dots}{for additional arguments}
}
\details{
Functions \code{Zadj.d}, \code{Zadj.p}, \code{Zadj.q} and \code{Zadj.r} allow to create the density function, distribution function, quantile function and random generation, respectively. Function \code{plotZadj} can be used to create a  plot for the distribution.

Alternatively, the function \code{gen.Zadj} creates the all the standard \code{d,p,q,r} functions plus the plotting function.}
\value{
 The function \code{gen.Zadj} returns  the \code{d},  \code{p}, \code{q} and \code{r} functions plus the plotting function. 
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gamlssZadj}}, ~~~
}
\examples{
# 1.
gen.Zadj("BCT")

plotBCTZadj()
plotBCTZadj(mu=3,sigma=0.35,xi0=0.5)

# 2.
gen.Family("SST", "log")
gen.Zadj("logSST") 
 plotlogSSTZadj()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
\keyword{regression}
