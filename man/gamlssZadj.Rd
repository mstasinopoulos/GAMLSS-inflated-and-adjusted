\name{gamlssZadj}
\alias{gamlssZadj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting  positive real line response variable with zeros.
}
\description{
Function \code{gamlssZadj()} allows to fit zero adjusted   gamlss models when the response variable distribution is defined on the  positive real line. 
The \code{gamlssZadj} model for adjusted  positive variables is a \code{gamlss} model provides one extra parameters for the mass point at zero. This is equivalent to fit two separate models, a gamlss model for the (0,Inf) part, and a logit model for zero part versus the non-zero  part. The function works similarly but provides one fitted object.
}
\usage{
gamlssZadj(y = NULL, mu.formula = ~1, sigma.formula = ~1, 
         nu.formula = ~1, tau.formula = ~1, 
         xi0.formula = ~1, data = NULL, 
         family = GA, 
         weights = rep(1, length(Y_)), trace = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the response variable}
  \item{mu.formula}{a model formula for \code{mu}}
  \item{sigma.formula}{ a model formula for \code{sigma}}
  \item{nu.formula}{ a model formula for \code{nu}}
  \item{tau.formula}{a model formula for \code{tau}}
  \item{xi0.formula}{a model formula for \code{xi0}} 
  \item{data}{a data frame containing the variables occurring in the formula.}
  \item{family}{any \code{gamlss} distribution family defined on the rael line}
  \item{weights}{a vector of weights as in \link{gamlss} }
  \item{trace}{logical, if TRUE information on model estimation will be printed during the fitting}
  \item{\dots}{for extra arguments to pass to \code{gamlss}}
}  
\details{
The default family is a gamma distribution (GA), but other  distributions on the positive rael line can be used, e.g.  those generated from existing continuous \code{gamlss.family} distributions  using say \code{\link{gen.Family}()} with "log"
or \code{gen.trun()} from package \pkg{gamlss.tr}}
\value{.
Returns a \code{gamlssZadj} object which has its own methods 
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

\author{Mikis Stasinopoulos, Robert Rigby and Marco Enea}


\seealso{
\code{\link{gamlss.family}}, \code{\link{ZAGA}}, \code{\link{ZAIG}}}

\examples{
  y0 <- rZAGA(1000, mu=.3, sigma=.4, nu=.15)# p0=0.13
  g0 <- gamlss(y0~1, family=ZAGA)
 t0 <- gamlssZadj(y=y0, mu.formula=~1, family=GA, trace=TRUE)
AIC(g0,t0, k=0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}

