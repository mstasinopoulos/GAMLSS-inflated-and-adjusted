\name{centiles.Inf0to1}
\alias{centiles.Inf0to1}
\alias{centiles.Zadj}

\title{Plotting centile curves for a \code{gamlssInf0to1} and  \code{gamlssZadj} object }
\description{Function \code{centiles.Inf0to1()} plots centile curves for distributions belonging to the GAMLSS family of distributions defined in the intervals (0,1],[0,1) and [0,1].  
The function also tabulates the sample percentages below each centile curve (for comparison with the model percentages given by the argument \code{cent}).
A restriction of the function is that it applies to models with one explanatory variable only.    
}
\usage{
centiles.Inf0to1(obj, xvar = NULL, cent = c(0.4, 2, 10, 25, 50, 75, 
    90, 98, 99.6), legend = TRUE, ylab = "y", xlab = "x", main = NULL, 
    main.gsub = "@", xleg = min(xvar), yleg = max(obj$y), xlim = range(xvar), 
    ylim = range(obj$y), save = FALSE, plot = TRUE, points = TRUE, 
    pch = 15, cex = 0.5, col = gray(0.7), col.centiles = 1:length(cent) + 
        2, lty.centiles = 1, lwd.centiles = 1, ...)
centiles.Zadj(obj, xvar = NULL, cent = c(0.4, 2, 10, 25, 50, 75, 
   90, 98, 99.6), legend = TRUE, ylab = "y", xlab = "x", main = NULL, 
   main.gsub = "@", xleg = min(xvar), yleg = max(obj$y), xlim = range(xvar), 
   ylim = range(obj$y), save = FALSE, plot = TRUE, points = TRUE, 
   pch = 15, cex = 0.5, col = gray(0.7), col.centiles = 1:length(cent) + 
   2, lty.centiles = 1, lwd.centiles = 1, ...)
}

\arguments{
  \item{obj}{a fitted gamlss object from fitting a gamlss distribution}
  \item{xvar}{the unique explanatory variable}
  \item{cent}{a vector with elements the \% centile values for which the centile curves have to be evaluated}
  \item{legend}{whether a legend is required in the plot or not, the default is \code{legent=TRUE}  }
  \item{ylab}{the y-variable label}
  \item{xlab}{the x-variable label}
  \item{main}{the main title here as character. 
              If NULL the default title "centile curves using  NO" (or the relevant distributions name) is shown}
  \item{main.gsub}{if the \code{main.gsub} (with default "@") appears in the \code{main} 
                 title then it is substituted with the default title.}
  \item{xleg}{position of the legend in the x-axis}
  \item{yleg}{position of the legend in the y-axis}
  \item{xlim}{the limits of the x-axis}
  \item{ylim}{the limits of the y-axis}
  \item{save}{whether to save the sample percentages or not with default equal to \code{FALSE}. 
              In this case the sample percentages are printed but are not saved }
  \item{plot}{whether to plot the centiles}
  \item{points}{whether the data points should be plotted, default is \code{TRUE}} 
  \item{pch}{the character to be used as the default in plotting points see \code{par}}
   \item{cex}{size of  character see \code{par}}
  \item{col}{plotting colour see \code{par}}
  \item{col.centiles}{Plotting colours  for the centile curves}
  \item{lty.centiles}{line type for the centile curves}
  \item{lwd.centiles}{The line width for the centile curves}
  \item{\dots}{for extra arguments}          
}
\details{
 Centiles are calculated using the fitted values in \code{obj} and \code{xvar} must
correspond exactly to the predictor in \code{obj} to plot correctly.

\code{col.centiles}, \code{lty.centiles} and \code{lwd.centiles} may be vector arguments
and are recycled to the length \code{cent} if necessary.
}
\value{
A centile plot is produced and the sample centiles below each centile curve are printed (or saved)
}
\references{Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Houssain, A., Stasinopoulos, M., Rigby, R. and Enea, M. (2015). Centile estimation for a proportion response variable. Accepted for publication on \emph{Statistics in Medicine}.  \cr

}

\author{Mikis Stasinopoulos \email{mikis.stasinopoulos@gamlss.org}, Bob Rigby \email{r.rigby@londonmet.ac.uk} with contribution from Steve Ellison}


\section{Warning }{This function is appropriate only when one continuous explanatory variable is fitted in the model}  
\seealso{  \code{\link{gamlssInf0to1}}, \code{\link{gamlss}}, \code{\link{centiles.split}} , \code{\link{centiles.com}}}
\examples{
\dontrun{
gen.Family("SST", "logit")
gen.Inf0to1("logitSST","One") 
set.seed(30)
x <- seq(0,1,l=2000)
dat <- data.frame(x)
dat$Y <- rlogitSSTInf1(2000,mu=-3+10*x-0.7*x^2,sigma=0.9,nu=0.5,
                       tau=5,xi1=plogis(-0.5*ifelse(x>0.7,-1,13)),log=FALSE)

m1 <- gamlssInf0to1(y=Y,mu.formula=~pb(x), sigma.formula=~pb(x),
                    nu.formula=~pb(x), tau.formula=~pb(x),
                    xi1.formula=~pb(x),
                    data=dat, family=logitSST)

centiles.Inf0to1(m1,xvar= dat$x, cent=c(2,10,25,50,75,90,98),
col.centiles=c(1,7:2),
ylab="proportion", xlab="x",legend=FALSE,main="(c) Inf. logitSST")
}
}
\keyword{regression}% 
