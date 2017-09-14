\name{predict.gamlssZadj}
\alias{predict.gamlssZadj}

\title{Extract Predictor Values and Standard Errors For New Data in a gamlssZadj Model }
\description{
 \code{predict.gamlssZadj} is the specific method which produce predictors for a new data set   
for a specified parameter from a gamlssZadj objects.  
The \code{predict.gamlssZadj} can be used to extract the linear predictors, fitted values and specific  terms in the model at new 
data values in the same way that the \code{predict.lm()} and \code{predict.glm()} functions can be used for 
\code{lm} or \code{glm} objects. Note that linear predictors, fitted values and specific terms in the model at the current  
data values can also be extracted using the function \code{lpred()} (which is called from predict if new data is NULL).
}\usage{
\method{predict}{gamlssZadj}(object, parameter = c("mu", "sigma", "nu", "tau", "xi0"), 
 newdata = NULL, type = c("link", "response", "terms"),                           
   terms = NULL, se.fit = FALSE,data = NULL, ...)                                         
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a gamlssZadj fitted model}
  \item{parameter}{ which distribution (or inflation) parameter is required, default \code{parameter="mu"}}
  \item{newdata}{a data frame containing new values for the explanatory variables used in the model}
  \item{type}{the default, gets the linear predictor for the specified distribution (or inflation) parameter. 
             \code{type="response"} gets the fitted values for the parameter while \code{type="terms"} gets the fitted terms contribution}
  \item{terms}{if \code{type="terms"}, which terms to be selected (default is all terms)}
   \item{se.fit}{if TRUE the approximate standard errors of the appropriate type are extracted if exist}
  \item{data}{ the data frame used in the original fit if is not defined in the call}
  \item{\dots}{ for extra arguments }
}
\details{The predict function  assumes that the object given in \code{newdata} is a data frame containing the right x-variables
used in the model. This could possible cause problems if transformed variables are used in the fitting of the original model.
For example, let us assume that a transformation of age is needed in the model i.e. \code{nage<-age^.5}. This could be fitted as 
\code{mod<-gamlss(y~cs(age^.5),data=mydata)} or as  \code{nage<-age^.5;  mod<-gamlss(y~cs(nage), data=mydata)}.
The later could more efficient if the data are in thousands rather in hundreds. In the first case,   
the code \code{predict(mod,newdata=data.frame(age=c(34,56)))} would produce the right results.  
In the second case a new data frame has to be created containing the old data plus any new transform data. This data frame has to 
be declared in the \code{data} option. The option \code{newdata} should 
contain a data.frame with the new names and the transformed values in which prediction is required, (see the last example). }
 
\value{
  A vector or a matrix  depending on the options.
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

plot(predict(g0p,"nu",type="link"),
     predict(t0p$binom,"mu",type="link")) 

plot(predict(g0p,"nu",type="link"),
     predict(t0p,"xi0",type="link")) 
}
\keyword{regression}
