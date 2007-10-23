\name{monoproc}
\alias{monoproc}
\alias{monoproc,locfit,numeric,missing,missing-method}
\alias{monoproc,locfit,numeric,missing,character-method}
\alias{monoproc,locfit,numeric,numeric,missing-method}
\alias{monoproc,locfit,numeric,list,character-method}
\alias{monoproc,fit,missing,ANY,ANY-method}
\alias{monoproc,loess,numeric,missing,missing-method}
\alias{monoproc,loess,numeric,numeric,missing-method}
\alias{monoproc,locpoly,numeric,missing,missing-method}
\alias{monoproc,locpoly,numeric,numeric,missing-method}
\alias{monoproc,ksmooth,numeric,missing,missing-method}
\alias{monoproc,ksmooth,numeric,numeric,missing-method}
\alias{monoproc,list,numeric,missing,missing-method}
\alias{monoproc,list,numeric,missing,character-method}
\alias{monoproc,list,numeric,numeric,missing-method}
\alias{monoproc,list,numeric,list,character-method}
\title{ Monotonizing procedure }
\description{
 this procedure monotonizes a given fit in either one variable or two
}
\usage{
monoproc(fit, bandwidth, xx, dir,...)
}
\arguments{
  \item{fit}{so far an object either of class \code{"list"}, \code{"loess"}, \code{"ksmooth"},\code{"locpoly"}, or \code{"locfit"} }
  \item{bandwidth}{ a single smoothing paramter of class \code{"numeric"}. Missing values are not accepted.}
  \item{xx}{ points where the monotone fit is to be evaluated. If missing the independent variables of \code{fit} are used instead. }
\item{dir}{gives the order, direction, and the variable for twodimensional problems. This variable has to be missing for onedimensional problems otherwise there is an error. Possible values for \code{dir} are \code{"x"}, \code{"y"}, \code{"xy"}, and \code{"yx"}. }
\item{\dots}{further parameters about the monotonization procedure. See also in Details.}
}
\details{ If fit is of class \code{"list"}, the independent variable should be equidistant. The length of a list object can be either 2 (one independent variable as first element) or 3 (two independent variables as first two elements).

In twodimensional problems, the value \code{"x"} for \code{dir} refers to a single monotonization with respect to the first variable and \code{"y"} correspondingly to the second variable. With the values \code{"xy"} and \code{"yx"}, the order of monotonization can be determined. The difference of the resulting monotone fits may be marginal, but in some cases of interest.

The kernel of the monotonizing procedure can be specified as well. Currently, only the Epanechnikov kernel is implemented which is the default value of \code{kernel}.

With the arguments \code{mono1} and \code{mono2}, the kind of monotonization can be specified - either \code{"increasing"} or \code{"decreasing"}. \code{mono1} refers to the first independent variable and correspondingly \code{mono2} to the second one. For an onedimensional problem only \code{mono1} can be specified. The default value for \code{mono1} and \code{mono2} is \code{"increasing"}. 
 
The gridsize can be further spezified, if the fit is not of class \code{"list"}. However, do not use too large gridsizes. A gridsize larger than 50 is for twodimensional problems not recommended.
}
\value{
An object of class \code{\link{monoproc.1d-class}} or  \code{\link{monoproc.2d-class}} according to the number of independent variables
}

\references{Dette, H., Neumeyer, N., and Pilz, K. (2004) \emph{A simple nonparametric estimator of a monotone regression function}.

Dette, H. and Scheder, R. (2005) \emph{Striclty monotone and smooth nonparametric regression for two or more variable}.}

\author{ Regine Scheder \email{Regine.Scheder@rub.de}}

\seealso{ \code{\link[=monoproc.1d-class]{monoproc.1d}} and \code{\link[=monoproc.2d-class]{monoproc.2d}} }
\examples{
#Spencer's Mortality Dataset
if(require(locfit)){
data(spencer)
attach(spencer)
fit<-locfit.raw(age,mortality, alpha=0.3, kern="epan")
fitmono<-monoproc(fit, bandwidth=0.0003, gridsize=30)
plot(age,mortality)
lines(fit)
lines(fitmono, col=2)}

#Fat Data to predict Bodyfat (two independent variables)
if(require(UsingR)&&require(locfit)){
data(fat)
fat<-fat[-39,] ##two extreme observations
fat<-fat[-41,] ##are deleted
attach(fat)
fit<-locfit.raw(cbind(weight, height),body.fat.siri, alpha=0.3, deg=1, kern="epan")
fitmono<-monoproc(fit,bandwidth=1, mono1= "increasing", mono2="decreasing", dir="xy", gridsize=30)
nf<- layout(matrix(c(1,1,1,2,2,3,3,3,4,4), 2, 5, byrow = TRUE))
layout.show(nf)
plot(fit, type="persp", theta = 135, phi = 30,col="lightblue",
 	cex=0.7,main="unconstraint Bodyfat estimate")
plot(fit)
plot(fitmono,theta = 135, phi = 30, col="lightblue",cex=0.7, 
	main="monotone Bodyfat estimate")
plot(fitmono, type="contour")
t<-cv(fitmono)
CV<-sum((t[,2]-body.fat.siri)^2)/250	#Cross Validation for the unconstraint estimator
CV2<-sum((t[,1]-body.fat.siri)^2)/250	#Cross Validation for the monotone estimator
}


#Two Examples about the cars Data
#first example somehow trivial since the loess-fit is already monotone increasing
cars.lo <- loess(dist ~ speed, cars,
control = loess.control(surface = "direct"))
predict<-predict(cars.lo, data.frame(speed = seq(5, 30, 1)))
plot(cars.lo, xlab="speed", ylab="dist")
lines(seq(5, 30, 1),predict)
monofit<-monoproc(cars.lo, bandwidth=0.3, gridsize=40)
lines(monofit,col=2)


data(cars)
speed<-cars$speed
dist<-cars$dist
fit1<-ksmooth(speed, dist, "normal", bandwidth=2)
##computes the Nadaraya-Watson estimate
fit2<-monoproc(fit1,bandwidth=0.7)
##calculates the monotone estimates
         plot(speed, dist)
         lines(fit1, col=2)
         lines(fit2, col=3)

#Comparison with R-function isoreg
fit1<-ksmooth(speed, dist, "normal", bandwidth=2.5)
fit2<-monoproc(fit1,bandwidth=0.7)
fit3<-isoreg(speed,dist)

plot(fit3,plot.type="single", main="monotone regression", xlab="speed", ylab="distance")
  lines(fit1, col=3, lwd=1.5)
lines(fit2, lwd=1.5)
legend(5,100, c("isoreg","ksmooth","monoproc"),col=c(2,3,1), lty=c(1,1,1))

}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}