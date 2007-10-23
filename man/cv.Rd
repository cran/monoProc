\name{cv}
\alias{cv}
\alias{cv,monoproc.1d-method}
\alias{cv,monoproc.2d-method}
\alias{cv,monoproclocfit.1d-method}
\alias{cv,monoproclocfit.2d-method}
\title{Cross validation function }
\description{
computes the cross validation value (leave one out method) for the monotonized and the unconstraint fit
}
\usage{
cv(fit, ...)
}

\arguments{
  \item{fit}{an object of class \code{"monoproclocfit.1d"} or \code{"monoproclocfit.2d"}}
  \item{\dots}{ currently not in use}
}
\details{
so far, this function can only be used for objects of class \code{"monoproc"} if the originally fit came from the \code{locfit}-function. This function is currently not  really computational efficient. 
}
\value{
  returns a matrix where column 1 and 2 represents the values for the new \code{fit} and the \code{fitold}, respectively, under leave one out cross valdiation. The rows correspond to the observation number.
}

\seealso{ \code{\link{monoproc}}}
\examples{
if(require(UsingR)&&require(locfit)){
	data(fat)
	fat<-fat[-39,] ##two extreme observations
	fat<-fat[-41,] ##are deleted
	attach(fat)
	x<-as.matrix(cbind(weight, height))
	fit<-locfit.raw(x,body.fat.siri, alpha=0.3, deg=1, kern="epan")
	fitmono<-monoproc(fit,bandwidth=1,mono1="increasing", mono2="decreasing", dir="xy", gridsize=30)
	nf<- layout(matrix(c(1,1,1,2,2,3,3,3,4,4), 2, 5, byrow = TRUE))
	layout.show(nf)
	plot(fit, type="persp", theta = 135, phi = 30,col="lightblue",
 	cex=0.7,main="unconstraint Bodyfat estimate")
	plot(fit)
	plot(fitmono,theta = 135, phi = 30, col="lightblue",cex=0.7, 
	main="monotone Bodyfat estimate")
	plot(fitmono, type="contour")
	t<-cv(fitmono)
	#Cross Validation for the unconstraint estimator
	sum((t[,2]-body.fat.siri)^2)/250    
	#Cross Validation for the monotone estimator
	sum((t[,1]-body.fat.siri)^2)/250   
	plot(seq(1:250),rep(1,250), type="l",col=2, xlab="observation index", 
	ylab="Ratio of CV-unconstraint over CV-monotone")
		points((t[,2]-body.fat.siri)^2/(t[,1]-body.fat.siri)^2)
	}
}
\keyword{regression}
\keyword{nonparametric}