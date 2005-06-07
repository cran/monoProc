\name{mono.2d}
\alias{mono.2d}
\title{strictly monotone and smooth function}
\description{
  this function applies a kernel smoothing method to monotonize a given fit with two independent variable
}
\usage{
mono.2d(fit, bandwidth, xx, kernel = "epanech", dir = "xy", mono1, mono2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{ a list containing x-values and y-values as independent variables and their corresponding z-values. The length of \code{fit[[1]]} determines degree of acuteness of the monotonizing procedure }
  \item{bandwidth}{ a single number which represents the kernel bandwidth smoothing parameter. Missing values are not accepted.}
  \item{xx}{an additional list of x-values and y-values where the monotonizing procedure is to be evaluated. If missing \code{list(x=fit[[1]],y=fit[[2]])} is used instead. }
  \item{kernel}{  "`epanech"' - the Epanechnikov kernel }
  \item{dir}{ with respect to which variable and in which order the function is to be monotonized. Possible is "x", "y", "xy", and "yx" }
  \item{mono1}{either "increasing" or "decreasing"  for x-variable}
  \item{mono2}{ either "increasing" or "decreasing" for y-variable}
}

\details{
 this function is used within "monoproc"
}
\value{
  returns an object of class "monofit"}
\references{Dette, H., Neumeyer, N., and Pilz, K. (2004) \emph{A simple nonparametric estimator of a monotone regression function}.

Dette, H. and Scheder, R. (2005) \emph{Striclty monotone and smooth nonparametric regression for two or more variable}.}

\author{ Regine Scheder \email{Regine.Scheder@rub.de}}


\seealso{ \code{\link{mono.1d}} and \code{\link{monoproc}} }

\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}