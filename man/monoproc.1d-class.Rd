\name{monoproc.1d-class}
\alias{monoproc.1d-class}
\alias{monoproclocfit.1d-class}
\alias{lines,monoproc.1d-method}
\alias{plot,monoproc.1d,missing-method}
\alias{print,monoproc.1d-method}
\alias{show,monoproc.1d-method}
\alias{summary,monoproc.1d-method}

\title{Class "monoproc.1d"}
\description{A representation of a monotonized fit with information about the original fit and more}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("monoproc.1d", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"monofit"} where fit@z is null. This is the new monotonized fit.}
    \item{\code{fitold}:}{Object of class \code{"fit"}. This slot represents the old fit and can be of class \code{"list"}, \code{"locfit"}, \code{"ksmooth"},\code{"locpoly"}, or \code{"loess"}. }
    \item{\code{gridsize}:}{Object of class \code{"numeric"} }
    \item{\code{bandwidth}:}{Object of class \code{"numeric"}. The used bandwidth for the monotonizing procedure.}
    \item{\code{kernel}:}{Object of class \code{"character"}. The used kernel. Currently only \code{"epanech"}.}
    \item{\code{mono}:}{Object of class \code{"character"}. Either \code{"increasing"} or \code{"decreasing"} }
    \item{\code{name}:}{Object of class \code{"character"}. Names of the variables.}
\item{\code{call}:}{Call of this object.}
  }}

\details{the class \code{"monoproclocfit.1d"} is a \code{"monoproc.1d"}-object where the original fit is of class \code{"locfit"}. Both \code{"monoproc.1d"} as well as \code{"monoproclocfit.1d"} belong to the class \code{"monoproc"}}

\section{Methods}{
  \describe{
    \item{lines}{\code{signature(x = "monoproc.1d")}: adds the line of the monotone fit to a given plot. }
    \item{plot}{\code{signature(x = "monoproc.1d", y = "missing")}: plots a monotone fit.}
    \item{print}{\code{signature(x = "monoproc.1d")}: prints a monotone fit.}
    \item{summary}{\code{signature(object = "monoproc.1d")}: gives a summary of the monotone fit.}
\item{cv}{\code{signature(fit = "monoproclocfit.1d")}: calculates the cross validation value of the monotone and the unconstraint fit.}
  }
}

\seealso{ \code{\link[=monoproc.2d-class]{monoproc.2d}} 
}



\keyword{classes}
