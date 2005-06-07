\name{monoproc.2d-class}
\docType{class}
\alias{monoproc.2d-class}
\alias{monoproclocfit.2d-class}
\alias{plot,monoproc.2d,missing-method}
\alias{print,monoproc.2d-method}
\alias{summary,monoproc.2d-method}

\title{Class "monoproc.2d" }
\description{A representation of a monotonized fit with information about the original fit and more}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("monoproc.2d", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"monofit"} . This is the new monotonized fit .}
    \item{\code{fitold}:}{Object of class \code{"fit"} . This slot represents the old fit and can be of class \code{"list"}, \code{"locfit"}, \code{"ksmooth" }, \code{"locpoly"}, or \code{"loess"}. }
    \item{\code{gridsize}:}{Object of class \code{"numeric"} .}
    \item{\code{bandwidth}:}{Object of class \code{"numeric"}. The used bandwidth for the monotonizing procedure.}
    \item{\code{kernel}:}{Object of class \code{"character"}. The used kernel. Currently only \code{"epanech"}.}
    \item{\code{mono}:}{Object of class \code{"character"}.  }
\item{\code{dir}:}{Directions and order of the monotonization.  Possible values are "xy", "yx", "x", or "y".} 
    \item{\code{name}:}{Object of class \code{"character"}. Names of the variables. }
\item{\code{call}:}{Call of this object.}
  }
}

\details{the class \code{"monoproclocfit.2d"} is a \code{"monoproc.2d"}-object where the original fit is of class \code{"locfit"}. Both \code{"monoproc.2d"} as well as \code{"monoproclocfit.2d"} belong to the class \code{"monoproc"}
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "monoproc.2d", y = "missing")}:plots the monotone fit. Three types are available: \code{"persp"}, \code{"contour"}, \code{"image"}. Default is \code{"persp"}.}
    \item{print}{\code{signature(x = "monoproc.2d")}: prints the monotone fit.}
    \item{summary}{\code{signature(object = "monoproc.2d")}: gives a summary of the monotone fit.}
\item{cv}{\code{signature(fit= "monoproclocfit.2d")}: calculates the cross validation value for the monotone and the unconstraint fit.}
  }
}

\seealso{ \code{\link{monoproc.1d-class}} 
}


\keyword{classes}
