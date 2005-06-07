\name{monofit-class}
\docType{class}
\alias{monofit-class}
\alias{lines,monofit-method}
\alias{plot,monofit,missing-method}

\title{Class "monofit" }
\description{A object of class \code{monofit} is the ouput of either \code{\link{mono.1d}} or \code{\link{mono.2d}}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("monofit", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{x-variable is of class \code{"numeric"}}
    \item{\code{y}:}{y-variable is of class \code{"numeric"}}
    \item{\code{z}:}{z-variable is of class \code{"matrix"}}
  }
}
\section{Methods}{
  \describe{
    \item{lines}{\code{signature(x = "monofit")}: ... }
    \item{plot}{\code{signature(x = "monofit", y = "missing")}: ... }
  }
}

\details{This Class is usually used within the function \code{monoproc}. It is the value of the functions \code{mono.1d} and \code{mono.2d}, respectively, which are called in the \code{monoproc} function. For one independent variable in the regression function, the matrix in slot \code{z} is a 0 x 0 matrix and slot \code{y} corresponds to the values of the regression function.  
}

\seealso{
  \code{\link{mono.1d}}, \code{\link{mono.2d}} ,\code{\link{monoproc.1d-class}} , or \code{\link{monoproc.2d-class}}
}

\keyword{classes}
