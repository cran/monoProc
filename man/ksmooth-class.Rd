\name{ksmooth-class}
\alias{ksmooth-class}
\alias{lines,ksmooth-method}
\alias{plot,ksmooth,missing-method}

\title{Class "ksmooth"  }
\description{S4-Class for the value of the \code{ksmooth} function}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ksmooth", ...)}. This class is the value of the function \code{ksmooth} transfered to a S4-Class.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"numeric"}. x-values of the smoothed fit in increasing order.}
    \item{\code{y}:}{Object of class \code{"numeric"}. The corresponding fitted y-values of the Nadaraya-Watson estimator.}
    \item{\code{call}:}{Object of class \code{"call"}. The call of the function \code{ksmooth} for this object.}
  }
}
\section{Methods}{
  \describe{
    \item{lines}{\code{signature(x = "ksmooth")}: line method to plot this class.}
    \item{plot}{\code{signature(x = "ksmooth", y = "missing")}: plot method for this class.}
  }
}



\seealso{
  \code{\link{ksmooth}}
}

\keyword{classes}
