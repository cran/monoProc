\name{monoproc-class}
\alias{monoproc-class}

\title{Class "monoproc" }
\description{A representation of a monotonized fit with information about the original fit and more}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("monoproc", ...)}.
   
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"monofit"} }
    \item{\code{fitold}:}{Object of class \code{"fit"}}
    \item{\code{gridsize}:}{Object of class \code{"numeric"} }
    \item{\code{bandwidth}:}{Object of class \code{"numeric"}  }
    \item{\code{kernel}:}{Object of class \code{"character"} }
    \item{\code{mono}:}{Object of class \code{"character"}  }
    \item{\code{name}:}{Object of class \code{"character"}  }
    \item{\code{call}:}{Object of class \code{"call"}  }
  }
}


\details{ For further details please see also \code{\link[=monoproc.1d-class]{monoproc.1d}} and \code{\link[=monoproc.2d-class]{monoproc.2d} }}


\keyword{classes}
