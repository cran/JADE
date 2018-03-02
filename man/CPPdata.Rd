\name{CPPdata}
\alias{CPPdata}
\docType{data}
\title{
Cocktail Party Problem Data
}
\description{
This data set is a toy example for the so called cocktail party problem. In this case three sounds are mixed together with one noise source using four microphones.
}
\usage{data("CPPdata")}
\format{
  A data frame with 50000 observations on the following 4 variables.
  \describe{
    \item{\code{Mic1}}{the mixture recorded by the first microphone.}
    \item{\code{Mic2}}{the mixture recorded by the second microphone.}
    \item{\code{Mic3}}{the mixture recorded by the third microphone.}
    \item{\code{Mic4}}{the mixture recorded by the fourth microphone.}
  }
}
\details{
The three original source files were kindly provided by Ella Bingham and are also available online at the following locations:
\url{research.ics.aalto.fi/ica/cocktail/source5.wav}, \url{research.ics.aalto.fi/ica/cocktail/source7.wav} and \url{research.ics.aalto.fi/ica/cocktail/source9.wav}.

Note that the original sound files are included in the package's subfolder datafiles. In the example section we illustrate how the CPPdata was created. 
An example analysis of the data is given in Miettinen et al. (2017).
}
\source{
\url{http://research.ics.aalto.fi/ica/cocktail/cocktail_en.cgi}
}
\references{
\cite{Miettinen, J., Nordhausen, K. and Taskinen, S. (2017), Blind Source Separation Based on Joint Diagonalization in R: The Packages JADE and BSSasymp, \emph{Journal of Statistical Software}, \bold{76}, 1--31, <doi:10.18637/jss.v076.i02>.}
}
\examples{
\dontrun{
# the data was created as follows:
library("tuneR")
S1 <- readWave(system.file("datafiles/source5.wav", package = "JADE"))
S2 <- readWave(system.file("datafiles/source7.wav", package = "JADE"))
S3 <- readWave(system.file("datafiles/source9.wav", package = "JADE"))

set.seed(321)
NOISE <- noise("white", duration = 50000)
S <- cbind(S1@left, S2@left, S3@left, NOISE@left)
S <- scale(S, center = FALSE, scale = apply(S, 2, sd))
St <- ts(S, start = 0, frequency = 8000)
p <- 4
A <- matrix(runif(p^2, 0, 1), p, p)
A

X <- tcrossprod(St, A)
Xt <- as.ts(X)

colnames(X) <- c("Mic1", "Mic2", "Mic3", "Mic4")
CPPdata <- as.data.frame(X)
}
}
\keyword{datasets}
