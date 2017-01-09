\name{FOBI}
\alias{FOBI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to perform FOBI for ICA
}
\description{
The FOBI method for independent component analysis (ICA). We assume that all components have different kurtosis values.
}
\usage{
FOBI(X, na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{a numeric matrix.}
\item{na.action}{A function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}


\value{
A list with class 'bss' containing the following components:
  \item{W}{estimated unmixing matrix.}
  \item{EV}{eigenvectors of autocovariance matrix.}
  \item{Xmu}{the original mean of the data.}
  \item{S}{estimated sources as time series objected standardized to have mean 0 and unit variances.}
}
\references{
\cite{Cardoso, J.-F. (1989), Source separation using higher order moments, in Proceedings of {IEEE} International Conference
                       on Accoustics, Speech and Signal Processing, 2109--2112.}

\cite{Miettinen, J., Taskinen S., Nordhausen, K. and Oja, H. (2015), Fourth Moments and Independent
Component Analysis, \emph{Statistical Science}, 30, 372--390.}

\cite{Miettinen, J., Nordhausen, K. and Taskinen, S. (2017), Blind Source Separation Based on Joint Diagonalization in R: The Packages JADE and BSSasymp, \emph{Journal of Statistical Software}, \bold{76}, 1--31, <doi:10.18637/jss.v076.i02>.}
}
\author{Klaus Nordhausen}

\note{
More general is the function \code{\link[ICS]{ics}} in the \pkg{ICS} package. 
}


\seealso{
\code{\link[ICS]{ics}}
}


\examples{
# 3 source and 3 signals

S <- cbind(rt(1000, 4), rnorm(1000), runif(1000))
A <- matrix(rnorm(9), ncol = 3)
X <- S \%*\% t(A)
res<-FOBI(X)
MD(coef(res),A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }

