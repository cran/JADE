\name{SIR}
\alias{SIR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Signal to Interference Ratio }
\description{
  Computes the signal to interference ratio between true and estimated signals 
}
\usage{
SIR(S, S.hat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{Matrix or dataframe with the true numeric signals.}
  \item{S.hat}{Matrix or dataframe with the estimated numeric signals.}
}
\details{
  The signal to interference ratio is measured in dB and values over 20 are thought to be good.
  It is scale and permutation invariant and can be seen as measuring the correlation between the matched true and estimated signals. 
}
\value{
  The value of the signal to interference ratio.
}
\references{
   \cite{Eriksson, J., Karvanen, J. and Koivunen, V. (2000), Source distribution adaptive maximum likelihood estimation in ICA model, \emph{Proceedings of
   the second international workshop on independent component analysis and blind source separation (ICA 2000)},  227--232. }  
}

\author{Klaus Nordhausen}

\seealso{\code{\link{amari.error}}, \code{\link{ComonGAP}}}
\examples{
S <- cbind(rt(1000, 4), rnorm(1000), runif(1000))
A <- matrix(rnorm(9), ncol = 3)
X <- S \%*\% t(A)

S.hat <- JADE(X, 3)$S
SIR(S, S.hat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }



