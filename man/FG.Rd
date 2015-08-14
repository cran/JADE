\name{FG}
\alias{FG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Joint Diagonalization of Real Positive-definite Matrices
}
\description{
This is a slightly modified version of Flury's FG algorithm for the joint diagonalization of k positive-definite matrices.
The underlying function is written in C.
}
\usage{
FG(X, weight = NULL, init = NULL, maxiter = 100, eps = 1e-06, na.action = na.fail)
}
\arguments{
  \item{X}{A matrix of k stacked pxp matrices with dimension c(kp,p) or an array with dimension c(p,p,k).}
  \item{weight}{A vector of length k to give weight to the different matrices, if NULL, all matrices have equal weight}
  \item{init}{???, if NULL, ???}
  \item{maxiter}{ Maximum number of iterations.}
  \item{eps}{ Convergence tolerance.}
  \item{na.action}{A function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\value{
A list with the components
  \item{V }{An orthogonal matrix.}
  \item{D }{A stacked matrix with the diagonal matrices or an array with the diagonal matrices. The form of the output
            depends on the form of the input.}
  \item{iter}{The Fortran function returns also the number of iterations.}
}
\references{
\cite{Flury, B. D. (1998), Common principal components and related models, Wiley, New York.}
}

\author{Jari Miettinen}


\seealso{
 \code{\link{rjd}}, \code{\link{rjd.fortran}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }

