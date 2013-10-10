# Method SOBI

SOBI <- function(X,...) UseMethod("SOBI")

# main function for SOBI
#
# input:
#  x = data matrix
#  k = either a intereger value, or a vector with integers. If a, then the lagged ACOVs 1:k are used otherwise those mentioned in k. 
#  method = method used for joint diagonalization, either rjd or djd
#  eps = eps for the JD
#  maxiter = maxiter for the JD

# output
#
# list of class "bss" with components
#   W = unmixing matrix
#   EV = Eigenvalues
#   k = lag used
#   S = sources as a time series object

SOBI.default <- function(X, k=12, method="rjd.fortran", eps = 1e-06, maxiter = 100, ...)
    {
    if (length(k)==1) k <- 1:k 
    nk <- length(k)
    method <- match.arg(method, c("rjd.fortran","rjd", "djd"))
    
    MEAN <- colMeans(X)
    COV <- cov(X)
    EVD <- eigen(COV, symmetric = TRUE)
    COV.sqrt.i <- EVD$vectors %*% (diag(EVD$values^(-0.5))) %*% t(EVD$vectors)
    X.C <- sweep(X,2,MEAN,"-")
    Y <- X.C %*% COV.sqrt.i
    p <- ncol(X)
    R <- array(0, dim=c(p,p,nk))
    n <- nrow(X) 
    
    for (i in 1:nk){
        Yt <- Y[1:(n-k[i]),]
        Yti <- Y[(1+k[i]):n,]
        Ri <- crossprod(Yt,Yti)/nrow(Yt)
        R[,,i] <- (Ri+t(Ri))/2
        }
    
    
    JD <- switch(method,
        "rjd.fortran"={
               rjd.fortran(R, eps=eps, maxiter=maxiter)$V
               }
        ,
        "rjd"={
               rjd(R, eps=eps, maxiter=maxiter)$V
               }
        ,
        "djd"={
               djd(R, eps=eps, maxiter=maxiter,...)
               }
        )
    W <-crossprod(JD, COV.sqrt.i)
    
    S <- tcrossprod(X.C,W)
    S <- ts(S, names=paste("Series",1:p))
    
    RES <- list(W=W, k=k, method=method, S=S)
    class(RES) <- "bss"
    RES
    }

SOBI.ts <- function(X, ...)
    {
    x <- as.matrix(X)
    RES <- SOBI.default(x,...)
    S <- RES$S
    attr(S, "tsp") <- attr(X, "tsp")
    RES$S <- S
    RES
    }
