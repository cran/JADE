`FOBI` <-
function(X, na.action = na.fail)
    {
    X <- na.action(X)
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    Col.center <- colMeans(X)
    X <- sweep(X, 2, Col.center, "-")
    
    COV <- crossprod(X)/n
    
    EVD.COV <- eigen(COV, symmetric=TRUE)
    
    COV.inv.sqrt <-  EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*% t(EVD.COV$vectors)
    Y <- tcrossprod(X, COV.inv.sqrt)
    r <- sqrt(rowSums(Y^2))
    Y <- r * Y

    COV4 <- crossprod(Y)
    
    EVD.COV4  <- eigen(COV4, symmetric=TRUE)
    W <- crossprod(EVD.COV4$vectors, COV.inv.sqrt)

    S <- tcrossprod(X, W)
    colnames(S) <- paste("IC.", 1:p, sep="")
    res<-list(W=W, EV=EVD.COV4$values, Xmu=Col.center, S=S)
    class(res) <- "bss"
    return(res)
    }
