## Least Square via QR factorization
## Basic from : ||Ax - b||^2 = 0
## Find x to make least square be approximate to 0

## QR factorization function (through Gram-Schmidt Algorithm)
GSQR <- function(x)
{
  # Get the number of rows and columns of A
  n <- ncol(x)
  m <- nrow(x)
  
  # Initialize Q and R
  Q <- matrix(0, m, n)
  R <- matrix(0, n, n)
  
  for (j in 1 : n)
  {
    # 1st term : v1 = a1
    v = x[, j]
    
    # 2nd term ~ :
    if (j > 1)
    {
      for (i in 1 : (j - 1))
      {
        # inner product
        R[i, j] <- t(Q[, i]) %*% x[, j]
        
        # make all elements orthogonal
        v <- v - R[i, j] * Q[, i] 
      }      
    }
    
  # Norm of the j-th diagonal of R
  R[j, j] <- sqrt(sum(v ^ 2))
  
  # Finalize elements of Q
  Q[, j] <- v / R[j, j]
}
  
  # Return Q and R
  QR <- list('Q' = Q, 'R' = R)
  return(QR)
}

## Processing b function (before backward substitution)
b_ <- function(Q, b)
{
  QT <- t(Q)
  b_ <- QT %*% b
  
  return(b_)
}
## Back Susbstitution function
BackSubstitution <- function(R, b)
{ 
  n = dim(R)[1]
  
  for (j in seq(n, 1, -1))
  { 
    # 1st term : x_n,1 = b_n,1 / R_n,n 
    b[j, 1] = b[j, 1] / R[j, j]
    
    # 2nd term ~ : x_n-1,1 = (b_n-1,1 - Rn-1,n * x_n,1) / R_n-1,n-1 
    if((j - 1) > 0)
    {
      b[1 : (j - 1), 1] = b[1 : (j - 1), 1] - (b[j, 1] * R[1 : (j - 1), j])
    }
  }  
  return(b)
}

## Input matrix : examples 12.4 (p. 228)
A <- rbind(c(0.97, 1.86, 0.41), c(1.23, 2.18, 0.53), c(0.80, 1.24, 0.62), 
           c(1.29, 0.98, 0.51), c(1.10, 1.23, 0.69), c(0.67, 0.34, 0.54),
           c(0.87, 0.26, 0.62), c(1.10, 0.16, 0.48), c(1.92, 0.22, 0.71),
           c(1.29, 0.12, 0.62))

b <- cbind(rep(1e3, 10))

## QR decomposition via Gram-schmidt algorithm
QR <- GSQR(A)
R <- QR$R
Q <- QR$Q
QR

## Rx = (Q^T) * b
b <- b_(Q, b)

## Back Substitution
x <- BackSubstitution(R, b)
x

