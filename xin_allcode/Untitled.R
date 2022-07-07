 library("bigalgebra")
 library("irlba")
 # Define an efficent matrix/transpose product:  
 matmul <- function(A, x, transpose=FALSE)
  {
     if(transpose)
       return(t( t(x) %*% A)) # i.e., t(A) %*% x + return (A %*% x)
    }
# Compute a small example and compare with other methods:  set.seed(1)
a = matrix(rnorm(1200),30)
b = matrix(rnorm(1200),40)
a = as.big.matrix(a)
b = as.big.matrix(b)
c = a%*%b
head(c)
