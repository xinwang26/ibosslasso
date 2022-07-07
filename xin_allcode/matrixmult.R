n = 500
p = 1000
y_tall<-matrix(rnorm(6000),1000,6) #y tall
x_tall<-matrix(rnorm(5000),5,1000)
y_wide<-t(y_tall)
x_wide<-t(x_tall)
y_tall_ff<-as.ff(y_tall) #y tall and ff
x_tall_ff<-as.ff(x_tall)
y_wide_ff<-as.ff(y_wide) #y tall and ff
x_wide_ff<-as.ff(x_wide)

dim(ffmatrixmult(x_tall_ff,y_tall_ff))



ffmatrixmult
function (x, y = NULL, xt = FALSE, yt = FALSE, ram.output = FALSE, 
          override.big.error = FALSE, ...) 
{
  {
    i1 <- NULL
    i2 <- NULL
  }
  dimx <- dim(x)
  if (!is.null(y)) 
    dimy <- dim(y)
  if (is.null(y)) 
    dimy <- dimx
  p <- max(c(dimx, dimy))
  n <- max(min(dimx), min(dimy))
  outDim <- inDim <- rep(NA, 2)
  outDim[1] <- dimx[xt + 1]
  outDim[2] <- dimy[2 - yt]
  inDim[1] <- dimx[2 - xt]
  inDim[2] <- dimy[yt + 1]
  if (inDim[1] != inDim[2]) 
    stop("non-conformable arguments")
  if (all(outDim > n) & (!override.big.error)) 
    stop("Returned value is at risk of being extremely large. Both dimensions of output will be fairly large.")
  if (xt & yt) 
    stop("For ff matrix algebra, set only one of xt or yt to TRUE")
  if (all(outDim == n) | (!"ff" %in% c(class(x), class(y))) | 
      ram.output) {
    out <- matrix(0, outDim[1], outDim[2])
  }
  else {
    out <- ff(0, dim = outDim, ...)
  }
  if (all(outDim == n)) {
    if ((xt) & (!yt)) 
      ffapply({
        out <- out + crossprod(x[i1:i2, ], y[i1:i2, ])
      }, X = x, MARGIN = 1)
    if ((!xt) & (yt)) 
      ffapply({
        out <- out + tcrossprod(x[, i1:i2], y[, i1:i2])
      }, X = x, MARGIN = 2)
    if ((!xt) & (!yt)) 
      ffapply({
        out <- out + x[, i1:i2] %*% y[i1:i2, ]
      }, X = x, MARGIN = 2)
  }
  if (outDim[1] > outDim[2] | (outDim[1] == p & outDim[2] == 
                               p)) {
    if ((xt) & (!yt)) 
      ffapply({
        out[i1:i2, ] <- crossprod(x[, i1:i2], y)
      }, X = x, MARGIN = 2)
    if ((!xt) & (yt)) 
      ffapply({
        out[i1:i2, ] <- tcrossprod(x[i1:i2, ], y)
      }, X = x, MARGIN = 1)
    if ((!xt) & (!yt)) 
      ffapply({
        out[i1:i2, ] <- x[i1:i2, ] %*% y
      }, X = x, MARGIN = 1)
  }
  if (outDim[1] < outDim[2]) {
    if ((xt) & (!yt)) 
      ffapply({
        out[, i1:i2] <- crossprod(x, y[, i1:i2])
      }, X = y, MARGIN = 2)
    if ((!xt) & (yt)) 
      ffapply({
        out[, i1:i2] <- tcrossprod(x, y[i1:i2, ])
      }, X = y, MARGIN = 1)
    if ((!xt) & (!yt)) 
      ffapply({
        out[, i1:i2] <- x %*% y[, i1:i2]
      }, X = y, MARGIN = 2)
  }
  return(out)
}
<environment: namespace:bootSVD>