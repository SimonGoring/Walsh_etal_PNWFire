get.pos.min <- function (x, nx = 128, pointplot = TRUE, harm = 0, expand = 0, 
                         cycle = FALSE, ...) 
{
  pcafd <- x
  if (!(inherits(pcafd, "pca.fd"))) 
    stop("Argument 'x' is not a pca.fd object.")
  harmfd <- pcafd[[1]]
  basisfd <- harmfd$basis
  rangex <- basisfd$rangeval
{
  if (length(nx) > 1) {
    x <- nx
    nx <- length(x)
  }
  else x <- seq(rangex[1], rangex[2], length = nx)
}
  fdmat <- eval.fd(x, harmfd)
  meanmat <- eval.fd(x, pcafd$meanfd)
  dimfd <- dim(fdmat)
  nharm <- dimfd[2]
  plotsPerPg <- sum(par("mfrow"))
  harm <- as.vector(harm)
  if (harm[1] == 0)
    harm <- (1:nharm)
  if (length(dimfd) == 2) {
    
    pcmat <- list()
    
    for (jharm in 1:length(harm)) {
      if (jharm == 2) {
        op <- par(ask = TRUE)
        on.exit(par(op))
      }
      iharm <- harm[jharm]
      if (expand == 0) {
        fac <- sqrt(pcafd$values[iharm])
      }
      else {
        fac <- expand
      }
      vecharm <- fdmat[, iharm]
      pcmat[[jharm]] <- cbind(meanmat + fac * vecharm, meanmat - 
        fac * vecharm)
    }
  }
  else {
    if (cycle && dimfd[3] == 2) {
      meanmat <- drop(meanmat)
      for (jharm in 1:length(harm)) {
        if (jharm == 2) {
          op <- par(ask = TRUE)
          on.exit(par(op))
        }
        iharm <- harm[jharm]
{
  if (expand == 0) 
    fac <- 2 * sqrt(pcafd$values[iharm])
  else fac <- expand
}
        matharm <- fdmat[, iharm, ]
        mat1 <- meanmat + fac * matharm
        mat2 <- meanmat - fac * matharm
      }
    }
    else {
      for (jharm in 1:length(harm)) {
        if (jharm == 2) {
          op <- par(ask = TRUE)
          on.exit(par(op))
        }
        iharm <- harm[jharm]
        fac <- {
          if (expand == 0) 
            sqrt(pcafd$values[iharm])
          else expand
        }
        meanmat <- drop(meanmat)
        matharm <- fdmat[, iharm, ]
        nvar <- dim(matharm)[2]
        for (jvar in 1:nvar) {
          pcmat[[jharm]] <- cbind(meanmat[, jvar] + fac * matharm[, 
                                                         jvar], meanmat[, jvar] - fac * matharm[, 
                                                                                                jvar])
       
        }
      }
    }
  }
  
  outputs <- lapply(1:length(pcmat), function(jharm) data.frame(x = x, high = pcmat[[jharm]][,1], mean = meanmat, low = pcmat[[jharm]][,2]))
  
  return(list(sets = outputs, yrange = lapply(pcmat, range)))
}
