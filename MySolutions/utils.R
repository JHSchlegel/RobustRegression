library(MASS)

f.tau.huber <- function(dat, loc, tc=1.345){
  ## Purpose: Korrekturfaktor Tau f�r die Varianz des Huber-M-Sch�tzern
  ## -------------------------------------------------------------------------
  ## Arguments: dat =  Daten
  ##            loc =  Lokations-Punkt
  ##            tc  =  Parameter der Huber Psi-Funktion
  ## -------------------------------------------------------------------------
  ## Author: Rene Locher
  ## Update: R. Frisullo 23.4.02 (Kommentar eingef�hrt)
  resid <- (dat - loc)/mad(dat)
  psi <- ifelse(abs(resid) <= tc, resid, sign(resid)*tc)
  psiP <- ifelse(abs(resid) <= tc, 1, 0)
  return(length(dat) * sum(psi^2) / sum(psiP)^2)
}



p.Serie2Aufg2c <- function (x, cex = 1){
  ## Purpose: Plot a profile.nls Object
  ## -------------------------------------------------------------------------
  ## Arguments: x =  profile.nls Object
  ##
  ## Description:  Modification of the R-function p.profileTraces:
  ##               Displays a series of plots of the profile t
  ##               function and the likelihood profile traces for the
  ##               parameters - without the first -
  ##               in a nonlinear regression model that has
  ##               been fitted with `nls' and profiled with
  ##               `profile.nls'
  ##
  ## -------------------------------------------------------------------------
  ## Author of p.profileTraces: Andreas Ruckstuhl
  ##                            modified by Roberto Frisullo
  ##
  nx <- names(x)
  np <- length(x)
  opar <- par(oma = c(2, 2, 1.5, 0), mfrow = c(np-1, np-1), mar = c(2,
                                                                    4, 0, 0) + 0.2)
  on.exit(par(opar))
  for (i in 2:np) {
    for (j in 2:i) {
      if (i == j) {
        if (!is.null(this.comp <- x[[i]])) {
          xx <- this.comp$par[, nx[i]]
          tau <- this.comp[[1]]
          plot(spline(xx, tau), xlab = "", ylab = "",
               type = "l", las = 1, mgp = c(3, 0.8, 0),
               cex = 0.5 * cex)
          points(xx[tau == 0], 0, pch = 3)
          pusr <- par("usr")
          if (is.R()) {
            mtext(side = 1, line = 0.8, at = -1/(2 *
                                                   (np-1)) + (i-1)/(np-1), text = nx[j], outer = TRUE,
                  cex = cex)
            mtext(side = 2, line = 0.8, at = 1 + 1/(2 *
                                                      (np-1)) - (i-1)/(np-1), text = nx[i], outer = TRUE,
                  cex = cex)
          }
          else {
            mtext(side = 1, line = 0.8, at = mean(pusr[1:2]),
                  text = nx[j], outer = TRUE, cex = cex)
            mtext(side = 2, line = 0.8, at = mean(pusr[3:4]),
                  text = nx[i], outer = TRUE, cex = cex)
          }
        }
      }
      else {
        if ((!is.null(x.comp <- x[[j]])) & (!is.null(y.comp <- x[[i]]))) {
          xx <- x.comp$par[, nx[j]]
          xy <- x.comp$par[, nx[i]]
          yx <- y.comp$par[, nx[j]]
          yy <- y.comp$par[, nx[i]]
          plot(xx, xy, xlab = "", ylab = "", las = 1,
               mgp = c(3, 0.8, 0), type = "n", xlim = range(c(xx,
                                                              yx)), ylim = range(c(xy, yy)), cex = 0.5 *
                 cex)
          lines(xx, xy, col = 2)
          lines(yx, yy, col = 3)
        }
      }
    }
    if (i < np)
      for (k in 1:(np - i + if (is.R())
        0
        else 1)) frame()
  }
  mtext(side = 3, line = 0.2, text = "t-Profil-Plot und Profilspuren",
        outer = TRUE, cex = 1.2 * cex)
}


f.calib <- function(cal.x, cal.pred, cal.limit, cal.y){
  ## Purpose: Bestimmung des Kalibrationsintervalls
  ## -------------------------------------------------------------------------
  ## Arguments: cal.x:      x-Werte der Kalibrationskurve
  ##            cal.pred    predict-Objekt fuer Kalibrationskurve
  ##            cal.limit   Halbe Breite des Vorhersageintervalls
  ##            cal.y:      beobachteter y-Wert, zu dem der x-Wert gesucht wird
  ##
  ## -------------------------------------------------------------------------
  ## Author: Andreas Ruckstuhl, Roberto Frisullo
  
  I1 <- uniroot(function(x) cal.y -
                  approx(cal.x,cal.pred-cal.limit,xout=x)$y,
                range(cal.x))$root
  Ic <- uniroot(function(x) cal.y -
                  approx(cal.x,cal.pred,xout=x)$y,
                range(cal.x))$root
  I2 <- uniroot(function(x) cal.y -
                  approx(cal.x,cal.pred+cal.limit,xout=x)$y,
                range(cal.x))$root
  
  return(data.frame(lower=min(I1,I2), center=Ic, upper=max(I1,I2)))
}


p.ldv <- function(ldafit, data=D[,2:3], group=D$Klasse){
  ## Rks, Juni 02
  ldaDV <- predict(ldafit, data, dimen=2)$x
  eqscplot(ldaDV[,1], ldaDV[,2], type="n",
           xlab="first discriminant variable",
           ylab="second discriminant variable", main="", las=1)
  if(is.factor(group)) NoL <- length(levels(group))
  else NoL <- length(unique(group))
  for(i in 1:NoL){
    x.c <- LETTERS[1:NoL][i]; ok <- group==x.c
    text(ldaDV[ok,1], ldaDV[ok,2], labels=rep(x.c, sum(ok)),
         cex=0.6, col=i+1)
  }
  ## Add centre of groups
  x <- apply(diag(ldafit$prior) %*% ldafit$means, 2, sum)
  z <- scale(ldafit$means, x, F)      ## remove overall means
  points(z %*% ldafit$scaling , col=2:5, pch=20)
  points(z %*% ldafit$scaling, pch=1)
}


p.predplot <- function(ldafit, data=D[,2:4], group=D$Klasse, len=100){
  ## Rks, Juni 02
  ## From V&R, 2nd ed
  g <- as.factor(group)
  lev <- levels(g)
  plot(data[,1], data[,2], type="n", las=1)
  text(data[,1], data[,2], labels=as.character(g), col=unclass(g)+1, cex=0.6)
  h <- range(data[,1]); xp <- seq(h[1], h[2], length=len)
  h <- range(data[,2]); yp <- seq(h[1], h[2], length=len)
  xyp <- expand.grid(x1=xp, x2=yp)
  Z <- predict(ldafit, newdata=xyp)
  zp <- Z$post[,2] - pmax(Z$post[,1],Z$post[,3],Z$post[,4])
  contour(xp, yp, matrix(zp, len), add=T, levels=0, drawlabels=F)
  zp <- Z$post[,3] - pmax(Z$post[,1],Z$post[,2],Z$post[,4])
  contour(xp, yp, matrix(zp, len), add=T, levels=0, drawlabels=F)
  zp <- Z$post[,4] - pmax(Z$post[,1],Z$post[,2],Z$post[,3])
  contour(xp, yp, matrix(zp, len), add=T, levels=0, drawlabels=F)
  invisible()
}

rlda <- function(x, grouping, prior=proportions,
                 tol=1e-04, nu=5, ...){
  ## Rks, Juni 02
  ## modified version of lda.default: Estimation with cov.mcd and
  ## robust estimation of location
  if (is.null(dim(x)))
    stop("x is not a matrix")
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  if (n != length(grouping))
    stop("nrow(x) and length(grouping) are different")
  g <- as.factor(grouping)
  lev <- lev1 <- levels(g)
  counts <- as.vector(table(g))
  if (any(counts == 0)) {
    warning(paste("group(s)", paste(lev[counts == 0], collapse = " "),
                  "are empty"))
    lev1 <- lev[counts > 0]
    g <- factor(g, levels = lev1)
    counts <- as.vector(table(g))
  }
  proportions <- counts/n
  ng <- length(proportions)
  if (any(prior < 0) || round(sum(prior), 5) != 1)
    stop("invalid prior")
  if (length(prior) != ng)
    stop("prior is of incorrect length")
  names(prior) <- names(counts) <- lev1
  ## robust estiamtion of group location
  group.means <- matrix(NA, ncol=ncol(x), nrow=length(lev1),
                        dimnames=list(lev1,names(x)))
  for(i in lev1)
    group.means[i,] <- cov.rob(x[g==i,], method="mcd")$center
  f1 <- sqrt(diag(var(x - group.means[g, ])))
  if (any(f1 < tol))
    stop(paste("variable(s)", paste(format((1:p)[f1 < tol]),
                                    collapse = " "), "appear to be constant within groups"))
  scaling <- diag(1/f1, , p)
  cov <- n/(n - ng) * cov.rob((x - group.means[g, ]) %*% scaling,
                              method="mcd")$cov
  sX <- svd(cov, nu = 0)
  rank <- sum(sX$d > tol^2)
  if (rank < p)
    warning("variables are collinear")
  scaling <- scaling %*% sX$v[, 1:rank] %*% diag(sqrt(1/sX$d[1:rank]),
                                                 , rank)
  ##
  xbar <- apply(prior %*% group.means, 2, sum)
  fac <- 1/(ng - 1)
  X <- sqrt((n * prior) * fac) * scale(group.means, center = xbar,
                                       scale = FALSE) %*% scaling
  X.s <- svd(X, nu = 0)
  rank <- sum(X.s$d > tol * X.s$d[1])
  scaling <- scaling %*% X.s$v[, 1:rank]
  if (is.null(dimnames(x)))
    dimnames(scaling) <- list(NULL, paste("LD", 1:rank, sep = ""))
  else {
    dimnames(scaling) <- list(colnames(x), paste("LD", 1:rank,
                                                 sep = ""))
    dimnames(group.means)[[2]] <- colnames(x)
  }
  structure(list(prior = prior, counts = counts, means = group.means,
                 scaling = scaling, lev = lev, svd = X.s$d[1:rank], N = n,
                 call = match.call()), class = "lda")
}