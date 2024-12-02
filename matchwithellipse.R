library("hadron")

## an ellipse is given by the fixed paramters of the ratio of the axes, its center points, its angle, and the variable length of the major axis.
## a line is given by its slope and intercept.
## We take the line as a polar, and calculate its pole with respect to the ellipse.
## If the polar is a tangent to the ellipse, the pole will be on the ellipse and fulfill the ellipse equation.
## Otherwise, the equation is not fulfilled.
## With this, we can reduce the problem of finding the touching point between a line and an ellipse with variable diamter to a root-finding procedure.
## We return the pole, and its value when inserted in the ellipse equation, and the distance between center and pole.
tangent_general <- function(a, m, n, r, xc, yc, theta = pi / 4) {
  ## slope and intersect of hamiltonian parametrisation, ellipse radius in x-direction, ratio of a/b, central points of ellipse
  ## angle by which the ellipse is tilted, is always 45° here.
  ## source formula: wikipedia (ellipse, pole and polar)
  ## a/b = r
  res <- data.frame(
    x0 = c(), y0 = c(),
    onellipse = c(),
    distancex = c(), distancey = c(),
    radx = c(), rady = c()
  )
  for (i in seq_along(a)) {
    s <- sin(theta)
    c <- cos(theta)
    b <- a[i] / r
    Axx <- a[i]^2 * s^2 + b^2 * c^2
    Axy <- (b^2 - a[i]^2) * s * c
    Ayy <- a[i]^2 * c^2 + b^2 * s^2
    Bx <- (-Axx * xc - Axy * yc)
    By <- (-Axy * xc - Ayy * yc)
    C <- Axx * xc^2 + 2 * Axy * xc * yc + Ayy * yc^2 - a[i]^2 * b^2
    D <- m
    E <- -1
    Fpar <- n

    amat <- matrix(c(Axx, Axy, Bx, Axy, Ayy, By, Bx, By, C), nrow = 3, byrow = T)
    xval <- solve(a = amat, b = c(D, E, Fpar))
    x0 <- xval[1] / xval[3]
    y0 <- xval[2] / xval[3]
    newline <- data.frame(
      x0 = x0, y0 = y0,
      onellipse = (Axx * x0^2 + 2 * Axy * x0 * y0 + Ayy * y0^2 + 2 * Bx * x0 + 2 * By * y0 + C),
      distancex = (x0 - xc), distancey = (y0 - yc), radx = a[i], rady = b
    )
    # set calculated points into ellipse to see if the tangent matches the radius
    res <- rbind(res, newline)
  }
  return(res)
}

## Use the secant method to find the diameter of an ellipse that makes the provided line a tangent (of the upper half).
## return the major axis length needed for this.
findrootellipse <- function(m, n, r, xc, yc, theta = pi / 4, tol = 1e-12, max = 1000) {
  myf <- function(x) {
    myt <- tangent_general(a = x, m = m, n = n, r = r, xc = xc, yc = yc, theta = theta)
    return(myt$onellipse * 1e+6)
  }
  root <- cmna::secant(myf, 10, tol = tol, m = max)
  # root <- uniroot(myf, interval=c(0.01, 100), tol = tol, maxiter=max)
  return(root)
}

## draw ellipse with correlation coefficient estimated from bootstraps
## x, y: bootstrapsamples
## mean, sd: can be given separately
## std: how large the ellipse should be in terms of standard deviations
## choose whether error bars, rectanfle around ellipse should be drawn
## shifted by 45°
## Draw confidence ellipse based on https://stackoverflow.com/a/41821484 and  https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html, which is based on https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html.

draw_ellipse <- function(x, y,
                         meanx = mean(x), meany = mean(y),
                         sdx = sd(x), sdy = sd(y),
                         points = 10000, verbose = FALSE,
                         nstd = 1, rectangle = FALSE,
                         bars = FALSE, ...) {
  ## if one bootstrap sample is NA, we set the corresponding other one to NA as well
  x[is.na(y)] <- NA
  y[is.na(x)] <- NA
  x <- na.omit(x)
  y <- na.omit(y)
  corcof <- cor(x, y) ## pearson correlation coefficient
  ## set major axes of elippsis
  ell_rad_x <- sqrt(1 + corcof)
  ell_rad_y <- sqrt(1 - corcof)
  # angle of major axis with x axis phi or tau - in correlation ellipses, this is always 45°
  phi <- pi / 4
  if (verbose) {
    print(corcof)
  }

  ## ellipsis from parametric equation
  ## rescale ellipsis by standard deviation, number of standard deviations that should be displayd
  t <- seq(0, 2 * pi, len = points)
  xpoints <- (ell_rad_x * cos(t) * cos(phi) - ell_rad_y * sin(t) * sin(phi)) * sdx * nstd
  ypoints <- (ell_rad_x * cos(t) * sin(phi) + ell_rad_y * sin(t) * cos(phi)) * sdy * nstd
  # ~   xpoints <- (ell_rad_x*cos(t)*cos(phi)*sdx - ell_rad_y*sdy*sin(t)*sin(phi)) * nstd
  # ~   ypoints <- (ell_rad_x*cos(t)*sin(phi)*sdx + ell_rad_y*sdy*sin(t)*cos(phi)) * nstd

  ## plot shifted ellipsis, point with errors, rectangle made up of errors
  plotwitherror(x = meanx + xpoints, meany + ypoints, ...)
  if (bars) plotwitherror(x = meanx, y = meany, dx = sdx, dy = sdx, rep = TRUE)
  if (rectangle) {
    polygon(
      x = c(meanx - sdx, meanx + sdx, meanx + sdx, meanx - sdx),
      y = c(meany - sdy, meany - sdy, meany + sdy, meany + sdy)
    )
  }
}



draw_ellipse_general <- function(meanx, meany, radx, rady, phi = pi / 4,
                                 points = 10000, verbose = FALSE,
                                 nstd = 1, rectangle = FALSE,
                                 bars = FALSE, ...) {
  t <- seq(0, 2 * pi, len = points)
  xpoints <- (radx * cos(t) * cos(phi) - rady * sin(t) * sin(phi)) * nstd
  ypoints <- (radx * cos(t) * sin(phi) + rady * sin(t) * cos(phi)) * nstd
  # ~   xpoints <- (ell_rad_x*cos(t)*cos(phi)*sdx - ell_rad_y*sdy*sin(t)*sin(phi)) * nstd
  # ~   ypoints <- (ell_rad_x*cos(t)*sin(phi)*sdx + ell_rad_y*sdy*sin(t)*cos(phi)) * nstd

  ## plot shifted ellipsis, point with errors, rectangle made up of errors
  plotwitherror(x = meanx + xpoints, meany + ypoints, ...)
}

draw_ellipse_general_lines <- function(meanx, meany, radx, rady, phi = pi / 4,
                                       points = 10000, verbose = FALSE,
                                       nstd = 1, rectangle = FALSE,
                                       bars = FALSE, ...) {
  if (length(nstd) > 1) stop("this plot will be too messsy, use draw_ellipse_general instead")
  t <- seq(0, 2 * pi, len = points)
  xpoints <- (radx * cos(t) * cos(phi) - rady * sin(t) * sin(phi)) * nstd
  ypoints <- (radx * cos(t) * sin(phi) + rady * sin(t) * cos(phi)) * nstd
  # ~   xpoints <- (ell_rad_x*cos(t)*cos(phi)*sdx - ell_rad_y*sdy*sin(t)*sin(phi)) * nstd
  # ~   ypoints <- (ell_rad_x*cos(t)*sin(phi)*sdx + ell_rad_y*sdy*sin(t)*cos(phi)) * nstd

  ## plot shifted ellipsis, point with errors, rectangle made up of errors
  lines(x = meanx + xpoints, meany + ypoints, ...)
}


## get matching beta, P to Hamiltonian
getmatchingbeta <- function(finalresult, interpolation) {
  interpolatedbeta <- c()
  for (i in seq(1, length(finalresult$betaiso), 1)) {
    Vs <- finalresult$p3contlim[i]
    index <- NA
    # determine in which interval of the V_t V_s is sitting by selecting
    # lowest possible interval for which the upper limit is bigger than V_s
    for (j in seq(1, length(hamresults$beta))) {
      if (length(interpolation$upper[interpolation$index == j]) > 0) {
        if (Vs < interpolation$upper[interpolation$index == j]) {
          index <- j
          break
        }
      }
    }
    # highest V_s could be bigger than highest V_t, so there is not an index for every y
    if (!is.na(index)) {
      t <- (Vs - interpolation$b[index]) / interpolation$a[index]
      interpolatedbeta[i] <- t
    }
  }
  return(interpolatedbeta)
}


## columns/content:
# data: betacontlim, dbetacontlim, plaq3contlim, dplaq3contlim, betaiso
# bsdata: beta, plaq3
# hamres: a, b, lowerx, upperx, with a and b giving the interpolation a*x+b in that interval
# indices: uprange, lowrange, bsindex, resindex

# Ellipse angle taken from https://demonstrations.wolfram.com/EllipseRepresentingTheConfidenceRegionOfACovarianceMatrix/, https://web2.ba.infn.it/~pompili//teaching/data_analysis_lab/Glen_Cowan_Statistical_Data_Analysis_1998.pdf.
getmatchingellipse <- function(data, bsdata, hamres, indices, verbose = F) {
  res <- data.frame(betaiso = c(), devstd = c(), root = c(), x0 = c(), y0 = c(), segment = c(), radx = c(), rady = c(), theta = c(), cor = c(), index = c())

  for (i in seq_along(indices$uprange)) {
    if (verbose) print(data[i, ])

    resbs$beta[, indices$bsindex[i]][is.na(resbs$plaq3[, indices$bsindex[i]])] <- NA
    resbs$plaq3[, indices$bsindex[i]][is.na(resbs$beta[, indices$bsindex[i]])] <- NA
    ## calculate covariance matrix to set ellipse parameters
    cor <- cor(
      x = na.omit(resbs$beta[, indices$bsindex[i]]),
      y = na.omit(resbs$plaq3[, indices$bsindex[i]])
    )
    xrad <- data$dbetacontlim[indices$resindex[i]]
    yrad <- data$dplaq3contlim[indices$resindex[i]]
    theta <- 0.5 * atan(2 * cor * xrad * yrad / (xrad^2 - yrad^2))
    ## calculate minimum distance for several line segmants, but only use the line segment if the touching point is on it.
    if (xrad < yrad) theta <- theta - pi / 2
    for (j in seq(indices$lowrange[i], indices$uprange[i])) {
      root <- findrootellipse(
        m = hamres$a[j], n = hamres$b[j],
        r = xrad / yrad,
        xc = data$betacontlim[indices$resindex[i]], yc = data$plaq3contlim[indices$resindex[i]],
        theta = theta, tol = 1e-6
      )
      dat <- tangent_general(
        m = hamres$a[j], n = hamres$b[j],
        r = xrad / yrad,
        xc = data$betacontlim[indices$resindex[i]], yc = data$plaq3contlim[indices$resindex[i]],
        a = root, theta = theta
      )
      if (dat$x0 > hamres$lowerx[j] && dat$x0 < hamres$upperx[j]) {
        newline <- data.frame(
          betaiso = data$betaiso[indices$resindex[i]], devstd = root / xrad, root = root,
          x0 = dat$x0, y0 = dat$y0, segment = j, radx = dat$radx, rady = dat$rady, theta = theta, cor = cor, index = i
        )
        if (verbose) print(i)
        if (verbose) print(newline)
        res <- rbind(res, newline)
      }
    }
    # check correlation
    # check relation between radii and area/sigma
    # If there is a match with two line segments, choose the one with the smaller distance
    if (length(res$index[res$index == i]) > 1) {
      minindex <- which.min(res$devstd[res$index == i])
      doubleentries <- rep(T, length(res$index[res$index == i]))
      doubleentries[minindex] <- F
      if (verbose) print(paste("double results, keep index", minindex))
      res[res$index == i, ][doubleentries, ] <- NA
      res <- na.omit(res)
    } else if (length(res$index[res$index == i]) == 0) {
      print(paste("no result found for index", i))
      newline <- data.frame(
        betaiso = data$betaiso[indices$resindex[i]], devstd = "NA", root = "NA",
        x0 = "NA", y0 = "NA", segment = "NA", radx = "NA", rady = "NA", theta = theta, cor = cor, index = i
      )
      res <- rbind(res, newline)
    }
  }
  na.omit(res)
}
