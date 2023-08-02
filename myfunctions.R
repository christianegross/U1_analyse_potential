library(hadron)

# All readloopfile type of functions read in data
# from the measurement into a cf object. Later these
# objects can be used to determine the potential or
# other ratios of the observables.
# All work by selecting the appropriate columns of data.
# The data are in the format
# one line per configuration, one column per observable.
# This is different to the usual structure used by the hadron-package
# but makes it possible to measure the normal and the
# sideways potential with one file.
# file is the filename, path the path to the file
# the number of rows that are read in can be limited with the maxrows parameter,
# a negative number means all rows are read

readloopfilecfrotated <- function (file, path  = "",
                                    skip = 0, Nsmax,
                                    yt, zerooffset = 0, every = 1, maxrows = -1,  maxt) {
    # reads in the Wilson loop data for the sideways( = rotated) potential,
    # used to calculate V(yt), into a cf-container
    # Read in W(x, yt) for all possible x
    # Nsmax: maximum extent measured in the spacial x-direction,
    # which is kept constant,
    # typically half of spatial extent of lattice
    # yt: read columns with R = sqrt(x^2 + yt^2) for x starting at zero
    # set zerooffset to 1 if W(x = 0, y = 0) and so on has also been measured
    # skip: number of configurations that are stripped from the beginnig
    # every: step between read lines
    # read in data, select columns that are needed, give columns
    # and additional info to cf object
  if (missing(maxt))  maxt <- Nsmax
  tmp <- as.matrix(read.table(paste(path, file, sep = ""), skip = skip, nrows = maxrows))
  Dm <- dim(tmp)
  confno <- tmp[[Dm[2]]]
  nom <- length(tmp[, 1])
  ii <- seq((yt + zerooffset - 1) * (Nsmax + zerooffset) + 1 + zerooffset,
         (yt + zerooffset) * (Nsmax + zerooffset), 1)
  ret <- cf_meta(nrObs = 1, Time =  maxt, nrStypes = 1)
  ret <- cf_orig(ret, cf = tmp[seq(1, nom, every), ii])
  ret$conf.index <- confno

  return(invisible(ret))
}

readloopfilecfsub <- function (file, path = "", skip = 0,
                            Nytmax, x, Nsmax,
                            zerooffset = 0, every = 1, maxrows = -1, maxt) {
    # reads in the Wilson loop data for the normal potential,
    # from which the anisotropy could be determined by subtracting the points,
    # used to calculate V(x), into a cf-container
    # Nytmax: maximum extent measured in the direction that is kept constant,
    # typically half of spatial extent of lattice
    # skip: number of configurations that are stripped from the beginnig
    # every: step between read lines
    # Read in W(x, yt) for all possible yt
    # x: read columns with x = x
    # set zerooffset to 1 if W(x = 0, y = 0) and so on has also been measured
    # read in data, select columns that are needed, give columns
    # and additional info to cf object

  if(missing(maxt))  maxt <- Nytmax
  tmp <- as.matrix(read.table(paste(path, file, sep = ""), skip = skip, nrows = maxrows))
  Dm <- dim(tmp)
  confno <- tmp[[Dm[2]]]
  nom <- length(tmp[, 1])
  ii <- seq((Nsmax + zerooffset) * zerooffset + zerooffset + x,
            (Nsmax + zerooffset) * (Nytmax + zerooffset), Nsmax + zerooffset)
  ret <- cf_meta(nrObs = 1, Time =  maxt, nrStypes = 1)
  ret <- cf_orig(ret, cf = tmp[seq(1, nom, every), ii])
  ret$conf.index <- confno

  return(invisible(ret))
}

readloopfilecfsmall <- function (file, path = "", skip = 0,
                                Nsmax, x, y, Ntmax, start = 0, maxrows = -1,  maxt) {
    # Reads in nonplanar Wilson loops, measured only for small distances.
    # Nsmax: maximum extent measured in the spacial x- and y-direction,
    # typically min(4, lattice size, those directions are kept fixed here.
    # Ntmax: maximum extent measured in the spatial direction, t is changed here.
    # x, y: read columns with W(t, x, y) for all possible t
    # Ntmax: temporal extent of the lattice
    # skip: number of configurations that are stripped from the beginnig
    # start: smallest t that was measured
    # read in data, select columns that are needed, give columns
    # and additional info to cf object
  if (missing(maxt))  maxt <- Ntmax
  tmp <- as.matrix(read.table(paste(path, file, sep = ""), skip = skip, nrows = maxrows))
  Dm <- dim(tmp)
  confno <- tmp[[Dm[2]]]
  ii <- seq(((x) * (Nsmax + 1) + y + 1) + (start * (Nsmax + 1) * (Nsmax + 1)),
            (Nsmax + 1) * (Nsmax + 1) * (Ntmax + 1), (Nsmax + 1) * (Nsmax + 1))
  ret <- cf_meta(nrObs = 1, Time =  maxt, nrStypes = 1)
  ret <- cf_orig(ret, cf = tmp[, ii])
  ret$conf.index <- confno

  return(invisible(ret))
}

readloopfilecfplaquette <- function (file, path = "", skip = 0,
                                    Nsmax, Ntmax, every = 1, maxrows = -1) {
    # reads W(x = 1, t = 1, y = 0) and W(x = 1, y = 1, t = 0)
    # into one cf-container from the file with nonplanar loops
    # Nsmax: maximum extent measured in the spacial x-direction,
    # typically min(4, lattice size)
    # Ntmax: temporal extent of the lattice
    # skip: number of configurations that are stripped from the beginnig
    # every: step between read lines
  tmp <- as.matrix(read.table(paste(path, file, sep = ""), skip = skip, nrows = maxrows))
  Dm <- dim(tmp)
  confno <- tmp[[Dm[2]]]
  Time <- 2
  nom <- length(tmp[, 1])
  ii <- c((Nsmax + 1) + 2,
        2 * (Nsmax + 1) * (Nsmax + 1), (Nsmax + 1) * (Nsmax + 1))
#~   print(ii)
  ret <- cf_meta(nrObs = 1, Time = Time, nrStypes = 1)
  ret <- cf_orig(ret, cf = tmp[seq(1, nom, every), ii])
  ret$conf.index <- confno

  return(invisible(ret))
}

readloopfilecfonecolumn <- function (file, path = "",
                                    skip = 0, column = 1, every = 1, maxrows = -1) {
    # reads column into a cf container
    # column has to be read twice so the expectation value can be
    # determined with bootstrap.cf
    # column: index of the column that is to be read
    # skip: number of configurations that are stripped from the beginnig
    # every: step between read lines
    # used instead of uwerr to get the expectation value
    # to get bootstrap samples
  tmp <- as.matrix(read.table(paste(path, file, sep = ""), skip = skip, nrows = maxrows))
  Dm <- dim(tmp)
  confno <- tmp[[Dm[2]]]
  Time <- 1
  ii <- c(column, column)
  ret <- cf_meta(nrObs = 1, Time = Time, nrStypes = 1)
  ret <- cf_orig(ret, cf = tmp[, ii])
  ret$conf.index <- confno

  return(invisible(ret))
}

# reads in the data for the given parameters,
# determines the correlators = expectation values
# of the Wilson loop via bootstrap
# then plots the correlators on a log-scale
# to the file that is currently open
# If at least one correlator is negative,
# the correlators starting from x = 2 are also plotted on a linear scale
# bootsamples: number of bootstrap samples that are generated
# title: title for the plot, e.g. details about loop extent.


calcplotWloopsideways <- function (file, skip, Ns, yt, bootsamples,
                                 title, path = "", nsave = 100,
                                 zerooffset = 0, every = 1, l = 2, fraction = 0.5, maxrows = -1) {
    # plots loops for the sideways potential
    # x is changed, yt is kept constant
    # nsave: number of steps MC-sweeps between saving configurations
    # set zerooffset to 1 if W(x = 0, y = 0) and so on has also been measured
    # skip: number of configurations that are stripped from the beginnig
    # every: step between read lines
    # l: median blocking size for bootstrap
    WL <- readloopfilecfrotated(file = file, path = path, skip = skip,
            Nsmax = Ns * fraction, yt = yt, zerooffset = zerooffset, every = every, maxrows = maxrows)
    WL <- bootstrap.cf(WL, boot.R = bootsamples, boot.l = l)
    title <- sprintf("%s, %d configs\n used every %d, boot.l = %d\n",
            title, (length(WL$cf[, 1]) + skip) * nsave, nsave * every, l)
    plot(WL, log = "y", xlab = "x/a_s", ylab = "C(x)",
            main = sprintf("%s, logscale", title))
    if (-1 %in% sign(WL$cf0)) {
      plot(WL, xlab = "x/a_s", ylab = "C(x)",
            main = sprintf("%s, linear scale", title),
            ylim = c(-WL$cf0[3], WL$cf0[3]))
      points(x = seq(-1, Ns * fraction + 1, 1),
            y = rep(0, times = length(seq(-1, Ns * fraction + 1, 1))),
            col = 2, type = "l")
      message(WL$cf0, " ")
    }
    return(WL)
}


calcplotWloopnormalspatial <- function (file, skip, Ns, x, bootsamples,
                                    title, path = "", zerooffset = 0,
                                    every = 1, nsave = 100, l = 2, fraction = 0.5, maxrows = -1) {
    # plots loops for the spatial normal potential
    # y is cahnged, x is kept constant
    # nsave: number of steps MC-sweeps between saving configurations
    # set zerooffset to 1 if W(x = 0, y = 0) and so on has also been measured
    # skip: number of configurations that are stripped from the beginnig
    # every: step between read lines
    # l: median blocking size for bootstrap
    WL <- readloopfilecfsub(file = file, skip = skip, Nytmax = Ns * fraction, x = x,
                Nsmax = Ns * fraction, zerooffset = zerooffset, every = every, maxrows = maxrows)
    WL <- bootstrap.cf(WL, boot.R = bootsamples, boot.l = l)
    title <- sprintf("%s, %d configs\n used every %d, boot.l = %d\n", title,
                    (length(WL$cf[, 1]) + skip) * nsave, nsave * every, l)

    plot(WL, log = "y", xlab = "y/a_s", ylab = "C(y)",
            main = sprintf("%s, logscale", title))
    if (-1 %in% sign(WL$cf0)) {
      plot(WL, xlab = "y/a_s", ylab = "C(y)",
            main = sprintf("%s, linear scale", title),
            ylim = c(-WL$cf0[3], WL$cf0[3]))
      points(x = seq(-1, Ns * fraction + 1, 1),
            y = rep(0, times = length(seq(-1, Ns * fraction + 1, 1))),
                col = 2, type = "l")
   }
    return(WL)
}

calcplotWloopnormaltemporal <- function (file, skip, Ns, Nt, x, bootsamples,
                                        title, path = "", zerooffset = 0,
                                        every = 1, nsave = 100, l = 2, fraction = 0.5, maxrows = -1) {
    # plots loops for the spatial temporal potential
    # t is changed, x is kept constant
    # nsave: number of steps MC-sweeps between saving configurations
    # set zerooffset to 1 if W(x = 0, y = 0) and so on has also been measured
    # skip: number of configurations that are stripped from the beginnig
    # every: step between read lines
    # l: median blocking size for bootstrap
    WL <- readloopfilecfsub(file = file, skip = skip, Nytmax = Nt * fraction, x = x,
                    Nsmax = Ns * fraction, zerooffset = zerooffset, every = every, maxrows = maxrows)
    WL <- bootstrap.cf(WL, boot.R = bootsamples, boot.l = l)
    title <- sprintf("%s, %d configs\n used every %d, boot.l = %d\n", title,
                    (length(WL$cf[, 1]) + skip) * nsave, nsave * every, l)
   plot(WL, log = "y", xlab = "t/a_s", ylab = "C(t)",
            main = sprintf("%s, logscale", title))
    if (-1 %in% sign(WL$cf0)) {
      plot(WL, xlab = "t/a_s", ylab = "C(t)",
            main = sprintf("%s, linear scale", title),
            ylim = c(-WL$cf0[3], WL$cf0[3]))
      points(x = seq(-1, Nt * fraction + 1, 1),
            y = rep(0, times = length(seq(-1, Nt * fraction + 1, 1))),
                col = 2, type = "l")
    }
    return(WL)
}

calcplotWloopsmall <- function (filename, skip, Nsmax, Ntmax,
                                Nt, x, y, start, bootsamples,
                                title, nsave = 100, l = 2, maxrows = -1) {
    # plots loops for the non-integer potential
    # t is changed, x and y are kept constant
    # nsave: number of steps MC-sweeps between saving configurations
    # skip: number of configurations that are stripped from the beginnig
    # every: step between read lines
    # l: median blocking size for bootstrap
    # additionally plots lines signifying half the lattice extent
    WL <- readloopfilecfsmall(file = filename, skip = skip,
                    Nsmax = Nsmax, x = x, y = y, Ntmax = Ntmax, start = start, maxrows = maxrows)
    WL <- bootstrap.cf(WL, boot.R = bootsamples, boot.l = l)
    title <- sprintf("%s%d configs, skipped %d", title,
            (length(WL$cf[, 1]) + skip - 1) * nsave, skip * nsave)
    try(plot(WL, log = "y", xlab = "t/a_t", ylab = "C(x)",
            main = sprintf("%s, logscale", title)))
    arrows((Nt - start) / 2, 10, (Nt - start) / 2, 1e-08,
            angle = 90, length = 0.1, code = 0, col = 2)
    if (-1 %in% sign(WL$cf0)) {
      try(plot(WL, xlab = "t/a_t", ylab = "C(x)",
            main = sprintf("%s, linear scale", title),
            ylim = c(-WL$cf0[3], WL$cf0[3])))
      points(x = seq(-1, Ntmax + 1, 1),
            y = rep(0, times = length(seq(-1, Ntmax + 1, 1))),
            col = 2, type = "l")
      arrows((Nt - start) / 2, 10, (Nt - start) / 2, -10,
            angle = 90, length = 0.1, code = 0, col = 2)
      message(WL$cf0, " ")
    }
    return(WL)
}

# the deteffmass-functions use the determined expectation values
# of the Wilson loops to calculate the ratio of the neighbouring loops,
# the estimator for the effective masses, and then use the
# provided boundaries t1 and t2 to fit the effective mass.

deteffmass <- function (WL, yt, potential, t1, t2, isspatial, type = "log") {
    # effective masses for the integer potentials
    # WL: bootstrapped data for the Wilson loops
    # yt: the effective mass is V(yt)
    # potential: dataframe used to collect results for all possible extents
    # t1, t2: lower and upper boundary for effective mass
    # isspatial: bool used for the summarising of results:
    # temporal or spatial potential
    # returns effective masses, including fit,
    # plus ylimits to be able to zoom in on regions of the plot with
    # larger loop extent = closer to plateau
    # the fit is only done if effective masses were determined,
    # if fit was successful, summary and limits are made, else
    # effective masses are recalculated
    ylim <- c(0, 0)
    WL.effmass <- bootstrap.effectivemass(WL, type = type)
    if (length(na.omit(as.vector(WL.effmass$effMass))) > 1) {
        WL.effmass <- try(fit.effectivemass(WL.effmass, t1 = t1,
                t2 = t2, useCov = FALSE), silent = FALSE)
        ## usecov=FALSE
        if (!inherits(WL.effmass, "try-error")) {
            newline <- data.frame(R = yt, m = WL.effmass$effmassfit$t0[1],
                    dm = WL.effmass$effmassfit$se[1], space = isspatial,
                    p = WL.effmass$effmassfit$Qval,
                    chi = WL.effmass$effmassfit$chisqr / WL.effmass$effmassfit$dof)
            potential <- rbind(potential, newline)
            regionplotzoom <- seq(ceiling(0.4 * length(WL.effmass$effMass)) + 1,
                    length(WL.effmass$effMass) - 1)
            if (length(regionplotzoom) > 0) {
            if (length(regionplotzoom) > 8) regionplotzoom <- regionplotzoom[1:8]
                # max <- min(max(na.omit(WL.effmass$effMass[regionplotzoom])), 2)
                max <- min(max(na.omit(WL.effmass$effMass[regionplotzoom])), 2)
                min <- max(min(na.omit(WL.effmass$effMass[regionplotzoom])), 0)
                ylim <- c(min, max)
                # print(ylim)
            }
        } else {
            WL.effmass <- bootstrap.effectivemass(WL, type = type)
        }
    }
    return(list(WL.effmass, ylim))
}

deteffmasssmall <- function (WL, x, y, t1, t2, type = "log") {
    # effective masses for the non-integer potentials
    # WL: bootstrapped data for the Wilson loops
    # t1, t2: lower and upper boundary for effective mass
    # x, y: the potential is calculated for sqrt(x^2 + y^2)
    # returns effective masses, including fit,
    # plus ylimits to be able to zoom in on regions of the plot with
    # larger loop extent = closer to plateau,
    # plus a line for a dataframe to summarise the result
    # the fit is only done if effective masses were determined,
    # if fit was successful, summary and limits are made, else
    # effective masses are recalculated
    ylim <- c(0, 0)
    listresults <- list(FALSE, FALSE, FALSE)
    newline <- data.frame(R = NA, m = NA, dm = NA)
    WL.effmass <- bootstrap.effectivemass(WL, type = type)
    if (length(na.omit(as.vector(WL.effmass$effMass))) > 1) {
        WL.effmass <- try(fit.effectivemass(WL.effmass, t1 = t1,
                    t2 = t2, useCov = FALSE), silent = FALSE)
        listresults <- list(WL.effmass, x, y)
        if (inherits(WL.effmass, "try-error")) {
            #have to recalculate so the correct effective masses are stored instead of error
            WL.effmass <- bootstrap.effectivemass(WL, type = type)
            listresults <- list(FALSE, FALSE, FALSE)
            newline <- data.frame(R = sqrt(x^2 + y^2), m = NA, dm = NA)
        } else {
            regionplotzoom <- seq(floor(length(WL.effmass$effMass)) * 0.25,
                    min(floor(length(WL.effmass$effMass)) * 0.25, length(WL.effmass$effMass)))
            if (length(regionplotzoom > 0)) {
            if (length(regionplotzoom) > 8) regionplotzoom <- regionplotzoom[1:8]
                max <- min(max(na.omit(WL.effmass$effMass[regionplotzoom])), 2)
                min <- max(min(na.omit(WL.effmass$effMass[regionplotzoom])), 0)
                ylim <- c(min, max)
            }
            newline <- data.frame(R = sqrt(x^2 + y^2),
                        m = WL.effmass$effmassfit$t0[1],
                        dm = WL.effmass$effmassfit$se[1])
        }
    }
    return(list(newline, ylim, listresults))
}

# several plotting functions,
# TODO: replace everywhere with plotwitherror from hadron package

plotyerr <- function (x, y, dy, arlength = 0.02, log = "", ...) {
    xlo <- x
    xhi <- x
    ylo <- y - dy
    yhi <- y + dy
    dots <- list(...)
    if (!"ylim" %in% names(dots)) {
        ylim <- c(min(na.omit(ylo)), max(na.omit(yhi)))
        plot(x = x, y = y, ylim = ylim, log = log, ...)
    } else {
        plot(x = x, y = y, log = log, ...)
    }
    arrows(xhi, yhi, xlo, ylo, angle = 90, length = arlength, code = 3, ...)
}

pointsyerr <- function (x, y, dy, arlength = 0.02, ...) {
    xlo <- x
    xhi <- x
    ylo <- y - dy
    yhi <- y + dy
    points(x = x, y = y, ...)
    arrows(xhi, yhi, xlo, ylo, angle = 90, length = arlength, code = 3, ...)
}

plotxyerr <- function (x, y, dx, dy, arlength = 0.02, ...) {
    xlo <- x - dx
    xhi <- x + dx
    ylo <- y
    yhi <- y
    plot(x = x, y = y, ...)
    arrows(xhi, yhi, xlo, ylo, angle = 90, length = arlength, code = 0, ...)
    xlo <- x
    xhi <- x
    ylo <- y - dy
    yhi <- y + dy
    arrows(xhi, yhi, xlo, ylo, angle = 90, length = arlength, code = 0, ...)
}

pointsxyerr <- function (x, y, dx, dy, arlength = 0.02, ...) {
    xlo <- x - dx
    xhi <- x + dx
    ylo <- y
    yhi <- y
    points(x = x, y = y, ...)
    arrows(xhi, yhi, xlo, ylo, angle = 90, length = arlength, code = 0, ...)
    xlo <- x
    xhi <- x
    ylo <- y - dy
    yhi <- y + dy
    arrows(xhi, yhi, xlo, ylo, angle = 90, length = arlength, code = 0, ...)
}

errorpolygon <- function (X, fitresult, col.p, col.band = "gray",
                        polygon = TRUE, arlength = 0.1, pch = 1, ...) {
    # like plot of bootstrapfit, but with rep=TRUE
    prediction <- predict(fitresult, X)
    if (missing(pch)) {
        p.pch <- col.p
    } else {
        p.pch <- pch
    }
    if (fitresult$errormodel == "yerrors") {
      limits <- pointsyerr(x = fitresult$x, y = fitresult$y,
            dy = fitresult$dy, col = col.p, pch = p.pch,
            arlength = arlength, ...)
    } else {
      limits <- pointsxyerr(x = fitresult$x, y = fitresult$y,
            dy = fitresult$dy, dx = fitresult$dx, col = col.p,
            pch = p.pch, arlength = arlength, ...)
    }
   ylim <- limits$ylim

    if (polygon) {
        polyval <- c(prediction$val + prediction$err,
                rev(prediction$val - prediction$err))
        if (any(polyval < ylim[1]) || any(polyval > ylim[2])) {
          polyval[polyval < ylim[1]] <- ylim[1]
          polyval[polyval > ylim[2]] <- ylim[2]
        }
        col.band <- "gray"
        pcol <- col2rgb(col.band, alpha = TRUE) / 255
        pcol[4] <- 0.65
        pcol <- rgb(red = pcol[1], green = pcol[2],
                    blue = pcol[3], alpha = pcol[4])
        polygon(x = c(X, rev(X)), y = polyval, col = pcol,
                lty = 0, lwd = 0.001, border = pcol)
    }
    lines(x = X, y = prediction$val, col = col.p, ...)
}

points.effectivemass <- function (x, ..., ref.value, col, col.fitline) {
    # plot.effectivemass with rep=TRUE
  effMass <- x
  if (missing(col)) {
    col <- c("black", rainbow(n = (effMass$nrObs - 1)))
  }
  if (missing(col.fitline)) {
    col.fitline <- col[1]
  }
  # BaKo: is this also valid for acosh type effective masses?
  t <- effMass$t.idx
  suppressWarnings(pointsyerr(x = t - 1, y = effMass$effMass[t],
                dy = effMass$deffMass[t], col = col[1], ...))
  if (effMass$nrObs > 1) {
    for (i in 1:(effMass$nrObs - 1)) {
      suppressWarnings(pointsyerr(x = t - 1, y = effMass$t0[t + i * length(t)],
            dy = effMass$se[t + i * length(t)], rep = TRUE,
            col = col[i + 1], ...))
    }
  }
  if (!missing(ref.value)) {
    abline(h = ref.value, col = c("darkgreen"), lwd = c(3))
  }
  if (!is.null(effMass$effmassfit)) {
    lines(x = c(effMass$t1, effMass$t2),
          y = c(effMass$effmassfit$t0[1], effMass$effmassfit$t0[1]),
          col = col.fitline,
          lwd = 1.3)
      pcol <- col2rgb(col.fitline, alpha = TRUE) / 255
      pcol[4] <- 0.65
      pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
      rect(xleft = effMass$t1, xright = effMass$t2,
            ybottom = effMass$effmassfit$t0[1] - effMass$effmassfit$se[1],
            ytop = effMass$effmassfit$t0[1] + effMass$effmassfit$se[1],
            col = pcol, border = NA)
  }
}

determinerzero <- function (fit.result, bootsamples, c = -1.65, xi = 1) {
    # Determine r0 as solution of equation -r^2 d / dr V(r) = c,
    # solve for each bootstrapsample, r0 = mean pm sd
    # V(r) = a + sigma * r + b * ln(r)
    # -d/dr V(r) = -sigma -b / r
    # c = -sigma * r^2 - b * r, use p-q-formula for p = b / sigma, q = c / sigma
    # if a_t V is determined, c has to be rescaled so effectively a_s V is used
    # return mean, sd, bootstrapsamples
    c <- c * xi
    rzerolist <- c()
    for (bs in seq(1, bootsamples)) {
    #~     a <- fit.resultscaled$t[bs, 1]
        sigma <- fit.result$t[bs, 2]
        b <- fit.result$t[bs, 3]
        r1 <- -b / (2 * sigma) + sqrt((b / (2 * sigma))^2 - c / sigma)
        rzerolist[bs] <- r1
    }
    rzerolist <- na.omit(rzerolist)
    rzero <- mean(rzerolist)
    drzero <- sd(rzerolist)
    return(list(rzero, drzero, rzerolist))
}

drawallticks <- function (all = FALSE, inward = FALSE) {
    # want to draw tick labels on every side, but do not always want tick labels
    if (!inward) {
        tick <- -0.02
    } else {
        tick <- 0.02
    }
    axis(1, labels = TRUE, tck = tick)
    axis(2, labels = TRUE, tck = tick)
    axis(3, labels = all, tck = tick)
    axis(4, labels = all, tck = tick)
}

find_bootsamples <- function (Ns, Nt, beta, xi, plotpath = "",
                            smearing = FALSE, nape = 0, alpha = 1) {
    # for a configuration given by the parameters,
    # determine if there are any results for the non-integer potential.
    # if yes, determine how many bootstrap samples there are.
    basestring <- sprintf("%slistmeffsmallNs%dNt%dbeta%fxi%fbsamples*",
                    plotpath, Ns, Nt, beta, xi)
    if (smearing) {
        basestring <- sprintf("%slistmeffsmallNs%dNt%dbeta%fxi%fnape%dalpha%fbsamples*",
                plotpath, Ns, Nt, beta, xi, nape, alpha)
    }
    len <- nchar(basestring)
    files <- list.files(pattern = basestring)
    if (length(files) == 0) {
        stop("There are no data for this set of parameters")
    }
    numbers <- c()
    for (i in seq(1, length(files))) {
        numbers[i] <- substring(files[i], len, len + 4)
    }
    numbers <- suppressWarnings(as.integer(numbers))
    # check if number are of the right format, else strip characters until this is true
    for (i in seq(1, length(numbers))) {
        for (j in seq(3, 1)) {
            if (!is.na(numbers[i])) {
                break
            } else {
                numbers[i] <- suppressWarnings(as.integer(substring(files[i], len, len + j)))
            }
        }
    }
    print("simulations with the following amount of bootstrapsamples present:")
    print(numbers)
    return(max(numbers))
}

fitplotfunctions <- function (fun1, fun2, x1, y1, dy1, x2, y2, dy2,
                            par1, par2, title = "", boot.R = 500, bs = FALSE,
                            bsamples1 = array(), bsamples2 = array(),
                            verbose = FALSE, polygon = TRUE, ...) {
    # fits two different functions to two sets of data,
    # plots the functions together with limits fitting to both functions
    #generate bootstrapsamples, do the fit
    if (!bs) {
        bsamples1 <- parametric.bootstrap(boot.R, c(y1), c(dy1))
        bsamples2 <- parametric.bootstrap(boot.R, c(y2), c(dy2))
    }
    fit.result1 <- try(bootstrap.nlsfit(fun1, par1, y1, x1, bsamples1))
    fit.result2 <- try(bootstrap.nlsfit(fun2, par2, y2, x2, bsamples2))


    #determine proper limits
    dots <- list(...)
    my.ylim <- FALSE
    my.xlim <- FALSE
    if (!"ylim" %in% names(dots)) {
        ylimnew <- c(min(na.omit(y1 - dy1), na.omit(y2 - dy2)),
                max(na.omit(y1 + dy1), na.omit(y2 + dy2)))
        my.ylim <- TRUE
    }
    if (!"xlim" %in% names(dots)) {
        xlimnew <- c(min(na.omit(x1), na.omit(x2)),
                max(na.omit(x1), na.omit(x2)))
        my.xlim <- TRUE
        xmin <- min(xlimnew)
        xmax <- max(xlimnew)
    } else {
        xmin <- min(dots$xlim)
        xmax <- max(dots$xlim)
    }
    if (verbose) print(c(xmin, xmax))
    if (verbose) print("dy1")
    if (verbose) print(dy1)
    if (verbose) print("dy2")
    if (verbose) print(dy2)

    #points for plotting the result
    xvalues <- seq(from = xmin - 1, to = xmax + 1, length.out = 200)
    #plot everything together
    if (my.ylim && my.xlim) {
        plotyerr(x = x1, y = y1, dy = dy1, main = title,
                col = 1, pch = 1, xlim = xlimnew, ylim = ylimnew, ...)
        if (verbose) print("no limits given")
    }
    if (!my.ylim && my.xlim) {
        plotyerr(x = x1, y = y1, dy = dy1, main = title,
                col = 1, pch = 1, xlim = xlimnew, ...)
        if (verbose) print("y-limits given")
    }
    if (my.ylim && !my.xlim) {
        plotyerr(x = x1, y = y1, dy = dy1, main = title,
                col = 1, pch = 1, ylim = ylimnew, ...)
        if (verbose) print("x-limits given")
    }
    if (!my.ylim && !my.xlim) {
        plotyerr(x = x1, y = y1, dy = dy1, main = title, col = 1, pch = 1, ...)
        if (verbose) print("both limits given")
    }
    pointsyerr(x = x2, y = y2, dy = dy2, col = 2, pch = 2, ...)

    if (!inherits(fit.result1, "try-error")) {
        errorpolygon(xvalues, fit.result1, col.p = 1,
                    col.band = "gray", polygon = polygon)
        if (verbose) print("res1")
        if (verbose) print(names(fit.result1))
        if (verbose) print("bsamples")
        if (verbose) print(apply(bsamples1, 2, mean))
        if (verbose) print(apply(bsamples1, 2, sd))
        if (verbose) print("fit")
        if (verbose) print(fit.result1$x)
        if (verbose) print(fit.result1$y)
        if (verbose) print(fit.result1$dy)
    }
    if (!inherits(fit.result2, "try-error")) {
        errorpolygon(xvalues, fit.result2, col.p = 2,
                col.band = "red", polygon = polygon)
        if (verbose) print("res2")
        if (verbose) print("bsamples")
        if (verbose) print(apply(bsamples2, 2, mean))
        if (verbose) print(apply(bsamples2, 2, sd))
        if (verbose) print("fit")
        if (verbose) print(fit.result2$x)
        if (verbose) print(fit.result2$y)
        if (verbose) print(fit.result2$dy)
    }
    return(list(fit.result1, fit.result2))
}

readinbootstrapsamples <- function (beta, Ns, Nt, xi, bootsamples = 500,
                            path = "", names, columns, filename = "resultspot") {
    # reads in bootstrapsamples that were saved in a list to use
    # them for further analysis
    # the filename is constructed from the base filename and the parameters
    # the names and columns show which data are to be read in and how they are divided
    # bootsamples show how many samples there are
    # the name is constructed and the data are read in
    # an empty array for the results is made and filled
    # if only one column is available for a name, it is read directly from a list
    # if several columns are available, they are read from a sublist
    # returns an array containing all bootstrapsamples in the order of the given names
    if (length(columns) != length(names)) {
        stop("The columns per name and names do not have the same length! check input")
    }
    filename <- sprintf("%s%sNs%dNt%dbeta%fxi%fbs%d.RData",
                path, filename, Ns, Nt, beta, xi, bootsamples)
    listres <- readRDS(filename)
    bsamples <- array(c(rep(NA, sum(columns) * bootsamples)),
                dim = c(bootsamples, sum(columns)))
    index <- 1
    for (i in seq(1, length(names))) {
        name <- names[i]
        # print(length(listres[[name]]))
        for (j in seq(1, columns[i])) {
            if (columns[i] > 1) {
                bsamples[, index] <- listres[[name]][[j]]
            } else {
                bsamples[, index] <- listres[[name]]
            }
            index <- index + 1
        }
    }
    return(bsamples)
}

# handy for plotting
mmtoinches <- function (size) return(size / 25.4)

# determine the intercept between a linear function and a constant
#c = a*x + b
#c-b = a*x
#(c-b)/a = x
getintercept <- function(fitresult, rzeroone, bootsamples = 500) {
    intercepts <- c()
    for (bs in seq(1, bootsamples)){
        intercepts[bs] <- (rzeroone[bs] - fitresult$t[bs, 1]) / fitresult$t[bs, 2]
    }
    return(intercepts)
}

# determine the intercept between a linear function and a constant
#c = a*x + b
#c-b = a*x
#(c-b)/a = x
getinterceptfromparams <- function(bsa, bsb, rzeroone, bootsamples = 500) {
    intercepts <- c()
    for (bs in seq(1, bootsamples)){
        intercepts[bs] <- (rzeroone[bs] - bsb[bs]) / bsa[bs]
        # print(paste("a=", bsa[bs], "b=", bsb[bs], "rzero=", rzeroone[bs], "intercept=", intercepts[bs]))
    }
    return(intercepts)
}

# polynomials of first to fifth order, of a form that can be used by bootstrap.nlsfit
fnlin <- function (par, x, boot.r, ...) {
    return (par[1]  +  par[2] * x)
}
fnpar <- function (par, x, boot.r, ...) {
    return (par[1]  +  par[2] * x  +  par[3] * x^2)
}
fncub <- function (par, x, boot.r, ...) {
    return (par[1]  +  par[2] * x  +  par[3] * x^2  +  par[4] * x^3)
}
fnqar <- function (par, x, boot.r, ...) {
    return (par[1]  +  par[2] * x  +  par[3] * x^2  +  par[4] * x^3  +  par[5] * x^4)
}
fnqin <- function (par, x, boot.r, ...) {
    return (par[1]  +  par[2] * x  +  par[3] * x^2  +  par[4] * x^3  +  par[5] * x^4  +  par[6] * x^5)
}

# if a fit was done excluding NA values, but the order of bootstrapsamples is important,
# fill the result in another error, in which the same indices as in the original array are NA
# insert NA at thze given indices
fillexceptna <- function (indices, resultarray) {
    result <- resultarray
    for (i in seq_len(length(indices))) {
        result <- append(result, NA, after = indices - 1)
    }
    return (result)
}

printgitcommit <- function(pathtogit) {
    #get git commit hash of myfunctions.R, should be the same as of any script
    cwd <- setwd(pathtogit)
    githash <- try(system("git rev-parse --short HEAD", intern = TRUE))
    setwd(cwd)
    print(paste("## on git commit", githash))
    return(githash)
}

#' Predict values for bootstrapfit
#' based on hadron-function predict, but here the x-values differ for each bootstrapsample as well
#'
#' @param object Object of type bootstrapfit.
#' @param x Numeric vector with independent variable.
#' @param error Function to compute error from samples.
#' @param ... additional parameters to be passed on to the prediction function.
#'
#' @return
#' List with independent variable `x`, predicted central value `val`, error
#' estimate `err` and sample matrix `boot`.
#'
#' @export
#' @family NLS fit functions
predictwithxerror.bootstrapfit <- function (object, x, error = object$error.function, ...) {
  ## to include additional parameter to x$fn originally given as ... to
  ## bootstrap.nlsfit requires some pull-ups
  npar <- length(object$par.guess)
  val <- do.call(object$fn, c(list(par = object$t0[1:npar], x = mean(x), boot.r = 0), object$tofn))

  prediction <- list(x = x, val = val)
  stopifnot(length(x) == length(object$t[, 1]))

  if (!is.null(object$t)) {
    ## error band
    ## define a dummy function to be used in apply
    prediction_boot_fn <- function (boot.r) {
      par <- object$t[boot.r, 1:npar, drop = FALSE]
      do.call(object$fn, c(list(par = par, x = x[boot.r], boot.r = boot.r), object$tofn))
    }
    prediction_boot <- do.call(rbind, lapply(1:nrow(object$t), prediction_boot_fn))
    prediction$boot <- prediction_boot
  }

  return (prediction)
}

## routine to block a vector or a two dimensional array in blocks of length l
## discards anything that doesn't fit into exact multiples of the block length
## this is much faster than the blocking
## of tsboot
## directly copied from main repository
blockts <- function(data, l=2) {
  if(l == 1) {
    return(invisible(data))
  }
  if(is.vector(data)) {
    N <- floor(length(data)/l)*l
    return( apply(array(data, dim=c(l, N/l)), 2, mean))
  }
  if(length(dim(data))!=2) {
    stop("block.ts currently only implemented for vectors of 2-dim arrays\n")
  }
  N <- floor(length(data[,1])/l)*l
  ncf <- array(0, dim=c(N/l,length(data[1,])))
  j <- 1
  for ( i in seq(1,N,l)) {
    ncf[j,] <- apply(data[i:(i+l-1),], 2, mean)
    j <- j+1
  }
  return(invisible(ncf))
}

## define a function that uses the actual bootstrapsamples
## based on function from fit_xiyey
#' fit boostrap by boostrap (generated from y and dy)
#' x = matrix of size n \times N_pts
#' returns a list with bootstraps of parameters, their mean, stderr, chi^2, reduced chi^2
usebootstrap.fit_xiyey <-function(ansatz, x, y, ey, y_bts, guess, maxiter=10000, method="BFGS", N_bts = 500){
    N_pts <- length(y) # number of points
    N_par <- length(guess)
    N_dof <- N_pts - N_par

    par_bts <- matrix(data=NA, N_bts, N_par)
    ch2_bts <- matrix(data=NA, N_bts, N_par)

    ## First do the fit on the mean data:
    mini <- fit_xiyey(
                ansatz = ansatz, x=x, y=y, ey=ey,
                guess = guess,
                maxiter=maxiter, method=method)

    par_val <- mini[["par"]]

    ch2_val <- mini[["ch2"]]

    for (i in 1:N_bts){
        mini <- fit_xiyey(
            ansatz = ansatz, x=x, y=y_bts[i,], ey=ey,
            guess = guess,
            maxiter=maxiter, method=method)

        par_bts[i,] <- mini[["par"]]
    }

    par_sd <- apply(par_bts, 2, sd)

    res <- list(
        ansatz=ansatz,
        N_bts=N_bts, N_pts = N_pts, N_par = N_par, N_dof = N_dof,
        x=x,
        par=list(bts=par_bts, val=par_val, dval=par_sd),
        ch2=ch2_val,
        y = list(y = y, dy = ey, y_bts = y_bts), boot=TRUE
        )

    attr(res, "class") <- c("multifit", "list")
    return(res)
}


#' Predict values for multifit, based on predict.bootstrapfit
#'
#' @param object result of multifit.
#' @param x Numeric vector with independent variable.
#' @param error Function to compute error from samples.
#' @param ... additional parameters to be passed on to the prediction function.
#'
#' @return
#' List with independent variable `x`, predicted central value `val`, error
#' estimate `err` and sample matrix `boot`.
#'
#' @export
#' @family NLS fit functions
predict.anyfunction <- function (fn, par, par_bts, x, error = sd, ...) {
  ## to include additional parameter to x$fn originally given as ... to
  ## bootstrap.nlsfit requires some pull-ups
  if(length(par)!=length(par_bts[1, ])) stop("length of parameters and bootstrapsamples do not match")
  getresvector <- function(fn, par, x, ...) {
    res <- c()
    for (i in seq(1, length(x[1, ]))) {
        tmp <- fn(x = x[, i], par = par, ...)
        res[i] <- tmp
    }
    return(res)
  }
  npar <- length(par)
  val <- getresvector(fn = fn, par = par, x = x, ...)

  prediction <- list(x = x, val = val)

  if (!missing(par_bts)) {
    ## error band
    ## define a dummy function to be used in apply
    prediction_boot_fn <- function (boot.r) {
      par <- par_bts[boot.r, 1:npar, drop = FALSE]
      res <- getresvector(fn = fn, par = par, x = x, ...)
      return(res)
    }
    prediction_boot <- do.call(rbind, lapply(1:nrow(par_bts), prediction_boot_fn))
    prediction$boot <- prediction_boot

    err <- apply(prediction_boot, 2, error, na.rm = TRUE)
    stopifnot(length(err) == length(x[1, ]))
    prediction$err <- err
  }

  return (prediction)
}


## plot result of multifit with rep=TRUE
## plotindex: which row of matrix is used as the independent variable
errorpolygonmultifit <- function (X, fitresult, col.p, col.band = "gray",
                        polygon = TRUE, arlength = 0.1, pch = 1, plotindex = 1, mask, ylim, ...) {
    # like plot of bootstrapfit, but with rep=TRUE
    prediction <- predict.anyfunction(fn = fitresult$ansatz, par = fitresult$par$val, par_bts = fitresult$par$bts, X)
    if (missing(pch)) {
        p.pch <- col.p
    } else {
        p.pch <- pch
    }
    if(missing(mask)) mask <- seq(1:length(X[1, ]))
    limits <- plotwitherror(x = fitresult$x[plotindex, mask], y = fitresult$y$y[mask],
            dy = fitresult$y$dy[mask], col = col.p, pch = p.pch, rep=TRUE, ...)
    if(missing(ylim)) ylim <- limits$ylim

    if (polygon) {
        polyval <- c(prediction$val + prediction$err,
                rev(prediction$val - prediction$err))
        if (any(polyval < ylim[1]) || any(polyval > ylim[2])) {
          polyval[polyval < ylim[1]] <- ylim[1]
          polyval[polyval > ylim[2]] <- ylim[2]
        }
        col.band <- "gray"
        pcol <- col2rgb(col.band, alpha = TRUE) / 255
        pcol[4] <- 0.65
        pcol <- rgb(red = pcol[1], green = pcol[2],
                    blue = pcol[3], alpha = pcol[4])
        polygon(x = c(X[plotindex, ], rev(X[plotindex, ])), y = polyval, col = pcol,
                lty = 0, lwd = 0.001, border = pcol)
    }
    lines(x = X[plotindex, ], y = prediction$val, col = col.p, ...)
}


## TODO: make more pretty
## taken from summary.bootstrapnlsfit
summary.multifit <- function(object, ..., digits = 2, title="") {
    # print(object)
    y <- object$y$y
    dy <- object$y$dy
    npar <- object$N_par
    values <- object$par$val
    errors <- object$par$dval
    tmp <- apply(X=array(c(values, errors), dim=c(length(values), 2)), MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, with.cdot=FALSE)
    cat(title)
    cat("\nresult of multifit with bootstrap\n\n")
    if (object$boot) {
        bias <- values-apply(X=object$par$bts, MARGIN=2, FUN=mean, na.rm=TRUE)
        dim(bias) <- c(length(bias), 1)
        bias <- apply(X=bias, MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, with.cdot=FALSE)
        ci16 <- apply(X=object$par$bts, MARGIN=2, FUN=quantile, probs=c(0.16), drop=FALSE, na.rm=TRUE)
        dim(ci16) <- c(length(ci16), 1)
        ci16 <- apply(X=ci16, MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, with.cdot=FALSE)
        ci84 <- apply(X=object$par$bts, MARGIN=2, FUN=quantile, probs=c(0.84), drop=FALSE, na.rm=TRUE)
        dim(ci84) <- c(length(ci84), 1)
        ci84 <- apply(X=ci84, MARGIN=1, FUN=tex.catwitherror, digits=digits, with.dollar=FALSE, with.cdot=FALSE)
        cat("    best fit parameters with errors, bootstrap bias and 68% confidence interval\n\n")
        print(data.frame(par=tmp[1:npar], bias=bias[1:npar], ci16=ci16[1:npar], ci84=ci84[1:npar]))
    } else {
        cat("best fit parameters with error\n\n")
        print(data.frame(par=tmp[1:npar]))
    }
    cat("\n   chi^2 and fit quality\n")
    cat("chisqr / dof =", object$ch2, "/", object$N_dof, "=", object$ch2 / object$N_dof, "\n\n\n")


}