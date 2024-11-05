library("hadron")

mmtoinches <- function (size) return(size / 25.4)

pointsyerr <- function (x, y, dy, arlength = 0.02, ...) {
    try(print(pch))
    xlo <- x
    xhi <- x
    ylo <- y - dy
    yhi <- y + dy
    points(x = x, y = y, ...)
    arrows(xhi, yhi, xlo, ylo, angle = 90, length = arlength, code = 3, ...)
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
    # like plot of bootstrapfit, but with rep = TRUE
    prediction <- predict(fitresult, X)
    if (missing(pch)) {
        p.pch <- col.p
    } else {
        p.pch <- pch
    }
    if (fitresult$errormodel == "yerrors") {
      pointsyerr(x = fitresult$x, y = fitresult$y,
            dy = fitresult$dy, col = col.p, pch = p.pch,
            ...)
    } else {
      pointsxyerr(x = fitresult$x, y = fitresult$y,
            dy = fitresult$dy, dx = fitresult$dx, col = col.p,
            pch = p.pch, ...)
    }
   ylim <- c(fitresult$fn(fitresult$t0, X[1], 0), fitresult$fn(fitresult$t0, X[length(X)], 0))

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

drawcircle <- function(x, y, radius, points = 100, aspectratio = 1, ...){
    stopifnot(length(x) == length(y))
    stopifnot(length(x) > 0)
    # calculate circle points
    theta <- seq(0, 2*pi, length = points)
    # draw circle around each point in x, y
    for (i in 1:length(x)){
        lines(x = x[i] + radius*cos(theta), y = y[i] + aspectratio*radius*sin(theta), ...)
    }
}

# plot all results of the fits of r_0(beta) together
# if outputformat is none, it is assumed a resource for plotting is already opened
# betalist has to have all betas, also xi = 1
plotresultwithcircles <- function(file, fitlim = 0.3, betalist, radius = 0.005,
                                    outputname = "", outputformat = "pdf", title = "", 
                                    legendpos = "topleft", tikzsize = 400){
    stopifnot(outputformat == "pdf" || outputformat == "png" || outputformat == "none")
    # set graphical parameters
    cols <- c(1, 3, 4, 5, 6, 9, 10, 8, 7)
    pchs <- cols
    fontsize <- 2.6
    distance <- 3.5
    linewidth <- 2 * fontsize
    par(lwd = linewidth)
    # read in file
    data <- readRDS(file)
    # print(sum(data$nalist))
    data <- data$fitsrzero
    nxi <- length(data)
    stopifnot (nxi == length(betalist))
    # get points for determining plotlims
    xpoints <- c(unlist(sapply(data[2:nxi], getElement, "x")), betalist)
    xlim <- c(1.0*min(xpoints), max(xpoints))
    ylim <- c(data[[1]]$rzero - fitlim, data[[1]]$rzero + fitlim)
    aspectratio <- (ylim[2] - ylim[1]) / (xlim[2] - xlim[1])
    # open output file according to format
    if (outputformat == "pdf") {
        defaultmargin <- par(c("mai"))
        par(mai = c(defaultmargin[1] * max(1, 0.8 * fontsize),
            defaultmargin[2] * max(1, 1 * fontsize), 0.1, 0.1))
        packages <- c("\\usepackage{tikz}",
                "\\usepackage[active,tightpage]{preview}",
                "\\PreviewEnvironment{pgfpicture}",
                "\\setlength\\PreviewBorder{0pt}",
                "\\usepackage{amsmath}")
        tikzfile <- tikz.init(outputname, width = mmtoinches(tikzsize),
                    height = mmtoinches(tikzsize * 0.65), packages = packages)
    }
    if (outputformat == "png") png(paste0(outputname, ".png"), width = 500, height = 500*0.65, units = "mm", res = 200)
    # plot line for r_0(xi = 1)
    plot(NA, xlim = xlim, ylim = ylim, xlab = "$\\beta$", ylab = "$r_0 / a_s$", 
        cex.lab = fontsize, cex.axis = fontsize, cex = fontsize, lwd = linewidth, main = title)
    plotwitherror(x = betalist[1], y = data[[1]]$rzero, dy = data[[1]]$drzero, rep = TRUE, cex = fontsize, lwd = linewidth)
    lines(x = c(0.9, 1.1)*xlim, y = rep(data[[1]]$rzero, 2), lty = 2, cex = fontsize, lwd = linewidth)
    lines(x = c(0.9, 1.1)*xlim, y = rep(data[[1]]$rzero + data[[1]]$drzero, 2), lty = 3, cex = fontsize, lwd = linewidth)
    lines(x = c(0.9, 1.1)*xlim, y = rep(data[[1]]$rzero - data[[1]]$drzero, 2), lty = 3, cex = fontsize, lwd = linewidth)
    # plot lines for different xi
    # plot circles
    for (i in seq(2, length(betalist))) {
        errorpolygon(X = seq(0.9*xlim[1], 1.1*xlim[2], length = 1000), fitresult = data[[i]],
                col.p = cols[i], pch = pchs[i], col.band="grey", cex = fontsize, lwd = linewidth)
        drawcircle(x = betalist[i], y = data[[i]]$y[abs(data[[i]]$x - betalist[i]) < 1e-4][1], 
        radius = radius, aspectratio = aspectratio, col = cols[i], cex = fontsize, lwd = linewidth)
    }
    legend(x = legendpos, legend = c("1.00", "0.80", "0.67", "0.50", "0.40", "0.33", "0.25", "0.20", "0.17"),
            col = cols[1:9], pch = cols[1:9], title = "$\\xi_\\text{in} = $", 
            cex = fontsize, lwd = linewidth, ncol=2)
    if (outputformat == "pdf") {
        tikz.finalize(tikzfile, margins = "0")
        par(mai = defaultmargin)
    }
    if (outputformat == "png") dev.off()
}

# pdf("testplotall.pdf", width=mmtoinches(500), height = mmtoinches(500*0.65))
# plotresultwithcircles(file = "~/Documents/masterthesis/more_measurements/limitaic/contlim/listresultsrenormalizationnormalbeta1.650000omit0xiconstllxi1llr00fl0.30aicscaletauintetp0.RData", 
#     fitlim = 0.3, betalist = c(1.65, 1.6, 1.56, 1.49, 1.47, 1.44, 1.42), outputformat = "none")
# plotresultwithcircles(file = "~/Documents/masterthesis/more_measurements/limitaic/contlim/listresultsrenormalizationnormalbeta1.700000omit0xiconstllxi1llr00fl0.30aicscaletauintetp0.RData", 
#     fitlim = 0.3, betalist = c(1.70, 1.6, 1.56, 1.49, 1.47, 1.44, 1.42), outputformat = "none")
# plotresultwithcircles(file = "~/Documents/masterthesis/more_measurements/limitaic/contlim/listresultsrenormalizationsidewaysbeta1.650000omit0xiconstllxi1llr00fl0.30aicscaletauintetp0.RData", 
#     fitlim = 0.3, betalist = c(1.65, 1.6, 1.56, 1.49, 1.47, 1.44, 1.42), outputformat = "none")
# plotresultwithcircles(file = "~/Documents/masterthesis/more_measurements/limitaic/contlim/listresultsrenormalizationsidewaysbeta1.700000omit0xiconstllxi1llr00fl0.30aicscaletauintetp0.RData", 
#     fitlim = 0.3, betalist = c(1.70, 1.6, 1.56, 1.49, 1.47, 1.44, 1.42), outputformat = "none")
# dev.off()
# try(plotresultwithcircles(file = sprintf("~/Documents/masterthesis/more_measurements/limitaic/contlim/listresultsrenormalizationnormalbeta1.650000omit%dxiconstllxi1llr00fl0.30aicscaletauintetp%d.RData", 0, 0), 
#     fitlim = 0.3, betalist = c(1.65, 1.6, 1.56, 1.49, 1.48, 1.44, 1.42), outputformat = "pdf", 
#     title = "", outputname="normalb1.65omit0etp0", legendpos = "bottomright", tikzsize = 300))
       
if(TRUE){
stringmode <- FALSE
fitlim <- 0.3
if(stringmode) fitlim <- 0.03

if(!stringmode) pdf("plotallpointslist_tauint_more_stat_xi0.17.pdf", width=mmtoinches(500), height = mmtoinches(500*0.65), title="")
if(stringmode) pdf("plotallpointslist_strings.pdf", width=mmtoinches(500), height = mmtoinches(500*0.65), title="")
# sideways, beta=1.65
betalist <- list(
    c(1.65, 1.6, 1.56, 1.49, 1.47, 1.44, 1.42, 1.41, 1.4),
    c(1.65, 1.6, 1.56, 1.49, 1.47, 1.44, 1.42, 1.41, 1.4),
    c(1.65, 1.6, 1.56, 1.49, 1.47, 1.44, 1.42, 1.41, 1.4),
    c(1.65, 1.6, 1.56, 1.49, 1.47, 1.44, 1.42, 1.41, 1.4)
)
i <- 1
for (etp in c(0, 1)){
    for (omit in c(0, 1)) {
        # if(!stringmode) resfile <- sprintf("~/Documents/masterthesis/more_measurements/limitaic/contlim/listresultsrenormalizationsidewaysbeta1.650000omit%dxiconstllxi1llr00fl0.30aicscaletauintetp%d.RData", omit, etp)
        if(!stringmode) resfile <- sprintf("~/mnt_qbig/listresultsrenormalizationsidewaysbeta1.650000omit%dxiconstllxi1llr00fl0.30aicscaletauintetp%d.RData", omit, etp)
        if(stringmode) resfile <- sprintf("~/mnt_qbig/listresultsrenormalizationsidewaysstringbeta1.650000omit%dxiconstllxi1llr00fl0.03aicscaletauintetp%d.RData", omit, etp)
        try(plotresultwithcircles(file = resfile, fitlim = fitlim, betalist = betalist[[i]], outputformat = "none", title = paste("sideways beta 1.65 omit", omit, "total error used", etp)))
        i <- i+1
        # print(i)
    }
}
# normal, beta=1.65
betalist <- list(
    c(1.65, 1.6, 1.56, 1.49, 1.48, 1.44, 1.42, 1.41, 1.39),
    c(1.65, 1.6, 1.56, 1.49, 1.48, 1.44, 1.42, 1.41, 1.41),
    c(1.65, 1.6, 1.565, 1.49, 1.48, 1.45, 1.42, 1.41, 1.41),
    c(1.65, 1.6, 1.565, 1.49, 1.48, 1.45, 1.42, 1.41, 1.4)
)
i <- 1
for (etp in c(0, 1)){
    for (omit in c(0, 1)) {
        if(!stringmode) resfile <- sprintf("~/mnt_qbig/listresultsrenormalizationnormalbeta1.650000omit%dxiconstllxi1llr00fl0.30aicscaletauintetp%d.RData", omit, etp)
        if(stringmode) resfile <- sprintf("~/mnt_qbig/listresultsrenormalizationnormalstringbeta1.650000omit%dxiconstllxi1llr00fl0.03aicscaletauintetp%d.RData", omit, etp)
        try(plotresultwithcircles(file = resfile, fitlim = fitlim, betalist = betalist[[i]], outputformat = "none", title = paste("normal beta 1.65 omit", omit, "total error used", etp)))
        i <- i+1
        # print(i)
    }
}
# i <- 2
# print(par("cex"))
# if(!stringmode) resfile <- sprintf("~/mnt_qbig/listresultsrenormalizationnormalbeta1.650000omit%dxiconstllxi1llr00fl0.30aicscaletauintetp%d.RData", omit, etp)
# try(plotresultwithcircles(file = resfile, fitlim = fitlim, betalist = betalist[[i]], outputformat = "pdf", title = "", outputname="plotr0ofbetapapernormal1.65o1e0"))
        

# sideways, beta=1.70
betalist <- list(
    c(1.7, 1.6521, 1.6, 1.555, 1.525, 1.48, 1.49, 1.46, 1.47),
    c(1.7, 1.6500, 1.6, 1.54, 1.525, 1.48, 1.49, 1.46, 1.47),
    c(1.7, 1.6521, 1.595, 1.55, 1.520, 1.4814, 1.49, 1.46, 1.47),
    c(1.7, 1.6521, 1.6, 1.55, 1.520, 1.48, 1.49, 1.46, 1.47)
)
i <- 1
for (etp in c(0, 1)){
    for (omit in c(0, 1)) {
        if(!stringmode) resfile <- sprintf("~/mnt_qbig/listresultsrenormalizationsidewaysbeta1.700000omit%dxiconstllxi1llr00fl0.30aicscaletauintetp%d.RData", omit, etp)
        if(stringmode) resfile <- sprintf("~/mnt_qbig/listresultsrenormalizationsidewaysstringbeta1.700000omit%dxiconstllxi1llr00fl0.03aicscaletauintetp%d.RData", omit, etp)
        try(plotresultwithcircles(file = resfile, fitlim = fitlim, betalist = betalist[[i]], outputformat = "none", title = paste("sideways beta 1.70 omit", omit, "total error used", etp)))
        i <- i+1
        # print(i)
    }
}
# normal, beta=1.70
betalist <- list(
    c(1.7, 1.6521, 1.6075, 1.555, 1.525, 1.50, 1.49, 1.46, 1.45),
    c(1.7, 1.6521, 1.6075, 1.540, 1.525, 1.50, 1.49, 1.46, 1.45),
    c(1.7, 1.6400, 1.6075, 1.550, 1.525, 1.50, 1.49, 1.46, 1.45),
    c(1.7, 1.6400, 1.5950, 1.540, 1.520, 1.50, 1.49, 1.46, 1.45)
)
i <- 1
for (etp in c(0, 1)){
    for (omit in c(0, 1)) {
        if(!stringmode) resfile <- sprintf("~/mnt_qbig/listresultsrenormalizationnormalbeta1.700000omit%dxiconstllxi1llr00fl0.30aicscaletauintetp%d.RData", omit, etp)
        if(stringmode) resfile <- sprintf("~/mnt_qbig/listresultsrenormalizationnormalstringbeta1.700000omit%dxiconstllxi1llr00fl0.03aicscaletauintetp%d.RData", omit, etp)
        try(plotresultwithcircles(file = resfile, fitlim = fitlim, betalist = betalist[[i]], outputformat = "none", title = paste("normal beta 1.70 omit", omit, "total error used", etp)))
        # try(plotresultwithcircles(file = resfile, fitlim = fitlim, betalist = betalist[[i]], outputformat = "pdf", title = "", 
            # outputname=sprintf("plotr0ofbetapapernormal1.70o%de%d", omit, etp)))
        i <- i+1
        # print(i)
    }
}  
dev.off()
}
