library("hadron")
source("myfunctions.R")

# This script is used to determine the renormalised beta and
# the continuum limit after each ensemble has been analysed

# first, the data are read in, then the optimal beta, xi_ren and P are determined,
# and then the continuum limit is determined with different polynomial fits

# functions needed for the programm
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

# set whether the results of the sideways or the normal potential should be considered
sideways <- FALSE
omit <- 0

# set filenames, read in results, set up containers for bootstrapsamples
if (sideways) {
data <- read.table("resultsummary2p1dsidewaysb1.700Ns16.csv", header = TRUE, sep = " ")
filenameres <- "resultsrotated"
side <- 2
} else {
data <- read.table("resultsummary2p1dnormalb1.700Ns16.csv", header = TRUE, sep = " ")
filenameres <- "resultssubtracted"
side <- 2
}
data <- na.omit(data)
if (sideways) {
data <- data[data$omit == omit, ]
}

nom <- length(data$beta)

# read in bootstrapsamples
bootsamples <- 500
arrayrzero <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))
arrayp <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))
arrayxi <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))
intercepts <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))

for (i in seq(1, nom)) {
    string <- sprintf("i = %d, beta = %f, Ns = %d, Nt = %d, xi = %f",
                    i, data$beta[i], data$Ns[i], data$Nt[i], data$xi[i])
    print(string)
    result <- readinbootstrapsamples(beta = data$beta[i], Ns = data$Ns[i],
                    Nt = data$Nt[i], xi = data$xi[i], columns = c(1, 1, 1),
                    names = c("bsrzeros", "bsp", "bsxicalc"), filename = filenameres)
    arrayrzero[, i] <- result[, 1]
    arrayp[, i] <- result[, 2]
    arrayxi[, i] <- result[, 3]
}

# set up parameters for plotting, initialize output file
fontsize <- 2.6
distance <- 3.5
linewidth <- 2 * fontsize
path <- "tikz"
packages <- c("\\usepackage{tikz}",
                "\\usepackage[active,tightpage]{preview}",
                "\\PreviewEnvironment{pgfpicture}",
                "\\setlength\\PreviewBorder{0pt}",
                "\\usepackage{amsmath}")

for (size in c(0.65, 0.75, 1.33)) {


nameplot <- sprintf("%srenormr0ratio%.2f", path, size)
if (sideways) {
    nameplot <- sprintf("%srenormr0sidewaysratio%.2fomit%d", path, size, omit)
}

if (size <= 1) {
    tikzfile <- tikz.init(nameplot, width = mmtoinches(400),
                    height = mmtoinches(400 * size), packages = packages)
    }
if (size > 1) {
    tikzfile <- tikz.init(nameplot, width = mmtoinches(200),
                    height = mmtoinches(200 * size), packages = packages)
    }
defaultmargin <- par(c("mai"))
par(mai = c(defaultmargin[1] * max(1, 0.8 * fontsize),
            defaultmargin[2] * max(1, 0.8 * fontsize), 0.1, 0.1))
if (sideways) {
    par(mai = c(defaultmargin[1] * max(1, 0.8 * fontsize),
                0.01, 0.1, defaultmargin[2] * max(1, 0.8 * fontsize)))
}
print(nameplot)

# set up limits, masks for selecting the right points
ylim <- c(min(data$r0 - data$dr0), max(data$r0 + data$dr0))
if (size > 1) { ylim <- c(3, 5.5)}
mask <- data$beta == 1.7 & data$xi == 1 & data$c == -1.65
maskone <- mask
rzeroone <- data$r0[mask]
xlim <- c(min(data$beta), max(data$beta))
if (size > 1) { xlim <- c(1.45, 1.75)}

# set input anisotropies that were considered, container for results
xis <- c(1, 0.8, 2/3, 0.5, 0.4, 1/3, 0.25)
cols <- c(1, 3, 4, 5, 6, 9, 10, 8)
fitsrzero <- list(NULL)
fitsplaquette <- list(NULL)
fitsxi <- list(NULL)
intercepts <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
plaqren <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
xiphys <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))

legendtext <- c("1.000")
xvalues <- seq(1.4, 1.8, 0.01)

par(lwd = linewidth)

# plot result of xi_input = 1 as a line and as a point
plot(x = seq(1.4, 1.8, by = 0.01), y = rep(data$r0[mask], 41),
            cex.lab = fontsize, cex.axis = fontsize, cex = fontsize, ylab = "",
            xlab = "", type = "l", ylim = ylim, xlim = xlim)
points(x = seq(1.4, 1.8, by = 0.01),
        y = rep(data$r0[mask] - data$dr0[mask], 41), type = "l", lty = 2)
points(x = seq(1.4, 1.8, by = 0.01),
        y = rep(data$r0[mask] + data$dr0[mask], 41), type = "l", lty = 2)
points(x = seq(1.4, 1.8, by = 0.01),
        y = rep(data$r0[mask] - 0.6, 41), type = "l", lty = 3, col = cols[length(xis) + 1])
points(x = seq(1.4, 1.8, by = 0.01),
        y = rep(data$r0[mask] + 0.6, 41), type = "l", lty = 3, col = cols[length(xis) + 1])
plotwitherror(x = data$beta[maskone], y = data$r0[maskone],
        dy = data$dr0[maskone], col = 1, pch = 1, cex = fontsize, rep=TRUE)


intercepts[, 1] <- rep(1.7, bootsamples)
plaqren[, 1] <- arrayp[, mask]
xiphys[, 1] <- arrayxi[, mask]

interceptsimple <- c(1.7)

fitresults <- data.frame(xiin = NA, r0slope = NA, r0intercept = NA, chir0 = NA, pr0 = NA,
                        plaqslope = NA, plaqintercept = NA, chiplaq = NA, pplaq = NA,
                        xislope = NA, xiintercept = NA, chixi = NA, pxi = NA)

# for each xi input != 1:
# select the points, make linear fits to r_zero, plaquette, xi_ren as a function of beta
# plot r_zero(beta) including error regions
# r_zero was determined from the combination of the spatial and rescaled temporal potential
# determine beta_ren: r_zero(beta_ren, xi_input) = r_zero(xi_input=1)
# determine P(beta_ren), xi(beta_ren)
# put everything in one dataframe, nicely formatted for easy printing
for (i in seq(2, length(xis))) {
    mask <- abs(data$xi - xis[i]) < 0.01 & data$c == -1.65 & abs(data$r0 - data$r0[maskone]) < 0.6
    maskplot <- abs(data$xi - xis[i]) < 0.01 & data$c == -1.65
    fitsrzero[[i]] <- try(bootstrap.nlsfit(fnlin, c(1, 1), data$r0[mask],
                            data$beta[mask], arrayrzero[, mask]))
    fitsplaquette[[i]] <- try(bootstrap.nlsfit(fnlin, c(1, 1), data$p[mask],
                            data$beta[mask], arrayp[, mask]))
    fitsxi[[i]] <- try(bootstrap.nlsfit(fnlin, c(1, 1), data$xicalc[mask],
                            data$beta[mask], arrayxi[, mask]))
    xvaluesplot <- seq((data$r0[maskone] - 0.6 - fitsrzero[[i]]$t0[1]) / fitsrzero[[i]]$t0[2],
                        (data$r0[maskone] + 0.6 - fitsrzero[[i]]$t0[1]) / fitsrzero[[i]]$t0[2], by = 0.01)

    if (!inherits(fitsrzero[[i]], "try-error")) {
        try(errorpolygon(X = xvalues, fitsrzero[[i]], col.p = cols[i],
                col.band = cols[i], cex = fontsize, arlength = 0.05 * fontsize))
        try(plotwitherror(x = data$beta[maskplot], y = data$r0[maskplot],
                dy = data$dr0[maskplot], col = cols[i], pch = cols[i], cex = fontsize, rep = TRUE))
        interceptsimple[i] <- (rzeroone - fitsrzero[[i]]$t0[1]) / fitsrzero[[i]]$t0[2]
        intercepts[, i] <- getintercept(fitsrzero[[i]], arrayrzero[, maskone], bootsamples)
        prediction <- predict(fitsplaquette[[i]], intercepts[, i])
        plaqren[, i] <- prediction$val
        prediction <- predict(fitsxi[[i]], intercepts[, i])
        xiphys[, i] <- prediction$val
        newline <- data.frame(xiin = xis[i],
            r0slope = tex.catwitherror(fitsrzero[[i]]$t0[2], fitsrzero[[i]]$se[2], with.dollar = FALSE, digits = 2),
            r0intercept = tex.catwitherror(fitsrzero[[i]]$t0[1], fitsrzero[[i]]$se[1], with.dollar = FALSE, digits = 2),
            chir0 = fitsrzero[[i]]$chisqr / fitsrzero[[i]]$dof, pr0 = fitsrzero[[i]]$Qval,
            plaqslope = tex.catwitherror(fitsplaquette[[i]]$t0[2], fitsplaquette[[i]]$se[2], with.dollar = FALSE, digits = 2),
            plaqintercept = tex.catwitherror(fitsplaquette[[i]]$t0[1], fitsplaquette[[i]]$se[1], with.dollar = FALSE, digits = 2),
            chiplaq = fitsplaquette[[i]]$chisqr / fitsplaquette[[i]]$dof, pplaq = fitsplaquette[[i]]$Qval,
            xislope = tex.catwitherror(fitsxi[[i]]$t0[2], fitsxi[[i]]$se[2], with.dollar = FALSE, digits = 2),
            xiintercept = tex.catwitherror(fitsxi[[i]]$t0[1], fitsxi[[i]]$se[1], with.dollar = FALSE, digits = 2),
            chixi = fitsxi[[i]]$chisqr / fitsxi[[i]]$dof, pxi = fitsxi[[i]]$Qval)
        fitresults <- rbind(fitresults, newline)
    }

    legendtext[i] <- sprintf("%.3f", xis[i])
}

# make legend for plot
legendtext[length(xis) + 1] <- "bounds fit"

fitresults <- fitresults[-1, ]
legend(legend = legendtext, x = "topleft", title = "$\\xi_\\text{input} = $",
        col = c(cols), pch = c(cols), cex = fontsize)
mtext("$r_0 / a_s$", side = side, line = distance, cex = fontsize) #ylab
mtext("$\\beta$", side = 1, line = distance, cex = fontsize) #xlab
if (sideways) {
    if (size == 1.33) {
        legend(legend = "sideways", x = "bottomright", cex = fontsize, bty = "n")
    }
} else {
    if (size == 1.33) {
        legend(legend = "normal", x = "bottomright", cex = fontsize, bty = "n")
    }
}
box(bty = "o")

tikz.finalize(tikzfile)

# make a dataframe of results
result <- data.frame(xiin = xis, beta = apply(intercepts, 2, mean),
                    dbeta = apply(intercepts, 2, sd),
                    xiphys = apply(xiphys, 2, mean),
                    dxiphys = apply(xiphys, 2, sd),
                    p = apply(plaqren, 2, mean), dp = apply(plaqren, 2, sd))
result$xiphys[result$xiin == 1] <- data$xicalc[data$xi == 1]
result$dxiphys[result$xiin == 1] <- data$dxicalc[data$xi == 1]
result <- cbind(result, data.frame(betasimple = interceptsimple))
par(mai = defaultmargin)
}

# add more data to results for easy printing
newframe <- data.frame(xiopt = tex.catwitherror(result$xiphys[seq(2, length(xis))], result$dxiphys[seq(2, length(xis))], digits = 2, with.dollar = FALSE),
betaopt = tex.catwitherror(result$beta[seq(2, length(xis))], result$dbeta[seq(2, length(xis))], digits = 2, with.dollar = FALSE),
popt = tex.catwitherror(result$p[seq(2, length(xis))], result$dp[seq(2, length(xis))], digits = 2, with.dollar = FALSE))

fitresults <- cbind(fitresults, newframe)

#determine P(0) renormalized
# first: fit with both xi and beta renormalized to a cubic polynomial

bsamplescontlimit <- array(rep(NA, bootsamples * (2 * length(xis))), dim = c(bootsamples, 2 * length(xis)))
bsamplescontlimit[, seq(1, length(xis))] <- plaqren
bsamplescontlimit[, seq(length(xis) + 1, 2 * length(xis))] <- xiphys^2
fitplaq <- try(bootstrap.nlsfit(fncub, c(0.7, 1, 1, 1), x = result$xiphys^2, y = result$p, bsamples = bsamplescontlimit))

# plot results to pdf
if (sideways) pdf(sprintf("tikzplotallfitssidewaysomit%d.pdf", omit), title = "")
if (!sideways) pdf("tikzplotallfits.pdf", title = "")

# result of all linear fits
for (i in seq(1, length(xis))){
    try(plot(fitsrzero[[i]], main = sprintf("rzero, xiin = %f, chi = %f, p = %f", xis[i], fitsrzero[[i]]$chi / fitsrzero[[i]]$dof, fitsrzero[[i]]$Qval)))
    try(plot(fitsplaquette[[i]], main = sprintf("P, xiin = %f, chi = %f, p = %f", xis[i], fitsplaquette[[i]]$chi / fitsplaquette[[i]]$dof, fitsplaquette[[i]]$Qval)))
    try(plot(fitsxi[[i]], main = sprintf("xi, xiin = %f, chi = %f, p = %f", xis[i], fitsxi[[i]]$chi / fitsxi[[i]]$dof, fitsxi[[i]]$Qval)))
}
# cubic continuum limit
try(plot(fitplaq, main = sprintf("continuum limit plaquette: %f + /-%f, chi = %f, p = %f", fitplaq$t0[1], fitplaq$se[1], fitplaq$chi / fitplaq$dof, fitplaq$Qval), plot.range = c(-0.2, 1.2), ylim = c(0.98 * (fitplaq$t0[1] - fitplaq$se[1]), 1.02 * max(result$p)), xaxs = "i", xlim = c(0, max(xiphys^2))))
try(plotwitherror(x = c(0), y = c(fitplaq$t0[1]), dy = (fitplaq$se[1]), col = 2, pch = 2, rep = TRUE))



if (TRUE) {
# fits also with xi_input with only one measurement, empty containers for result
xis <- c(1, 0.8, 2/3, 0.5, 0.4, 1/3, 2/7, 0.25)
bsamplescontlimitnaive <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
bsamplescontlimitnaivexiren <- array(rep(NA, 2 * bootsamples * (length(xis))), dim = c(bootsamples, 2 * length(xis)))
bsamplescontlimitbeta <- array(rep(NA, 2 * bootsamples * (length(xis))), dim = c(bootsamples, 2 * length(xis)))
pnaive <- c()
xirennaive <- c()

# fill up containers: naive indicates no renormalized values (of beta) are used for fit
# naive: no renormalization at all
# naivexiren: use renormalized xi, but leave beta fix
for (i in seq(1, length(xis))) {
    row <- data$beta == 1.7 & abs(data$xi - xis[i]) < 0.01
    pnaive[i] <- data$p[row]
    xirennaive[i] <- data$xicalc[row]
    bsamplescontlimitnaive[, i] <- arrayp[, row]
    bsamplescontlimitnaivexiren[, i] <- arrayp[, row]
    bsamplescontlimitnaivexiren[, length(xis) + i] <- arrayxi[, row]^2
    bsamplescontlimitbeta[, i] <- intercepts[, i]
    bsamplescontlimitbeta[, length(xis) + i] <- arrayxi[, row]^2
}

# cubic fits with nothing renormalized and only xi renormalized
fitplaqnaive <- try(bootstrap.nlsfit(fncub, c(0.71, 1, 1, 1),
            x = xis^2, y = pnaive, bsamples = bsamplescontlimitnaive))


fitplaqnaivexiren <- try(bootstrap.nlsfit(fncub, c(0.71, 1, 1, 1),
            x = xirennaive^2, y = pnaive, bsamples = bsamplescontlimitnaivexiren))


# save and plot results for naive limits
resultslist <- list(intercepts = intercepts, xiphys = xiphys, plaqren = plaqren,
                    fitsrzero = fitsrzero, fitsxi = fitsxi,
                    fitsp = fitsplaquette, fitplaq = fitplaq,
                    fitplaqnaive = fitplaqnaive, fitplaqnaivexiren = fitplaqnaivexiren)



defaultusr <- par("usr")
par("usr" = c(0, 1, defaultusr[3], defaultusr[4]))
try(plot(fitplaqnaive, plot.range = c(-0.2, 1.2),
        main = sprintf("continuum limit plaquette: %f + /-%f, chi = %f",
        fitplaqnaive$t0[1], fitplaqnaive$se[1], fitplaqnaive$chi / fitplaqnaive$dof),
        ylim = c((fitplaqnaive$t0[1] - fitplaqnaive$se[1]), max(result$p)),
        xlab = "", xaxs = "i", xlim = c(0, 1)))
try(plotwitherror(x = c(0), y = c(fitplaqnaive$t0[1]),
        dy = (fitplaqnaive$se[1]), col = 2, pch = 2, rep = TRUE))
par("usr" = defaultusr)

defaultusr <- par("usr")
par("usr" = c(0, 1, defaultusr[3], defaultusr[4]))
try(plot(fitplaqnaivexiren, plot.range = c(-0.2, 1.2),
        main = sprintf("continuum limit plaquette: %f + /-%f, chi = %f",
        fitplaqnaivexiren$t0[1], fitplaqnaivexiren$se[1],
        fitplaqnaivexiren$chi / fitplaqnaivexiren$dof),
        ylim = c((fitplaqnaivexiren$t0[1] - fitplaqnaivexiren$se[1]), max(result$p)),
        xlab = "", xaxs = "i", xlim = c(0, 1)))
try(mtext(fitplaqnaivexiren$x, side = 1))
print(fitplaqnaivexiren$x)
try(plotwitherror(x = c(0), y = c(fitplaqnaivexiren$t0[1]),
        dy = (fitplaqnaivexiren$se[1]), col = 2, pch = 2, rep = TRUE))
par("usr" = defaultusr)

# for each polynomial:
# fitplaqnaive: xi=xi_input, p=p_meas -> neither beta nor xi reenormalized
# fitplaqnaivexiren: xi=xi_ren, p=p_meas -> beta not renormalized
# fitplaq: xi=xi_ren, p=p_interpolated -> xi and beta renormalized
# for each: do fit to continuum limit, plot, add results to list and table
# list: region 1-5 fitplaqnaive, region 6-10 fitplaqnaivexiren, region 11-15 fitplaq

fitspolynomial <- list()
resultspolynomial <- data.frame(degree = NA, lim = NA, chi = NA,
                    p = NA, type = NA, limplot = NA, dlimplot = NA)
i <- 1
for (fun in c(fnlin, fnpar, fncub, fnqar, fnqin)){
    fitplaqnaive <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
                x = xis^2, y = pnaive, bsamples = bsamplescontlimitnaive))
    fitspolynomial[[i]] <- fitplaqnaive
    print(summary(fitplaqnaive))
    plot(fitplaqnaive, main = sprintf("continuum limit plaquette: %f + /-%f, chi = %f, p = %f,\ndegree of polynomial:%d",
            fitplaqnaive$t0[1], fitplaqnaive$se[1],
            fitplaqnaive$chi / fitplaqnaive$dof, fitplaqnaive$Qval, i),
            plot.range = c(-0.2, 1.2), xlab = "", xaxs = "i", xlim = c(0, 1),
            ylim = c((fitplaqnaive$t0[1] - fitplaqnaive$se[1]), max(result$p)))
    resultspolynomial <- rbind(resultspolynomial,
            data.frame(degree = i, lim = tex.catwitherror(fitplaqnaive$t0[1],
            fitplaqnaive$se[1], digits = 2, with.dollar = FALSE),
            chi = fitplaqnaive$chi / fitplaqnaive$dof, p = fitplaqnaive$Qval,
            type = "naive", limplot = fitplaqnaive$t0[1],
            dlimplot = fitplaqnaive$se[1]))

    fitplaqnaivexiren <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
            x = xirennaive^2, y = pnaive,
            bsamples = bsamplescontlimitnaivexiren))
    fitspolynomial[[5 + i]] <- fitplaqnaivexiren
    plot(fitplaqnaivexiren, main = sprintf("continuum limit plaquette: %f + /-%f, chi = %f, p = %f,\ndegree of polynomial:%d",
            fitplaqnaivexiren$t0[1], fitplaqnaivexiren$se[1],
            fitplaqnaivexiren$chi / fitplaqnaivexiren$dof, fitplaqnaivexiren$Qval, i),
            plot.range = c(-0.2, 1.2), xlab = "", xaxs = "i", xlim = c(0, ),
            ylim = c((fitplaqnaivexiren$t0[1] - fitplaqnaivexiren$se[1]), max(result$p)))
    resultspolynomial <- rbind(resultspolynomial,
            data.frame(degree = i, lim = tex.catwitherror(fitplaqnaivexiren$t0[1],
            fitplaqnaivexiren$se[1], digits = 2, with.dollar = FALSE),
            chi = fitplaqnaivexiren$chi / fitplaqnaivexiren$dof,
            p = fitplaqnaivexiren$Qval, type = "naivexiren",
            limplot = fitplaqnaivexiren$t0[1], dlimplot = fitplaqnaivexiren$se[1]))

    fitplaq <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
                x = result$xiphys^2, y = result$p, bsamples = bsamplescontlimit))
    fitspolynomial[[10 + i]] <- fitplaq
    plot(fitplaq, main = sprintf("continuum limit plaquette: %f + /-%f, chi = %f, p = %f,\ndegree of polynomial:%d",
            fitplaq$t0[1], fitplaq$se[1],
            fitplaq$chi / fitplaq$dof, fitplaq$Qval, i),
            plot.range = c(-0.2, 1.2),
            ylim = c((fitplaq$t0[1] - fitplaq$se[1]), max(result$p)),
            xlab = "", xaxs = "i", xlim = c(0, 1))
    resultspolynomial <- rbind(resultspolynomial,
            data.frame(degree = i, lim = tex.catwitherror(fitplaq$t0[1],
            fitplaq$se[1], digits = 2, with.dollar = FALSE),
            chi = fitplaq$chi / fitplaq$dof, p = fitplaq$Qval,
            type = "plaq", limplot = fitplaq$t0[1], dlimplot = fitplaq$se[1]))

    fitbeta <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
                x = result$xiphys^2, y = result$beta, bsamples = bsamplescontlimitbeta))
    fitspolynomial[[15 + i]] <- fitbeta
    plot(fitplaq, main = sprintf("continuum limit beta: %f + /-%f, chi = %f, p = %f,\ndegree of polynomial:%d",
            fitbeta$t0[1], fitbeta$se[1],
            fitbeta$chi / fitbeta$dof, fitbeta$Qval, i),
            plot.range = c(-0.2, 1.2),
            ylim = c((fitbeta$t0[1] - fitbeta$se[1]), max(result$p)),
            xlab = "", xaxs = "i", xlim = c(0, 1))
    resultspolynomial <- rbind(resultspolynomial,
            data.frame(degree = i, lim = tex.catwitherror(fitbeta$t0[1],
            fitbeta$se[1], digits = 2, with.dollar = FALSE),
            chi = fitbeta$chi / fitbeta$dof, p = fitbeta$Qval,
            type = "beta", limplot = fitbeta$t0[1], dlimplot = fitbeta$se[1]))
    i <- i + 1
}
resultspolynomial <- resultspolynomial[-1, ]
namepol <- "plotstikz/polynomialnormal.csv"
if (sideways) namepol <- "plotstikz/polynomialsideways.csv"
# write out result
write.table(resultspolynomial, namepol, col.names = TRUE, row.names = FALSE)

}

# write out results
# contains:
# result: data frame with xiin, beta, dbeta(renormalized), xiphys, dxiphys, p, dp, betasimple (renormed beta from intercept from mean values)
# fitresult: data frame with xiin, r0slope, r0intercept,  chir0, pr0,
#                        plaqslope, plaqintercept, chiplaq, pplaq,
#                        xislope, xiintercept, chixi, pxi,
#                       xiopt, betaopt, pop
# written out with tex.catwitherror
# resultspolynomial: dataframe with degree of polynomial, result of cont. lim. with catwitherror, chi and p of lim, type of lim, result of lim as two doubles
# resultlist: intercepts = beta_ren, xiphys, plaqren , fitsrzero, fitsxi, fitsp (results of linear interpolations), fitplaq, fitplaqnaive, fitplaqnaivexiren(cubic polynomial cont limits)
# fitspolynomial: bootstrapnlsfit results of cont limit, region 1-5 fitplaqnaive, region 6-10 fitplaqnaivexiren, region 11-15 fitplaq
print(fitresults)
print(result)
if (!sideways) {
write.table(result, "plotstikz/resultsrenormalization.csv",
            col.names = TRUE, row.names = FALSE, append = FALSE)
write.table(fitresults, "plotstikz/fitresultsrenormalization.csv",
            col.names = TRUE, row.names = FALSE, append = FALSE)
saveRDS(resultslist, "plotstikz/listresultsrenormalization.RData")
saveRDS(fitspolynomial, "plotstikz/listpolynomialrenormalization.RData")
}

if (sideways) {
write.table(result, sprintf("plotstikz/resultsrenormalizationsidewaysomit%d.csv", omit),
            col.names = TRUE, row.names = FALSE, append = FALSE)
write.table(fitresults, sprintf("plotstikz/fitresultsrenormalizationsidewaysomit%d.csv", omit),
            col.names = TRUE, row.names = FALSE, append = FALSE)
saveRDS(resultslist, sprintf("plotstikz/listresultsrenormalizationsidewaysomit%d.RData", omit))
saveRDS(fitspolynomial, sprintf("plotstikz/listpolynomialrenormalizationsidewaysomit%d.RData", omit))
}
# move all plots into subfolder
system("mv -v tikz* plotstikz/")
