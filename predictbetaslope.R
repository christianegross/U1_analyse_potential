library("hadron")
# source("myfunctions.R")
# args <- commandArgs(trailingOnly = TRUE)
# source(args[1])


library(optparse)


if (TRUE) {
option_list <- list(
    make_option(c("-s", "--bootsamples"), type = "integer", default = 500,
    help = "how many bootstrapsamples should be drawn [default %default]"),
    make_option(c("-b", "--betaone"), type = "double", default = 1.7,
    help = "beta at xi=1 [default %default]"),
    make_option(c("-L", "--length"), type = "integer", default = 16,
    help = "spatial extent of lattice at xi=1 [default %default]"),

    make_option(c("-T", "--timeextent"), type = "integer", default = 16,
    help = "time extent of lattice at xi=1 [default %default]"),
    make_option(c("--myfunctions"), type = "character",
        default = "/hiskp4/gross/masterthesis/su2/build/debug/analysisscripts/myfunctions.R",
#~     make_option(c("--myfunctions"), type = "character", default = "myfunctions.R",
    help = "path to where additional functions are stored [default %default]"),
    make_option(c("--respath"), type = "character", default = "plotstikz/",
    help = "path to where the resulting plots and data are stored [default %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
}

# This script is used to determine the renormalised beta and
# the continuum limit after each ensemble has been analysed
# The renormalisation is done with the slope of the lowest two
# potential points, V(1a_s) and V(sqrt(2)a_s), scaled to units
# of a_s with the input anisotropy

# first, the data are read in, then the optimal beta, xi_ren and P are determined,
# and then the continuum limit is determined with different polynomial fits

# functions needed for the programm
# determine the intercept between a linear function and a constant
#c = a*x + b
#c-b = a*x
#(c-b)/a = x
getintercept <- function(fitresult, slopeone, bootsamples = 500) {
    intercepts <- c()
    for (bs in seq(1, bootsamples)){
        intercepts[bs] <- (slopeone[bs] - fitresult$t[bs, 1]) / fitresult$t[bs, 2]
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


# set filenames, read in results, set up containers for bootstrapsamples
data <- read.table(sprintf("summarysmallbetaone%f.csv", opt$b), header = TRUE, sep = " ")
datanormal <- read.table(sprintf("resultsummary2p1dnormalb%.3fNs%d.csv", opt$b, opt$L), header = TRUE, sep = " ")
filenameres <- "resultspot"
side <- 2
# data <- na.omit(data)
nom <- length(data$beta)
print(nom)

# read in bootstrapsamples, scale slope to a_s
bootsamples <- opt$bootsamples
arrayslope <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))
arrayp <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))
intercepts <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))


for (i in seq(1, nom)) {
    string <- sprintf("i = %d, beta = %f, Ns = %d, Nt = %d, xi = %f",
                    i, data$beta[i], data$Ns[i], data$Nt[i], data$xiin[i])
    print(string)
    correspondnormal <- abs(datanormal$beta - data$beta[i]) < 0.001 & abs(datanormal$xi - data$xiin[i]) < 0.001
    print(datanormal$xi[correspondnormal])
    print(datanormal$xicalc[correspondnormal])
    result <- readinbootstrapsamples(beta = data$beta[i], Ns = data$Ns[i],
                    Nt = data$Nt[i], xi = data$xiin[i], columns = c(1, 1),
                    names = c("bsslope", "bsp"), filename = filenameres)
    arrayp[, i] <- result[, 2]
    # arrayslope[, i] <- result[, 1] / data$xiin[i]
    # data$slopescale[i] <- data$slope[i] / data$xiin[i]
    # data$dslopescale[i] <- data$dslope[i] / data$xiin[i]
    arrayslope[, i] <- result[, 1] / datanormal$xicalc[correspondnormal]
    data$slopescale[i] <- data$slope[i] / datanormal$xicalc[correspondnormal]
    data$dslopescale[i] <- data$dslope[i] / datanormal$xicalc[correspondnormal]
}

# print(data$slope/data$xiin)

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

for (size in c(0.65)) {


nameplot <- sprintf("%srenormsloperatio%.2f", path, size)

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
# print(nameplot)

# set up limits, masks for selecting the right points
ylim <- c(min(data$slopescale - data$dslopescale), max(data$slopescale + data$dslopescale))
# if (size > 1) { ylim <- c(3, 5.5)}
mask <- data$beta == 1.7 & data$xiin == 1 & data$c == -1.65
maskone <- mask
slopeone <- data$slopescale[mask]
xlim <- c(min(data$beta), max(data$beta))
# if (size > 1) { xlim <- c(1.45, 1.75)}

# set input anisotropies that were considered, container for results
xis <- c(1, 0.8, 2/3, 0.5, 0.4, 1/3, 0.25)
cols <- c(1, 3, 4, 5, 6, 9, 10, 8) # colors for prettier plots
fitsslope <- list(NULL)
fitsplaquette <- list(NULL)
intercepts <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
plaqren <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))

legendtext <- c("1.000")
xvalues <- seq(1.4, 1.8, 0.01)

par(lwd = linewidth)

# plot result of xi_input = 1 as a line and as a point
plot(x = seq(1.4, 1.8, by = 0.01), y = rep(data$slopescale[mask], 41),
            cex.lab = fontsize, cex.axis = fontsize, cex = fontsize, ylab = "",
            xlab = "", type = "l", ylim = ylim, xlim = xlim)
points(x = seq(1.4, 1.8, by = 0.01),
        y = rep(data$slopescale[mask] - data$dslopescale[mask], 41), type = "l", lty = 2)
points(x = seq(1.4, 1.8, by = 0.01),
        y = rep(data$slopescale[mask] + data$dslopescale[mask], 41), type = "l", lty = 2)
points(x = seq(1.4, 1.8, by = 0.01),
        y = rep(data$slopescale[mask] - 0.6, 41), type = "l", lty = 3, col = cols[length(xis) + 1])
points(x = seq(1.4, 1.8, by = 0.01),
        y = rep(data$slopescale[mask] + 0.6, 41), type = "l", lty = 3, col = cols[length(xis) + 1])
plotwitherror(x = data$beta[maskone], y = data$slopescale[maskone],
        dy = data$dslopescale[maskone], col = 1, pch = 1, cex = fontsize, rep = TRUE)


intercepts[, 1] <- rep(opt$b, bootsamples)
plaqren[, 1] <- arrayp[, mask]

interceptsimple <- c(opt$b)

fitresults <- data.frame(xiin = NA, slopeslope = NA, slopeintercept = NA, chislope = NA, pslope = NA,
                        plaqslope = NA, plaqintercept = NA, chiplaq = NA, pplaq = NA)

# for each xi input != 1:
# select the points, make linear fits to slope, plaquette, xi_ren as a function of beta
# plot slope(beta) including error regions
# slope was determined from the combination of the spatial and rescaled temporal potential
# determine beta_ren: slope(beta_ren, xi_input) = slope(xi_input=1)
# determine P(beta_ren), xi(beta_ren)
# put everything in one dataframe, nicely formatted for easy printing
for (i in seq(2, length(xis))) {
    # print(i)
    mask <- abs(data$xiin - xis[i]) < 0.01 & data$c == -1.65 #& abs(data$slopescale - data$slopescale[maskone]) < 0.6
    maskplot <- abs(data$xiin - xis[i]) < 0.01 & data$c == -1.65
    fitsslope[[i]] <- try(bootstrap.nlsfit(fnlin, c(0.3, 0.3), data$slopescale[mask],
                            data$beta[mask], na.omit(arrayslope[, mask])))
    fitsplaquette[[i]] <- try(bootstrap.nlsfit(fnlin, c(1, 1), data$p[mask],
                            data$beta[mask], na.omit(arrayp[, mask])))
    # print(fitsslope[i])
    # print(fitsplaquette[i])
    xvaluesplot <- seq(1.4, 1.8, by = 0.01)
    # xvaluesplot <- seq((data$slopescale[maskone] - 0.6 - fitsslope[[i]]$t0[1]) / fitsslope[[i]]$t0[2],
    #                     (data$slopescale[maskone] + 0.6 - fitsslope[[i]]$t0[1]) / fitsslope[[i]]$t0[2], by = 0.01)

    if (!inherits(fitsslope[[i]], "try-error")) {
        try(errorpolygon(X = xvalues, fitsslope[[i]], col.p = cols[i],
                col.band = cols[i], cex = fontsize, arlength = 0.05 * fontsize))
        try(plotwitherror(x = data$beta[maskplot], y = data$slopescale[maskplot],
                dy = data$dslopescale[maskplot], col = cols[i], pch = cols[i], cex = fontsize, rep = TRUE))
        interceptsimple[i] <- (slopeone - fitsslope[[i]]$t0[1]) / fitsslope[[i]]$t0[2]
        interceptintermediate <- getintercept(fitsslope[[i]], arrayslope[, maskone], bootsamples = length(fitsslope[[i]]$t[, 1]))
        intercepts[, i] <- c(interceptintermediate, rep(NA, (500 - length(interceptintermediate))))
        prediction <- predict(fitsplaquette[[i]], intercepts[, i])
        plaqren[, i] <- prediction$val
        newline <- data.frame(xiin = xis[i],
            slopeslope = tex.catwitherror(fitsslope[[i]]$t0[2], fitsslope[[i]]$se[2], with.dollar = FALSE, digits = 2),
            slopeintercept = tex.catwitherror(fitsslope[[i]]$t0[1], fitsslope[[i]]$se[1], with.dollar = FALSE, digits = 2),
            chislope = fitsslope[[i]]$chisqr / fitsslope[[i]]$dof, pslope = fitsslope[[i]]$Qval,
            plaqslope = tex.catwitherror(fitsplaquette[[i]]$t0[2], fitsplaquette[[i]]$se[2], with.dollar = FALSE, digits = 2),
            plaqintercept = tex.catwitherror(fitsplaquette[[i]]$t0[1], fitsplaquette[[i]]$se[1], with.dollar = FALSE, digits = 2),
            chiplaq = fitsplaquette[[i]]$chisqr / fitsplaquette[[i]]$dof, pplaq = fitsplaquette[[i]]$Qval)
        fitresults <- rbind(fitresults, newline)
    } else {
        try(plotwitherror(x = data$beta[maskplot], y = data$slopescale[maskplot],
                dy = data$dslopescale[maskplot], col = cols[i], pch = cols[i], cex = fontsize, rep = TRUE))
        newline <- data.frame(xiin = xis[i], slopeslope = NA, slopeintercept = NA, chislope = NA, pslope = NA,
                        plaqslope = NA, plaqintercept = NA, chiplaq = NA, pplaq = NA)
        fitresults <- rbind(fitresults, newline)
        interceptsimple[i] <- NA
    }

    legendtext[i] <- sprintf("%.3f", xis[i])
}

# make legend for plot
legendtext[length(xis) + 1] <- "bounds fit"

fitresults <- fitresults[-1, ]
legend(legend = legendtext, x = "bottomleft", title = "$\\xi_\\text{input} = $",
        col = c(cols), pch = c(cols), cex = fontsize)
mtext("$slope / a_s$", side = side, line = distance, cex = fontsize) #ylab
mtext("$\\beta$", side = 1, line = distance, cex = fontsize) #xlab
box(bty = "o")

tikz.finalize(tikzfile)

# make a dataframe of results
result <- data.frame(xiin = xis, beta = apply(intercepts, 2, mean, na.rm = T),
                    dbeta = apply(intercepts, 2, sd, na.rm = T),
                    p = apply(plaqren, 2, mean, na.rm = T), dp = apply(plaqren, 2, sd, na.rm = T))
result <- cbind(result, data.frame(betasimple = interceptsimple))
par(mai = defaultmargin)
}
print(result)

# add more data to results for easy printing
# newframe <- data.frame(betaopt = tex.catwitherror(result$beta[seq(2, length(xis))], result$dbeta[seq(2, length(xis))], digits = 2, with.dollar = FALSE),
# popt = tex.catwitherror(result$p[seq(2, length(xis))], result$dp[seq(2, length(xis))], digits = 2, with.dollar = FALSE))

# fitresults <- cbind(fitresults, newframe)


#determine P(0) renormalized
# first: fit with beta renormalized to a cubic polynomial

bsamplescontlimit <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
bsamplescontlimit[, seq(1, length(xis))] <- plaqren
fitplaq <- try(bootstrap.nlsfit(fncub, c(0.7, 1, 1, 1), x = xis^2, y = result$p, bsamples = bsamplescontlimit))

# plot results to pdf
pdf("tikzplotallfitsnaive.pdf", title = "")

# result of all linear fits
for (i in seq(1, length(xis))) {
    try(plot(fitsslope[[i]], main = sprintf("slope, xiin = %f, chi = %f, p = %f", xis[i], fitsslope[[i]]$chi / fitsslope[[i]]$dof, fitsslope[[i]]$Qval)))
    try(plot(fitsplaquette[[i]], main = sprintf("P, xiin = %f, chi = %f, p = %f", xis[i], fitsplaquette[[i]]$chi / fitsplaquette[[i]]$dof, fitsplaquette[[i]]$Qval)))
}
# cubic continuum limit
try(plot(fitplaq, main = sprintf("continuum limit plaquette: %f+/-%f, chi = %f, p = %f", fitplaq$t0[1], fitplaq$se[1], fitplaq$chi / fitplaq$dof, fitplaq$Qval), plot.range = c(-0.2, 1.2), ylim = c(0.98 * (fitplaq$t0[1] - fitplaq$se[1]), 1.02 * max(result$p)), xaxs = "i", xlim = c(0, max(xis^2))))
try(plotwitherror(x = c(0), y = c(fitplaq$t0[1]), dy = (fitplaq$se[1]), col = 2, pch = 2, rep = TRUE))


if (TRUE) {
# fits also with xi_input with only one measurement, empty containers for result
xis <- c(1, 0.8, 2/3, 0.5, 0.4, 1/3, 2/7, 0.25)
bsamplescontlimitnaive <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))

# fill up containers: naive indicates no renormalized values (of beta) are used for fit
# naive: no renormalization at all

pnaive <- c()
for (i in seq(1, length(xis))) {
    row <- data$beta == opt$b & abs(data$xiin - xis[i]) < 0.01
    pnaive[i] <- data$p[row]
    bsamplescontlimitnaive[, i] <- arrayp[, row]
}
print(length(pnaive))

# cubic fits with nothing renormalized
fitplaqnaive <- try(bootstrap.nlsfit(fncub, c(0.71, 1, 1, 1), x = xis^2, y = pnaive, bsamples = bsamplescontlimitnaive))

# save and plot results for naive limits
resultslist <- list(intercepts = intercepts, plaqren = plaqren,
                    fitsslope = fitsslope, fitsp = fitsplaquette,
                    fitplaqnaive = fitplaqnaive, fitplaq = fitplaq)


defaultusr <- par("usr")
par("usr" = c(0, 1, defaultusr[3], defaultusr[4]))
try(plot(fitplaqnaive, plot.range = c(-0.2, 1.2),
        main = sprintf("continuum limit plaquette: %f+/-%f, chi = %f",
        fitplaqnaive$t0[1], fitplaqnaive$se[1], fitplaqnaive$chi / fitplaqnaive$dof),
        ylim = c((fitplaqnaive$t0[1] - fitplaqnaive$se[1]), max(result$p)),
        xlab = "", xaxs = "i", xlim = c(0, 1)))
try(plotwitherror(x = c(0), y = c(fitplaqnaive$t0[1]), dy = (fitplaqnaive$se[1]),
        col = 2, pch = 2, rep = TRUE))
par("usr" = defaultusr)

# for each polynomial:
# fitplaq: xi=xi_input, p=p_interpolated -> only beta renormalized
# for each: do fit to continuum limit, plot, add results to list and table
# list: region 1-5 fitplaq

fitspolynomial <- list()
i <- 1
resultspolynomial <- data.frame(degree = NA, lim = NA, chi = NA,
                    p = NA, type = NA, limplot = NA, dlimplot = NA)
for (fun in c(fnlin, fnpar, fncub, fnqar, fnqin)) {
    print(i)

    fitplaq <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
            x = xis[c(1, 2, 3, 4, 5, 6, 8)]^2, y = result$p,
            bsamples = na.omit(bsamplescontlimit)))
    fitspolynomial[[i]] <- fitplaq
    if (!inherits(fitplaq, "try-error")) {
    plot(fitplaq, main = sprintf("cont. lim. plaq.: %f+/-%f, chi = %f, p = %f,\ndegree of polynomial:%d",
            fitplaq$t0[1], fitplaq$se[1],
            fitplaq$chi / fitplaq$dof, fitplaq$Qval, i),
            plot.range = c(-0.2, 1.2), xlab = "", xaxs = "i")
    resultspolynomial <- rbind(resultspolynomial,
            data.frame(degree = i, lim = tex.catwitherror(fitplaq$t0[1],
            fitplaq$se[1], digits = 2, with.dollar = FALSE),
            chi = fitplaq$chi / fitplaq$dof, p = fitplaq$Qval,
            type = "plaq", limplot = fitplaq$t0[1], dlimplot = fitplaq$se[1]))
    } else {
        plotwitherror(x = xis[c(1, 2, 3, 4, 5, 6, 8)]^2, y = result$p, dy = result$dp,
            main = sprintf("degree of polynomial:%d", i),
            plot.range = c(-0.2, 1.2), xlab = "", xaxs = "i")
            print("fit not successful")
    }
    i <- i + 1
}
resultspolynomial <- resultspolynomial[-1, ]
namepol <- "plotstikz/polynomialslope.csv"
# write out result
write.table(resultspolynomial, namepol, col.names = TRUE, row.names = FALSE)

}

# write out results
print(result)
print(fitresults)
print(resultspolynomial)

write.table(result, sprintf("%fresultsrenormalizationslope.csv", opt$respath),
             col.names = TRUE, row.names = FALSE, append = FALSE)
write.table(fitresults, sprintf("%ffitresultsrenormalizationslope.csv", opt$respath),
            col.names = TRUE, row.names = FALSE, append = FALSE)
saveRDS(resultslist, sprintf("%flistresultsrenormalizationslope.RData", opt$respath))
saveRDS(fitspolynomial, sprintf("%flistpolynomialrenormalizationslope.RData", opt$respath))

# move all plots into subfolder
system(sprintf("mv -v tikz* %f", opt$respath))
