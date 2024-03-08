library("hadron")
library(optparse)


if (TRUE) {
option_list <- list(
    make_option(c("-s", "--bootsamples"), type = "integer", default = 500,
    help = "how many bootstrapsamples should be drawn [default %default]"),
    make_option(c("-b", "--beta"), type = "double", default = 1.7,
    help = "beta at xi=1 [default %default]"),
    make_option(c("-L", "--length"), type = "integer", default = 16,
    help = "spatial extent of lattice at xi=1 [default %default]"),

    make_option(c("-T", "--timeextent"), type = "integer", default = 16,
    help = "time extent of lattice at xi=1 [default %default]"),
    make_option(c("--crzero"), type = "double", default = -1.65,
    help = "c used for determining r_0 [default %default]"),
    make_option(c("--myfunctions"), type = "character",
        default = "/hiskp4/gross/masterthesis/analyse/code/U1_analyse_potential/",
#~     make_option(c("--myfunctions"), type = "character", default = "myfunctions.R",
    help = "path to where additional functions are stored [default %default]"),

    make_option(c("--respath"), type = "character", default = "plotstikz/",
    help = "path to where the resulting plots and data are stored [default %default]"),
    make_option(c("--summaryname"), type = "character", default = "summaryfile",
    help = "name of the summary file of the data [default %default]"),
    make_option(c("--datapath"), type = "character", default = "./",
    help = "path to where the data for the analyzed configs are stored [default %default]"),
    make_option(c("--type"), type = "character", default = "normal",
    help = "type of ensembles that shoule be analysed, one of normal, sideways and slope [default %default]"),
    make_option(c("-o", "--omit"), type = "integer", default = -1,
    help = "omission of highest points from the potential, -1=any [default %default]"),
    make_option(c("--indexfitcontlim"), type = "integer", default = 1,
    help = "for the fit to the contlim, only xi with index larger than this are used [default %default]"),
    make_option(c("--lowlimxi"), type = "integer", default = 1,
    help = "points lower than this are not used to determine xi [default %default]"),
    make_option(c("--lowlimpot"), type = "integer", default = -1,
    help = "points lower than this are not used to determine the potential and r0. If negative, the same as lowlimxi [default %default]"),
    make_option(c("--fitlim"), type = "double", default = 0.3,
    help = "how much may the value of r0 deviate from the xi=1 value to still be considered? [default %default]"),

    make_option(c("--naive"), action = "store_true", default = FALSE,
    help = "if true, continuum limit with beta and xi unrenormalized is calculated [default %default]"),
    make_option(c("--xiconst"), action = "store_true", default = FALSE,
    help = "if true, xi(beta) is assumed to be constant [default %default]"),

    make_option(c("--aic"), action = "store_true", default = FALSE,
    help = "if true, effective masses are determined with weighting
        according to the akaike information criterion [default %default]"),
    make_option(c("--scaletauint"), action = "store_true", default = FALSE,
    help = "if true, errors and bootstrapsamples for the correlator are rescaled
    such that the effects of autocorrelation are taken into account [default %default]"),
    make_option(c("--errortotpot"), action = "store_true", default = FALSE,
    help = "if true, potential bootstrap samples are rescaled before analysis
        to take total error into account [default %default]"),
    make_option(c("--xisingle"), action = "store_true", default = FALSE,
    help = "if true, continuum limit is fit in powers of xi instead of xi^2 [default %default]"),
    make_option(c("--xifix"), action = "store_true", default = FALSE,
    help = "if true, use a list with more xis for the extrapolation [default %default]")

)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
source(paste(opt$myfunctions, "myfunctions.R", sep = ""))

githash <- printgitcommit(opt$myfunctions)

factorrzero <- 1

if (opt$lowlimpot < 0) opt$lowlimpot <- opt$lowlimxi

type <- opt$type
if (! (type == "normal" || type == "sideways" || type == "slope" || type == "sidewaysstring" || type == "normalstring")) {
    stop(paste("type", type, "not allowed"))
}

## sanity check for fit mask:
if (opt$indexfitcontlim < 1) {
    print("indexfitcontlim has to be larger than zero, setting it to 1")
    opt$indexfitcontlim <- 1
}

}

# This script is used to determine the renormalised beta and
# the continuum limit after each ensemble has been analysed

# first, the data are read in, then the optimal beta, xi_ren and P are determined,
# and then the continuum limit is determined with different polynomial fits


# set filenames, read in results, set up containers for bootstrapsamples
if (type == "sideways" || type == "sidewaysstring") {
dataname <- sprintf("%s/resultsummary2p1dsidewaysb%.3fNs%d.csv", opt$datapath, opt$beta, opt$length)
filenameres <- sprintf("%s/resultssideways", opt$datapath)
side <- 2
}
if (type == "normal" || type == "normalstring") {
dataname <- sprintf("%s/resultsummary2p1dnormalb%.3fNs%d.csv", opt$datapath, opt$beta, opt$length)
filenameres <- sprintf("%s/resultsnormal", opt$datapath)
side <- 2
}
if (type == "slope") {
dataname <- sprintf("%s/summarysmallbetaone%fL%d.csv", opt$datapath, opt$beta, opt$length)
filenameres <- sprintf("%s/resultsmallscaled", opt$datapath)
side <- 2
end <- ""
}

## also possibility to give a custom name for the summaryfile
## if not default, assume custom name is given
if (opt$summaryname != "summaryfile") {
    dataname <- opt$summaryname
    paste("assuming custom filename for summary:", dataname)
}

if (!file.exists(dataname)) {
    stop(paste("file for results does not exist, check for", dataname))
}
data <- read.table(dataname, header = TRUE, sep = " ")
# data <- na.omit(data)
if (type == "sideways" || type == "sidewaysstring") {
data <- data[data$c == opt$crzero, ]
data <- data[data$lowlim == opt$lowlimxi, ]
data <- data[data$lowlimpot == opt$lowlimpot, ]
}

if (type == "normal" || type == "normalstring") {
data <- data[data$c == opt$crzero, ]
data <- data[data$lowlim == opt$lowlimxi, ]
data <- data[data$lowlimpot == opt$lowlimpot, ]
}

if (opt$omit >= 0) {
data <- data[data$omit == opt$omit, ]
}

# print(data)
end2 <- ""


## filter for AIC and scaled errors
if (opt$aic) {
data <- data[data$aic == TRUE, ]
end2 <- sprintf("%saic", end2)
} else {
data <- data[data$aic == FALSE, ]
}

if (opt$scaletauint) {
data <- data[data$scaletauint == TRUE, ]
end2 <- sprintf("%sscaletauintetp%d", end2, opt$errortotpot)
} else {
data <- data[data$scaletauint == FALSE, ]
}

## filter if results with total or onyl statistical uncertainty should be used,
## only aplicable if AIC was used

if (opt$errortotpot && opt$aic) {
    data <- data[data$errortotpot == TRUE, ]
} else {
    data <- data[data$errortotpot == FALSE, ]
}


data$r0 <- data$r0 * factorrzero

nom <- length(data$beta)

# read in bootstrapsamples
bootsamples <- opt$bootsamples
arrayrzero <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))
arrayp <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))
arrayxi <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))
intercepts <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))

nalist <- rep(FALSE, opt$bootsamples)
for (i in seq(1, nom)) {
    string <- sprintf("i = %d, beta = %f, Ns = %d, Nt = %d, xi = %f",
                    i, data$beta[i], data$Ns[i], data$Nt[i], data$xi[i])
    print(string)
    if (type == "normal" || type == "sideways") {
    if (opt$crzero == -1.65) end <- sprintf("omit%dllxi%dllr0%d%s", data$omit[i], data$lowlim[i], data$lowlimpot[i], end2)
    if (opt$crzero != -1.65) end <- sprintf("omit%dllxi%dllr0%dc%.2f%s", data$omit[i], data$lowlim[i], data$lowlimpot[i], opt$crzero, end2)
    result <- readinbootstrapsamples(beta = data$beta[i], Ns = data$Ns[i],
                    Nt = data$Nt[i], xi = data$xi[i], columns = c(1, 1, 1),
                    names = c("bsrzeros", "bsp", "bsxicalc"), filename = filenameres, end = end)
    }    
    if (type == "normalstring" || type == "sidewaysstring") {
    if (opt$crzero == -1.65) end <- sprintf("omit%dllxi%dllr0%d%s", data$omit[i], data$lowlim[i], data$lowlimpot[i], end2)
    if (opt$crzero != -1.65) end <- sprintf("omit%dllxi%dllr0%dc%.2f%s", data$omit[i], data$lowlim[i], data$lowlimpot[i], opt$crzero, end2)
    result <- readinbootstrapsamples(beta = data$beta[i], Ns = data$Ns[i],
                    Nt = data$Nt[i], xi = data$xi[i], columns = c(1, 1, 1),
                    names = c("bsst", "bsp", "bsxicalc"), filename = filenameres, end = end)
    }
    if (type == "slope") {
    result <- readinbootstrapsamples(beta = data$beta[i], Ns = data$Ns[i],
                    Nt = data$Nt[i], xi = data$xi[i], columns = c(1, 1, 1),
                    names = c("bsslopescaled", "bsp", "bsxiren"), filename = filenameres)
    }
    # print(head(result))
    arrayrzero[, i] <- result[, 1] * factorrzero
    arrayp[, i] <- result[, 2]
    arrayxi[, i] <- result[, 3]
    nalist <- nalist | is.na(result[, 1]) | is.na(result[, 2]) | is.na(result[, 3])
    if (sum(is.na(result)) > 0) print(sum(is.na(result)))
}
print(sum(nalist))
# stop if there are too many NAs in the data
stopifnot(sum(nalist) < 0.1*opt$bootsamples)
# remove all NAs from dataset, save this list to be able to remove the same bootstrapsample from other observables
arrayrzero <- arrayrzero[!nalist, ]
arrayxi <- arrayxi[!nalist, ]
arrayp <- arrayp[!nalist, ]
bootsamples <- bootsamples - sum(nalist)

arraypst <- parametric.bootstrap(boot.R = bootsamples, x = data$puwst, dx = data$dpuwst, seed = 123456)

if (type == "slope") {
    data$xicalc <- apply(arrayxi, 2, mean, na.rm = F)
    data$dxicalc <- apply(arrayxi, 2, sd, na.rm = F)
    data$r0 <- apply(arrayrzero, 2, mean, na.rm = F)
    data$dr0 <- apply(arrayrzero, 2, sd, na.rm = F)
    data$xi <- data$xiin
    try(data[c("xiin", "ratio", "dratio", "st", "dst", "chipot", "icslope", "dicslope", "icpot", "dicpot", "logpot", "dlogpot", "ratioslope", "dratioslope", "ratiopot", "dratiopot", "rzero", "drzero", "puw", "dpuw", "job")] <- NULL)
    # print(data)
}

##  we want to use the string tension to fix a_s, but still use the same code, so we fix the name by renaming
if (type == "sidewaysstring" || type == "normalstring") {
    data$r0 <- data$st
    data$dr0 <- data$dst
    
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

xiconststr <- ""
if (opt$xiconst) xiconststr <- "xiconst"
endname <- sprintf("%sbeta%fomit%d%sllxi%dllr0%dfl%.2f", type, opt$beta, opt$omit, xiconststr, opt$lowlimxi, opt$lowlimpot, opt$fitlim)
if(opt$crzero != -1.65) endname <- sprintf("%sc%.2f", endname, opt$crzero)
if (opt$aic) endname <- sprintf("%saic", endname)
if (opt$scaletauint) endname <- sprintf("%sscaletauintetp%d", endname, opt$errortotpot)

for (size in c(0.65)) {

nameplot <- sprintf("%srenorm%s", path, endname)

if (size <= 1) {
    tikzfile <- tikz.init(nameplot, width = mmtoinches(400),
                    height = mmtoinches(400 * size), packages = packages)
    }
if (size > 1) {
    tikzfile <- tikz.init(nameplot, width = mmtoinches(200),
                    height = mmtoinches(200 * size), packages = packages)
    }
defaultmargin <- par(c("mai"))
if (type == "normal" || type == "sideways" || type == "normalstring" || type == "sidewaysstring") {
    par(mai = c(defaultmargin[1] * max(1, 0.8 * fontsize),
            defaultmargin[2] * max(1, 0.8 * fontsize), 0.1, 0.1))
}
# if (type == "sideways") {
#     par(mai = c(defaultmargin[1] * max(1, 0.8 * fontsize),
#                 0.01, 0.1, defaultmargin[2] * max(1, 0.8 * fontsize)))
# }
print(nameplot)

# set up limits, masks for selecting the right points
mask <- data$beta == opt$beta & data$xi == 1 #& data$c == -1.65
if(opt$omit == -1) mask <- mask & data$omit == 0
maskone <- mask
rzeroone <- data$r0[mask]
xlim <- c(min(data$beta, na.rm = T), max(data$beta, na.rm = T))
ylim <- c(max(rzeroone - 2 * opt$fitlim, min(data$r0 - data$dr0, na.rm = T)), min(rzeroone + 2 * opt$fitlim, max(data$r0 + data$dr0, na.rm = T)))
if (size > 1) { xlim <- c(1.45, 1.75)}

# set input anisotropies that were considered, container for results
if(opt$xifix) xis <- c(1, 0.8, 2/3, 0.5, 0.4, 1/3, 0.25, 0.2)

xis <- sort(unique(as.numeric(format(data$xi, digits=6))), decreasing=TRUE)
if(opt$indexfitcontlim >= length(xis)) stop(paste("indexfitcontlim is too large! Maximum", length(xis)-1))

cols <- c(1, 3, 4, 5, 6, 9, 10, 11, 12, 8)
fitsrzero <- list(NULL)
fitsplaquette <- list(NULL)
fitsplaquettest <- list(NULL)
fitsxi <- list(NULL)
intercepts <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
plaqren <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
plaqrenst <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
xiphys <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
interceptsimple <- c()
xisimple <- c()
psimple <- c()
pstsimple <- c()

legendtext <- c("1.000")
xvalues <- seq(opt$beta-0.4, opt$beta+0.2, length.out=41)

par(lwd = linewidth)

# plot result of xi_input = 1 as a line and as a point
plot(x = xvalues, y = rep(data$r0[mask], 41),
            cex.lab = fontsize, cex.axis = fontsize, cex = fontsize, ylab = "",
            xlab = "", type = "l", ylim = ylim, xlim = xlim)
points(x = xvalues,
        y = rep(data$r0[mask] - data$dr0[mask], 41), type = "l", lty = 2)
points(x = xvalues,
        y = rep(data$r0[mask] + data$dr0[mask], 41), type = "l", lty = 2)
points(x = xvalues,
        y = rep(data$r0[mask] - opt$fitlim, 41), type = "l", lty = 3, col = cols[length(xis) + 1])
points(x = xvalues,
        y = rep(data$r0[mask] + opt$fitlim, 41), type = "l", lty = 3, col = cols[length(xis) + 1])
plotwitherror(x = data$beta[maskone], y = data$r0[maskone],
        dy = data$dr0[maskone], col = 1, pch = 1, cex = fontsize, rep = TRUE)


if (opt$xiconst) {
    fitresults <- data.frame(xiin = NA, r0slope = NA, r0intercept = NA, chir0 = NA, pr0 = NA,
                        plaqslope = NA, plaqintercept = NA, chiplaq = NA, pplaq = NA,
                        plaqstslope = NA, plaqstintercept = NA, chiplaqst = NA, pplaqst = NA,
                        xiintercept = NA, chixi = NA, pxi = NA)
} else {
    fitresults <- data.frame(xiin = NA, r0slope = NA, r0intercept = NA, chir0 = NA, pr0 = NA,
                        plaqslope = NA, plaqintercept = NA, chiplaq = NA, pplaq = NA,
                        plaqstslope = NA, plaqstintercept = NA, chiplaqst = NA, pplaqst = NA,
                        xislope = NA, xiintercept = NA, chixi = NA, pxi = NA)
}

# for each xi input != 1:
# select the points, make linear fits to r_zero, plaquette, xi_ren as a function of beta
# plot r_zero(beta) including error regions
# r_zero was determined from the combination of the spatial and rescaled temporal potential
# determine beta_ren: r_zero(beta_ren, xi_input) = r_zero(xi_input=1)
# determine P(beta_ren), xi(beta_ren)
# put everything in one dataframe, nicely formatted for easy printing
start <- 1
maskfitlim <- abs(data$r0 - data$r0[maskone]) < opt$fitlim
if (length(data$xi[data$xi == 1 & maskfitlim] < 3)) start <- 2

## for the independent fits, it is not necessary to reomve the unneeded anisotropied from the fit,
## because they do not influence the degrees of freedom or the results for the other fits.
# start <- max(start, opt$indexfitcontlim)

if (start != 1) {
    intercepts[, 1] <- rep(opt$beta, bootsamples)
    plaqren[, 1] <- arrayp[, maskone]
    plaqrenst[, 1] <- arraypst[, maskone]
    xiphys[, 1] <- arrayxi[, maskone]

    interceptsimple <- c(opt$beta)
    xisimple <- data$xicalc[maskone]
    psimple <- data$p[maskone]
    pstsimple <- data$puwst[maskone]

    fitsrzero[[1]] <- list(rzero = data$r0[maskone], drzero = data$dr0[maskone], fitlim = opt$fitlim)
}

print(paste("start=", start))
for (i in seq(start, length(xis))) {
    mask <- abs(data$xi - xis[i]) < 0.01 & abs(data$r0 - data$r0[maskone]) < opt$fitlim #&data$c == -1.65
    maskplot <- abs(data$xi - xis[i]) < 0.01 #& data$c == -1.65
    fitsrzero[[i]] <- try(bootstrap.nlsfit(fnlin, c(1, 1), data$r0[mask],
                            data$beta[mask], na.omit(arrayrzero[, mask])))
    # cat("\n\nxi=", xis[i], "rzero\n")
    # summary(fitsrzero[[i]])
    removed <- attributes(na.omit(arrayrzero[, mask]))$na.action
    # print(removed)
    fitsplaquette[[i]] <- try(bootstrap.nlsfit(fnlin, c(1, 1), data$p[mask],
                            data$beta[mask], na.omit(arrayp[, mask])))
    fitsplaquettest[[i]] <- try(bootstrap.nlsfit(fnlin, c(1, 1), data$puwst[mask],
                            data$beta[mask], na.omit(arraypst[, mask])))
    # cat("\n\nxi=", xis[i], "plaq\n")
    # summary(fitsplaquette[[i]])
    if (opt$xiconst) {
        fitsxi[[i]] <- try(bootstrap.nlsfit(fncon, c(1), data$xicalc[mask],
                            data$beta[mask], na.omit(arrayxi[, mask]), success.infos = 1:4))

    } else {
        fitsxi[[i]] <- try(bootstrap.nlsfit(fnlin, c(1, 1), data$xicalc[mask],
                            data$beta[mask], na.omit(arrayxi[, mask]), success.infos = 1:4))
    }
    # cat("\n\nxi=", xis[i], "xiphys\n")
    # summary(fitsxi[[i]])

    if (!inherits(fitsrzero[[i]], "try-error") && !inherits(fitsplaquette[[i]], "try-error") && !inherits(fitsxi[[i]], "try-error")) {
        try(errorpolygon(X = xvalues, fitsrzero[[i]], col.p = cols[i],
                col.band = cols[i], cex = fontsize, arlength = 0.05 * fontsize))
        try(plotwitherror(x = data$beta[maskplot], y = data$r0[maskplot],
                dy = data$dr0[maskplot], col = cols[i], pch = cols[i], cex = fontsize, rep = TRUE))
        interceptsimple[i] <- (rzeroone - fitsrzero[[i]]$t0[1]) / fitsrzero[[i]]$t0[2]
        interceptintermediate <- getintercept(fitsrzero[[i]], arrayrzero[, maskone], bootsamples = length(fitsrzero[[i]]$t[, 1]))
        intercepts[, i] <- fillexceptna(removed, interceptintermediate)

        # cannot use predict, because taking the sd of all means underestimates the error, for each intercept have to take the appropriate boot
        prediction <- predictwithxerror.bootstrapfit(fitsplaquette[[i]], interceptintermediate)
        plaqren[, i] <- fillexceptna(removed, prediction$boot)
        prediction <- predictwithxerror.bootstrapfit(fitsplaquettest[[i]], interceptintermediate)
        plaqrenst[, i] <- fillexceptna(removed, prediction$boot)
        prediction <- predictwithxerror.bootstrapfit(fitsxi[[i]], interceptintermediate)

        if (opt$xiconst) {
            xiphys[, i] <- fitsxi[[i]]$t[, 1]
        } else {
            xiphys[, i] <- fillexceptna(removed, prediction$boot)
        }

        xisimple[i] <- fitsxi[[i]]$fn(par = fitsxi[[i]]$t0, x = interceptsimple[i], boot.r = 0)
        psimple[i] <- fitsplaquette[[i]]$fn(par = fitsplaquette[[i]]$t0, x = interceptsimple[i], boot.r = 0)
        pstsimple[i] <- fitsplaquettest[[i]]$fn(par = fitsplaquettest[[i]]$t0, x = interceptsimple[i], boot.r = 0)

        newline <- data.frame(xiin = xis[i],
            r0slope = tex.catwitherror(fitsrzero[[i]]$t0[2], fitsrzero[[i]]$se[2], with.dollar = FALSE, digits = 2),
            r0intercept = tex.catwitherror(fitsrzero[[i]]$t0[1], fitsrzero[[i]]$se[1], with.dollar = FALSE, digits = 2),
            chir0 = fitsrzero[[i]]$chisqr / fitsrzero[[i]]$dof, pr0 = fitsrzero[[i]]$Qval,

            plaqslope = tex.catwitherror(fitsplaquette[[i]]$t0[2], fitsplaquette[[i]]$se[2], with.dollar = FALSE, digits = 2),
            plaqintercept = tex.catwitherror(fitsplaquette[[i]]$t0[1], fitsplaquette[[i]]$se[1], with.dollar = FALSE, digits = 2),
            chiplaq = fitsplaquette[[i]]$chisqr / fitsplaquette[[i]]$dof, pplaq = fitsplaquette[[i]]$Qval,

            plaqstslope = tex.catwitherror(fitsplaquettest[[i]]$t0[2], fitsplaquettest[[i]]$se[2], with.dollar = FALSE, digits = 2),
            plaqstintercept = tex.catwitherror(fitsplaquettest[[i]]$t0[1], fitsplaquettest[[i]]$se[1], with.dollar = FALSE, digits = 2),
            chiplaqst = fitsplaquettest[[i]]$chisqr / fitsplaquettest[[i]]$dof, pplaqst = fitsplaquettest[[i]]$Qval,
            xiintercept = tex.catwitherror(fitsxi[[i]]$t0[1], fitsxi[[i]]$se[1], with.dollar = FALSE, digits = 2),
            chixi = fitsxi[[i]]$chisqr / fitsxi[[i]]$dof, pxi = fitsxi[[i]]$Qval)

        if(!opt$xiconst) {
            newline2 <- data.frame(
            xislope = tex.catwitherror(fitsxi[[i]]$t0[2], fitsxi[[i]]$se[2], with.dollar = FALSE, digits = 2))
            newline <- cbind(newline, newline2)
        }
        fitresults <- rbind(fitresults, newline)
    } else {
        try(plotwitherror(x = data$beta[maskplot], y = data$r0[maskplot],
                dy = data$dr0[maskplot], col = cols[i], pch = cols[i], cex = fontsize, rep = TRUE))
        newline <- data.frame(xiin = xis[i], r0slope = NA, r0intercept = NA, chir0 = NA, pr0 = NA,
                        plaqslope = NA, plaqintercept = NA, chiplaq = NA, pplaq = NA,
                        xiintercept = NA, chixi = NA, pxi = NA)
        if(!opt$xiconst) {
            newline2 <- data.frame(xislope = NA)
            newline <- cbind(newline, newline2)
        }
        fitresults <- rbind(fitresults, newline)
        interceptsimple[i] <- NA
        xisimple[i] <- NA
        psimple[i] <- NA
        pstsimple[i] <- NA
    }

    legendtext[i] <- sprintf("%.3f", xis[i])
}


# arrows(x0=1.42, x1=1.42, y0=0, y1=10)
# make legend for plot
legendtext[length(xis) + 1] <- "bounds fit"

fitresults <- fitresults[-1, ]
if (type == "normal" || type == "sideways") {
mtext("$r_0 / a_s$", side = side, line = distance, cex = fontsize) #ylab
legend(legend = legendtext, x = "topleft", title = "$\\xi_\\text{input} = $",
        col = c(cols), pch = c(cols), cex = fontsize)
}
if (type == "normalstring" || type == "sidewaysstring") {
mtext("$a_s \\sigma$", side = side, line = distance, cex = fontsize) #ylab
legend(legend = legendtext, x = "bottomleft", title = "$\\xi_\\text{input} = $",
        col = c(cols), pch = c(cols), cex = fontsize)
}
if (type == "slope") {
mtext("$a_s \\Delta V$", side = side, line = distance, cex = fontsize) #ylab
legend(legend = legendtext, x = "bottomleft", title = "$\\xi_\\text{input} = $",
        col = c(cols), pch = c(cols), cex = fontsize)
}
mtext("$\\beta$", side = 1, line = distance, cex = fontsize) #xlab
if (type == "sideways") {
    if (size  ==  1.33) {
        legend(legend = "sideways", x = "bottomright", cex = fontsize, bty = "n")
    }
}
if (type == "normal") {
    if (size == 1.33) {
        legend(legend = "normal", x = "bottomright", cex = fontsize, bty = "n")
    }
}
box(bty = "o")

tikz.finalize(tikzfile)

# make a dataframe of results
result <- data.frame(xiin = xis, beta = apply(intercepts, 2, mean, na.rm = F),
                    dbeta = apply(intercepts, 2, sd, na.rm = F),
                    xiphys = apply(xiphys, 2, mean, na.rm = F),
                    dxiphys = apply(xiphys, 2, sd, na.rm = F),
                    p = apply(plaqren, 2, mean, na.rm = F), dp = apply(plaqren, 2, sd, na.rm = F),
                    pst = apply(plaqrenst, 2, mean, na.rm = F), dpst = apply(plaqrenst, 2, sd, na.rm = F))
# result$xiphys[result$xiin == 1] <- data$xicalc[data$xi == 1]
# result$dxiphys[result$xiin == 1] <- data$dxicalc[data$xi == 1]
result <- cbind(result, data.frame(betasimple = interceptsimple, xisimple = xisimple, psimple = psimple, pstsimple = pstsimple))
par(mai = defaultmargin)
print(result)
}

## add more data to results for easy printing
# newframe <- data.frame(xiopt = tex.catwitherror(result$xiphys[seq(start, length(xis))], result$dxiphys[seq(2, length(xis))], digits = 2, with.dollar = FALSE),
# betaopt = tex.catwitherror(result$beta[seq(start, length(xis))], result$dbeta[seq(2, length(xis))], digits = 2, with.dollar = FALSE),
# popt = tex.catwitherror(result$p[seq(start, length(xis))], result$dp[seq(2, length(xis))], digits = 2, with.dollar = FALSE))


# plot results to pdf
pdf(sprintf("tikzplotallfits%sranges.pdf", endname), title="")

# result of all linear fits
for (i in seq(start, length(xis))){
    try(plot(fitsrzero[[i]], main = sprintf("rzero, xiin = %f, chi = %f, p = %f",
            xis[i], fitsrzero[[i]]$chi / fitsrzero[[i]]$dof, fitsrzero[[i]]$Qval), xlab="beta", ylab="r0/as"))
    try(plot(fitsplaquette[[i]], main = sprintf("P, xiin = %f, chi = %f, p = %f",
            xis[i], fitsplaquette[[i]]$chi / fitsplaquette[[i]]$dof, fitsplaquette[[i]]$Qval), xlab="beta", ylab="P"))
    try(plot(fitsplaquettest[[i]], main = sprintf("P_st, xiin = %f, chi = %f, p = %f",
            xis[i], fitsplaquettest[[i]]$chi / fitsplaquettest[[i]]$dof, fitsplaquettest[[i]]$Qval), xlab="beta", ylab="P_st"))
    try(plot(fitsxi[[i]], main = sprintf("xi, xiin = %f, chi = %f, p = %f",
            xis[i], fitsxi[[i]]$chi / fitsxi[[i]]$dof, fitsxi[[i]]$Qval), xlab="beta", ylab="xi_ren"))
   }



# save and plot results for naive limits
resultslist <- list(intercepts = intercepts, xiphys = xiphys, plaqren = plaqren, plaqrenst = plaqrenst,
                    fitsrzero = fitsrzero, fitsxi = fitsxi,
                    fitsp = fitsplaquette, fitspst = fitsplaquettest, githash = githash, nalist = nalist)




namesave <- sprintf("%s/resultsrenormalization%s.csv", opt$respath, endname)
write.table(result, file=namesave, col.names = TRUE, row.names = FALSE, append = FALSE)

namesave <- sprintf("%s/fitresultsrenormalization%s.csv", opt$respath, endname)
write.table(fitresults, file=namesave, col.names = TRUE, row.names = FALSE, append = FALSE)

namesave <- sprintf("%s/listresultsrenormalization%s.RData", opt$respath, endname)
saveRDS(resultslist, file=namesave)


if (TRUE) {
# empty containers for result
bsamplescontlimit <- array(rep(NA, bootsamples * (2 * length(xis))), dim = c(bootsamples, 2 * length(xis)))
bsamplescontlimitnaive <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
bsamplescontlimitnaivexiren <- array(rep(NA, 2 * bootsamples * (length(xis))), dim = c(bootsamples, 2 * length(xis)))
bsamplescontlimitbeta <- array(rep(NA, 2 * bootsamples * (length(xis))), dim = c(bootsamples, 2 * length(xis)))
bsamplescontlimitst <- array(rep(NA, bootsamples * (2 * length(xis))), dim = c(bootsamples, 2 * length(xis)))

pnaive <- c()
xirennaive <- c()

# fill up containers: naive indicates no renormalized values (of beta) are used for fit
# naive: no renormalization at all
# naivexiren: use renormalized xi, but leave beta fixed
# beta: take continuum limit for beta(xi) as well
bsamplescontlimit[, seq(1, length(xis))] <- plaqren
bsamplescontlimit[, seq(length(xis) + 1, 2 * length(xis))] <- xiphys^2
bsamplescontlimitst[, seq(1, length(xis))] <- plaqrenst
bsamplescontlimitst[, seq(length(xis) + 1, 2 * length(xis))] <- xiphys^2
bsamplescontlimitbeta[, seq(length(xis) + 1, 2 * length(xis))] <- xiphys^2
if(opt$xisingle) {
bsamplescontlimit[, seq(length(xis) + 1, 2 * length(xis))] <- xiphys
bsamplescontlimitbeta[, seq(length(xis) + 1, 2 * length(xis))] <- xiphys
}
for (i in seq(1, length(xis))) {
    if (opt$naive) {
    row <- data$beta == opt$beta & abs(data$xi - xis[i]) < 0.01
    pnaive[i] <- data$p[row]
    xirennaive[i] <- data$xicalc[row]
    bsamplescontlimitnaive[, i] <- arrayp[, row]
    bsamplescontlimitnaivexiren[, i] <- arrayp[, row]
    bsamplescontlimitnaivexiren[, length(xis) + i] <- arrayxi[, row]^2
    if (opt$xisingle) bsamplescontlimitnaivexiren[, length(xis) + i] <- arrayxi[, row]
    }
    bsamplescontlimitbeta[, i] <- intercepts[, i]
    # bsamplescontlimitbeta[, length(xis) + i] <- arrayxi[, row]^2
}

if (start != 1) {
    bsamplescontlimitbeta[, 1] <- parametric.bootstrap(bootsamples, c(opt$beta), c(1e-4))
}






# for each polynomial:
# fitplaqnaive: xi=xi_input, p=p_meas -> neither beta nor xi reenormalized
# fitplaqnaivexiren: xi=xi_ren, p=p_meas -> beta not renormalized
# fitplaq: xi=xi_ren, p=p_interpolated -> xi and beta renormalized
# for each: do fit to continuum limit, plot, add results to list and table
# list: region 1-5 fitplaqnaive, region 6-10 fitplaqnaivexiren, region 11-15 fitplaq, region 16-20 beta
# na.omit removes entire row containing na
# do this for all possible fits in the region lower=xis[1], upper=length(xis)-2

        listfits <- list(fits = list(), lowlimfit = c(), uplimfit = c())
        resultspolynomial <- data.frame(
            degree = NA, lim = NA, chi = NA,
            p = NA, type = NA, limplot = NA, dlimplot = NA, 
            lowlimfit = NA, uplimfit = NA
        )
indexlim <- 1

listnames <- paste0(rep(c("naive", "xirennaive", "plaq", "beta", "p_st"), each=5), seq(1, 5))

for (lowlimfit in seq(opt$indexfitcontlim, length(xis) - 2)) {
    for (uplimfit in seq(max(length(xis) - 2, lowlimfit + 2), length(xis))) {
        print(paste("lowlim", lowlimfit, "uplim", uplimfit))
        fitspolynomial <- list()

        xiinfit <- xis^2
        xinaivefit <- xirennaive^2
        xirenfit <- result$xisimple^2
        if (opt$xisingle) {
            xiinfit <- xis
            xinaivefit <- xirennaive
            xirenfit <- result$xisimple
        }

        maskfitcontlim <- seq(lowlimfit, uplimfit)
        i <- 1
        for (fun in c(fnlin, fnpar, fncub, fnqar, fnqin)) {
            if (uplimfit - lowlimfit > i) {
            print(paste("Doing fits of polynomial degree", i))
            if (opt$naive) {
                # naive xi and beta
                fitplaqnaive <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
                    x = xiinfit, y = pnaive, bsamples = bsamplescontlimitnaive,
                    mask = maskfitcontlim
                ))
                fitspolynomial[[listnames[i]]] <- fitplaqnaive

                try(plot(fitplaqnaive,
                    main = sprintf(
                        "continuum limit plaquette: %4f + /-%4f, chi = %4f, p = %4f,\ndegree of polynomial:%d\nlowlim fit: %d, uplim fit: %d",
                        fitplaqnaive$t0[1], fitplaqnaive$se[1],
                        fitplaqnaive$chi / fitplaqnaive$dof, fitplaqnaive$Qval, i, lowlimfit, uplimfit
                    ),
                    plot.range = c(-0.2, 1.2), xaxs = "i", xlim = c(0, 1.05),
                    ylim = c((fitplaqnaive$t0[1] - fitplaqnaive$se[1]), max(na.omit(result$p))), xlab = "xi_in^2", ylab = "P"
                ))
                try(plotwitherror(
                    rep = TRUE, col = "red",
                    x = fitplaqnaive$x[fitplaqnaive$mask],
                    dx = fitplaqnaive$dx[fitplaqnaive$mask],
                    y = fitplaqnaive$y[fitplaqnaive$mask],
                    dy = fitplaqnaive$dy[fitplaqnaive$mask]
                ))

                try(resultspolynomial <- rbind(
                    resultspolynomial,
                    data.frame(
                        degree = i, lim = tex.catwitherror(fitplaqnaive$t0[1],
                            fitplaqnaive$se[1],
                            digits = 2, with.dollar = FALSE
                        ),
                        chi = fitplaqnaive$chi / fitplaqnaive$dof, p = fitplaqnaive$Qval,
                        type = "naive", limplot = fitplaqnaive$t0[1],
                        dlimplot = fitplaqnaive$se[1], lowlimfit = lowlimfit, uplimfit = uplimfit
                    )
                ))

                # naive beta, xi renorm
                fitplaqnaivexiren <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
                    x = xinaivefit, y = pnaive,
                    bsamples = bsamplescontlimitnaivexiren,
                    mask = maskfitcontlim
                ))
                fitspolynomial[[listnames[5 + i]]] <- fitplaqnaivexiren

                try(plot(fitplaqnaivexiren,
                    main = sprintf(
                        "continuum limit plaquette: %4f + /-%4f, chi = %4f, p = %4f,\ndegree of polynomial:%d\nlowlim fit: %d, uplim fit: %d",
                        fitplaqnaivexiren$t0[1], fitplaqnaivexiren$se[1],
                        fitplaqnaivexiren$chi / fitplaqnaivexiren$dof, fitplaqnaivexiren$Qval, i, lowlimfit, uplimfit
                    ),
                    plot.range = c(-0.2, 1.2), xaxs = "i", xlim = c(0, 1.05),
                    ylim = c((fitplaqnaivexiren$t0[1] - fitplaqnaivexiren$se[1]), max(na.omit(result$p))), xlab = "xi_ren^2", ylab = "P"
                ))

                try(plotwitherror(
                    rep = TRUE, col = "red",
                    x = fitplaqnaivexiren$x[fitplaqnaivexiren$mask],
                    dx = fitplaqnaivexiren$dx[fitplaqnaivexiren$mask],
                    y = fitplaqnaivexiren$y[fitplaqnaivexiren$mask],
                    dy = fitplaqnaivexiren$dy[fitplaqnaivexiren$mask]
                ))

                try(resultspolynomial <- rbind(
                    resultspolynomial,
                    data.frame(
                        degree = i, lim = tex.catwitherror(fitplaqnaivexiren$t0[1],
                            fitplaqnaivexiren$se[1],
                            digits = 2, with.dollar = FALSE
                        ),
                        chi = fitplaqnaivexiren$chi / fitplaqnaivexiren$dof,
                        p = fitplaqnaivexiren$Qval, type = "naivexiren",
                        limplot = fitplaqnaivexiren$t0[1], dlimplot = fitplaqnaivexiren$se[1],
                        lowlimfit = lowlimfit, uplimfit = uplimfit
                    )
                ))

                ## print if there was an error in fitting
                if (inherits(fitplaqnaive, "try-error")) {
                    print(paste("There was an error with fitplaqnaive for polynomial degree", i))
                }
                if (inherits(fitplaqnaivexiren, "try-error")) {
                    print(paste("There was an error with fitplaqnaivexiren for polynomial degree", i))
                }
            }

            # xi and beta renorm
            # print(attributes(na.omit(bsamplescontlimit))$na.action)
            fitplaq <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
                x = xirenfit, y = result$psimple, bsamples = bsamplescontlimit,
                mask = maskfitcontlim
            ))

            if (!inherits(fitplaq, "try-error")) {
                fitspolynomial[[listnames[10 + i]]] <- fitplaq

                try(plot(fitplaq,
                    main = sprintf(
                        "continuum limit plaquette: %4f + /-%4f, chi = %4f, p = %4f,\ndegree of polynomial:%d\nlowlim fit: %d, uplim fit: %d",
                        fitplaq$t0[1], fitplaq$se[1],
                        fitplaq$chi / fitplaq$dof, fitplaq$Qval, i, lowlimfit, uplimfit
                    ),
                    plot.range = c(-0.2, 1.2),
                    ylim = c((fitplaq$t0[1] - fitplaq$se[1]), max(na.omit(result$p))),
                    xaxs = "i", xlim = c(0, 1), xlab = "xi_ren^2", ylab = "P"
                ))
                try(plotwitherror(
                    rep = TRUE, col = "red",
                    x = fitplaq$x[fitplaq$mask],
                    dx = fitplaq$dx[fitplaq$mask],
                    y = fitplaq$y[fitplaq$mask],
                    dy = fitplaq$dy[fitplaq$mask]
                ))
                try(resultspolynomial <- rbind(
                    resultspolynomial,
                    data.frame(
                        degree = i, lim = tex.catwitherror(fitplaq$t0[1],
                            fitplaq$se[1],
                            digits = 2, with.dollar = FALSE
                        ),
                        chi = fitplaq$chi / fitplaq$dof, p = fitplaq$Qval,
                        type = "plaq", limplot = fitplaq$t0[1], dlimplot = fitplaq$se[1],
                        lowlimfit = lowlimfit, uplimfit = uplimfit
                    ))
                )
            }

            # beta cont limit
            fitbeta <- try(bootstrap.nlsfit(fun, rep(0.1, i + 1),
                x = xirenfit, y = result$betasimple, bsamples = bsamplescontlimitbeta,
                mask = maskfitcontlim
            ))
            fitspolynomial[[listnames[15 + i]]] <- fitbeta

            try(plot(fitbeta,
                main = sprintf(
                    "continuum limit beta: %4f + /-%4f, chi = %4f, p = %4f,\ndegree of polynomial:%d\nlowlim fit: %d, uplim fit: %d",
                    fitbeta$t0[1], fitbeta$se[1],
                    fitbeta$chi / fitbeta$dof, fitbeta$Qval, i, lowlimfit, uplimfit
                ),
                plot.range = c(-0.2, 1.2),
                ylim = c(1.3, 1.75),
                xlab = "xi_ren^2", ylab = "beta_ren", xaxs = "i", xlim = c(0, 1.05)
            ))
                try(plotwitherror(
                    rep = TRUE, col = "red",
                    x = fitbeta$x[fitbeta$mask],
                    dx = fitbeta$dx[fitbeta$mask],
                    y = fitbeta$y[fitbeta$mask],
                    dy = fitbeta$dy[fitbeta$mask]
                ))

            try(resultspolynomial <- rbind(
                resultspolynomial,
                data.frame(
                    degree = i, lim = tex.catwitherror(fitbeta$t0[1],
                        fitbeta$se[1],
                        digits = 2, with.dollar = FALSE
                    ),
                    chi = fitbeta$chi / fitbeta$dof, p = fitbeta$Qval,
                    type = "beta", limplot = fitbeta$t0[1], dlimplot = fitbeta$se[1],
                    lowlimfit = lowlimfit, uplimfit = uplimfit
                )
            ))

            # xi and beta renorm, spatial-temporal plaquette
            # print(attributes(na.omit(bsamplescontlimit))$na.action)
            fitplaqst <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
                x = xirenfit, y = result$pstsimple, bsamples = bsamplescontlimitst,
                mask = maskfitcontlim
            ))

            if (!inherits(fitplaqst, "try-error")) {
                fitspolynomial[[listnames[20 + i]]] <- fitplaqst

                try(plot(fitplaqst,
                    main = sprintf(
                        "continuum limit temporal plaquette: %4f + /-%4f, chi = %4f, p = %4f,\ndegree of polynomial:%d\nlowlim fit: %d, uplim fit: %d",
                        fitplaqst$t0[1], fitplaqst$se[1],
                        fitplaqst$chi / fitplaqst$dof, fitplaqst$Qval, i, lowlimfit, uplimfit
                    ),
                    plot.range = c(-0.2, 1.2),
                    ylim = c(min(na.omit(result$pst)), (fitplaqst$t0[1] - fitplaqst$se[1])),
                    xaxs = "i", xlim = c(0, 1), xlab = "xi_ren^2", ylab = "P_st"
                ))
                try(plotwitherror(
                    rep = TRUE, col = "red",
                    x = fitplaqst$x[fitplaqst$mask],
                    dx = fitplaqst$dx[fitplaqst$mask],
                    y = fitplaqst$y[fitplaqst$mask],
                    dy = fitplaqst$dy[fitplaqst$mask]
                ))

                try(resultspolynomial <- rbind(
                    resultspolynomial,
                    data.frame(
                        degree = i, lim = tex.catwitherror(fitplaqst$t0[1],
                            fitplaqst$se[1],
                            digits = 2, with.dollar = FALSE
                        ),
                        chi = fitplaqst$chi / fitplaqst$dof, p = fitplaqst$Qval,
                        type = "p_st", limplot = fitplaqst$t0[1], dlimplot = fitplaqst$se[1],
                        lowlimfit = lowlimfit, uplimfit = uplimfit
                    )
                ))
            }

            ## print if there were any errors in fitting
            if (inherits(fitplaq, "try-error")) {
                print(paste("There was an error with fitplaq for polynomial degree", i))
            }
            if (inherits(fitbeta, "try-error")) {
                print(paste("There was an error with fitbeta for polynomial degree", i))
            }
            if (inherits(fitplaqst, "try-error")) {
                print(paste("There was an error with fitplaqst for polynomial degree", i))
            }
            i <- i + 1
            }
        }
        listfits$fits[[indexlim]] <- fitspolynomial
        listfits$lowlimfit[indexlim] <- lowlimfit
        listfits$uplimfit[indexlim] <- uplimfit
        indexlim <- indexlim + 1
    }
}
}
        resultspolynomial <- resultspolynomial[-1, ]
        # write out result
        print(resultspolynomial)


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
        # print(fitresults)
        # print(result)



listfits$githash <- githash
listfits$nalist <- nalist
if (opt$xisingle) endname <- paste0(endname, "xisingle")

namepol <- sprintf("%s/polynomial%s.csv", opt$respath, endname)
write.table(resultspolynomial, namepol, col.names = TRUE, row.names = FALSE)

namesave <- sprintf("%s/listpolynomialrenormalization%s.RData", opt$respath, endname)
saveRDS(listfits, file = namesave)
# move all plots into subfolder
# system(sprintf("mv -v tikz* %s", opt$respath))
