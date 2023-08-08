library("hadron")
library(optparse)


if (TRUE) {
option_list <- list(
    make_option(c("-s", "--bootsamples"), type = "integer", default = 500,
    help = "how many bootstrapsamples should be drawn [default %default]"),
    make_option(c("-b", "--beta"), type = "double", default = 1.7,
    help = "beta at xi = 1 [default %default]"),
    make_option(c("-L", "--length"), type = "integer", default = 16,
    help = "spatial extent of lattice at xi = 1 [default %default]"),

    make_option(c("-T", "--timeextent"), type = "integer", default = 16,
    help = "time extent of lattice at xi = 1 [default %default]"),
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
    help = "omission of highest points from the potential, -1 = any [default %default]"),
    make_option(c("--indexfitcontlim"), type = "integer", default = 1,
    help = "for the fit to the contlim, only xi with index larger than this are used [default %default]"),
    make_option(c("--fitlim"), type = "double", default = 0.3,
    help = "how much may the value of r0 deviate from the xi = 1 value to still be considered? [default %default]"),

    make_option(c("--naive"), action = "store_true", default = FALSE,
    help = "if true, continuum limit with beta and xi unrenormalized is calculated [default %default]"),
    make_option(c("--xiconst"), action = "store_true", default = FALSE,
    help = "if true, xi(beta) is assumed to be constant [default %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
source(paste(opt$myfunctions, "myfunctions.R", sep = ""))
source(paste(opt$myfunctions, "fit_xiyey.R", sep = ""))

githash <- printgitcommit(opt$myfunctions)

type <- opt$type
if (! (type == "normal" || type == "sideways" || type == "slope")) {
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
if (type == "sideways") {
dataname <- sprintf("%s/resultsummary2p1dsidewaysb%.3fNs%d.csv", opt$datapath, opt$beta, opt$length)
filenameres <- sprintf("%s/resultsrotated", opt$datapath)
side <- 2
}
if (type == "normal") {
dataname <- sprintf("%s/resultsummary2p1dnormalb%.3fNs%d.csv", opt$datapath, opt$beta, opt$length)
filenameres <- sprintf("%s/resultssubtracted", opt$datapath)
side <- 2
}
if (type == "slope") {
dataname <- sprintf("%s/summarysmallbetaone%fL%d.csv", opt$datapath, opt$beta, opt$length)
filenameres <- sprintf("%s/resultsmallscaled", opt$datapath)
side <- 2
}

## also possibility to give a custom name for the summaryfile
## if not default, assume custom name is given
if (opt$summaryname !=  "summaryfile") {
    dataname <- opt$summaryname
    paste("assuming custom filename for summary:", dataname)
}

if (!file.exists(dataname)) {
    stop(paste("file for results does not exist, check for", dataname))
}
data <- read.table(dataname, header = TRUE, sep = " ")
# data <- na.omit(data)
if (type == "sideways") {
data <- data[data$c == opt$crzero, ]
}

if (type == "normal") {
data <- data[data$c == opt$crzero, ]
}

if (opt$omit >=  0) {
data <- data[data$omit == opt$omit, ]
}

nom <- length(data$beta)

# read in bootstrapsamples
bootsamples <- opt$bootsamples
arrayrzero <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))
arrayp <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))
arrayxi <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))
intercepts <- array(rep(NA, bootsamples * nom), dim = c(bootsamples, nom))


for (i in seq(1, nom)) {
    string <- sprintf("i = %d, beta = %f, Ns = %d, Nt = %d, xi = %f",
                    i, data$beta[i], data$Ns[i], data$Nt[i], data$xi[i])
    print(string)
    if (type == "normal" || type == "sideways") {
    end <- sprintf("omit%d", data$omit[i])
    result <- readinbootstrapsamples(beta = data$beta[i], Ns = data$Ns[i],
                    Nt = data$Nt[i], xi = data$xi[i], columns = c(1, 1, 1),
                    names = c("bsrzeros", "bsp", "bsxicalc"), filename = filenameres, end = end)
    }
    if (type == "slope") {
    result <- readinbootstrapsamples(beta = data$beta[i], Ns = data$Ns[i],
                    Nt = data$Nt[i], xi = data$xi[i], columns = c(1, 1, 1),
                    names = c("bsslopescaled", "bsp", "bsxiren"), filename = filenameres)
    }
    # print(head(result))
    arrayrzero[, i] <- result[1:opt$bootsamples, 1]
    arrayp[, i] <- result[1:opt$bootsamples, 2]
    arrayxi[, i] <- result[1:opt$bootsamples, 3]
    if (sum(is.na(result)) > 0) print(sum(is.na(result)))
}

if (type == "slope") {
    data$xicalc <- apply(arrayxi, 2, mean, na.rm = T)
    data$dxicalc <- apply(arrayxi, 2, sd, na.rm = T)
    data$r0 <- apply(arrayrzero, 2, mean, na.rm = T)
    data$dr0 <- apply(arrayrzero, 2, sd, na.rm = T)
    data$xi <- data$xiin
    try(data[c("xiin", "ratio", "dratio", "st", "dst", "chipot", "icslope", "dicslope", "icpot", "dicpot", "logpot", "dlogpot", "ratioslope", "dratioslope", "ratiopot", "dratiopot", "rzero", "drzero", "puw", "dpuw", "job")] <- NULL)
    # print(data)
}

## define ansatz for multifit
## goal: fit obs(beta, xi) with one single fit.
## x is a matrix containing all x-values x_1, x_2, ...
ansatzmulti <- function(x, par) {
  beta <- x[1]
  xi <- x[2]
  a <- par[1]
  if (abs(xi-1.0) < 1e-3) {
    b <- par[2]
  } else if (abs(xi - 0.8) < 1e-3) {
    b <- par[3]
  } else if (abs(xi - 2./3.) < 1e-3) {
    b <- par[4]
  } else if (abs(xi - 0.5) < 1e-3) {
    b <- par[5]
  } else if (abs(xi - 0.4) < 1e-3) {
    b <- par[6]
  } else if (abs(xi - 1./3.) < 1e-3) {
    b <- par[7]
  } else if (abs(xi - 0.25) < 1e-3) {
    b <- par[8]
  } else {
    stop(paste("invalid xi in fit ansatz, you gave", xi))
  }
  return (a * beta + b)
}

## ansatz for multifit
## goal: fit obs(beta, xi) with one single fit.
## x is a matrix containing all x-values x_1, x_2, ...
ansatzmulticonst <- function(x, par) {
  xi <- x[2]
  if (abs(xi-1.0) < 1e-3) {
    b <- par[1]
  } else if (abs(xi - 0.8) < 1e-3) {
    b <- par[2]
  } else if (abs(xi - 2./3.) < 1e-3) {
    b <- par[3]
  } else if (abs(xi - 0.5) < 1e-3) {
    b <- par[4]
  } else if (abs(xi - 0.4) < 1e-3) {
    b <- par[5]
  } else if (abs(xi - 1./3.) < 1e-3) {
    b <- par[6]
  } else if (abs(xi - 0.25) < 1e-3) {
    b <- par[7]
  } else {
    stop(paste("invalid xi in fit ansatz, you gave", xi))
  }
  return (b)
}




xis <- c(1, 0.8, 2/3, 0.5, 0.4, 1/3, 0.25)
intercepts <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
plaqren <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
xiphys <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
intercepts[, 1] <- rep(opt$beta, bootsamples)

maskomit <- data$beta == opt$beta & data$xi == 1
if (opt$omit == -1) maskomit <- maskomit & data$omit == 0
plaqren[, 1] <- arrayp[, maskomit]
xiphys[, 1] <- arrayxi[, maskomit]


## do fit - similar to predictbeta, but hopefully less messy
## for r0, P, xi

rzerozero <- data$r0[data$beta == opt$beta & data$xi == 1]
if(opt$omit == -1) rzerozero <- data$r0[data$beta == opt$beta & data$xi == 1 & data$omit == 0]
mask <- abs(data$r0 - rzerozero) < opt$fitlim
x <- matrix(data = c(data$beta[mask], data$xi[mask]), byrow = TRUE, nrow = 2)

start <- 1
if (length(data$xi[data$xi == 1 & mask] < 3)) start <- 2

## r0

y <- data$r0[mask]
dy <- data$dr0[mask]
y_bts <- arrayrzero[, mask]
resultrzero <- usebootstrap.fit_xiyey(ansatz = ansatzmulti, x = x, y = y, ey = dy, y_bts = y_bts, guess = c(5, rep(-6, 7)), N_bts = opt$bootsamples)
summary.multifit(resultrzero)

for (i in seq(start, length(xis))) {
    intercepts[, i] <- getinterceptfromparams(bsa = resultrzero$par$bts[, 1], bsb = resultrzero$par$bts[, 1 + i], rzeroone = arrayrzero[, data$beta == opt$beta & data$xi == 1], bootsamples = opt$bootsamples)
}

## xi

y <- data$xicalc[mask]
dy <- data$dxicalc[mask]
y_bts <- arrayxi[, mask]

print("xiphys")
if(!opt$xiconst) {
resultxi <- usebootstrap.fit_xiyey(ansatz = ansatzmulti, x = x, y = y, ey = dy, y_bts = y_bts, guess = c(5, rep(-6, 7)), N_bts = opt$bootsamples)
summary.multifit(resultxi)
for (i in seq(start, length(xis))) {
    xren <- matrix(data = c(intercepts[, i], rep(xis[i], opt$bootsamples)), byrow = TRUE, nrow = 2)
    for (bs in seq(1, opt$bootsamples)) {
        xiphys[bs, i] <- ansatzmulti(x = xren[, bs], par = resultxi$par$bts[bs, ])
    }
}
}
if(opt$xiconst) {
resultxi <- usebootstrap.fit_xiyey(ansatz = ansatzmulticonst, x = x, y = y, ey = dy, y_bts = y_bts, guess = xis, N_bts = opt$bootsamples)
summary.multifit(resultxi)
xiphys <- resultxi$par$bts
}




## p

y <- data$p[mask]
dy <- data$dp[mask]
y_bts <- arrayp[, mask]

print("plaq")
resultp <- usebootstrap.fit_xiyey(ansatz = ansatzmulti, x = x, y = y, ey = dy, y_bts = y_bts, guess = c(5, rep(-6, 7)), N_bts = opt$bootsamples)
summary.multifit(resultp)

for (i in seq(start, length(xis))) {
    xren <- matrix(data = c(intercepts[, i], rep(xis[i], opt$bootsamples)), byrow = TRUE, nrow = 2)
    for (bs in seq(1, opt$bootsamples)) {
        plaqren[bs, i] <- ansatzmulti(x = xren[, bs], par = resultp$par$bts[bs, ])
    }
}

# make a dataframe of results
result <- data.frame(xiin = xis, beta = apply(intercepts, 2, mean, na.rm = T),
                    dbeta = apply(intercepts, 2, sd, na.rm = T),
                    xiphys = apply(xiphys, 2, mean, na.rm = T),
                    dxiphys = apply(xiphys, 2, sd, na.rm = T),
                    p = apply(plaqren, 2, mean, na.rm = T), dp = apply(plaqren, 2, sd, na.rm = T))
print(result)

resultslist <- list(resultp = resultp, plaqren = plaqren, resultxi = resultxi,
                    xiphys = xiphys, resultrzero = resultrzero, intercepts = intercepts)


## contlimit: like in predictbeta

title <- paste("multifit", opt$type, "b", opt$beta, "omit", opt$omit, "cont", opt$indexfitcontlim, sep = "")
if(opt$xiconst) title <- paste(title, "xiconst", sep="")
pdf(paste(title, ".pdf", sep=""), title = "")

if (TRUE) {

bsamplescontlimit <- array(rep(NA, bootsamples * (2 * length(xis))), dim = c(bootsamples, 2 * length(xis)))
bsamplescontlimit[, seq(1, length(xis))] <- plaqren
bsamplescontlimit[, seq(length(xis) + 1, 2 * length(xis))] <- xiphys^2


# empty containers for result
bsamplescontlimit <- array(rep(NA, bootsamples * (2 * length(xis))), dim = c(bootsamples, 2 * length(xis)))
bsamplescontlimitnaive <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
bsamplescontlimitnaivexiren <- array(rep(NA, 2 * bootsamples * (length(xis))), dim = c(bootsamples, 2 * length(xis)))
bsamplescontlimitbeta <- array(rep(NA, 2 * bootsamples * (length(xis))), dim = c(bootsamples, 2 * length(xis)))
pnaive <- c()
xirennaive <- c()

# fill up containers: naive indicates no renormalized values (of beta) are used for fit
# naive: no renormalization at all
# naivexiren: use renormalized xi, but leave beta fixed
# beta: take continuum limit for beta(xi) as well
bsamplescontlimit[, seq(1, length(xis))] <- plaqren
bsamplescontlimit[, seq(length(xis) + 1, 2 * length(xis))] <- xiphys^2
bsamplescontlimitbeta[, seq(length(xis) + 1, 2 * length(xis))] <- xiphys^2
for (i in seq(1, length(xis))) {
    bsamplescontlimitbeta[, i] <- intercepts[, i]
    # bsamplescontlimitbeta[, length(xis) + i] <- arrayxi[, row]^2
}

if (length(data$xi[data$xi == 1 & mask] < 3)) {
    bsamplescontlimitbeta[, 1] <- parametric.bootstrap(bootsamples, c(opt$beta), c(1e-4))
}






# for each polynomial:
# fitplaq: xi = xi_ren, p = p_interpolated -> xi and beta renormalized
# for each: do fit to continuum limit, plot, add results to list and table
# list: region 1-5 fitplaqnaive, region 6-10 fitplaqnaivexiren, region 11-15 fitplaq, region 16-20 beta
# na.omit removes entire row containing na

fitspolynomial <- list()
resultspolynomial <- data.frame(degree = NA, lim = NA, chi = NA,
                    p = NA, type = NA, limplot = NA, dlimplot = NA)
maskfitcontlim <- seq(opt$indexfitcontlim, length(xis))
i <- 1
for (fun in c(fnlin, fnpar, fncub, fnqar, fnqin)){
    print(paste("Doing fits of polynomial degree", i))
# xi and beta renorm
# print(attributes(na.omit(bsamplescontlimit))$na.action)
    fitplaq <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
                x = result$xiphys^2, y = result$p, bsamples = na.omit(bsamplescontlimit),
                mask = maskfitcontlim))

    if (!inherits(fitplaq, "try-error")) {
    fitspolynomial[[i]] <- fitplaq
    plot(fitplaq, main = sprintf("continuum limit plaquette: %f + /-%f, chi = %f, p = %f,\ndegree of polynomial:%d",
            fitplaq$t0[1], fitplaq$se[1],
            fitplaq$chi / fitplaq$dof, fitplaq$Qval, i),
            plot.range = c(-0.2, 1.2),
            ylim = c((fitplaq$t0[1] - fitplaq$se[1]), max(result$p)),
            xaxs = "i", xlim = c(0, 1), xlab = "xi_ren^2", ylab = "P")
    resultspolynomial <- rbind(resultspolynomial,
            data.frame(degree = i, lim = tex.catwitherror(fitplaq$t0[1],
            fitplaq$se[1], digits = 2, with.dollar = FALSE),
            chi = fitplaq$chi / fitplaq$dof, p = fitplaq$Qval,
            type = "plaq", limplot = fitplaq$t0[1], dlimplot = fitplaq$se[1]))
    }

# beta cont limit
    fitbeta <- try(bootstrap.nlsfit(fun, rep(0.1, i + 1),
                x = result$xiphys^2, y = result$beta, bsamples = na.omit(bsamplescontlimitbeta),
                mask = maskfitcontlim))
    fitspolynomial[[5 + i]] <- fitbeta
    try(plot(fitbeta, main = sprintf("continuum limit beta: %f + /-%f, chi = %f, p = %f,\ndegree of polynomial:%d",
            fitbeta$t0[1], fitbeta$se[1],
            fitbeta$chi / fitbeta$dof, fitbeta$Qval, i),
            plot.range = c(-0.2, 1.2),
            ylim = c(1.3, 1.75),
            xlab = "xi_ren^2", ylab = "beta_ren", xaxs = "i", xlim = c(0, 1.05)))
    try(resultspolynomial <- rbind(resultspolynomial,
            data.frame(degree = i, lim = tex.catwitherror(fitbeta$t0[1],
            fitbeta$se[1], digits = 2, with.dollar = FALSE),
            chi = fitbeta$chi / fitbeta$dof, p = fitbeta$Qval,
            type = "beta", limplot = fitbeta$t0[1], dlimplot = fitbeta$se[1])))

## print if there were any errors in fitting
    if (inherits(fitplaq, "try-error")) {
        print(paste("There was an error with fitplaq for polynomial degree", i))
    }
    if (inherits(fitbeta, "try-error")) {
        print(paste("There was an error with fitbeta for polynomial degree", i))
    }
    i <- i + 1
}
# plotwitherror(x = result$xiphys^2, y = result$beta, dy = apply(bsamplescontlimitbeta[, seq(1, length(xis))], 2, sd), dx = apply(bsamplescontlimitbeta[, seq(length(xis)+1, 2 * length(xis))], 2, sd))
# plotwitherror(x = result$xiphys^2, y = result$p, dy = result$dp, dx = apply(bsamplescontlimitbeta[, seq(length(xis)+1, 2 * length(xis))], 2, sd))
print(resultspolynomial)

}


## plot
supports <- 1000
cols <- c(1, 3, 4, 5, 6, 9, 10, 8)

xlim <- c(0.95 * min(data$beta[mask]), 1.05 * max(data$beta[mask]))
ylim <- c(0.95 * min(data$r0[mask]), 1.05 * max(data$r0[mask]))

plot(NA, xlim = xlim * c(1.05, 0.95), ylim = ylim * c(1.05, 0.95), xlab = "beta", ylab = "r_0/a_s", main = "r_0/a_s(beta, xi) with multifit")
legendtext <- c()

# if (start !=  1) {
#     plotwitherror(x = data$beta[mask & data$xi == 1], y = data$r0[mask & data$xi == 1],
#         dy = data$dr0[mask & data$xi == 1], col = cols[1], pch = cols[1], rep = TRUE)
#     append(legendtext, "1.00")
# }


maskomit <- data$beta == opt$beta & data$xi == 1
if (opt$omit == -1) maskomit <- maskomit & data$omit == 0
drzerozero <- data$dr0[maskomit]
print(rzerozero)
print(drzerozero)
lines(x = seq(xlim[1], xlim[2], len = 100), y = rep(rzerozero, 100), lty = 1)
lines(x = seq(xlim[1], xlim[2], len = 100), y = rep(rzerozero + drzerozero, 100), lty = 2)
lines(x = seq(xlim[1], xlim[2], len = 100), y = rep(rzerozero - drzerozero, 100), lty = 2)

for (i in seq(1, length(xis))) {
    print("xi")
    print(xis[i])
    xplot <- matrix(c(seq(xlim[1], xlim[2], len = supports), rep(xis[i], supports)), byrow = TRUE, nrow = 2)
    errorpolygonmultifit(X = xplot, fitresult = resultrzero, col.p = cols[i], pch = cols[i],
        mask = abs(resultrzero$x[2, ] - xis[i]) < 1e-3, ylim = ylim)
    legendtext[i] <-  sprintf("%.2f", xis[i])
}

legend(x = "topleft", legend = legendtext, col = cols, pch = cols)



xlim <- c(0.95 * min(data$beta[mask]), 1.05 * max(data$beta[mask]))
ylim <- c(0.95 * min(data$xicalc[mask]), 1.05 * max(data$xicalc[mask]))

plot(NA, xlim = xlim * c(1.05, 0.95), ylim = ylim * c(1.05, 0.95), xlab = "beta", ylab = "xicalc", main = "xicalc(beta, xi) with multifit")
legendtext <- c()

# if (start !=  1) {
#     plotwitherror(x = data$beta[mask & data$xi == 1], y = data$xicalc[mask & data$xi == 1],
#         dy = data$dxicalc[mask & data$xi == 1], col = cols[1], pch = cols[1], rep = TRUE)
#     append(legendtext, "1.00")
# }

for (i in seq(1, length(xis))) {
    print("xi")
    print(xis[i])
    xplot <- matrix(c(seq(xlim[1], xlim[2], len = supports), rep(xis[i], supports)), byrow = TRUE, nrow = 2)
    errorpolygonmultifit(X = xplot, fitresult = resultxi, col.p = cols[i], pch = cols[i],
        mask = abs(resultrzero$x[2, ] - xis[i]) < 1e-3, ylim = ylim)
    legendtext[i] <-  sprintf("%.2f", xis[i])
}

legend(x = "topleft", legend = legendtext, col = cols, pch = cols)



xlim <- c(0.95 * min(data$beta[mask]), 1.05 * max(data$beta[mask]))
ylim <- c(0.95 * min(data$p[mask]), 1.05 * max(data$p[mask]))

plot(NA, xlim = xlim * c(1.05, 0.95), ylim = ylim * c(1.05, 0.95), xlab = "beta", ylab = "p", main = "p(beta, xi) with multifit")
legendtext <- c()

# if (start !=  1) {
#     plotwitherror(x = data$beta[mask & data$xi == 1], y = data$p[mask & data$xi == 1],
#         dy = data$dp[mask & data$xi == 1], col = cols[1], pch = cols[1], rep = TRUE)
#     append(legendtext, "1.00")
# }

for (i in seq(1, length(xis))) {
    print("xi")
    print(xis[i])
    xplot <- matrix(c(seq(xlim[1], xlim[2], len = supports), rep(xis[i], supports)), byrow = TRUE, nrow = 2)
    errorpolygonmultifit(X = xplot, fitresult = resultp, col.p = cols[i], pch = cols[i],
        mask = abs(resultrzero$x[2, ] - xis[i]) < 1e-3, ylim = ylim)
    legendtext[i] <-  sprintf("%.2f", xis[i])
}

legend(x = "topleft", legend = legendtext, col = cols, pch = cols)


# write out results like in predictbeta
# contains:
# result: data frame with xiin, beta, dbeta(renormalized), xiphys, dxiphys, p, dp, betasimple (renormed beta from intercept from mean values)
# resultspolynomial: dataframe with degree of polynomial, result of cont. lim. with catwitherror, chi and p of lim, type of lim, result of lim as two doubles
# resultlist: intercepts = beta_ren, xiphys, plaqren , fitsrzero, fitsxi, fitsp (results of linear interpolations), fitplaq, fitplaqnaive, fitplaqnaivexiren(cubic polynomial cont limits)
# fitspolynomial: bootstrapnlsfit results of cont limit, region 1-5 fitplaqnaive, region 6-10 fitplaqnaivexiren, region 11-15 fitplaq
# print(fitresults)
# print(result)


resultspolynomial <- resultspolynomial[-1, ]
if (type == "normal") namepol <- sprintf("%s/polynomialnormalmultibeta%fomit%dcont%d", opt$respath, opt$beta, opt$omit, opt$indexfitcontlim)
if (type == "sideways") namepol <- sprintf("%s/polynomialsidewaysmultibeta%fomit%dcont%d", opt$respath, opt$beta, opt$omit, opt$indexfitcontlim)
if (type == "slope") namepol <- sprintf("%s/polynomialslopemultibeta%fomit%dcont%d", opt$respath, opt$beta, opt$omit, opt$indexfitcontlim)
if(opt$xiconst) namepol <- paste(namepol, "xiconst", sep="")
# write out result
write.table(resultspolynomial, paste(namepol, ".csv", sep=""), col.names = TRUE, row.names = FALSE)

xiconststr <- ""
if(opt$xiconst) xiconststr <- "xiconst"

if (type == "normal") {
write.table(result, sprintf("%s/resultsrenormalizationmultibeta%fomit%d%s.csv", opt$respath, opt$beta, opt$omit, xiconststr),
            col.names = TRUE, row.names = FALSE, append = FALSE)
saveRDS(resultslist, sprintf("%s/listresultsrenormalizationmultibeta%fomit%d%s.RData", opt$respath, opt$beta, opt$omit, xiconststr))
saveRDS(fitspolynomial, sprintf("%s/listpolynomialrenormalizationmultibeta%fomit%dcont%d%s.RData",
opt$respath, opt$beta, opt$omit, opt$indexfitcontlim, xiconststr))
}

if (type == "sideways") {
write.table(result, sprintf("%s/resultsrenormalizationsidewaysmultibeta%fomit%d%s.csv", opt$respath, opt$beta, opt$omit, xiconststr),
            col.names = TRUE, row.names = FALSE, append = FALSE)
saveRDS(resultslist, sprintf("%s/listresultsrenormalizationsidewaysmultibeta%fomit%d%s.RData", opt$respath, opt$beta, opt$omit, xiconststr))
saveRDS(fitspolynomial, sprintf("%s/listpolynomialrenormalizationsidewaysmultibeta%fomit%dcont%d%s.RData",
opt$respath, opt$beta, opt$omit, opt$indexfitcontlim, xiconststr))
}

if (type == "slope") {
write.table(result, sprintf("%s/resultsrenormalizationslopemultibeta%fomit%d%s.csv", opt$respath, opt$beta, opt$omit, xiconststr),
            col.names = TRUE, row.names = FALSE, append = FALSE)
saveRDS(resultslist, sprintf("%s/listresultsrenormalizationslopemultibeta%fomit%d%s.RData", opt$respath, opt$beta, opt$omit, xiconststr))
saveRDS(fitspolynomial, sprintf("%s/listpolynomialrenormalizationslopemultibeta%fomit%dcont%d%s.RData",
opt$respath, opt$beta, opt$omit, opt$indexfitcontlim, xiconststr))
}
