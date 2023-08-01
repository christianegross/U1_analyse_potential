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
    make_option(c("--fitlim"), type = "double", default = 0.3,
    help = "how much may the value of r0 deviate from the xi=1 value to still be considered? [default %default]"),

    make_option(c("--naive"), action = "store_true", default = FALSE,
    help = "if true, continuum limit with beta and xi unrenormalized is calculated [default %default]")
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
if (opt$summaryname != "summaryfile") {
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

if (opt$omit >=0) {
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
    result <- readinbootstrapsamples(beta = data$beta[i], Ns = data$Ns[i],
                    Nt = data$Nt[i], xi = data$xi[i], columns = c(1, 1, 1),
                    names = c("bsrzeros", "bsp", "bsxicalc"), filename = filenameres)
    }
    if (type == "slope") {
    result <- readinbootstrapsamples(beta = data$beta[i], Ns = data$Ns[i],
                    Nt = data$Nt[i], xi = data$xi[i], columns = c(1, 1, 1),
                    names = c("bsslopescaled", "bsp", "bsxiren"), filename = filenameres)
    }
    # print(head(result))
    arrayrzero[, i] <- result[, 1]
    arrayp[, i] <- result[, 2]
    arrayxi[, i] <- result[, 3]
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
ansatz <- function(x, par){
  beta <- x[1]
  xi <- x[2]
  a <- par[1]
  if (abs(xi-1.0) < 1e-3){
    b <- par[2]
  } else if(abs(xi-0.8) < 1e-3){
    b <- par[3]
  } else if(abs(xi-2./3.) < 1e-3){
    b <- par[4]
  } else if (abs(xi-0.5) < 1e-3){
    b <- par[5]
  } else if(abs(xi-0.4) < 1e-3){
    b <- par[6]
  } else if(abs(xi-1./3.) < 1e-3){
    b <- par[7]
  } else if(abs(xi-0.25) < 1e-3){
    b <- par[8]
  } else {
    stop("invalid xi in fit ansatz")
  }
  return (a*beta + b)
}


## define a function that uses the actual bootstrapsamples
## based on function from fit_xiyey
#' fit boostrap by boostrap (generated from y and dy)
#' x = matrix of size n \times N_pts
#' returns a list with bootstraps of parameters, their mean, stderr, chi^2, reduced chi^2
usebootstrap.fit_xiyey <-function(ansatz, x, y, ey, y_bts, guess, maxiter=10000, method="BFGS", N_bts = 500){
    N_pts <- length(y) # number of points
    N_par <- length(guess)
    N_dof = N_pts - N_par

    par_bts <- matrix(data=NA, N_bts, N_par)
    ch2_bts <- matrix(data=NA, N_bts, N_par)

    # y_bts <- matrix(data=rnorm(N_bts*N_pts, y, ey), ncol=N_pts, byrow=TRUE)

    ansatz <- NA
    for (i in 1:N_bts){
        mini <- fit_xiyey(
            ansatz = ansatz, x=x, y=y_bts[i,], ey=ey,
            guess = guess,
            maxiter=maxiter, method=method)

        if (i==1){
            ansatz <- mini[["ansatz"]]
            N_dof <- mini[["N_dof"]]
        }

        par_bts[i,] <- mini[["par"]]
        ch2_bts[i,] <- mini[["ch2"]]
    }

    # mean and standard error on the bootstraps

    par_val <- apply(par_bts, 2, mean)
    par_sd <- apply(par_bts, 2, sd)

    ch2_val <- apply(ch2_bts, 2, mean)
    ch2_sd <- apply(ch2_bts, 2, sd)

    res <- list(
        ansatz=ansatz,
        N_bts=N_bts, N_pts = N_pts, N_par=N_par,
        x=x,
        par=list(bts=par_bts, val=par_val, dval=par_sd),
        ch2=list(bts=ch2_bts, val=ch2_val, dval=ch2_sd),
        ch2_dof=list(
            bts=ch2_bts/N_dof, val=ch2_val/N_dof,
            dval=ch2_sd/N_dof)
        )

    return(res)
}

xis <- c(1, 0.8, 2/3, 0.5, 0.4, 1/3, 0.25)
intercepts <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
plaqren <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
xiphys <- array(rep(NA, bootsamples * (length(xis))), dim = c(bootsamples, length(xis)))
intercepts[, 1] <- rep(opt$beta, bootsamples)
plaqren[, 1] <- arrayp[, data$beta==opt$beta & data$xi==1]
xiphys[, 1] <- arrayxi[, data$beta==opt$beta & data$xi==1]


## do fit - similar to predictbeta, but hopefully less messy
## for r0, P, xi

rzerozero <- data$r0[data$beta==opt$beta & data$xi==1]
mask <- abs(data$r0 - rzerozero) < opt$fitlim
x <- matrix(data=c(data$beta[mask], data$xi[mask]), byrow=TRUE, nrow=2)

start <- 1
if (length(data$xi[data$xi==1 & mask]<3)) start <- 2

## r0

y <- data$r0[mask]
dy <- data$dr0[mask]
y_bts <- arrayrzero[, mask]

resultrzero <- usebootstrap.fit_xiyey(ansatz=ansatz, x=x, y=y, ey=dy, y_bts = y_bts, guess=c(5, rep(-6, 7)), N_bts = opt$bootsamples)
print(resultrzero)

for(i in seq(start, length(xis))) {
    intercepts[, i] <- getinterceptfromparams(bsa = resultrzero$par$bts[, 1], bsb = resultrzero$par$bts[, 1 + i], rzeroone = arrayrzero[, data$beta==opt$beta & data$xi==1], bootsamples = opt$bootsamples)
}

## xi

y <- data$xicalc[mask]
dy <- data$dxicalc[mask]
y_bts <- arrayxi[, mask]

resultxi <- usebootstrap.fit_xiyey(ansatz=ansatz, x=x, y=y, ey=dy, y_bts = y_bts, guess=c(5, rep(-6, 7)), N_bts = opt$bootsamples)
# print(resultxi)
for(i in seq(start, length(xis))) {
    xren <- matrix(data=c(intercepts[, i], rep(xis[i], opt$bootsamples)), byrow=TRUE, nrow=2)
    for(bs in seq(1, opt$bootsamples)) {
        xiphys[bs, i] <- ansatz(x = xren[, bs], par = resultxi$par$bts[bs, ])
    }
}

## p

y <- data$p[mask]
dy <- data$dp[mask]
y_bts <- arrayp[, mask]

resultp <- usebootstrap.fit_xiyey(ansatz=ansatz, x=x, y=y, ey=dy, y_bts = y_bts, guess=c(5, rep(-6, 7)), N_bts = opt$bootsamples)
# print(resultp)

for(i in seq(start, length(xis))) {
    xren <- matrix(data=c(intercepts[, i], rep(xis[i], opt$bootsamples)), byrow=TRUE, nrow=2)
    for(bs in seq(1, opt$bootsamples)) {
        plaqren[bs, i] <- ansatz(x = xren[, bs], par = resultp$par$bts[bs, ])
    }
}

# make a dataframe of results
result <- data.frame(xiin = xis, beta = apply(intercepts, 2, mean, na.rm = T),
                    dbeta = apply(intercepts, 2, sd, na.rm = T),
                    xiphys = apply(xiphys, 2, mean, na.rm = T),
                    dxiphys = apply(xiphys, 2, sd, na.rm = T),
                    p = apply(plaqren, 2, mean, na.rm = T), dp = apply(plaqren, 2, sd, na.rm = T))
print(result)


## plot: write separate function that calculates polygon for every xi
## TODO: generalize predict.bootstrapfit to take function, x, bootstraps of parameters
## TODO: generalize errorpolygon to take arbitrary function, params


## contlimit: like in predictbeta

pdf(paste("multifitb", opt$beta, ".pdf"))


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

if (length(data$xi[data$xi==1 & mask]<3)) {
    bsamplescontlimitbeta[, 1] <- parametric.bootstrap(bootsamples, c(opt$beta), c(1e-4))
}






# for each polynomial:
# fitplaq: xi=xi_ren, p=p_interpolated -> xi and beta renormalized
# for each: do fit to continuum limit, plot, add results to list and table
# list: region 1-5 fitplaqnaive, region 6-10 fitplaqnaivexiren, region 11-15 fitplaq, region 16-20 beta
# na.omit removes entire row containing na

fitspolynomial <- list()
resultspolynomial <- data.frame(degree = NA, lim = NA, chi = NA,
                    p = NA, type = NA, limplot = NA, dlimplot = NA)
i <- 1
for (fun in c(fnlin, fnpar, fncub, fnqar, fnqin)){
    print(paste("Doing fits of polynomial degree", i))
# xi and beta renorm
# print(attributes(na.omit(bsamplescontlimit))$na.action)
    fitplaq <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
                x = result$xiphys^2, y = result$p, bsamples = na.omit(bsamplescontlimit)))

    if (!inherits(fitplaq, "try-error")) {
    fitspolynomial[[10 + i]] <- fitplaq
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
                x = result$xiphys^2, y = result$beta, bsamples = na.omit(bsamplescontlimitbeta)))
    fitspolynomial[[15 + i]] <- fitbeta
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
    if (inherits(fitplaq, "try-error")){
        print(paste("There was an error with fitplaq for polynomial degree", i))
    }
    if (inherits(fitbeta, "try-error")){
        print(paste("There was an error with fitbeta for polynomial degree", i))
    }
    i <- i + 1
}
# plotwitherror(x=result$xiphys^2, y=result$beta, dy=apply(bsamplescontlimitbeta[, seq(1, length(xis))], 2, sd), dx=apply(bsamplescontlimitbeta[, seq(length(xis)+1, 2*length(xis))], 2, sd))
# plotwitherror(x=result$xiphys^2, y=result$p, dy=result$dp, dx=apply(bsamplescontlimitbeta[, seq(length(xis)+1, 2*length(xis))], 2, sd))

resultspolynomial <- resultspolynomial[-1, ]
if (type == "normal") namepol <- sprintf("%s/polynomialnormalbeta%f.csv", opt$respath, opt$beta)
if (type == "sideways") namepol <- sprintf("%s/polynomialsidewaysbeta%f.csv", opt$respath, opt$beta)
if (type == "slope") namepol <- sprintf("%s/polynomialslopebeta%f.csv", opt$respath, opt$beta)
# write out result
write.table(resultspolynomial, namepol, col.names = TRUE, row.names = FALSE)
print(resultspolynomial)

}