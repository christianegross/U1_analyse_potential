library("hadron")
library(optparse)

## goal of this script:
## For each xi, we want to be anle to separately decide whether we take the interpolated result
## or the result of a single point
## to do this, we read in the results from either predictbeta or multifitpredictbeta
## and then replace the results for those xi where we want to take a single point


option_list <- list(
    make_option(c("-i", "--inputfile"), type = "character",
        default = "inputpredictcontlim.csv",
    help = "file from which the input parameters are read [default %default]"),
    make_option(c("--myfunctions"), type = "character",
        default = "/hiskp4/gross/masterthesis/analyse/code/U1_analyse_potential/",
    help = "path to where additional functions are stored [default %default]"),
    make_option(c("-l", "--line"), type = "integer", default = 1,
    help = "line of inputfile that is to be used for readout [default %default]")
    )
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
source(paste(opt$myfunctions, "myfunctions.R", sep = ""))

## read in options and select correct line
options <- read.table(opt$inputfile, sep=",", header=TRUE)
options <- options[opt$line, ]

opt <- append(opt, options)
opt$singlemulti <- as.character(opt$singlemulti)
opt$respath <- as.character(opt$respath)
opt$datapath <- as.character(opt$datapath)
opt$summaryname <- as.character(opt$summaryname)
opt$xi1beta <- opt$beta

githash <- printgitcommit(opt$myfunctions)

type <- as.character(opt$type)
if (! (type == "normal" || type == "sideways" || type == "slope")) {
    stop(paste("type", type, "not allowed"))
}

if (! (opt$singlemulti == "single" || opt$singlemulti == "multi")) {
    stop(paste("singlemulti", opt$singlemulti, "not allowed"))
}

if (!file.exists(opt$summaryname)) {
    stop(paste("file for results does not exist, check for", opt$summaryname))
}
data <- read.table(opt$summaryname, header = TRUE, sep = " ")
# data <- na.omit(data)

if (type == "sideways") {
data <- data[data$c == opt$crzero, ]
data <- data[data$lowlim == opt$lowlimfitpot, ]
filenameres <- sprintf("%s/resultssideways", opt$datapath)
end <- sprintf("omit%dlowlim%d", opt$omit, opt$lowlimfitpot)
}

if (type == "normal") {
data <- data[data$c == opt$crzero, ]
data <- data[data$lowlim == opt$lowlimfitpot, ]
filenameres <- sprintf("%s/resultsnormal", opt$datapath)
end <- sprintf("omit%dlowlim%d", opt$omit, opt$lowlimfitpot)
}

if (type == "slope") {
filenameres <- sprintf("%s/resultsmallscaled", opt$datapath)
end <- ""
}

if (opt$omit >=0) {
data <- data[data$omit == opt$omit, ]
}

nom <- length(data$beta)


xiconststr <- ""
if (opt$xiconst) xiconststr <- "xiconst"
print(opt$singlemulti)
if (opt$singlemulti == "single") endname <- sprintf("%sbeta%fomit%d%slowlim%d", type, opt$beta, opt$omit, xiconststr, opt$lowlimfitpot)
if (opt$singlemulti == "multi") endname <- sprintf("multi%sbeta%fomit%d%slowlim%d", type, opt$beta, opt$omit, xiconststr, opt$lowlimfitpot)

print(sprintf("%s/listresultsrenormalization%s.RData", as.character(opt$respath), endname))
resultslist <- readRDS(sprintf("%s/listresultsrenormalization%s.RData", as.character(opt$respath), endname))

namesave <- sprintf("%s/resultsrenormalization%s.csv", as.character(opt$respath), endname)
result <- read.table(file=namesave, header=TRUE)

xis <- result$xiin
result$betachosen <- rep(FALSE, length(xis))
print(xis)

## for each xi, check if there is a fixed beta that should be used.
## If yes, replace results for plaquette and xi_ren, and the mean for beta.
## assume the error for beta that was determined before, and draw paramteric bootstrapsamples from mean and sd
for (i in seq(1, length(xis))){
    fix <- paste("xi", i, "fix", sep="")
    if (fix %in% names(opt)) {
        betafix <- paste("xi", i, "beta", sep="")
        if (as.logical(opt[[fix]])) {
            j <- which(abs(data$xi - xis[i]) < 1e-3 & abs(data$beta - opt[[betafix]]) < 1e-3 )
            readin <- readinbootstrapsamples(beta = data$beta[j], Ns = data$Ns[j],
                        Nt = data$Nt[j], xi = data$xi[j], columns = c(1, 1),
                        names = c("bsp", "bsxicalc"), filename = filenameres, end = end)
            resultslist[["plaqren"]][, i] <- readin[, 1]
            resultslist[["xiphys"]][, i] <- readin[, 2]
            resultslist[["intercepts"]][, i] <- parametric.bootstrap(boot.R = opt$bootsamples, x=opt[[betafix]], dx=result$dbeta[i], seed=123456)
            result$beta[i] <- opt[[betafix]]
            result$betasimple[i] <- opt[[betafix]]
            result$p[i] <- data$p[j]
            result$dp[i] <- data$dp[j]
            result$xiphys[i] <- data$xicalc[j]
            result$dxiphys[i] <- data$dxicalc[j]
            result$betachosen[i] <- TRUE
        }
    }
}

print(options)
print(result)

if (TRUE) {
# empty containers for result
bsamplescontlimit <- array(rep(NA, opt$bootsamples * (2 * length(xis))), dim = c(opt$bootsamples, 2 * length(xis)))
bsamplescontlimitnaive <- array(rep(NA, opt$bootsamples * (length(xis))), dim = c(opt$bootsamples, length(xis)))
bsamplescontlimitnaivexiren <- array(rep(NA, 2 * opt$bootsamples * (length(xis))), dim = c(opt$bootsamples, 2 * length(xis)))
bsamplescontlimitbeta <- array(rep(NA, 2 * opt$bootsamples * (length(xis))), dim = c(opt$bootsamples, 2 * length(xis)))
pnaive <- c()
xirennaive <- c()

# fill up containers: naive indicates no renormalized values (of beta) are used for fit
# naive: no renormalization at all
# naivexiren: use renormalized xi, but leave beta fixed
# beta: take continuum limit for beta(xi) as well
bsamplescontlimit[, seq(1, length(xis))] <- resultslist[["plaqren"]]
bsamplescontlimit[, seq(length(xis) + 1, 2 * length(xis))] <- resultslist[["xiphys"]]^2
bsamplescontlimitbeta[, seq(length(xis) + 1, 2 * length(xis))] <- resultslist[["xiphys"]]^2
bsamplescontlimitbeta[, seq(1, length(xis))] <- resultslist[["intercepts"]][, seq(1, length(xis))]

if (result$dbeta[1] < 1e-8) {
    bsamplescontlimitbeta[, 1] <- parametric.bootstrap(opt$bootsamples, c(opt$beta), c(1e-4))
}


# for each polynomial:
# fitplaqnaive: xi=xi_input, p=p_meas -> neither beta nor xi reenormalized
# fitplaqnaivexiren: xi=xi_ren, p=p_meas -> beta not renormalized
# fitplaq: xi=xi_ren, p=p_interpolated -> xi and beta renormalized
# for each: do fit to continuum limit, plot, add results to list and table
# list: region 1-5 fitplaqnaive, region 6-10 fitplaqnaivexiren, region 11-15 fitplaq, region 16-20 beta
# na.omit removes entire row containing na

fitspolynomial <- list()
resultspolynomial <- data.frame(degree = NA, lim = NA, chi = NA,
                    p = NA, type = NA, limplot = NA, dlimplot = NA)

maskfitcontlim <- seq(opt$indexfitcontlim, length(xis))
i <- 1

pdf(sprintf("%s/plotsbetachosen%s.pdf", opt$respath, endname), title="")
for (fun in c(fnlin, fnpar, fncub, fnqar, fnqin)){
    print(paste("Doing fits of polynomial degree", i))

# xi and beta renorm
# print(attributes(na.omit(bsamplescontlimit))$na.action)
    fitplaq <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
                x = result$xiphys^2, y = result$p, bsamples = bsamplescontlimit,
                mask = maskfitcontlim))

    if (!inherits(fitplaq, "try-error")) {
    fitspolynomial[[i]] <- fitplaq

    plot(fitplaq, main = sprintf("continuum limit plaquette: %f + /-%f, chi = %f, p = %f,\ndegree of polynomial:%d",
            fitplaq$t0[1], fitplaq$se[1],
            fitplaq$chi / fitplaq$dof, fitplaq$Qval, i),
            plot.range = c(-0.2, 1.2),
            ylim = c((fitplaq$t0[1] - fitplaq$se[1]), max(na.omit(result$p))),
            xaxs = "i", xlim = c(0, 1), xlab = "xi_ren^2", ylab = "P")

    resultspolynomial <- rbind(resultspolynomial,
            data.frame(degree = i, lim = tex.catwitherror(fitplaq$t0[1],
            fitplaq$se[1], digits = 2, with.dollar = FALSE),
            chi = fitplaq$chi / fitplaq$dof, p = fitplaq$Qval,
            type = "plaq", limplot = fitplaq$t0[1], dlimplot = fitplaq$se[1]))
    }

# beta cont limit
    fitbeta <- try(bootstrap.nlsfit(fun, rep(0.1, i + 1),
                x = result$xiphys^2, y = result$beta, bsamples = bsamplescontlimitbeta,
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
resultspolynomial <- resultspolynomial[-1, ]
# write out result

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
# resultslist: intercepts = beta_ren, xiphys, plaqren , fitsrzero, fitsxi, fitsp (results of linear interpolations), fitplaq, fitplaqnaive, fitplaqnaivexiren(cubic polynomial cont limits)
# fitspolynomial: bootstrapnlsfit results of cont limit, region 1-5 fitplaqnaive, region 6-10 fitplaqnaivexiren, region 11-15 fitplaq
# print(fitresults)
# print(result)
namepol <- sprintf("%s/polynomialbetachosen%scont%d.csv", opt$respath, endname, opt$indexfitcontlim)
write.table(resultspolynomial, namepol, col.names = TRUE, row.names = FALSE)
print(resultspolynomial)

namesave <- sprintf("%s/resultsrenormalizationbetachosen%s.csv", opt$respath, endname)
write.table(result, file=namesave, col.names = TRUE, row.names = FALSE, append = FALSE)

namesave <- sprintf("%s/listresultsrenormalizationbetachosen%s.RData", opt$respath, endname)
saveRDS(resultslist, file=namesave)

namesave <- sprintf("%s/listpolynomialrenormalizationbetachosen%scont%d.RData", opt$respath, endname, opt$indexfitcontlim)
saveRDS(fitspolynomial, file=namesave)

