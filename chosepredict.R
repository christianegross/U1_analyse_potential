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

if (opt$omit < 0) stop("omit must be non-negative!")
data <- read.table(opt$summaryname, header = TRUE, sep = " ")
# data <- na.omit(data)

if (type == "sideways") {
data <- data[data$c == opt$crzero, ]
data <- data[data$lowlim == opt$lowlimxi, ]
data <- data[data$lowlimpot == opt$lowlimpot, ]
filenameres <- sprintf("%s/resultssideways", opt$datapath)
end <- sprintf("omit%dllxi%dllr0%d", opt$omit, opt$lowlimxi, opt$lowlimpot)
}

if (type == "normal") {
data <- data[data$c == opt$crzero, ]
data <- data[data$lowlim == opt$lowlimxi, ]
data <- data[data$lowlimpot == opt$lowlimpot, ]
filenameres <- sprintf("%s/resultsnormal", opt$datapath)
end <- sprintf("omit%dllxi%dllr0%d", opt$omit, opt$lowlimxi, opt$lowlimpot)
}

if (type == "slope") {
filenameres <- sprintf("%s/resultsmallscaled", opt$datapath)
end <- ""
}

if (opt$omit >=0) {
data <- data[data$omit == opt$omit, ]
}

if (opt$crzero != -1.65) end <- sprintf("%sc%.2f", end, opt$crzero)

## filter for AIC and scaled errors
if (opt$aic) {
data <- data[data$aic == TRUE, ]
end <- sprintf("%saic", end)
} else {
data <- data[data$aic == FALSE, ]
}

if (opt$scaletauint) {
data <- data[data$scaletauint == TRUE, ]
end <- sprintf("%sscaletauintetp%d", end, opt$errortotpot)
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

nom <- length(data$beta)


xiconststr <- ""
if (opt$xiconst) xiconststr <- "xiconst"
print(opt$singlemulti)
if (opt$singlemulti == "single") endname <- sprintf("%sbeta%fomit%d%sllxi%dllr0%dfl%.2f", type, opt$beta, opt$omit, xiconststr, opt$lowlimxi, opt$lowlimpot, opt$fitlim)
if (opt$singlemulti == "multi") endname <- sprintf("multi%sbeta%fomit%d%slowlim%d", type, opt$beta, opt$omit, xiconststr, opt$lowlimxi)
if (opt$crzero != -1.65) endname <- sprintf("%sc%.2f", endname, opt$crzero)
if (opt$aic) endname <- sprintf("%saic", endname)
if (opt$scaletauint) endname <- sprintf("%sscaletauintetp%d", endname, opt$errortotpot)
endnamewrite <- endname
if (opt$xisingle) endnamewrite <- paste0(endname, "xisingle")



print(sprintf("%s/listresultsrenormalization%s.RData", as.character(opt$respath), endname))
resultslist <- readRDS(sprintf("%s/listresultsrenormalization%s.RData", as.character(opt$respath), endname))
# opt$bootsamples <- opt$bootsamples ## - sum(resultslist$nalist)

namesave <- sprintf("%s/resultsrenormalization%s.csv", as.character(opt$respath), endname)
result <- read.table(file = namesave, header = TRUE)

xis <- result$xiin
result$betachosen <- rep(FALSE, length(xis))
print(xis)

# also read in data for spatial-temporal plaquette, later do extrapolation to continuum limit with parametric bootstrap
# pst is only calculated for type normal

## for each xi, check if there is a fixed beta that should be used.
## If yes, replace results for plaquette and xi_ren, and the mean for beta.
## assume the error for beta that was determined before, and draw paramteric bootstrapsamples from mean and sd
for (i in seq(1, length(xis))){
    fix <- paste("xi", i, "fix", sep="")
    if (fix %in% names(opt)) {
        betafix <- paste("xi", i, "beta", sep="")
        if (as.logical(opt[[fix]])) {
            j <- which(abs(data$xi - xis[i]) < 1e-5 & abs(data$beta - opt[[betafix]]) < 1e-5)
            readin <- readinbootstrapsamples(beta = data$beta[j], Ns = data$Ns[j],
                        Nt = data$Nt[j], xi = data$xi[j], columns = c(1, 1),
                        names = c("bsp", "bsxicalc"), filename = filenameres,
                        end = end)
            resultslist[["plaqren"]][, i] <- readin[, 1]
            resultslist[["xiphys"]][, i] <- readin[, 2]
            # resultslist[["intercepts"]][, i] <- parametric.bootstrap(boot.R = opt$bootsamples, x=opt[[betafix]], dx=result$dbeta[i], seed=123456)
            resultslist[["intercepts"]][, i] <- resultslist[["intercepts"]][, i] - result$betasimple[i] + opt[[betafix]]
            resultslist[["plaqrenst"]][, i] <- parametric.bootstrap(boot.R = opt$bootsamples, x=data$puwst[j], dx=data$dpuwst[j], seed=123456)
            result$beta[i] <- opt[[betafix]]
            result$betasimple[i] <- opt[[betafix]]
            result$p[i] <- data$p[j]
            result$psimple[i] <- data$p[j]
            result$dp[i] <- data$dp[j]
            result$xiphys[i] <- data$xicalc[j]
            result$xisimple[i] <- data$xicalc[j]
            result$dxiphys[i] <- data$dxicalc[j]
            result$betachosen[i] <- TRUE
            result$pst[i] <- data$puwst[j]
            result$pstsimple[i] <- data$puwst[j]
            result$dpst[i] <- data$dpuwst[j]
        }
    }
}
## remove samples that contain a NA anywhere
            resultslist[["plaqren"]][resultslist$nalist, ] <- NA
            resultslist[["xiphys"]][resultslist$nalist, ] <- NA
            resultslist[["intercepts"]][resultslist$nalist, ] <- NA
            resultslist[["plaqrenst"]][resultslist$nalist, ] <- NA
            

print(options)
print(result)

if (TRUE) {
# empty containers for result
bsamplescontlimit <- array(rep(NA, opt$bootsamples * (2 * length(xis))), dim = c(opt$bootsamples, 2 * length(xis)))
bsamplescontlimitbeta <- array(rep(NA, 2 * opt$bootsamples * (length(xis))), dim = c(opt$bootsamples, 2 * length(xis)))
bsamplespst <- array(rep(NA, opt$bootsamples * (2 * length(xis))), dim = c(opt$bootsamples, 2 * length(xis)))


# fill up containers: 
# beta: take continuum limit for beta(xi) as well
bsamplescontlimit[, seq(1, length(xis))] <- resultslist[["plaqren"]]
bsamplescontlimit[, seq(length(xis) + 1, 2 * length(xis))] <- resultslist[["xiphys"]]^2
bsamplescontlimitbeta[, seq(length(xis) + 1, 2 * length(xis))] <- resultslist[["xiphys"]]^2
bsamplescontlimitbeta[, seq(1, length(xis))] <- resultslist[["intercepts"]][, seq(1, length(xis))]
if (opt$xisingle) {
    bsamplescontlimit[, seq(length(xis) + 1, 2 * length(xis))] <- resultslist[["xiphys"]]
    bsamplescontlimitbeta[, seq(length(xis) + 1, 2 * length(xis))] <- resultslist[["xiphys"]]
}
bsamplespst[, seq(1, length(xis))] <- resultslist[["plaqrenst"]]
bsamplespst[, seq(length(xis) + 1, 2 * length(xis))] <- resultslist[["xiphys"]]^2


if (result$dbeta[1] < 1e-8) {
    bsamplescontlimitbeta[, 1] <- parametric.bootstrap(opt$bootsamples, c(opt$beta), c(1e-4))
}



namesave <- sprintf("%s/resultsrenormalizationbetachosen%s.csv", opt$respath, endnamewrite)
write.table(result, file=namesave, col.names = TRUE, row.names = FALSE, append = FALSE)

resultslist$githash <- githash
namesave <- sprintf("%s/listresultsrenormalizationbetachosen%s.RData", opt$respath, endnamewrite)
saveRDS(resultslist, file=namesave)



# for each polynomial:
# fitplaq: xi=xi_ren, p=p_interpolated -> xi and beta renormalized
# for each: do fit to continuum limit, plot, add results to list and table
# na.omit removes entire row containing na

        listfits <- list(fits = list(), lowlimfit = c(), uplimfit = c())
        resultspolynomial <- data.frame(
            degree = NA, lim = NA, chi = NA,
            p = NA, type = NA, limplot = NA, dlimplot = NA, 
            lowlimfit = NA, uplimfit = NA
        )

listnames <- paste0(rep(c("plaq", "beta", "p_st"), each=5), seq(1, 5))
indexlim <- 1
pdf(sprintf("%s/plotsbetachosen%scont%d.pdf", opt$respath, endnamewrite, opt$indexfitcontlim), title = "")
for (lowlimfit in seq(opt$indexfitcontlim, length(xis) - 2)) {
    for (uplimfit in seq(max(length(xis) - 2, lowlimfit + 2), length(xis))) {
        print(paste("lowlim", lowlimfit, "uplim", uplimfit))
        fitspolynomial <- list()

maskfitcontlim <- seq(lowlimfit, uplimfit)
i <- 1
xirenfit <- result$xisimple^2
if(opt$xisingle) {
xirenfit <- result$xisimple
}

for (fun in c(fnlin, fnpar, fncub, fnqar, fnqin)) {
            if (uplimfit - lowlimfit > i) {
    print(paste("Doing fits of polynomial degree", i))

    # xi and beta renorm
    # print(attributes(na.omit(bsamplescontlimit))$na.action)
    err <- apply(bsamplescontlimit, 2, sd, na.rm=T)
    fitplaq <- try(bootstrap.nlsfit(fun, rep(1, i + 1),
        x = xirenfit, y = result$psimple, bsamples = bsamplescontlimit,
        mask = maskfitcontlim, na.rm=TRUE, 
                            dy = err[1:length(xirenfit)], 
                            dx = err[(1:length(xirenfit)) + length(xirenfit)]
    ))

        fitspolynomial[[listnames[i]]] <- fitplaq
    if (!inherits(fitplaq, "try-error")) {

      
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

        resultspolynomial <- rbind(
            resultspolynomial,
            data.frame(
                degree = i, lim = tex.catwitherror(fitplaq$t0[1],
                    fitplaq$se[1],
                    digits = 2, with.dollar = FALSE
                ),
                chi = fitplaq$chi / fitplaq$dof, p = fitplaq$Qval,
                type = "plaq", limplot = fitplaq$t0[1], dlimplot = fitplaq$se[1],
                        lowlimfit = lowlimfit, uplimfit = uplimfit
            )
        )
    }

    # beta cont limit
    err <- apply(bsamplescontlimitbeta, 2, sd, na.rm=T)
    fitbeta <- try(bootstrap.nlsfit(fun, rep(0.1, i + 1),
        x = xirenfit, y = result$beta, bsamples = bsamplescontlimitbeta,
        mask = maskfitcontlim, na.rm=TRUE, 
                            dy = err[1:length(xirenfit)], 
                            dx = err[(1:length(xirenfit)) + length(xirenfit)]
    ))
    fitspolynomial[[listnames[5 + i]]] <- fitbeta

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

    # pst cont limit, only plot, do not save
    err <- apply(bsamplespst, 2, sd, na.rm=T)
    fitpst <- try(bootstrap.nlsfit(fun, rep(0.1, i + 1),
        x = xirenfit, y = result$pst, bsamples = bsamplespst,
        mask = maskfitcontlim, na.rm=TRUE, 
                            dy = err[1:length(xirenfit)], 
                            dx = err[(1:length(xirenfit)) + length(xirenfit)]
    ))

    try(plot(fitpst,
        main = sprintf(
            "continuum limit pst: %f + /-%f, chi = %f, p = %f,\ndegree of polynomial:%d",
            fitpst$t0[1], fitpst$se[1],
            fitpst$chi / fitpst$dof, fitpst$Qval, i
        ),
        plot.range = c(-0.2, 1.2),
        ylim = c(0.6, 1.1),
        xlab = "xi_ren^2", ylab = "P_st", xaxs = "i", xlim = c(0, 1.05)
    ))

    try(resultspolynomial <- rbind(
        resultspolynomial,
        data.frame(
            degree = i, lim = tex.catwitherror(fitpst$t0[1],
                fitpst$se[1],
                digits = 2, with.dollar = FALSE
            ),
            chi = fitpst$chi / fitpst$dof, p = fitpst$Qval,
            type = "p_st", limplot = fitpst$t0[1], dlimplot = fitpst$se[1],
                        lowlimfit = lowlimfit, uplimfit = uplimfit
        )
    ))

    fitspolynomial[[listnames[10 + i]]] <- fitpst


    ## print if there were any errors in fitting
    if (inherits(fitplaq, "try-error")) {
        print(paste("There was an error with fitplaq for polynomial degree", i))
    }
    if (inherits(fitbeta, "try-error")) {
        print(paste("There was an error with fitbeta for polynomial degree", i))
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
        resultspolynomial <- resultspolynomial[-1, ]
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

listfits$githash <- githash
listfits$nalist <- resultslist$nalist

namepol <- sprintf("%s/polynomialbetachosen%scont%d.csv", opt$respath, endnamewrite, opt$indexfitcontlim)
write.table(resultspolynomial, namepol, col.names = TRUE, row.names = FALSE)
print(resultspolynomial)

namesave <- sprintf("%s/listpolynomialrenormalizationbetachosen%scont%d.RData", opt$respath, endnamewrite, opt$indexfitcontlim)
saveRDS(listfits, file=namesave)


