library("hadron")
library(optparse)

## goal of this script:
## Combine several different continuum limits into one with an AIC criterion

# plot cdf

AICweight <- function(chi, npar, nmeas) exp(-0.5 * (chi + 1 / chi + 2 * npar - nmeas))

removeunphysical <- function(cf, reldevmax = 10, removelimit = 50, verbose = T, fitname = "") {
    reldev <- abs((cf$t[, 1] - cf$t0[1]) / cf$t0[1])
    weirdfitmask <- reldev > reldevmax
    # if (sum(weirdfitmask, na.rm = T) > removelimit) stop(paste("removing this many samples does not make sense, check your assumptions for fit", fitname))
    if (sum(weirdfitmask, na.rm = T) > removelimit) message("removing this many samples does not make sense, check your assumptions for fit", fitname)
    if (verbose & sum(weirdfitmask, na.rm = T) > 1) {
        print(paste(sum(weirdfitmask, na.rm = T), "removed in", fitname))
        cf$bsamples[weirdfitmask, ] <- NA
        cf$t[weirdfitmask, ] <- NA
        cf$info.boot[weirdfitmask] <- -2
        sum(is.na(cf$t[, 1]))
        cf$se <- apply(X = cf$t, MARGIN = 2, FUN = sd, na.rm = T)
    }
    return(invisible(cf))
}


option_list <- list(
    make_option(c("-t", "--type"),
        type = "character",
        default = "plaq",
        help = "one of plaq, plaqinter, beta, betainter, p_st, p_stinter, plaqsmall, ratio [default %default]"
    ),
    make_option(c("-m", "--mode"),
        type = "character",
        default = "all",
        help = "average over these degrees/fit ranges [default %default]"
    ),
    make_option(c("-p", "--path"),
        type = "character",
        default = "./plotstikz/",
        help = "path to resultfiles [default %default]"
    ),
    make_option(c("--reldevmax"),
        type = "integer",
        default = 10000,
        help = "maximum relative deviation to mean cont limit result in bootstrap samples [default %default]"
    ),
    make_option(c("-v", "--verbose"),
        action = "store_true", default = FALSE,
        help = "print more information during process [default %default]"
    ),
    make_option(c("--xiinter"),
        action = "store_true", default = FALSE,
        help = "assume the interpolated xi_ren was used also for the directly chosen measurements [default %default]"
    ),
    make_option(c("--myfunctions"),
        type = "character",
        default = "/hiskp4/gross/masterthesis/analyse/code/U1_analyse_potential/",
        help = "path to where additional functions are stored,
            relative to folder where script is executed [default %default]"
    )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
source(paste(opt$myfunctions, "myfunctions.R", sep = ""))
nalim <- 30
reldevmax <- opt$reldevmax
removelimit <- 30

if (opt$xiinter) stopifnot(!(opt$type == "plaqinter" || opt$type == "betainter" || opt$type == "p_stinter"))

basename <- ""
elementname <- ""
if (opt$type == "plaq") {
    basename <- "listpolynomialrenormalizationbetachosen"
    tablename <- "polynomialbetachosen"
    elementname <- "plaq"
} else if (opt$type == "plaqinter") {
    basename <- "listpolynomialrenormalization"
    tablename <- "polynomial"
    elementname <- "plaq"
} else if (opt$type == "beta") {
    basename <- "listpolynomialrenormalizationbetachosen"
    tablename <- "polynomialbetachosen"
    elementname <- "beta"
} else if (opt$type == "betainter") {
    basename <- "listpolynomialrenormalization"
    tablename <- "polynomial"
    elementname <- "beta"
} else if (opt$type == "p_st") {
    basename <- "listpolynomialrenormalizationbetachosen"
    tablename <- "polynomialbetachosen"
    elementname <- "p_st"
} else if (opt$type == "p_stinter") {
    basename <- "listpolynomialrenormalization"
    tablename <- "polynomial"
    elementname <- "p_st"
} else if (opt$type == "plaqsmall") {
    basename <- "listpolynomialrenormalizationbetachosenL3"
    tablename <- "polynomialbetachosen"
    elementname <- "plaqsmall"
} else if (opt$type == "ratio") {
    basename <- "listpolynomialrenormalizationbetachosenL3"
    tablename <- "polynomialbetachosen"
    elementname <- "ratio"
} else {
    stop(paste("incorrect type given, you gave", opt$type))
}


degrees <- c()
uplim <- c()
lowlim <- c()

if (opt$mode == "all") {
    degrees <- c(3, 1, 1, 2, 1, 1, 2, 1, 1, 1)
    uplim <- c(10, 8, 8, 8, 9, 9, 9, 10, 10, 10)
    lowlim <- c(1, 5, 4, 3, 6, 5, 4, 7, 6, 5)
} else if (opt$mode == "xi0.20") {
    degrees <- c(1, 2)
    uplim <- c(8, 8)
    lowlim <- c(4, 3)
} else if (opt$mode == "xi0.19") {
    degrees <- c(1, 2)
    uplim <- c(9, 9)
    lowlim <- c(6, 4)
} else if (opt$mode == "xi0.18") {
    degrees <- c(1, 1, 2)
    uplim <- c(10, 10, 10)
    lowlim <- c(7, 6, 5)
} else {
    stop(paste("incorrect mode given, you gave", opt$mode))
}

res <- list(etp0b1.70 = list(), etp1b1.70 = list(), etp0b1.65 = list(), etp1b1.65 = list(), list(opt = opt, degrees = degrees, uplim = uplim, lowlim = lowlim))

restable <- data.frame(
    beta = c(), etp = c(), omit = c(), pot = c(), degree = c(), lowlim = c(), uplim = c(), contlim = c(), dcontlim = c(),
    AICweight = c(), errorweight = c(), chi = c(), npar = c(), nmeas = c(), isna = c()
)
averagedrestable <- data.frame(
    beta = c(), etp = c(),
    contlimAIC = c(), dcontlimAIC = c(), dcontlimspreadAIC = c(),
    contlimerror = c(), dcontlimerror = c(), dcontlimspreaderror = c(),
    contlimunweighted = c(), dcontlimunweighted = c(), dcontlimspreadunweighted = c(),
    naratio = c(), nainbs = c()
)

filename <- sprintf("%s/nalistmode%s%s.RData", opt$path, opt$mode, ifelse(opt$xiinter, "xiinter", ""))
if (!file.exists(filename)) indexNAlist <- c()
if (file.exists(filename)) indexNAlist <- readRDS(sprintf("%s/nalistmode%s%s.RData", opt$path, opt$mode, ifelse(opt$xiinter, "xiinter", "")))
localnalist <- c()
problemnames <- c()

for (beta in c(1.70, 1.65)) {
    for (etp in c(0, 1)) {
        myweight <- c()
        myerrorweight <- c()
        bs <- array(NA, dim = c(500, 4 * length(degrees)))
        mymean <- c()
        mysd <- c()
        index <- 1
        mynas <- c()
        mychi <- c()
        mynmeas <- c()
        mynpar <- c()
        for (omit in c(0, 1)) {
            for (pot in c("normal", "sideways")) {
                filename <- sprintf("beta%fomit%dxiconstllxi1llr00fl0.30aicscaletauintetp%d%scont1maxiter300", beta, omit, etp, ifelse(opt$xiinter, "xiinter", ""))
                if (opt$type == "plaqinter" || opt$type == "betainter") filename <- sprintf("beta%fomit%dxiconstllxi1llr00fl0.30aicscaletauintetp%d", beta, omit, etp)
                data <- readRDS(sprintf("%s/%s%s%s.RData", opt$path, basename, pot, filename))
                table <- read.table(sprintf("%s/%s%s%s.csv", opt$path, tablename, pot, filename), header = T)
                table <- table[table$type == elementname, ]
                for (lim_i in seq_along(lowlim)) {
                    # print(data$lowlim==lowlim[lim_i] & data$uplim==uplim[lim_i])
                    fitmask <- which(data$lowlim == lowlim[lim_i] & data$uplim == uplim[lim_i])
                    fitname <- paste0(elementname, degrees[lim_i])
                    if (inherits(data$fits[[fitmask]][[fitname]], "try-error")) {
                        # print(names(data$fits[[fitmask]][[fitname]]))
                        # stop(paste("fit does not exist", beta, "etp", etp, "omit", omit, "pot", pot, lim_i))
                        print(paste("fit does not exist", beta, "etp", etp, "omit", omit, "pot", pot, lim_i))
                        bs[, index] <- rep(NA, 500)
                    } else {
                        if (data$fits[[fitmask]][[fitname]]$t0[1] < 0.1 || data$fits[[fitmask]][[fitname]]$t0[1] > 5) stop(paste("this fit has an unphysical mean result", beta, etp, omit, pot, lim_i))
                        data$fits[[fitmask]][[fitname]] <- removeunphysical(data$fits[[fitmask]][[fitname]],
                            reldevmax = reldevmax, removelimit = removelimit,
                            fitname = paste(omit, pot, etp, lim_i, elementname), verbose = opt$verbose
                        )
                        mymean[index] <- data$fits[[fitmask]][[fitname]]$t0[1]
                        mysd[index] <- data$fits[[fitmask]][[fitname]]$se[1]
                        myweight[index] <- AICweight(
                            chi = data$fits[[fitmask]][[fitname]]$chisqr,
                            npar = length(data$fits[[fitmask]][[fitname]]$par.guess),
                            nmeas = uplim[lim_i] - lowlim[lim_i]
                        )
                        myerrorweight[index] <- 1 / mysd[index]^2
                        bs[, index] <- data$fits[[fitmask]][[fitname]]$t[, 1]
                        mychi[index] <- data$fits[[fitmask]][[fitname]]$chisqr
                        mynpar[index] <- length(data$fits[[fitmask]][[fitname]]$par.guess)
                        mynmeas[index] <- uplim[lim_i] - lowlim[lim_i]
                    }
                    mynas[index] <- sum(is.na(bs[, index]))
                    if (sum(is.na(bs[, index])) > nalim || index %in% indexNAlist) {
                        if (opt$verbose) message(paste("remove fit", lim_i, "from", etp, pot, omit))
                        myweight[index] <- 0
                        myerrorweight[index] <- 0
                        mymean[index] <- 0
                        mysd[index] <- ifelse(is.na(mysd[index]), 0, mysd[index])
                        bs[, index] <- rep(0, 500)
                        localnalist <- append(localnalist, index)
                        problemnames <- append(problemnames, paste(beta, "etp", etp, "omit", omit, "pot", pot, "lowlim", lowlim[lim_i], "uplim", uplim[lim_i], "degree", degrees[lim_i]))
                    }
                    index <- index + 1
                }
            }
        }
        # AIC
        myweight <- myweight / sum(myweight)
        meanAIC <- sum(mymean * myweight)
        bsAIC <- apply(MARGIN = 1, FUN = sum, X = sweep(x = bs, MARGIN = 2, STATS = myweight, FUN = "*"))
        if (opt$verbose) print(sum(is.na(bsAIC)))
        sdAIC <- sd(bsAIC, na.rm = T)
        sdAICspread <- sqrt(sdAIC^2 + sum((mymean - meanAIC)^2 * ifelse(myweight == 0, 0, 1)) / length(mymean[mymean != 0]))

        # weighted with 1/sd^2
        myerrorweight <- myerrorweight / sum(myerrorweight) # normalized 1/sd^2
        meanerror <- sum(mymean * myerrorweight)
        bserror <- apply(MARGIN = 1, FUN = sum, X = sweep(x = bs, MARGIN = 2, STATS = myerrorweight, FUN = "*"))
        sderror <- sd(bserror, na.rm = T)
        sderrorspread <- sqrt(sderror^2 + sum((mymean - meanerror)^2 * ifelse(myweight == 0, 0, 1)) / length(mymean[mymean != 0]))

        # unweighted
        meanunweighted <- sum(mymean) / length(mymean[mymean != 0])
        bsunweighted <- apply(bs, MARGIN = 1, FUN = mean)
        sdunweighted <- sd(bsunweighted, na.rm = T)
        sdunweightedspread <- sqrt(sdunweighted^2 + sum((mymean - meanunweighted)^2 * ifelse(myweight == 0, 0, 1)) / length(mymean[mymean != 0]))

        # save results
        averagedrestable <- rbind(averagedrestable, data.frame(
            beta = beta, etp = etp,
            contlimAIC = meanAIC, dcontlimAIC = sdAIC, dcontlimspreadAIC = sdAICspread,
            contlimerror = meanerror, dcontlimerror = sderror, dcontlimspreaderror = sderrorspread,
            contlimunweighted = meanunweighted, dcontlimunweighted = sdunweighted, dcontlimspreadunweighted = sdunweightedspread,
            naratio = sum(mynas) / (4 * length(degrees) * 500), nainbs = sum(is.na(bsAIC))
        ))
        tmpframe <- data.frame(
            beta = beta, etp = etp,
            omit = rep(c(0, 1), each = 2 * length(degrees)),
            pot = rep(c("normal", "sideways", "normal", "sideways"), each = length(degrees)),
            degree = rep(degrees, 4),
            lowlim = rep(lowlim, 4),
            uplim = rep(uplim, 4),
            contlim = mymean,
            dcontlim = mysd,
            AICweight = myweight,
            errorweight = myerrorweight,
            chi = mychi,
            npar = mynpar,
            nmeas = mynmeas,
            isna = mynas
        )
        restable <- rbind(restable, tmpframe)
        res[[paste0("etp", etp, "b", beta)]] <- list(data <- tmpframe,
            bsAIC = bsAIC, bserror = bserror, bsunweighted = bsunweighted,
            bsAICspread = bsAIC + (bsAIC - meanAIC) * (sdAICspread / sdAIC -1),
            bserrorspread = bserror + (bserror - meanerror) * (sderrorspread / sderror -1),
            bsunweightedspread = bsunweighted + (bsunweighted - meanunweighted) * (sdunweightedspread / sdunweighted -1)
        )
    }
}

if (opt$verbose) print(restable)
averagedrestable

savename <- sprintf("%s/contlimittype%smode%s%s", opt$path, opt$type, opt$mode, ifelse(opt$xiinter, "xiinter", ""))

write.table(restable,
    sprintf("%snotaveraged.csv", savename),
    row.names = T
)
write.table(averagedrestable,
    sprintf("%saveraged.csv", savename),
    row.names = T
)
saveRDS(
    res,
    sprintf("%saveraged.RData", savename)
)

saveRDS(localnalist, sprintf("%s/nalistmode%s%s.RData", opt$path, opt$mode, ifelse(opt$xiinter, "xiinter", "")))


pdf(sprintf("%s/limits%s%s%s.pdf", opt$path, opt$type, opt$mode, ifelse(opt$xiinter, "xiinter", "")), title = "")


for (beta in c(1.70, 1.65)) {
    for (etp in c(0, 1)) {
        mymask <- restable$beta == beta & restable$etp == etp & restable$contlim != 0
        plotwitherror(
            x = seq(1, length(which(mymask))), y = restable$contlim[mymask], dy = c(restable$dcontlim[mymask]),
            xlab = "index", ylab = "contlim",
            main = sprintf(
                "etp %d beta %.2f\nremoved fits %d\nhighest AICweight %.3f", etp, beta,
                sum(restable$beta == beta & restable$etp == etp & restable$contlim == 0),
                max(restable$AICweight[mymask])
            ),
            ylim = range(restable$contlim[mymask] + c(restable$dcontlim[mymask]), restable$contlim[mymask] - c(restable$dcontlim[mymask]))
        )
        mymask <- averagedrestable$beta == beta & averagedrestable$etp == etp
        abline(h = averagedrestable$contlimAIC[mymask], col = 2)
        abline(h = averagedrestable$contlimAIC[mymask] + c(-1, 1) * averagedrestable$dcontlimAIC[mymask], lty = 2, col = 2)
        # abline(h = averagedrestable$contlimAIC[mymask] + c(-1, 1) * averagedrestable$dcontlimspreadAIC[mymask], lty = 3, col = 2, lwd = 3)
        abline(h = averagedrestable$contlimerror[mymask], col = "blue")
        abline(h = averagedrestable$contlimerror[mymask] + c(-1, 1) * averagedrestable$dcontlimerror[mymask], lty = 2, col = "blue")
        # abline(h = averagedrestable$contlimerror[mymask] + c(-1, 1) * averagedrestable$dcontlimspreaderror[mymask], lty = 3, col = "blue", lwd = 3)
        abline(h = averagedrestable$contlimunweighted[mymask], col = "darkgreen")
        abline(h = averagedrestable$contlimunweighted[mymask] + c(-1, 1) * averagedrestable$dcontlimunweighted[mymask], lty = 2, col = "darkgreen")
        # abline(h = averagedrestable$contlimunweighted[mymask] + c(-1, 1) * averagedrestable$dcontlimspreadunweighted[mymask], lty = 3, col = "darkgreen", lwd = 3)
        legend(
            x = "top", pch = c(1, NA, NA, NA, NA, NA, NA), col = c(1, 1, 1, 1, 2, "blue", "darkgreen"), lty = c(NA, 1, 2, 3, 1, 1, 1),
            legend = c("data", "mean", "sd", "sd + spread", "AIC", "error", "unweighted"), ncol = 4
        )
    }
}


print("warnings:")
print(warnings())
print("problematic fits")
print(problemnames)
