library(hadron)
library(optparse)

if (TRUE) {
    # set option list
option_list <- list(
    make_option(c("-b", "--beta"), type = "double", default = 1.7,
    help = "beta-Parameter of simulation [default %default]"),
    make_option(c("--skip"), type = "integer", default = 1000,
    help = "how many lines are skipped when reading in [default %default]"),

    make_option(c("-s", "--bootsamples"), type = "integer", default = 500,
    help = "how many bootstrapsamples should be drawn [default %default]"),
    make_option(c("--maxrows"), type = "integer", default = -1,
    help = "Maximum of configurations that are read in, -1=all [default %default]"),
    make_option(c("-r", "--Ns"), type = "integer", default = 16,
    help = "L of lattice (spatial extent) [default %default]"),

    make_option(c("-t", "--Nt"), type = "integer", default = 16,
    help = "T of lattice (temporal extent) [default %default]"),
    make_option(c("-n", "--nape"), type = "integer", default = 0,
    help = "Number of APE-smears that were done [default %default]"),

    make_option(c("-l", "--lowerboundmeff"), type = "integer", default = 0,
    help = "meff is determined in region
            [lowerboundmeff, L/2 - 2] [default %default]"),
    make_option(c("-j", "--job"), type = "integer", default = 0,
    help = "qbig job number, written into results [default %default]"),

    make_option(c("-o", "--omit"), type = "integer", default = 0,
    help = "how many of the highest points should be omitted
            in calculating the potentials? [default %default]"),
    make_option(c("-a", "--alpha"), type = "double", default = 1.0,
    help = "alpha used in APE-smearing [default %default]"),

    make_option(c("--betaone"), type = "double", default = 0,
    help = "input beta at corresponding xi = 1 [default %default]"),
    make_option(c("-x", "--xi"), type = "double", default = 0,
    help = "xi used in lattice, only used if xidiff = TRUE, else
            xi is assumed to be L/T [default %default]"),

    make_option(c("--fraction"), type = "double", default = 0.5,
    help = "Maximal extent of the Wilson-Loops as a fraction of the lattice size [default %default]"),
    make_option(c("--zerooffset"), type = "integer", default = 0,
    help = "Offset for selecting the correct loops
            if W(x = 0, y = 0) was measured, set to 1 [default %default]"),

    make_option(c("-e", "--every"), type = "integer", default = 1,
    help = "only reads every eth line [default %default]"),
    make_option(c("--bootl"), type = "integer", default = 2,
    help = "block length in bootstrapping configurations [default %default]"),

    make_option(c("--nsave"), type = "integer", default = 0,
    help = "steps between saved configs [default %default]"),
    make_option(c("--crzero"), type = "double", default = -1.65,
    help = "c used for determining r_0 [default %default]"),
    make_option(c("--lowlim"), type = "integer", default = 1,
    help = "points lower than this are not used
            for determining xi [default %default]"),
    make_option(c("--lowlimpot"), type = "integer", default = -1,
    help = "points lower than this are not used
            for determining r_0, if negative, the same as lowlim [default %default]"),

    make_option(c("--aic"), action = "store_true", default = FALSE,
    help = "if true, effective masses are determined with weighting
        according to the akaike information criterion [default %default]"),
    make_option(c("--scaletauint"), action = "store_true", default = FALSE,
    help = "if true, errors and bootstrapsamples for the correlator are rescaled
        such that the effects of autocorrelation are taken into account [default %default]"),
    make_option(c("--errortotpot"), action = "store_true", default = FALSE,
    help = "if true, potential bootstrap samples are rescaled before analysis
        to take total error into account [default %default]"),
    make_option(c("--analyse"), action = "store_true", default = FALSE,
    help = "if true, correlators and effective masses
            are determined from provided table [default %default]"),
    make_option(c("--drawbootstrap"), action = "store_true", default = FALSE,
    help = "Draw bootstrap samples and determine effective mass
            before doing effective mass fits [default %default]"),

    make_option(c("--dofit"), action = "store_true", default = FALSE,
    help = "if true, potentials and xi are calculated [default %default]"),
    make_option(c("--smearing"), action = "store_true", default = FALSE,
    help = "are the calculations done based on
            smeared lattices? (changes filenames) [default %default]"),

    make_option(c("--xidiff"), action = "store_true", default = FALSE,
    help = "Is xi different to L/T? [default %default]"),
    make_option(c("--usecov"), action = "store_true", default = FALSE,
    help = "Use covariance for determining effective mass [default %default]"),
    make_option(c("--plotonlymeff"), action = "store_true", default = FALSE,
    help = "The plots of the fitted meff are
            plotted seperately [default %default]"),

    make_option(c("--respath"), type = "character", default = "",
    help = "path to where the resultfiles are stored [default %default]"),
    make_option(c("--plotpath"), type = "character", default = "",
    help = "path to where the plots are stored [default %default]"),
    make_option(c("--effmasstype"), type = "character", default = "log",
    help = "type of effective mass [default %default]"),
    make_option(c("--extra"), type = "character", default = "",
    help = "extra string for ending of fit analysis [default %default]"),

    make_option(c("--myfunctions"), type = "character",
    default = "/hiskp4/gross/masterthesis/analyse/code/U1_analyse_potential/",
#~     make_option(c("--myfunctions"), type = "character", default = "myfunctions.R",
    help = "path to where additional functions are stored,
            relative to folder where script is executed [default %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
}

if (TRUE) {
# set some constants
source(paste(opt$myfunctions, "myfunctions.R", sep = ""))
githash <- printgitcommit(opt$myfunctions)
beta <- opt$beta
skip <- opt$skip
bootsamples <- opt$bootsamples
Ns <- opt$Ns
Nt <- opt$Nt
xi <- Ns / Nt
if (opt$xidiff) {
    xi <- opt$xi
}

nape <- opt$nape
alpha <- opt$alpha

analyse <- opt$analyse
dofit <- opt$dofit

if (opt$scaletauint) {
        opt$usecov <- TRUE
        opt$bootl <- 1
}

if (opt$lowlimpot < 0) opt$lowlimpot <- opt$lowlim

endinganalysis <- sprintf("Nt%dNs%dbeta%fxi%fbootl%dusecov%d", Nt, Ns, beta, xi, opt$bootl, opt$usecov)
if (opt$smearing) {
        endinganalysis <- sprintf("Nt%dNs%dbeta%fxi%fnape%dalpha%fbootl%dusecov%d",
        Nt, Ns, beta, xi, nape, alpha, opt$bootl, opt$usecov)
}

endingdofit <- sprintf("Nt%dNs%dbeta%fxi%fbs%domit%dllxi%dllr0%d",
        Nt, Ns, beta, xi, opt$bootsamples, opt$omit, opt$lowlim, opt$lowlimpot)
if (opt$smearing) {
        endingdofit <- sprintf("Nt%dNs%dbeta%fxi%fnape%dalpha%fbs%domit%dllxi%dllr0%d",
        Nt, Ns, beta, xi, nape, alpha, opt$bootsamples, opt$omit, opt$lowlim, opt$lowlimpot)
}

if(opt$crzero != -1.65) endingdofit <- sprintf("%sc%.2f", endingdofit, opt$crzero)

if(opt$aic){
        endinganalysis <- sprintf("%saic", endinganalysis)
        endingdofit <- sprintf("%saic", endingdofit)
}

if(opt$scaletauint){
        endinganalysis <- sprintf("%sscaletauint", endinganalysis)
        endingdofit <- sprintf("%sscaletauintetp%d", endingdofit, opt$errortotpot)
}

endingdofit <- paste0(endingdofit, opt$extra)

}

if (analyse) {
# set names for plot, tables, open, lists for saving results
filenameforplots <- sprintf("%splotsmeffsideways%s.pdf", opt$plotpath, endinganalysis)
if (!opt$drawbootstrap) filenameforplots <- sprintf("%splotsmeffsidewaysnc%s.pdf", opt$plotpath, endinganalysis)
pdf(file = filenameforplots, title = "")
listfits <- list()
listtauint <- list()
potential <- data.frame(R = NA, m = NA, dm = NA, space = NA, p = NA, chi = NA)
if (opt$aic) {
#     potential <- data.frame(R = NA, m = NA, dm = NA, space = NA, setot = NA, bootstat = NA, meanboot = NA, bias = NA, m16 = NA, m84 = NA, bootsys = NA)
    potential <- data.frame(R = NA, space = NA, median = NA, errtotmean = NA, m16 = NA, m84 = NA,
                        bootmedian = NA, booterrtot = NA, bootstat = NA, bootsys = NA, bias = NA)
}

namelistbounds <- sprintf(
            "%stableanalysisrotatedNs%dNt%dbeta%fxi%f.csv",
            opt$respath, Ns, Nt, beta, xi)
if (opt$smearing) {
    namelistbounds <- sprintf(
            "%stableanalysisrotatedNs%dNt%dbeta%fxi%fnape%dalpha%f.csv",
            opt$respath, Ns, Nt, beta, xi, nape, alpha)
}
if (file.exists(namelistbounds)) {
    listbounds <- read.table(namelistbounds, sep = ",", header = TRUE)
} else {
    # if no table with boundaries can be found, create a dummy table
    listbounds <- cbind(data.frame(yt = c(seq(1, Ns / 2), seq(1, Nt / 2))),
                    data.frame(spacial = c(rep(TRUE, Ns / 2), rep(FALSE, Nt / 2))),
                    data.frame(lower = rep(opt$lowerboundmeff, Ns / 2 + Nt / 2)),
                    data.frame(upper = rep(Ns / 2 - 2, Ns / 2 + Nt / 2)))
}
print(listbounds)
negatives <- c()

## if the correlators and effective masses are not determined in-place, read in previous results
if (!opt$drawbootstrap) {
    readinres <- readRDS(sprintf("%slistmeffsideways%s.RData", opt$plotpath, endinganalysis))
}

# determine Wilson Loop correlators and effective masses
# for all possible combinations Meff(x0/a_s) (use W(x0, t=0, y0+1)/W(x0, t=0, y0))
# and Meff(x0/a_t) (use W(x0, t0+1, y=0)/W(x0, t0, y=0))
# try method as in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.63.074501
# x here as z there, y here as x there
# measures a_s V(r / a_s)
# also determine number of correlators smaller than zero per distance
for (y in seq(1, Ns / 2, 1)) {
    # first spatial loops: determine expectation values of Loops,
    # plot, determine effective masses, set boundaries, determine
    # plateau and plot. Savev results of m_eff to list
    # also save uwerr = autocorrelation time of non-bootstrapped loop
    # list has to have an entry everywhere
    listresults <- list(FALSE, FALSE, FALSE, FALSE)
    uwerrresults <- FALSE
    if (y > 0) {

    # loops
    print(y)
    filename <- sprintf(
        "result2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fcoarsedistance",
        Nt, Ns, beta, xi, nape, alpha)

    title <- sprintf("beta = %.2f coarse y = %d skipped %d",
                beta, y, skip * opt$nsave)
    if (opt$drawbootstrap) {
        WL <- calcplotWloopsideways(file = filename,
                path = opt$respath, skip = skip, Ns = Ns, yt = y,
                title = title, bootsamples = bootsamples,
                zerooffset = opt$zerooffset, nsave = opt$nsave,
                every = opt$every, l = opt$bootl, fraction = opt$fraction,
                maxrows = opt$maxrows)

        uwerrresults <- uwerr.cf(WL)

#         print(names(uwerrresults$uwcf))
#         "value"   "dvalue"  "ddvalue" "tauint"  "dtauint" "t"

        if (opt$scaletauint) {
        ## multiply errors by 2*tauint
            WL$tsboot.se <- WL$tsboot.se * 2 * uwerrresults$uwcf$tauint
        ## if bootstrap below mean, subtract (2*tauint-1)*difference
        ## if bootstrap above mean, add (2*tauint-1)*difference
        ## add/subtract ist automatically done by sign of difference
        ## before multiplication: difference a
        ## after mutliplication: difference 2*tauint*a
            # ## How to deal with dtauint?
            WL$cf.tsboot$t <- t(apply(X=WL$cf.tsboot$t, MARGIN=1, FUN=function(x, tauint, mean) x + (x - mean) * (2 * tauint - 1),
                tauint = uwerrresults$uwcf$tauint, mean = WL$cf0))
        }

        WL <- bootstrap.effectivemass(WL, type = opt$effmasstype)
    } else {
        WL <- readinres[[y]][[1]]
        if(opt$aic) WL <- readinres[[y]][[1]]$effmass
    }

    # m_eff
    if(!opt$aic) {
    t1 <- listbounds$lower[listbounds$spacial == TRUE & listbounds$yt == y]
    t2 <- listbounds$upper[listbounds$spacial == TRUE & listbounds$yt == y]

    WL.effmasslist <- deteffmass(WL.effmass = WL, yt = y,
            t1 = t1, t2 = t2, isspatial = 1, usecov=opt$usecov)
    # plot results, save results off meff in
    # short and long form (with and without bootstrapsamples
    if (WL.effmasslist[[2]][[2]] != 0) {
        listresults <- list(WL.effmasslist[[1]], y, FALSE, uwerrresults)
        potential <- rbind(potential, WL.effmasslist[[3]])
    }
    # coarse distance = spatial potential
    names(listresults) <- c("effmass", "yt", "coarse", "uwerr")

    title <- sprintf("%s, %d configs\n",
            title, (length(WL$cf$cf[, 1]) + skip) * opt$nsave)
    try(plot(WL.effmasslist[[1]], xlab = "x / a_s", ylab = "Meff",
            main = sprintf("%s, t1 = %d, t2 = %d, p = %f",
            title, t1, t2, WL.effmasslist[[1]]$effmassfit$Qval)))
    try(plot(WL.effmasslist[[1]], xlab = "x / a_s", ylab = "Meff",
            main = sprintf("%s, t1 = %d, t2 = %d, p = %f, zoomed",
            title, t1, t2, WL.effmasslist[[1]]$effmassfit$Qval),
            ylim = WL.effmasslist[[2]]))
    }

    if (opt$aic) {
        effmass <- deteffmassaic(WL)
        plot(effmass$effmass, xlab="x/a_s", ylab="meff", main=title)
        plot(effmass$effmass, ylim=c(effmass$effmassfit$t0 - 5*effmass$effmassfit$se, effmass$effmassfit$t0 + 5*effmass$effmassfit$se),
                xlab="x/a_s", ylab="meff", main=title)
        lines(x=c(-10, 2*Nt), y=rep(effmass$effmassfit$t0, 2))
        lines(x=c(-10, 2*Nt), y=rep(effmass$effmassfit$t0 + effmass$effmassfit$se, 2))
        lines(x=c(-10, 2*Nt), y=rep(effmass$effmassfit$t0 - effmass$effmassfit$se, 2))
        plot(effmass$effmass, ylim=c(effmass$boot$mean - 5*effmass$boot$errtot, effmass$boot$mean + 5*effmass$boot$errtot),
                xlab="x/a_s", ylab="meff", main=title)
        lines(x=c(-10, 2*Nt), y=rep(effmass$boot$mean, 2))
        lines(x=c(-10, 2*Nt), y=rep(effmass$boot$mean + effmass$boot$errtot, 2))
        lines(x=c(-10, 2*Nt), y=rep(effmass$boot$mean - effmass$boot$errtot, 2))
        listresults <- list(effmass, y, FALSE, uwerrresults)

        potentialnew <- data.frame(R = y, space = TRUE, median = effmass$effmassfit$t0, errtotmean = effmass$effmassfit$se,
                m16 = effmass$effmassfit$m16, m84 = effmass$effmassfit$m84,
                bootmedian = effmass$boot$mean, booterrtot = effmass$boot$errtot, bootstat = effmass$boot$errstat,
                bootsys = sqrt(effmass$boot$errtot^2 - effmass$boot$errstat^2),
                bias = effmass$effmassfit$t0 - effmass$boot$mean)
        potential <- rbind(potential, potentialnew)
    }
}
    listfits[[y]] <- listresults
    if (opt$drawbootstrap) listtauint[[y]] <- uwerrresults
    if (opt$drawbootstrap) negatives[y] <- sum(WL$cf0 < 0)
}
# same procedure as for y loop, here for temporal potential
# measured a_s V(r/a_t)
for (t in seq(1, Nt / 2, 1)) {
    listresults <- list(FALSE, FALSE, FALSE, FALSE)
    uwerrresults <- FALSE
    if (t < Nt / 2 + 1) {

#~     message("\nt = ", t)
    print(t)
    filename <- sprintf(
            "result2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%ffinedistance",
            Nt, Ns, beta, xi, nape, alpha)

    # loops
    title <- sprintf("beta = %.2f fine t = %d skipped %d",
                beta, t, skip * opt$nsave)

    if (opt$drawbootstrap) {
        WL <- calcplotWloopsideways(file = filename,
                path = opt$respath, skip = skip, Ns = Ns, yt = t, title = title,
                bootsamples = bootsamples, zerooffset = opt$zerooffset,
                nsave = opt$nsave, every = opt$every, l = opt$bootl, fraction = opt$fraction,
                maxrows = opt$maxrows)

        uwerrresults <- uwerr.cf(WL)

        if (opt$scaletauint) {
                ## multiply errors by 2*tauint
            WL$tsboot.se <- WL$tsboot.se * 2 * uwerrresults$uwcf$tauint
            ## if bootstrap below mean, subtract (2*tauint-1)*difference
            ## if bootstrap above mean, add (2*tauint-1)*difference
            ## add/subtract ist automatically done by sign of difference
            ## before multiplication: difference a
            ## after mutliplication: difference 2*tauint*a
            ## How to deal with dtauint?
            WL$cf.tsboot$t <- t(apply(X=WL$cf.tsboot$t, MARGIN=1, FUN=function(x, tauint, mean) x + (x - mean) * (2 * tauint - 1),
                tauint = uwerrresults$uwcf$tauint, mean = WL$cf0))
        }

        WL <- bootstrap.effectivemass(WL, type = opt$effmasstype)
    } else {
        WL <- readinres[[Ns / 2 + t]][[1]]
        if (opt$aic) WL <- readinres[[Ns / 2 + t]][[1]]$effmass
    }


    # m_eff
    if(!opt$aic) {
    t1 <- listbounds$lower[listbounds$spacial == FALSE & listbounds$yt == t]
    t2 <- listbounds$upper[listbounds$spacial == FALSE & listbounds$yt == t]

    WL.effmasslist <- deteffmass(WL.effmass = WL, yt = t,
            t1 = t1, t2 = t2, isspatial = 0, usecov=opt$usecov)
    if (WL.effmasslist[[2]][[2]] != 0) {
        listresults <- list(WL.effmasslist[[1]], t, TRUE, uwerrresults)
        potential <- rbind(potential, WL.effmasslist[[3]])
    }
    names(listresults) <- c("effmass", "yt", "fine", "uwerr")

    title <- sprintf("%s, %d configs\n",
            title, (length(WL$cf$cf[, 1]) + skip) * opt$nsave)
    try(plot(WL.effmasslist[[1]], xlab = "x / a_s", ylab = "Meff",
            main = sprintf("%s, t1 = %d, t2 = %d, p = %f",
            title, t1, t2, WL.effmasslist[[1]]$effmassfit$Qval)))
    try(plot(WL.effmasslist[[1]], xlab = "x / a_s", ylab = "Meff",
            main = sprintf("%s, t1 = %d, t2 = %d, p = %f",
            title, t1, t2, WL.effmasslist[[1]]$effmassfit$Qval),
            ylim = WL.effmasslist[[2]]))
    }
    if (opt$aic) {
        effmass <- deteffmassaic(WL)
        plot(effmass$effmass, xlab="x/a_s", ylab="meff", main=title)
        plot(effmass$effmass, ylim=c(effmass$effmassfit$t0 - 5*effmass$effmassfit$se, effmass$effmassfit$t0 + 5*effmass$effmassfit$se),
                xlab="x/a_s", ylab="meff", main=title)
        lines(x=c(-10, 2*Nt), y=rep(effmass$effmassfit$t0, 2))
        lines(x=c(-10, 2*Nt), y=rep(effmass$effmassfit$t0 + effmass$effmassfit$se, 2))
        lines(x=c(-10, 2*Nt), y=rep(effmass$effmassfit$t0 - effmass$effmassfit$se, 2))
        plot(effmass$effmass, ylim=c(effmass$boot$mean - 5*effmass$boot$errtot, effmass$boot$mean + 5*effmass$boot$errtot),
                xlab="x/a_s", ylab="meff", main=title)
        lines(x=c(-10, 2*Nt), y=rep(effmass$boot$mean, 2))
        lines(x=c(-10, 2*Nt), y=rep(effmass$boot$mean + effmass$boot$errtot, 2))
        lines(x=c(-10, 2*Nt), y=rep(effmass$boot$mean - effmass$boot$errtot, 2))
        listresults <- list(effmass, t, TRUE, uwerrresults)
        potentialnew <- data.frame(R = t, space = FALSE, median = effmass$effmassfit$t0, errtotmean = effmass$effmassfit$se,
                m16 = effmass$effmassfit$m16, m84 = effmass$effmassfit$m84,
                bootmedian = effmass$boot$mean, booterrtot = effmass$boot$errtot, bootstat = effmass$boot$errstat,
                bootsys = sqrt(effmass$boot$errtot^2 - effmass$boot$errstat^2),
                bias = effmass$effmassfit$t0 - effmass$boot$mean)
        potential <- rbind(potential, potentialnew)
    }
}
    listfits[[Ns / 2 + t]] <- listresults
    if (opt$drawbootstrap) listtauint[[Ns / 2 + t]] <- uwerrresults
    if (opt$drawbootstrap) negatives[Ns / 2 + t] <- sum(WL$cf0 < 0)
}
listfits[[Ns/2 + Nt/2 + 1]] <- githash

#write out results
potential <- na.omit(potential)
filenamepotential <- sprintf("%spotmeffsideways%s.csv", opt$plotpath, endinganalysis)
filenamelist <- sprintf("%slistmeffsideways%s.RData", opt$plotpath, endinganalysis)
filenameuwerr <- sprintf("%slistuwerrtauintsideways%s.RData",  opt$plotpath, endinganalysis)
filenamenegatives <- sprintf("%snegativessideways%s.csv",  opt$plotpath, endinganalysis)


write.table(potential, filenamepotential, row.names = FALSE)
print(filenamelist)
saveRDS(listfits, file = filenamelist)
if (opt$drawbootstrap) {
        saveRDS(listtauint, file = filenameuwerr)
        write.table(data.frame(x = c(seq(1, Ns / 2), seq(1, Nt / 2)), neg = negatives),
        file = filenamenegatives, row.names = FALSE)
}
}

if (dofit) {
# use data from the potential to fit the expected form to the fine=temporal
# and coarse=spatial potential,
# then determine xi_ren by interpolation, like in
# https://journals.aps.org/prd/abstract/10.1103/PhysRevD.63.074501
# then determine r_0 / a_s, determine the force, plot everything
# save results and bootstrapsamples of fits
# also try some other ratios to determine xi

# define functions for potential, force, force = -r^2 * derivation of potential
fnpot <- function (par, x, boot.r, ...) par[1] + par[2] * x + par[3] * log(x)
fnforce <- function (par, x, boot.r, ...) (- 1.0) * par[2] * x^2 - par[3] * x
fnforceerr <- function (parerr, x, boot.r, correlation, ...) {
    return (sqrt(parerr[2]^2 * x^4 + parerr[3]^2 * x^2 + 2 * x^3 * correlation[2, 3] * parerr[2] * parerr[3]))
}

# read in data
filenamelist <- sprintf(
            "%slistmeffsideways%s.RData",
            opt$plotpath, endinganalysis)

listmeff <- readRDS(file = filenamelist)



filenameforplots <- sprintf("%spotentialsideways%s.pdf", opt$plotpath, endingdofit)
pdf(file = filenameforplots, title = "")


# calculate plaquette as W(x = 1, y = 1, t = 0)
# with uwerr and with cf, to get bootstrapsamples

filename <- sprintf(
        "%sresult2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fcoarsedistance",
        opt$respath, Nt, Ns, beta, xi, nape, alpha)

alldata <- read.table(filename)
plaquettecolumn <- alldata[, (opt$zerooffset + Ns / 2) * opt$zerooffset + 1 + opt$zerooffset]
plaquettedata <- uwerrprimary(plaquettecolumn[seq(skip + 1, length(plaquettecolumn), opt$every)], S = 6, pl = TRUE)
nom <- floor(length(plaquettecolumn) / opt$every)
column <- (opt$zerooffset + Ns / 2) * opt$zerooffset + 1 + opt$zerooffset

plaquettecf <- readloopfilecfonecolumn(file = filename,
                skip = skip, column = column, every = opt$every)
bootl <- max(opt$bootl, 4 * plaquettedata$tauint^2)
plaquettecf <- bootstrap.cf(plaquettecf,
                boot.R = bootsamples, boot.l = bootl)

# coarse potential
# read in data points, bootstrap samples:
# have to initialise empty vectors,
# because not every list element has a valid entry
xc <- rep(NA, Ns / 2)
yc <- rep(NA, Ns / 2)

#filter out points without valid entry in bootstrap.nlsfit, mask for combining vectors
maskc <- rep(FALSE, Ns / 2)
bsamplesc <- array(c(rep(NA, Ns / 2 * bootsamples)), dim = c(bootsamples, Ns / 2))

#also have vectors for combined potential
# finemask: needed later for combining potentials
x <- rep(NA, Ns / 2 + Nt / 2)
dx <- rep(1e-4, Ns / 2 + Nt / 2)
y <- rep(NA, Ns / 2 + Nt / 2)
mask <- rep(FALSE, Ns / 2 + Nt / 2)
finemask <- rep(FALSE, Ns / 2 + Nt / 2)
bsamples <- array(c(rep(NA, 2 * (Ns / 2 + Nt / 2) * bootsamples)), dim = c(bootsamples, 2 * (Ns / 2 + Nt / 2)))
#~ message("defined empty arrays")

# print(names(listmeff[[1]][[1]]))
# print(names(listmeff[[1]][[2]]))

for (i in seq(1, Ns / 2 - opt$omit, 1)) {
  if (listmeff[[i]][[2]]) {
      xc[i] <- listmeff[[i]][[2]]
      yc[i] <- listmeff[[i]][[1]]$effmassfit$t0[1]
      maskc[i] <- i > opt$lowlimpot
      if(!opt$aic) {
        bsamplesc[, i] <- listmeff[[i]][[1]]$massfit.tsboot[, 1]
        bsamples[, i] <- listmeff[[i]][[1]]$massfit.tsboot[, 1]
      }
      else if (opt$aic) {
        # include systematical errors
        # bootsample = bootsample - (mean - bootsample) * (errtotmean / bootstat - 1)
        if (opt$errortotpot)  bsamplesc[, i] <- listmeff[[i]][[1]]$boot$m50 - (listmeff[[i]][[1]]$effmassfit$t0 - listmeff[[i]][[1]]$boot$m50) * (listmeff[[i]][[1]]$effmassfit$se / listmeff[[i]][[1]]$boot$errstat - 1)
        # only include statistical errors
        if (!opt$errortotpot) bsamplesc[, i] <- listmeff[[i]][[1]]$boot$m50
        bsamples[, i] <- bsamplesc[, i]
      }

      x[i] <- listmeff[[i]][[2]]
      y[i] <- listmeff[[i]][[1]]$effmassfit$t0[1]
      mask[i] <- i > opt$lowlimpot
      finemask[i] <- FALSE
  }
}
xall <- x
yall <- y
bsamplesall <- bsamples[, 1:Ns / 2 + Nt / 2]



# fine potential
# read in data points, bootstrap samples:
# have to initialise empty vectors,
# because not every list element has a valid entry
xf <- rep(NA, Nt / 2)
yf <- rep(NA, Nt / 2)

#filter out point without valid entry in bootstrap
maskf <- rep(FALSE, Nt / 2)
bsamplesf <- array(c(rep(NA, Nt / 2 * bootsamples)), dim = c(bootsamples, Nt / 2))
for (i in seq(1, Nt / 2 - opt$omit / xi, 1)) {
  if (listmeff[[Ns / 2 + i]][[2]]) {
      xf[i]     <- listmeff[[Ns / 2 + i]][[2]]
      yf[i]     <- listmeff[[Ns / 2 + i]][[1]]$effmassfit$t0[1]
      maskf[i]      <- i > opt$lowlimpot/opt$xi
      if (!opt$aic) {
        bsamplesf[, i] <- listmeff[[Ns / 2 + i]][[1]]$massfit.tsboot[, 1]
        bsamples[, Ns / 2 + i] <- listmeff[[Ns / 2 + i]][[1]]$massfit.tsboot[, 1]
      }
      else if (opt$aic) {
        # include systematical errors
        # bootsample = bootsample - (mean - bootsample) * (errtotmean / bootstat - 1)
        if (opt$errortotpot) bsamplesf[, i] <- listmeff[[Ns / 2 + i]][[1]]$boot$m50 - (listmeff[[Ns / 2 + i]][[1]]$effmassfit$t0 - listmeff[[Ns / 2 + i]][[1]]$boot$m50) * (listmeff[[Ns / 2 + i]][[1]]$effmassfit$se / listmeff[[Ns / 2 + i]][[1]]$boot$errstat - 1)
        # only include statistical errors
        if (!opt$errortotpot) bsamplesf[, i] <- listmeff[[Ns / 2 + i]][[1]]$boot$m50
        bsamples[, Ns / 2 + i] <- bsamplesf[, i]
      }

      x[Ns / 2 + i] <- listmeff[[Ns / 2 + i]][[2]]
      y[Ns / 2 + i] <- listmeff[[Ns / 2 + i]][[1]]$effmassfit$t0[1]
      mask[Ns / 2 + i]  <- i > opt$lowlimpot/opt$xi
      finemask[Ns / 2 + i] <- TRUE

  }
}


## for the fine and coarse potentials, only omit is necessary, lowlim does not change the result
#determine parameters of potentials by bootstrap, save results
# print(cor(bsamplesc[, maskc]))

fit.resultcoarse <- bootstrap.nlsfit(fnpot, c(0.2, 0.2, 0.2), yc, xc,
                                        bsamplesc, mask = maskc)
filenamecoarse <- sprintf("%sfitresultcoarsesideways%s.RData", opt$plotpath, endingdofit)

saveRDS(fit.resultcoarse, file = filenamecoarse)

# print(cor(bsamplesf[, maskf]))
fit.resultfine <- bootstrap.nlsfit(fnpot, c(0.2, 0.2, 0.2), yf, xf,
                                    bsamplesf, mask = maskf)
filenamefine <- sprintf("%sfitresultfinesideways%s.RData", opt$plotpath, endingdofit)
saveRDS(fit.resultfine, file = filenamefine)

# plot potentials
try(plot(fit.resultfine, main = sprintf(
        "beta=%.3f, xi=%.3f, (2 + 1)D, Ns=%d, fine\n
        %d measurements of which skipped %d",
        beta, xi, Ns, nom * opt$nsave, skip * opt$nsave),
        ylab = "a_s V(t)", xlab = "t / a_t"))
try(plot(fit.resultcoarse, main = sprintf(
        "beta=%.3f, xi=%.3f, (2 + 1)D, Ns=%d, coarse\n
        %d measurements of which skipped %d",
        beta, xi, Ns, nom * opt$nsave, skip * opt$nsave),
        ylab = "a_s V(y)", xlab = "y/a_s"))

# calculate xi in different ways: determine ratio of linear coefficients of the potentials
singlexis <- data.frame(xi = NA)
for (bs in seq(1, bootsamples)) {
    val <- fit.resultfine$t[bs, 2] / fit.resultcoarse$t[bs, 2]
    singlexis <- rbind(singlexis, data.frame(xi = val))
}
singlexis <- na.omit(singlexis)
xisingle <- mean(singlexis$xi)
dxisingle <- sd(singlexis$xi)
strings <- sprintf("xi from single potentials: %f+/-%f", xisingle, dxisingle)


# determine xi by interpolating between different V(t) to get V(x) = V(t xi),
# like in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.63.074501

# Determine parameters of interpolation for V_t, V_t(t) = a * t + b
# defined by two points, lower and upper,
# with distance 1 a_t, V(i) = lower, V(i + 1) = upper
# for all bootstrapsamples: calculate interpolation, determine xi from interpolation
getxisideways <- function(ydata, Nstmp = Ns, Nttmp = Nt, omit = opt$omit, lowlim = opt$lowlim, xitmp = xi) {
    interpolation <- data.frame(lower = NA, upper = NA,
                a = NA, b = NA, lowerx = NA)
    for (i in seq(1, Nttmp / 2 - 1 - omit / xitmp, 1)) {
        if (listmeff[[Nstmp / 2 + i + 1]][[2]]) {
    #     lower     <- listmeff[[Nstmp / 2 + i]][[1]]$massfit.tsboot[bs, 1]
    #     upper     <- listmeff[[Nstmp / 2 + i + 1]][[1]]$massfit.tsboot[bs, 1]
        lower <- ydata[Nstmp / 2 + i]
        upper <- ydata[Nstmp / 2 + i + 1]
        a <- upper - lower
        b <- lower - i * a
        newline <- data.frame(lower = lower, upper = upper,
                    a = a, b = b, lowerx = i)
        interpolation <- rbind(interpolation, newline)
    }
    }
    
    interpolation <- interpolation[-1, ]
    
    # solve equation V_s(x) = V_t(xi t) = a * t + b -> t = (V_s-b) / a
    # -> xir = x / t = xiresult
    # for all possible V_s with x > =  2 (variable name xi is already taken)
    # potential not dominated by linear part for y > 2
    xiboots <- c()
    for (i in seq(1 + lowlim, Nstmp / 2 - omit, 1)) {
      if (listmeff[[i]][[2]]) {
        Vs <- ydata[i]
        index <- NA
    # determine in which interval of the V_t V_s is sitting by selecting
    # lowest possible interval for which the upper limit is bigger than V_s
        for (j in seq(1, Nttmp / 2 - 1 - omit / xitmp)) {
            if (length(interpolation$upper[interpolation$lowerx == j]) > 0) {
            if (Vs < interpolation$upper[interpolation$lowerx == j]) {
                index <- interpolation$lowerx == j
                break
            }
            }
        }
        #highest V_s could be bigger than highest V_t, so there is not an index for every y
        if (!is.na(index)) {
            t <- (Vs - interpolation$b[index]) / interpolation$a[index]
            xir <- i / t
            newline <- data.frame(xi = xir, t = t, y = i,
                        j = interpolation$lowerx[index])
            xiboots <- append(xiboots, xir)
        }
      }
    }
    return(list(mean = mean(xiboots), all = xiboots))
}
xis <- apply(bsamples, 1, getxisideways)


xis_xi <- unlist(sapply(xis, getElement, "mean"))
xibootsamples <- unlist(sapply(xis, getElement, "all"))

# calculate xi, xi^2 as mean over all values
# compare to values of bootstrapsamples
xicalc <- mean(xis_xi)
dxicalc <-  sd(xis_xi)
xisquared <- mean(xis_xi^2)
dxisquared <- sd(xis_xi^2)
strings <- sprintf("%e +/- %e %f",
        xicalc - mean(xibootsamples), dxicalc - sd(xibootsamples),
        mean(xibootsamples))
message("deviation different means: ", strings)
message(xicalc)
message(warnings())

}

if (dofit) {
    #rescale t so potential, r0 can be determined, all r are given in units of a_s

# generate bootsamples for fit
bsamplesx <- parametric.bootstrap(bootsamples, c(x), c(dx))

# For the coarse potential, no rescaling is necessary,
# but an error has to be given to the fit-function. could error be zero?
dx[!finemask] <- 1e-8

# generate bootsamples for fit
bsamplesx <- parametric.bootstrap(bootsamples, c(x), c(dx))

x[finemask] <- x[finemask] * xicalc
dx[finemask] <- x[finemask] * dxicalc

#join bootstrapsamples for x, y together
# for (i in seq(1, Ns / 2 - opt$omit)) {
#     bsamples[, i + Ns / 2 + Nt / 2] <- bsamplesx[, i]
# }
# for (i in seq(1, Nt / 2 - opt$omit / xi)) {
#     bsamples[, i + Ns + Nt / 2] <- i * array(xibootsamples, dim = c(bootsamples, 1))
# }
bsamples[, Ns/2 + Nt/2 + (1:(Ns/2 - opt$omit))] <- bsamplesx[, (1:(Ns/2 - opt$omit))]
bsamples[, Ns + Nt/2 + seq(1, Nt / 2 - opt$omit / xi)] <- bsamplesx[, Ns/2 + seq(1, Nt / 2 - opt$omit / xi)] * array(rep(xibootsamples, Nt / 2 - opt$omit / xi), dim = c(bootsamples, Nt / 2 - opt$omit / xi))


title <- sprintf("beta = %f, xi = %f +/-%f, (2+1)D, Ns = %d\n
        %d measurements of which skipped %d",
        beta, xicalc, dxicalc, Ns, nom * opt$nsave, skip * opt$nsave)

# fit and save overall potential
# print(cor(bsamples))
fit.resultscaled <- bootstrap.nlsfit(fnpot, c(0.1, 0.13, 0.05),
                    y, x, bsamples, mask = mask)
filenamescaled <- sprintf(
            "%sfitresultscaledsideways%s.RData",
            opt$plotpath, endingdofit)

saveRDS(fit.resultscaled, file = filenamescaled)

strings <- sprintf("\nxi = %f  +/- %f\neta= %f  +- %f\n\n",
            xicalc, dxicalc, xicalc / xi, dxicalc / xi)

#if wanted determine correlation matrix seperately, as determined by the summary.bootstrapfit function
if (!is.null(fit.resultscaled$cov)) {
    cov.to.cor <- diag(1  /  fit.resultscaled$se)
    correlation <- cov.to.cor %*% fit.resultscaled$cov %*% cov.to.cor
}else if (!is.null(fit.resultscaled$t)) {
    correlation <- cor(fit.resultscaled$t,
            fit.resultscaled$t, use = "na.or.complete")
}

#~ print(correlation)

#graphically represent force: for each bootstrap sample of parameters, determine value of force on a sequence,
#determine error as standard deviation of values on each point of the sequence
xx <- seq(-0.2, Ns / 2 + 0.2, 0.05)
bootsforce <- data.frame(matrix(rep(NA, length(xx)), ncol = length(xx), nrow = 1))

for (i in seq(1, bootsamples, 1)) {
    bootsforce <- rbind(bootsforce, fnforce(fit.resultscaled$t[i, ], xx, 0))
}
bootsforce <- na.omit(bootsforce)

# forceerrs <- c()

# for (i in seq(1, length(xx), 1)) {
#     forceerrs[i] <- sd(bootsforce[, i])
# }
forceerrs <- apply(bootsforce, 2, sd)

#Determine r0 as solution of equation -r^2 d / dr V(r) = c, solve for each bootstrapsample, r0 = mean pm sd
# V(r) = a + sigma * r + b * ln(r)
# -d / dr V(r) = -sigma -b / r
# c = -sigma * r^2 - b * r, use p-q-formula for p = b / sigma, q = c / sigma
rzeroofc <- data.frame(r0 = NA, dr0 = NA, c = NA)

for (c in list(opt$crzero)) {
# rzerolist <- data.frame(r1 = NA, r2 = NA)
# for (bs in seq(1, bootsamples)) {
#     sigma <- fit.resultscaled$t[bs, 2]
#     b <- fit.resultscaled$t[bs, 3]
#     r1 <- -b / (2.0 * sigma) + sqrt((b / (2.0 * sigma))^2 - c / sigma)
#     r2 <- -b / (2.0 * sigma) - sqrt((b / (2.0 * sigma))^2 - c / sigma)
#     rzerolist <- rbind(rzerolist, data.frame(r1 = r1, r2 = r2))
# }
rzerolist <- apply(fit.resultscaled$t, 1, getrzero)
rzero <- getrzero(fit.resultscaled$t0)
# print(rzero - mean(rzerolist))
drzero <- sd(rzerolist)
rzeroofc <- rbind(rzeroofc, data.frame(r0 = rzero, dr0 = drzero, c = c))
bsrzero <- rzerolist
}
rzeroofc <- rzeroofc[-1, ]

#plot results: overall potential
title <- sprintf("beta = %.3f, xi = %.3f +/- %.3f, (2+ 1)D, Ns = %d\n
                %d measurements of which skipped %d",
                beta, xicalc, dxicalc, Ns, nom * opt$nsave, skip * opt$nsave)

try(plot(fit.resultscaled, main = title, ylab = "V(r)", xlab = "r / a_s"))

#force with error and r0, zoomed in
## values for force +/- errors
polyval <- c(fnforce(fit.resultscaled$t0, xx, 0) + forceerrs,
            rev(fnforce(fit.resultscaled$t0, xx, 0) - forceerrs))
pcol <- col2rgb("gray", alpha = TRUE) / 255
pcol[4] <- 0.65
pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
xlim <- c(0.6 * min(rzeroofc$r0), 1.4 * max(rzeroofc$r0))
ylim <- c(0.6 * min(abs(rzeroofc$c)), 1.4 * max(abs(rzeroofc$c)))
ylim <- rev(-ylim)
lowylim <- -2 * max(abs(rzeroofc$c))
try(plot(xx, fnforce(fit.resultscaled$t0, xx, 0),
        type = "l", xlab = "r / a_s", ylab = "-r^2 d / dr V(r)",
        main = title, xlim = xlim, ylim = ylim))
try(polygon(x = c(xx, rev(xx)), y = polyval,
        col = pcol, lty = 0, lwd = 0.001, border = pcol))


for (i in seq(1, length(rzeroofc$c))) {
    c <- rzeroofc$c[i]
    rzero <- rzeroofc$r0[i]
    drzero <- rzeroofc$dr0[i]
#also make lines for position of r0, polygon for error of r0
try(arrows(rzero, fnforce(fit.resultscaled$t0, rzero, 0),
                rzero, lowylim, angle = 90, length = 0.1, code = 0))
try(arrows(rzero, c, 0, c, angle = 90, length = 0.1, code = 0))
r0polyval1 <- c(fnforce(fit.resultscaled$t0, rzero - drzero, 0)
                - fnforceerr(fit.resultscaled$se, rzero - drzero, 0, correlation),
                fnforce(fit.resultscaled$t0, rzero + drzero, 0)
                + fnforceerr(fit.resultscaled$se, rzero + drzero, 0, correlation))
r0polyval1 <- c(r0polyval1, lowylim, lowylim)
r0polyval2 <- c(rzero - drzero, rzero + drzero)
r0polyval2 <- c(r0polyval2, rev(r0polyval2))
try(polygon(x = r0polyval2, y = r0polyval1,
        col = pcol, lty = 0, lwd = 0.001, border = pcol))
}

#plot visible combination of temporal and spatial potential
try(plot(fit.resultscaled,
        main = title, ylab = "V(r)", xlab = "r / a_s", pch = 3, cex = 0.1))
try(pointsyerr(col = 1, pch = 1, x = fit.resultscaled$x[finemask],
        y = fit.resultscaled$y[finemask], dy = fit.resultscaled$dy[finemask]))
try(pointsyerr(col = 2, pch = 2, x = fit.resultscaled$x[!finemask],
        y = fit.resultscaled$y[!finemask], dy = fit.resultscaled$dy[!finemask]))
legend(legend = c("spatial", "temporal"), col = c(2, 1), pch = c(2, 1), x = "topleft")

#plot thermalisation
plot(plaquettecolumn, main = "Thermalisation",
        xlab = sprintf("MCMC-steps / %d", opt$nsave), ylab = "P")

# plot qqplots to see if xi and rzero are normally distributed
try(qqnorm(xibootsamples, main = "qqplot of xi_ren"))
try(qqline(xibootsamples))
try(qqnorm(bsrzero, main = "qqplot of r_0"))
try(qqline(bsrzero))


#save values of force
resultforce <- data.frame(R = NA, force = NA, forceerr = NA)
for (i in seq (1, length(xx), 1)) {
    resultforce <- rbind(resultforce,
                data.frame(R = xx[i], force = fnpot(fit.resultscaled$t0, xx[i], 0),
                forceerr = forceerrs[i]))
}
resultforce <- na.omit(resultforce)

# determine r_0 from single potentials
rzerocoarse <- determinerzero(fit.resultcoarse, bootsamples = bootsamples)
rzerofine <- determinerzero(fit.resultfine, bootsamples = bootsamples)

# determine xi as ratio of different r_0s, string tensions
differentxis <- data.frame(xirzero = NA, xist = NA)

for (bs in seq(1, bootsamples)) {
xirzero <- rzerocoarse[[3]][bs] / rzerofine[[3]][bs]
xist <- fit.resultfine$t[bs, 2] / fit.resultcoarse$t[bs, 2]
differentxis <- rbind(differentxis, data.frame(xirzero = xirzero, xist = xist))
}
differentxis <- differentxis[-1, ]
differentxiresults <- c(mean(differentxis$xirzero), sd(differentxis$xirzero),
        mean(differentxis$xist), sd(differentxis$xist))


}
if (dofit) {


resultssummary <- list(st = fit.resultscaled$t0[2], dst = fit.resultscaled$se[2], bsst = fit.resultscaled$t[, 2],
        rzeros = rzero, drzeros = drzero, crz = opt$crzero, bsrzeros = bsrzero,
        p = plaquettecf$cf.tsboot$t0[1], dp = plaquettecf$tsboot.se[1], bsp = plaquettecf$cf.tsboot$t[, 1],
        beta = beta, xiin = xi, Nt = Nt, Ns = Ns, bootsamples = bootsamples,
        nom = nom, skip = opt$skip,
        xicalc = xicalc, dxicalc = dxicalc, bsxicalc = xibootsamples,
        lowlim = opt$lowlim, lowlimpot=opt$lowlimpot, omit = opt$omit, etp = opt$errortotpot)
class(resultssummary) <- "resultssummary"


resultlist <- data.frame(xi = NA, beta = NA, xicalc = NA, dxicalc = NA,
        r0 = NA, dr0 = NA, st = NA, dst = NA, p = NA, dp = NA,
        chi = NA, chifine = NA, chicoarse = NA, c = NA, bs = NA, Ns = NA,
        Nt = NA, nape = NA, alpha = NA, omit = NA, nom = NA, skip = NA,
        xisingle = NA, dxisingle = NA,
        job = NA, hash = NA, every = NA, tauint = NA, dtauint = NA, bootl = NA,
        xirzero = NA, dxirzero = NA, xist = NA, dxist = NA, lowlim = NA, lowlimpot = NA, aic = NA,
        scaletauint = NA, puw = NA, dpuw = NA, errortotpot = NA, coulpart = NA, dcoulpart = NA)

for (i in seq(1, max(1, length(rzeroofc$c)))) {
newline <- data.frame(xi = xi, beta = beta, xicalc = xicalc, dxicalc = dxicalc,
        r0 = rzeroofc$r0[i], dr0 = rzeroofc$dr0[i],
        st = fit.resultscaled$t0[2], dst = fit.resultscaled$se[2],
        p = plaquettecf$cf.tsboot$t0[1], dp = plaquettecf$tsboot.se[1],
        chi = fit.resultscaled$chisqr / fit.resultscaled$dof,
        chifine = fit.resultfine$chisqr / fit.resultfine$dof,
        chicoarse = fit.resultcoarse$chisqr / fit.resultcoarse$dof,
        c = rzeroofc$c[i], bs = bootsamples, Ns = Ns, Nt = Nt,
        nape = nape, alpha = alpha, omit = opt$omit, nom = nom, skip = skip,
        xisingle = xisingle, dxisingle = dxisingle,
        job = opt$job, hash = githash, every = opt$every,
        tauint = plaquettedata$tauint, dtauint = plaquettedata$dtauint[1],
        bootl = opt$bootl, xirzero = differentxiresults[1], dxirzero = differentxiresults[2],
        xist = differentxiresults[3], dxist = differentxiresults[4],
        lowlim = opt$lowlim, lowlimpot=opt$lowlimpot, aic = opt$aic, scaletauint = opt$scaletauint,
        puw = plaquettedata$value, dpuw = plaquettedata$dvalue, errortotpot = opt$errortotpot, 
        coulpart = fit.resultscaled$t0[3], dcoulpart = fit.resultscaled$se[3])

resultlist <- rbind(resultlist, newline)
}

resultlist <- resultlist[-1, ]
print(resultlist)
filename <- sprintf("%sresultsummary2p1dsidewaysb%.3fNs%d%s.csv",
                    opt$plotpath, opt$betaone, opt$Ns, opt$extra)
columnnames <- FALSE
if (!file.exists(filename)) {
    columnnames <- TRUE
}
write.table(resultlist, filename,
        append = TRUE, row.names = FALSE, col.names = columnnames)
nameresults <- sprintf("%sresultssideways%s.RData",
        opt$plotpath, endingdofit)
saveRDS(resultssummary, file = nameresults)

}

message(warnings())
