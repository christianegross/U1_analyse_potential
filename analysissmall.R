library(hadron)
library(optparse)


if (TRUE) {
option_list <- list(
    make_option(c("-b", "--beta"), type = "double", default = 1.7,
    help = "beta-Parameter of simulation [default %default]"),
    make_option(c("--skip"), type = "integer", default = 1000,
    help = "how many lines are skipped when reading in [default %default]"),

    make_option(c("-s", "--bootsamples"), type = "integer", default = 500,
    help = "how many bootstrapsamples should be drawn [default %default]"),
    make_option(c("-r", "--Ns"), type = "integer", default = 3,
    help = "L of lattice [default %default]"),

    make_option(c("-t", "--Nt"), type = "integer", default = 3,
    help = "T of lattice [default %default]"),
    make_option(c("-n", "--nape"), type = "integer", default = 0,
    help = "Number of APE-smears that were done [default %default]"),

    make_option(c("-j", "--job"), type = "integer", default = 0,
    help = "qbig job number [default %default]"),
    make_option(c("--t1"), type = "integer", default = 0,
    help = "boundaries for determining meff if no table is given [default %default]"),

    make_option(c("-a", "--alpha"), type = "double", default = 1.0,
    help = "alpha used in APE-smearing [default %default]"),
    make_option(c("--betaone"), type = "double", default = 0,
    help = "input beta at corresponding xi = 1 [default %default]"),

    make_option(c("--fraction"), type = "double", default = 1,
    help = "fraction of lattice that was measured [default %default]"),
    make_option(c("-x", "--xi"), type = "double", default = 1,
    help = "xi used in lattice, only used if xidiff = TRUE [default %default]"),
    make_option(c("--bootl"), type = "double", default = 2,
    help = "boot-l for bootstrapping samples [default %default]"),

    make_option(c("--analyse"), action = "store_true", default = FALSE,
    help = "if true, potential values are determined with
            effective mass [default %default]"),
    make_option(c("--plaquette"), action = "store_true", default = FALSE,
    help = "if true, plaquette and anisotropy are determined [default %default]"),
    make_option(c("--amplitude"), action = "store_true", default = FALSE,
    help = "if true, amplitude of spatial 1,1 plaquette is calculated [default %default]"),


    make_option(c("--smearing"), action = "store_true", default = FALSE,
    help = "are the calculations done based on
            smeared lattices? (changes filenames) [default %default]"),
    make_option(c("--xidiff"), action = "store_true", default = FALSE,
    help = "Is xi different to L / T? [default %default]"),
    make_option(c("--nsave"), type = "integer", default = 100,
    help = "Step between saved configurations [default %default]"),

    make_option(c("--myfunctions"), type = "character",
        default = "/hiskp4/gross/masterthesis/su2/build/debug/analysisscripts/",
#~     make_option(c("--myfunctions"), type = "character", default = "myfunctions.R",
    help = "path to where additional functions are stored [default %default]"),
    make_option(c("--respath"), type = "character", default = "",
    help = "path to where the resultfiles are stored [default %default]"),
    make_option(c("--plotpath"), type = "character", default = "",
    help = "path to where the plots are stored [default %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
}

if (TRUE) {
# set constants
source(paste(opt$myfunctions, "myfunctions.R", sep = ""))
githash <- printgitcommit(opt$myfunctions)
beta <- opt$beta
skip <- opt$skip
bootsamples <- opt$bootsamples
Ns <- opt$Ns
Nt <- opt$Nt
Nsmax <- min(Ns, 4)
Ntmax <- floor(Nt * opt$fraction)
message("Ntmax = ", Ntmax, " Nsmax = ", Nsmax, "\n")
xi <- Ns / Nt
if (opt$xidiff) {
    xi <- opt$xi
}

nape <- opt$nape
alpha <- opt$alpha
}

# analyse: determine potential points vie bootstrapping Wilson
# loops, determining effective masses and fitting to that
if (opt$analyse) {

# either read in existing table or generate new dataframe
# for boundaries of the plateaus of the effective masses
namelistbounds <- sprintf(
            "%stableanalysissmallNs%dNt%dbeta%fxi%f.csv",
            opt$plotpath, Ns, Nt, beta, xi)
if (opt$smearing) {
    namelistbounds <- sprintf(
            "%stableanalysissmallNs%dNt%dbeta%fxi%fnape%dalpha%f.csv",
            opt$plotpath, Ns, Nt, beta, xi, nape, alpha)
}

if (file.exists(namelistbounds)) {
    listbounds <- read.table(namelistbounds, sep = ",", header = TRUE)
} else {
    listbounds <- cbind(
                data.frame(x = c(rep(0, 4), rep(1, 4), rep(2, 4), rep(3, 4))),
                data.frame(y = rep(c(0, 1, 2, 3), 4)),
                data.frame(lower = rep(opt$t1, 16)),
                data.frame(upper = rep(Ntmax - opt$t1 - 2, 16)))
}


print(listbounds)

# open pdf for plots, open lists for results
filenameforplots <- sprintf(
            "%splots2p1dsmallNs%dNt%dbeta%fxi%fbs%d.pdf",
            opt$plotpath, Ns, Nt, beta, xi, bootsamples)
if (opt$smearing) {
    filenameforplots <- sprintf(
            "%splots2p1dsmallNs%dNt%dbeta%fxi%fnape%dalpha%fbs%d.pdf",
            opt$plotpath, Ns, Nt, beta, xi, nape, alpha, bootsamples)
}
pdf(file = filenameforplots, title = "")
listfits <- list(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
listtauint <- list()
listpotshort <- data.frame(R = NA, m = NA, dm = NA)

# iterate over every possible combination of (x, y)
# determine the expectation values, the effective masses
# set the plateau, plot the loops and the effective masses.
# save everything in lists
for (x in seq(0, 3)) {
    for (y in seq(0, 3)) {
    if (x > 0 || y > 0) {
    #save results, list has to have an entry everywhere
        resulteffmass <- list(data.frame(R = NA, m = NA, dm = NA), c(0, 0), list(FALSE, FALSE, FALSE))
        uwerrresults <- FALSE
        if (x < 2 && y < 2) {
        message("\nx = ", x, " y = ", y)
        filename <- sprintf(
                "%sresult2p1d.u1potential.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fnonplanar",
                opt$respath, Nt, Ns, beta, xi, nape, alpha)

        if (y == 0 || x == 0) {
            start <- 1
            title <- sprintf(
                "beta = %.3f x = %d  y = %d W(t = 0) not used\n", beta, x, y)
        } else {
            start <- 0
            title <- sprintf(
                "beta = %.3f x = %d  y = %d W(t = 0) used\n", beta, x, y)
        }
        WL <- calcplotWloopsmall(filename, skip, Nsmax, Ntmax,
                Nt, x, y, start, bootsamples, title, nsave = opt$nsave, l = opt$boot.l)

        title <- sprintf("%s%d configs, skipped %d", title,
                (length(WL$cf[, 1]) + skip - 1) * opt$nsave, skip * opt$nsave)


        # get boundaries for masses from table
        t1 <- listbounds$lower[listbounds$y == y & listbounds$x == x]
        t2 <- listbounds$upper[listbounds$y == y & listbounds$x == x]

        for (type in c("log", "acosh")) {

        resulteffmass <- deteffmasssmall(WL, x, y, t1, t2, type = type)

        # plot, with additional lines indicating half the lattice
        # extent, save results of meff
        # in short and long form (with and without bootstrapsamples)
        names(resulteffmass[[3]]) <- c("effmass", "x", "y")
        print(paste("type=", type, "effmass=", resulteffmass[[1]][2], "pm", resulteffmass[[1]][3]))

        try(plot(resulteffmass[[3]][[1]], xlab = "t / a_t", ylab = "Meff",
                main = sprintf("%s, t1 = %d, t2 = %d, p = %f",
                title, t1, t2, resulteffmass[[3]][[1]]$effmassfit$Qval)))
        arrows((Nt - start - 1) / 2, 10, (Nt - start - 1) / 2, -10,
                angle = 90, length = 0.1, code = 0, col = 2)
        drawallticks()
        try(plot(resulteffmass[[3]][[1]], xlab = "t / a_t", ylab = "Meff",
                main = sprintf("%s, t1 = %d, t2 = %d, p = %f",
                title, t1, t2, resulteffmass[[3]][[1]]$effmassfit$Qval),
                ylim = resulteffmass[[2]]))
        arrows((Nt - start - 1) / 2, 10, (Nt - start - 1) / 2, -10,
                angle = 90, length = 0.1, code = 0, col = 2)
        drawallticks(all = FALSE, inward = TRUE)
        }


        uwerrresults <- uwerr.cf(WL)
        # plot autocorrelation times of non-bootstrapped results
        try(plot(x = uwerrresults$uwcf$t - 1, y = uwerrresults$uwcf$tauint,
                main = sprintf("%s, autocorellation times", title),
                xlab = "t", ylab = "tauint"))
        uwerrresults <- uwerr.cf(WL)
        # plot autocorrelation times of non-bootstrapped results
        try(plot(x = uwerrresults$uwcf$t - 1, y = uwerrresults$uwcf$tauint,
                main = sprintf("%s, autocorellation times", title),
                xlab = "t", ylab = "tauint"))
        # uwerrbootstrap <- cbind(as.data.frame(t(apply(X=WL$cf, MARGIN=2L, FUN=uw_wrapper))),
        #                  t=(1:ncol(WL$cf)))

    }
        listfits[[4 * x + y]] <- resulteffmass[[3]]
        listtauint[[4 * x + y]] <- uwerrresults
        listpotshort <- rbind(listpotshort, resulteffmass[[1]])
    }
    }
}
listfits[[16]] <- githash
#write out results
filenamelist <- sprintf(
            "%slistmeffsmallNs%dNt%dbeta%fxi%fbsamples%d.RData",
            opt$plotpath, Ns, Nt, beta, xi, bootsamples)
filenameuwerr <- sprintf(
            "%slistuwerrtauintsmallNs%dNt%dbeta%fxi%fbsamples%d.RData",
            opt$plotpath, Ns, Nt, beta, xi, bootsamples)
filenamepot <- sprintf(
            "%slistpotsmallNs%dNt%dbeta%fxi%fbsamples%d.csv",
            opt$plotpath, Ns, Nt, beta, xi, bootsamples)

if (opt$smearing) {
    filenamelist <- sprintf(
            "%slistmeffsmallNs%dNt%dbeta%fxi%fnape%dalpha%fbsamples%d.RData",
            opt$plotpath, Ns, Nt, beta, xi, nape, alpha, bootsamples)
    filenameuwerr <- sprintf(
            "%slistuwerrtauintsmallNs%dNt%dbeta%fxi%fnape%dalpha%fbsamples%d.RData",
            opt$plotpath, Ns, Nt, beta, xi, nape, alpha, bootsamples)
    filenamepot <- sprintf(
            "%slistpotsmallNs%dNt%dbeta%fxi%fnape%dalpha%fbsamples%d.csv",
            opt$plotpath, Ns, Nt, beta, xi, nape, alpha, bootsamples)
}
saveRDS(listfits, file = filenamelist)
saveRDS(listtauint, file = filenameuwerr)
write.table(listpotshort, file = filenamepot, row.names = FALSE)
}


# plaquette: determine the slope of the potential, r_0,
# display the potential graphically, ...
if (opt$plaquette) {
    #read in files
    filenameforplots <- sprintf(
            "%spotential2p1dsmallplaqNs%dNt%dbeta%fxi%fbs%d.pdf",
            opt$plotpath, Ns, Nt, beta, xi, opt$bootsamples)
    if (opt$smearing) {
        filenameforplots <- sprintf(
            "%spotential2p1dsmallplaqNs%dNt%dbeta%fxi%fnape%dalpha%f.pdf",
            opt$plotpath, Ns, Nt, beta, xi, nape, alpha)
    }
    pdf(filenameforplots, title = "")
    if (opt$bootsamples == 0) {
        bootsamples <- find_bootsamples(Nt = Nt, Ns = Ns,
        xi = xi, beta = beta, plotpath = opt$plotpath)
    }

    filename <- sprintf(
            "%sresult2p1d.u1potential.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fnonplanar",
            opt$respath, Nt, Ns, beta, xi, nape, alpha)

    filenamelist <- sprintf(
            "%slistmeffsmallNs%dNt%dbeta%fxi%fbsamples%d.RData",
            opt$plotpath, Ns, Nt, beta, xi, bootsamples)

    if (opt$smearing) {
        filenamelist <- sprintf(
            "%slistmeffsmallNs%dNt%dbeta%fxi%fnape%dalpha%fbsamples%d.RData",
            opt$plotpath, Ns, Nt, beta, xi, nape, alpha, bootsamples)
    }
    listmeff <- readRDS(file = filenamelist)
    bootsamples <- length(listmeff[[5]][[1]]$massfit.tsboot[, 1])
    # print(listmeff)

    if (TRUE) { #plaquette
        # determine the plaquette with uwerr
        # (takes autocorrelation into account)
        # anf with bootstrap
        # plaquettedata <- readloopfilecfplaquette(file = filename,
        #                 skip = skip, Nsmax = Nsmax, Ntmax = Ntmax)
        plaquettedata <- readloopfilecfonecolumn(file = filename, skip = skip,
                            column = Nsmax + 3, every = 1)
        p1 <- plaquettedata$cf[, 1]
        dpvec <- c()
        for (i in seq(1, 10)) {
            print(i)
        plaquettedata <- bootstrap.cf(plaquettedata, boot.R = bootsamples, boot.l=2**i)

        plaquette <- plaquettedata$cf0[1]
        dp <- plaquettedata$tsboot.se[1]
        dpvec[i] <- dp
        }
        plaquettedata <- bootstrap.cf(plaquettedata, boot.R = bootsamples, boot.l=opt$bootl)

        plaquette <- plaquettedata$cf0[1]
        dp <- plaquettedata$tsboot.se[1]
#         print(names(plaquettedata))
#         [1] "nrObs"             "Time"              "nrStypes"
#  [4] "symmetrised"       "cf"                "conf.index"
#  [7] "boot.R"            "boot.l"            "seed"
# [10] "sim"               "endcorr"           "cf.tsboot"
# [13] "error_fn"          "cov_fn"            "resampling_method"
# [16] "cf0"               "tsboot.se"         "boot.samples"
        measurements <- read.table(filename, header = FALSE, skip = 0)
        plaquettedatauwerr <- uwerrprimary(measurements[seq(opt$skip + 1, length(measurements[, Nsmax + 3])), Nsmax + 3])
        p2 <- measurements[seq(opt$skip, length(measurements[, Nsmax + 3])), Nsmax + 3]
        # print(head(p1))
        # print(head(p2))
        # print(tail(p1))
        # print(tail(p2))
        # write.table("differenceplaquettes.csv", data.frame(p1=p1, p2=p2), row.names=F, col.names=T)
        # pdiff <- head(p1)-head(p2)
        # print(sum(pdiff))
        plaquetteuwerr <- plaquettedatauwerr$value
        dpuwerr <- plaquettedatauwerr$dvalue
        tauint <- plaquettedatauwerr$tauint
        plot(x=seq(1, length(measurements[, 1])), y=measurements[ , Nsmax + 3], xlab=sprintf("MCMC-steps/%d", opt$nsave), ylab="P", main="Thermalization")
        arrows(x0=opt$skip, y0=-1, y1=2)

        # print(dpvec)
        # print(dpuwerr)
        print(dpvec/dpuwerr)
        print(plaquette/plaquetteuwerr)
        print(plaquetteuwerr)
        print(plaquette)

    }

    if(FALSE) { # slope, ratio

    # Read in V(x=1, y=0) and V(sqrt(2)), then determine ratio and
    # slope of these data. Try to determine slope with a linear fit,
    # else determine slope for each bootstrapsample.
    # If no points are found, results are set to NA


    fnlin <- function (par, x, boot.r, ...) par[1]  + par[2]  *  x

    x <- c(1, sqrt(2))
    pot <- c(NA, NA)
    dpot <- c(NA, NA)
    bsamplespot <- array(c(rep(NA, 2 * bootsamples)), dim = c(bootsamples, 2))
    if (listmeff[[4]][[2]] && listmeff[[5]][[2]]) {
        pot[1] <- listmeff[[4]][[1]]$effmassfit$t0[1]
        pot[2] <- listmeff[[5]][[1]]$effmassfit$t0[1]
        bsamplespot[, 1] <- listmeff[[4]][[1]]$massfit.tsboot[, 1]
        bsamplespot[, 2] <- listmeff[[5]][[1]]$massfit.tsboot[, 1]
        dpot[1] <- listmeff[[4]][[1]]$effmassfit$se[1]
        dpot[2] <- listmeff[[5]][[1]]$effmassfit$se[1]

        ratiobootsamples <- c()

        for (bs in seq(1, bootsamples)) {
            ratiobootsamples[bs] <- bsamplespot[bs, 1]  /  bsamplespot[bs, 2]
        }
        ratio <- mean(ratiobootsamples)
        dratio <- sd(ratiobootsamples)

        # fit.result <- try(bootstrap.nlsfit(fnlin,
        #             c(0.2, 0.2), pot, x, bsamplespot))
        # if (!inherits(fit.result, "try-error")) {
        #     print(fit.result)
        #     slope <- fit.result$t0[2]
        #     dslope <- fit.result$se[2]
        #     bsslope <- fit.result$t[, 2]
        #     chislope <- fit.result$chisqr
        #     icslope <- fit.result$t0[1]
        #     dicslope <- fit.result$se[1]
        #     bsicslope <- fit.result$t[, 1]

        # } else {
            bsslope <- c()
            bsicslope <- c()
            for (bs in seq(1, bootsamples)) {
                bsslope[bs] <- (bsamplespot[bs, 2] - bsamplespot[bs, 1]) / (sqrt(2) - 1)
                bsicslope[bs] <- bsamplespot[bs, 1] - bsslope[bs]
            }
            slope <- mean(bsslope)
            dslope <- sd(bsslope)
            icslope <- mean(bsicslope)
            dicslope <- sd(bsicslope)
            chislope <- NA
        # }
    } else {
        print("no points for xi found")
        fit.result <- NA
        slope <- NA
        dslope <- NA
        bsslope <- rep(NA, bootsamples)
        icslope <- NA
        dicslope <- NA
        bsicslope <- rep(NA, bootsamples)
        chislope <- NA
        ratio <- NA
        dratio <- NA
        ratiobootsamples <- rep(NA, bootsamples)
    }

    # plot slope either as the fit or as the two points
    title <- sprintf("beta = %f, xi_in = %f, Ns = %d, Nt = %d",
                beta, xi, Ns, Nt)
    error3 <- try(plot(fit.result, xlab = "r / a_s",
                    ylab = "a_tV(r)", main = title))
    if (inherits(error3, "try-error")) {
        try(plotyerr(x = x, y = pot, dy = dpot,
                xlab = "r / a_s", ylab = "a_tV(r)", main = title))
    }
    }

    if (FALSE) { #potential +string tension, r0
    # read in all points of the potential
        fnpot <- function (par, x, boot.r, ...) par[1]  + par[2]  *  x  + par[3]  *  log(x)
        fnforce <- function (par, x, boot.r, ...) (- 1.0)  *  par[2]  *  x^2 - par[3]  *  x

        mask <- c(rep(FALSE, 15))
        x <- c(rep(NA, 15))
        y <- c(rep(NA, 15))
        V <- c(rep(NA, 15))
        bsamples <- array(c(rep(NA, 15 * bootsamples)), dim = c(bootsamples, 15))
        for (i in seq(0, 3)) {
            for (j in seq(0, 3)) {
                if (i != 0 || j != 0) {
                if (listmeff[[4 * i + j]][[2]] || listmeff[[4 * i + j]][[3]]) {
                    mask[4 * i + j] <- TRUE
                    x[4 * i + j] <- i
                    y[4 * i + j] <- j
                    V[4 * i + j] <- listmeff[[4 * i + j]][[1]]$effmassfit$t0[1]
                    bsamples[, 4 * i + j] <- listmeff[[4 * i + j]][[1]]$massfit.tsboot[, 1]
                } else {
                    print(paste("no result found for x=", i, " y=", j))
                }
                }
            }
        }

        # try to find a fit to the expected form,
        # from that determine string tension, rzero,
        # the offset of the potential and the factor in front of the logarithm
        # if that fails, set results to NA
        fit.pot <- try(bootstrap.nlsfit(fnpot, c(1, 1, 1),
                y = V, x = sqrt(x^2 + y^2), bsamples, mask = mask))
        if (!inherits(fit.pot, "try-error") && listmeff[[2]][[2]]) {
            print(fit.pot)
            plot(fit.pot, ylab = "a_tV(r)", xlab = "r / a_s",
                    main = "Potential including nonplanar points",
                    plot.range = c(1e-5, sqrt(18)), xlim = c(0, sqrt(18)))
            st <- fit.pot$t0[2]
            dst <- fit.pot$se[2]
            bsst <- fit.pot$t[, 2]
            icpot <- fit.pot$t0[1]
            dicpot <- fit.pot$se[1]
            bsicpot <- fit.pot$t[, 1]
            logpot <- fit.pot$t0[3]
            dlogpot <- fit.pot$se[3]
            bslogpot <- fit.pot$t[, 3]
            chipot <- fit.pot$chisqr / fit.pot$dof
            rzeros <- c()
            drzeros <- c()
            bsrzeros <- list()
            crz <- c(-1.65)
            for (i in seq(1, length(crz))) {
                rzerovec <- determinerzero(fit.pot, c = crz[i],
                            xi = xi, bootsamples = bootsamples)
                rzeros[i] <- rzerovec[[1]]
                drzeros[i] <- rzerovec[[2]]
                bsrzeros[[i]] <- rzerovec[[3]]
            }
        } else {
            print("fit to potential not sucessful")
            plotyerr(x = sqrt(x^2 + y^2), y = V, dy = apply(bsamples, 2, sd),
                    ylab = "a_tV(r)", xlab = "r",
                    main = "Potential including nonplanar points")
            st <- NA
            dst <- NA
            bsst <- rep(NA, bootsamples)
            icpot <- NA
            dicpot <- NA
            bsicpot <- rep(NA, bootsamples)
            logpot <- NA
            dlogpot <- NA
            bslogpot <- rep(NA, bootsamples)
            chipot <- NA
            rzeros <- c(NA)
            drzeros <- c(NA)
            bsrzeros <- list(rep(NA, bootsamples))
            crz <- c(NA)
        }
    }

    if (FALSE) {
        # calculate ratio, but subtract either
        # the intercept of the slope or the additive constant
        bsratioslope <- c()
        bsratiopot <- c()
        for (bs in seq(1, bootsamples)) {
            bsratioslope[bs] <- (bsamplespot[bs, 1] - bsicslope[bs])  /  (bsamplespot[bs, 2] - bsicslope[bs])
            bsratiopot[bs] <- (bsamplespot[bs, 1] - bsicpot[bs])  /  (bsamplespot[bs, 2] - bsicpot[bs])
        }
        ratioslope <- mean(bsratioslope)
        dratioslope <- sd(bsratioslope)
        ratiopot <- mean(bsratiopot)
        dratiopot <- sd(bsratiopot)
    }

    if (FALSE) {

    # Make a list of all results plus bootstrapsamples,
    # a data frame of the results,
    # and save them. Only print columnnames if the file is empty.
    resultssummary <- list(slope = slope, dslope = dslope, bsslope = bsslope, chislope = chislope,
        ratio = ratio, dratio = dratio, bsratio = ratiobootsamples,
        st = st, dst = dst, bsst = bsst,
        rzeros = rzeros, drzeros = drzeros, crz = crz, bsrzeros = bsrzeros,
        icslope = icslope, dicslope = dicslope, bsicslope = bsicslope,
        icpot = icpot, dicpot = dicpot, bsicpot = bsicpot,
        logpot = icpot, dlogpot = dicpot, bslogpot = bslogpot,
        ratioslope = ratioslope, dratioslope = dratioslope, bsratioslope = bsratioslope,
        ratiopot = ratiopot, dratiopot = dratiopot, bsratiopot = bsratiopot,
        p = plaquette, dp = dp, bsp = plaquettedata$cf.tsboot$t[, 1],
        beta = opt$beta, xiin = xi, Nt = opt$Nt, Ns = opt$Ns, bootsamples = bootsamples,
        nom = length(measurements[, 1]) + opt$skip, skip = opt$skip, fit.pot = fit.pot)
    class(resultssummary) <- "resultssummary"

    results <- data.frame(beta = NA, xiin = NA, Nt = NA, Ns = NA,
     p = NA, dp = NA, slope = NA, dslope = NA, chislope = NA,
     ratio = NA, dratio = NA, st = NA, dst = NA, chipot = NA,
     icslope = NA, dicslope = NA, icpot = NA, dicpot = NA, logpot = NA, dlogpot = NA,
     ratioslope = NA, dratioslope = NA, ratiopot = NA, dratiopot = NA,
     rzero = NA, drzero = NA, c = NA, puw = NA, dpuw = NA, tauint = NA, job = NA, hash = NA,
     boot = NA, nom = NA, skip = NA)
     #beta, xi_input, N_t, N_s = lattice extents, p = plaquette from bootstrap.cf
     #slope between V(1) and V(sqrt(2))
     #string tension = linear part of potential (a +st * r +c * ln(r))
     #intercept at r = 0 determined from slope
     #intercept at r = 0 determined from potential
     #ratio with subtracted intercepts from slope
     #ratio with subtracted intercepts from potential
     #rzero: -r^2d / dr V(r) |_r = r_0  = c, with used c
     #puw plaquette determined from uwerr
     #job SLURM jobID
     #boot number of b ootsamples

    for (i in seq(1, length(crz))) {
#~         message(i)
#~ message(opt$beta, " a ", xi, " a ", opt$Nt, " a ", opt$Ns, " a ", plaquette,
#~          " a ", dp, " a ", slope, " a ", dslope, " a ", chislope,
#~          " a ", ratio, " a ", dratio, " a ", st, " a ", dst, " a ", chipot,
#~          " a ", rzeros[i], " a ", drzeros[i], " a ", crz[i],
#~          " a ", plaquetteuwerr, " a ", dpuwerr, " a ", opt$job, " a ", bootsamples)
        newline <- data.frame(beta = opt$beta, xiin = xi, Nt = opt$Nt,
         Ns = opt$Ns, p = plaquette, dp = dp, slope = slope,
         dslope = dslope, chislope = chislope, ratio = ratio, dratio = dratio,
         st = st, dst = dst, chipot = chipot, icslope = icslope,
         dicslope = dicslope, icpot = icpot, dicpot = dicpot,
         logpot = logpot, dlogpot = dlogpot, ratioslope = ratioslope,
         dratioslope = dratioslope, ratiopot = ratiopot, dratiopot = dratiopot,
         rzero = rzeros[i], drzero = drzeros[i], c = crz[i],
         puw = plaquetteuwerr, dpuw = dpuwerr, tauint = tauint, job = opt$job, hash = githash, boot = bootsamples,
         nom = length(measurements[, 1]) + opt$skip, skip = opt$skip)
        results <- rbind(results, newline)
    }
    results <- results[-1, ]
    print(results)
    filenamexi <- sprintf("%ssummarysmallbetaone%fL%d.csv",
            opt$plotpath, opt$betaone, opt$Ns)
    columnnames <- FALSE
    if (!file.exists(filenamexi)) {
        columnnames <- TRUE
    }
    write.table(results, filenamexi, append = TRUE,
            row.names = FALSE, col.names = columnnames)
    nameresults <- sprintf("%sresultspotNs%dNt%dbeta%fxi%fbs%d.RData",
            opt$plotpath, Ns, Nt, beta, xi, bootsamples)
    saveRDS(resultssummary, file = nameresults)
    }

if (FALSE) {
    # read in results of the temporal normal potential,
    # to compare ineger and non-integer distances,
    # plot both, once zoomed out and once zoomed in
    filenamefineinteger <- sprintf(
            "%sfitresultfinesubtractedb%fxi%fNs%dbsamples%domit0l0.RData",
            opt$plotpath, beta, xi, Ns, bootsamples)
    if (opt$smearing) {
        filenamefineinteger <- sprintf(
            "%sfitresultfinesubtractedb%fxi%fNs%dnape%dalpha%fbsamples%domit0l0.RData",
            opt$plotpath, beta, xi, Ns, nape, alpha, bootsamples)
    }
    resultfineinteger <- try(readRDS(filenamefineinteger))
    # zoomed out
    try(plot(fit.pot, ylab = "a_tV(r)", xlab = "r / a_s",
            main = "Potential including nonplanar points"))
    try(errorpolygon(X = seq(0, 7, by = 0.05),
            resultfineinteger, col.p = 2, pch = 2))
    try(legend(x = "topleft", legend = c("non-integer", "integer"),
            col = c(1, 2), pch = c(1, 2)))

    # zoomed in
    try(plot(resultfineinteger, ylab = "a_tV(r)", xlab = "r / a_s",
            main = "Potential including nonplanar points", col = 2, pch = 2))
    try(errorpolygon(X = seq(0, 8, by = 0.05), fit.pot, col.p = 1, pch = 1))
    try(legend(x = "topleft", legend = c("non-integer", "integer"),
            col = c(1, 2), pch = c(1, 2)))
}
}

## determine matrix element and potential from 1,1 plaquette
if (opt$amplitude) {

if (opt$fraction != 1) {
    stop("not all timeslices are present, correlators cannot be symmetrised and amplitude cannot be determined")
}

## read in boundaries for amplitude determination, take the same boundaries as for the effective mass
namelistbounds <- sprintf(
            "%stableanalysissmallNs%dNt%dbeta%fxi%f.csv",
            opt$plotpath, Ns, Nt, beta, xi)
if (opt$smearing) {
    namelistbounds <- sprintf(
            "%stableanalysissmallNs%dNt%dbeta%fxi%fnape%dalpha%f.csv",
            opt$plotpath, Ns, Nt, beta, xi, nape, alpha)
}

if (file.exists(namelistbounds)) {
    listbounds <- read.table(namelistbounds, sep = ",", header = TRUE)
} else {
    listbounds <- cbind(
                data.frame(x = c(rep(0, 4), rep(1, 4), rep(2, 4), rep(3, 4))),
                data.frame(y = rep(c(0, 1, 2, 3), 4)),
                data.frame(lower = rep(opt$t1, 16)),
                data.frame(upper = rep(Ntmax - opt$t1 - 2, 16)))
}
## open pdf for plot
filenameforplots <- sprintf(
            "%samplitudeNs%dNt%dbeta%fxi%fbs%d.pdf",
            opt$plotpath, Ns, Nt, beta, xi, bootsamples)
if (opt$smearing) {
    filenameforplots <- sprintf(
            "%samplitudeNs%dNt%dbeta%fxi%fnape%dalpha%fbs%d.pdf",
            opt$plotpath, Ns, Nt, beta, xi, nape, alpha, bootsamples)
}
pdf(file = filenameforplots, title = "")

## read in data into cf, symmetrise, determine amplitude
filename <- sprintf(
                "%sresult2p1d.u1potential.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fnonplanar",
                opt$respath, Nt, Ns, beta, xi, nape, alpha)
    WL <- readloopfilecfsmall(file = filename, skip = skip,
                    Nsmax = Nsmax, x = 1, y = 1, Ntmax = Ntmax, start = 1)
    WL <- bootstrap.cf(WL, boot.R = bootsamples, boot.l = 5)
    WL <- symmetrise.cf(WL)
    title <- sprintf("symmetrised correlators of the spatial 1,1 plaquette\n%d configs, skipped %d",
            (length(WL$cf[, 1]) + opt$skip - 1) * opt$nsave, opt$skip * opt$nsave)
    try(plot(WL, log = "y", xlab = "t/a_t", ylab = "C(x)",
            main = sprintf("%s, logscale", title)))
    WL.matrix <- matrixfit(WL, t1 = listbounds$lower[listbounds$x == 1 & listbounds$y == 1],
            t2 = max(listbounds$upper[listbounds$x == 1 & listbounds$y == 1], Ntmax/2), sym.vec = "cosh",
            parlist = array(c(1, 1), dim = c(2, 1)), neg.vec = c(1), useCov = TRUE)
    plot(WL.matrix)
    summary(WL.matrix)

## determine spatial-spatial and temporal-spatial plaquette with uwerr
        measurements <- read.table(filename, header = FALSE, skip = opt$skip)
        plaquettedatauwerr <- uwerrprimary(measurements[, Nsmax + 3])
        plaquetteuwerr <- plaquettedatauwerr$value
        dpuwerr <- plaquettedatauwerr$dvalue
        plaquettetempdatauwerr <- uwerrprimary(measurements[, (Nsmax + 1)^2 + 2])
        plaquettetempuwerr <- plaquettetempdatauwerr$value
        dptempuwerr <- plaquettetempdatauwerr$dvalue

    result <- data.frame(beta = beta, L = Ns, T = Nt, xiin = xi,
            E = WL.matrix$t0[1], dE = WL.matrix$se[1], A = WL.matrix$t0[2], dA = WL.matrix$se[2],
            chi = WL.matrix$chisqr / WL.matrix$dof, pval = WL.matrix$Qval,
            pss = plaquetteuwerr, dpss = dpuwerr, pst = plaquettetempuwerr, dpst = dptempuwerr,
            bs = bootsamples, hash = githash)
    filenameres <- sprintf("%smatrixfitb%fxi%fT%dL%d.RData", opt$plotpath, beta, xi, Nt, Ns)
    saveRDS(WL.matrix, filenameres)
    filenamexi <- sprintf("%samplitudesbetaone%f.csv",
            opt$plotpath, opt$betaone)
    columnnames <- FALSE
    if (!file.exists(filenamexi)) {
        columnnames <- TRUE
    }
    write.table(result, filenamexi, append = TRUE,
            row.names = FALSE, col.names = columnnames)
}
