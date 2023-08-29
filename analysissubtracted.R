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
    help = "Ns of lattice [default %default]"),

    make_option(c("-t", "--Nt"), type = "integer", default = 16,
    help = "Nt of lattice [default %default]"),
    make_option(c("-n", "--nape"), type = "integer", default = 0,
    help = "Number of APE-smears that were done [default %default]"),

    make_option(c("-l", "--lowerboundmeff"), type = "integer", default = 0,
    help = "meff is determined in region
            [lowerboundmeff, Ns/2-2] [default %default]"),
    make_option(c("--rmin"), type = "integer", default = 2,
    help = "lowest radius used for calculation of xi [default %default]"),

    make_option(c("-m", "--rmax"), type = "integer", default = 4,
    help = "highest radius used for calculation of xi [default %default]"),
    make_option(c("-j", "--job"), type = "integer", default = 0,
    help = "qbig job number [default %default]"),

    make_option(c("-o", "--omit"), type = "integer", default = 0,
    help = "how many points should be omitted
            in calculating the potentials? [default %default]"),
    make_option(c("--lowlim"), type = "integer", default = 1,
    help = "points lower than this are not used
            for determining xi [default %default]"),
    make_option(c("-a", "--alpha"), type = "double", default = 1.0,
    help = "alpha used in APE-smearing [default %default]"),

    make_option(c("--betaone"), type = "double", default = 0,
    help = "input beta at corresponding xi = 1 [default %default]"),
    make_option(c("-x", "--xi"), type = "double", default = 0,
    help = "xi used in lattice [default %default]"),
    make_option(c("--fraction"), type = "double", default = 0.5,
    help = "Maximal extent of the Wilson-Loops as a fraction of the lattice size [default %default]"),

    make_option(c("--aic"), action = "store_true", default = FALSE,
    help = "if true, effective masses are determined with weighting
        according to the akaike information criterion [default %default]"),
    make_option(c("--scaletauint"), action = "store_true", default = FALSE,
    help = "if true, errors and bootstrapsamples for the correlator are rescaled
    such that the effects of autocorrelation are taken into account [default %default]"),

    make_option(c("--analyse"), action = "store_true", default = FALSE,
    help = "if true, correlators and effective masses
            are determined [default %default]"),

    make_option(c("--drawbootstrap"), action = "store_true", default = FALSE,
    help = "Draw bootstrap samples and determine effective mass
            before doing effective mass fits [default %default]"),
    make_option(c("--dofit"), action = "store_true", default = FALSE,
    help = "if true, potentials and xi are calculated [default %default]"),

    make_option(c("--smearing"), action = "store_true", default = FALSE,
    help = "are the calculations done based on smeared lattices?
            (changes filenames) [default %default]"),
    make_option(c("--xidiff"), action = "store_true", default = FALSE,
    help = "Is xi different to Ns/Nt? [default %default]"),
    make_option(c("--zerooffset"), type = "integer", default = 0,
    help = "Offset for selecting the correct loops if W(x = 0, y = 0)
            was measured [default %default]"),

    make_option(c("-e", "--every"), type = "integer", default = 1,
    help = "only reads every eth line [default %default]"),
    make_option(c("--bootl"), type = "integer", default = 2,
    help = "block length in bootstrapping configurations [default %default]"),

    make_option(c("--plotonlymeff"), action = "store_true", default = FALSE,
    help = "The plots of the fitted meff are plotted
            seperately [default %default]"),
    make_option(c("--usecov"), action = "store_true", default = FALSE,
    help = "use covariance when fitting effective masses"),
    make_option(c("--respath"), type = "character", default = "",
    help = "path to where the resultfiles are stored [default %default]"),

    make_option(c("--plotpath"), type = "character", default = "",
    help = "path to where the plots are stored [default %default]"),
    make_option(c("--effmasstype"), type = "character", default = "log",
    help = "type of effective mass [default %default]"),

    make_option(c("--nsave"), type = "integer", default = 0,
    help = "steps between saved configs [default %default]"),
    make_option(c("--crzero"), type = "double", default = -1.65,
    help = "c used for determining r_0 [default %default]"),
    make_option(c("--myfunctions"), type = "character",
    default = "/hiskp4/gross/masterthesis/analyse/code/U1_analyse_potential/",
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
# boundaries for effective masses if no table can be found


if (opt$scaletauint) {
        opt$usecov <- TRUE
        opt$bootl <- 1
}

t2 <- Ns / 2 - 2
t1 <- opt$lowerboundmeff
rmin <- opt$rmin
rmax <- opt$rmax
endinganalysis <- sprintf("Nt%dNs%dbeta%fxi%fbootl%dusecov%d", Nt, Ns, beta, xi, opt$bootl, opt$usecov)
if (opt$smearing) endinganalysis <- sprintf("Nt%dNs%dbeta%fxi%fnape%dalpha%fbootl%dusecov%d", Nt, Ns, beta, xi, nape, alpha, opt$bootl, opt$usecov)

endingdofit <- sprintf("Nt%dNs%dbeta%fxi%fbs%domit%dlowlim%d",Nt, Ns, beta, xi, opt$bootsamples, opt$omit, opt$lowlim)
if (opt$smearing) endingdofit <- sprintf("Nt%dNs%dbeta%fxi%fnape%dalpha%fbs%domit%dlowlim%d",Nt, Ns, beta, xi, nape, alpha, opt$bootsamples, opt$omit, opt$lowlim)


if(opt$aic){
        endinganalysis <- sprintf("%saic", endinganalysis)
        endingdofit <- sprintf("%saic", endingdofit)
}

if(opt$scaletauint){
        endinganalysis <- sprintf("%sscaletauint", endinganalysis)
        endingdofit <- sprintf("%sscaletauint", endingdofit)
}

}


if (analyse) {
# set names for plot, tables, open, lists for saving results

filenameforplots <- sprintf("%splotsmeffnormal%s.pdf", opt$plotpath, endinganalysis)
if (!opt$drawbootstrap) filenameforplots <- sprintf("%splotsmeffnormalnc%s.pdf", opt$plotpath, endinganalysis)

pdf(file = filenameforplots, title = "")
listfits <- list()
listtauint <- list()
potential <- data.frame(R = NA, m = NA, dm = NA, space = NA, p = NA, chi = NA)
if (opt$aic) {
    potential <- data.frame(R = NA, space = NA, median = NA, errtotmean = NA, m16 = NA, m84 = NA,
                        bootmedian = NA, booterrtot = NA, bootstat = NA, bootsys = NA, bias = NA)
}

namelistbounds <- sprintf(
            "%stableanalysissubtractedNs%dNt%dbeta%fxi%f.csv",
            opt$respath, Ns, Nt, beta, xi)
if (opt$smearing) {
    namelistbounds <- sprintf(
            "%stableanalysissubtractedNs%dNt%dbeta%fxi%fnape%dalpha%f.csv",
            opt$respath, Ns, Nt, beta, xi, nape, alpha)
}
if (file.exists(namelistbounds)) {
    listbounds <- read.table(namelistbounds, sep = ",", header = TRUE)
} else {
    # if no table with boundaries can be found, create a dummy table
    listbounds <- cbind(data.frame(yt = c(seq(1, Ns / 2), seq(1, Ns / 2))),
            data.frame(spacial = c(rep(TRUE, Ns / 2), rep(FALSE, Ns / 2))),
            data.frame(lower = c(rep(opt$lowerboundmeff, Ns / 2),
                    rep(floor(opt$lowerboundmeff / xi), Ns / 2))),
            data.frame(upper = c(rep(Ns / 2 - 2, Ns / 2),
                    rep(Nt / 2 - 2, Ns / 2))))
}
print(namelistbounds)
print(listbounds)
negatives <- c()

filename <- sprintf(
            "%sresult2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fcoarsedistance",
            opt$respath, Nt, Ns, beta, xi, nape, alpha)
alldata <- read.table(filename)
nsave <- alldata[2, length(alldata[1, ])] - alldata[1, length(alldata[1, ])]

## if the correlators and effective masses are not determined in-place, read in previous results
if (!opt$drawbootstrap) {
    readinres <- readRDS(sprintf("%slistmeffnormal%s.RData", opt$plotpath, endinganalysis))
}

# determine Wilson Loop correlators and effective masses
# for all possible combinations Meff(y0) (use W(x0+1, t=0, y0)/W(x0, t=0, y0))
# and Meff(t0) (use W(x0+1, t0, y=0)/W(x0, t0, y=0))
# measures a_s V(r / a_s) and a_tV(r / a_s)
# also determine number of correlators smaller than zero per distance
# for each x: first spatial and then temporal potential
# for all: determine expectation values of Loops,
# plot, determine effective masses, set boundaries, determine
# plateau and plot. Savev results of m_eff to list
# also save uwerr = autocorrelation time of non-bootstrapped loop


for (x in seq(1, Ns / 2, 1)) {
    #save results, list has to have an entry everywhere
    listresults <- list(FALSE, FALSE, FALSE)
    uwerrresults <- FALSE
    #calculate W(x, y + 1)/W(x, y)
    if (TRUE) {
# spatial, loops
    filename <- sprintf(
            "result2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fcoarsedistance",
            Nt, Ns, beta, xi, nape, alpha)

    title <- sprintf("beta = %.2f coarse x = %d skipped %d",
                    beta, x, skip * opt$nsave)

    if (opt$drawbootstrap) {
        WL <- calcplotWloopnormalspatial(file = filename, skip = skip,
                Ns = Ns, x = x, bootsamples = bootsamples,
                title = title, path = opt$respath, zerooffset = opt$zerooffset,
                every = opt$every, nsave = opt$nsave, l = opt$bootl, fraction = opt$fraction,
                maxrows = opt$maxrows)
        uwerrresults <- uwerr.cf(WL)
        negatives[x] <- sum(WL$cf0 < 0)

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
        WL <- readinres[[x]][[1]]
        if(opt$aic) WL <- readinres[[x]][[1]]$effmass
    }


# m_eff
    if (!opt$aic) {
    t1 <- listbounds$lower[listbounds$spacial == TRUE & listbounds$yt == x]
    t2 <- listbounds$upper[listbounds$spacial == TRUE & listbounds$yt == x]
    WL.effmasslist <- deteffmass(WL.effmass = WL, t1 = t1, t2 = t2, yt = x,
        isspatial = 1, usecov = opt$usecov)
    #plot results, save results off meff in
    #short and long form (with and without bootstrapsamples)
    if (WL.effmasslist[[2]][[2]] != 0) {
        listresults <- list(WL.effmasslist[[1]], x, TRUE)
        potential <- rbind(potential, WL.effmasslist[[3]])
    }
    names(listresults) <- c("effmass", "x", "spacial")

    try(plot(WL.effmasslist[[1]], xlab = "y/a_s", ylab = "Meff",
            main = sprintf("%s, t1 = %d, t2 = %d, p = %f",
            title, t1, t2, WL.effmasslist[[1]]$effmassfit$Qval)))
    try(plot(WL.effmasslist[[1]], xlab = "y/a_s", ylab = "Meff",
            main = sprintf("%s, t1 = %d, t2 = %d, p = %f",
            title, t1, t2, WL.effmasslist[[1]]$effmassfit$Qval),
            ylim = WL.effmasslist[[2]]))
    }
    if (opt$aic) {
        effmass <- deteffmassaic(WL)
        plot(effmass$effmass, xlab="y/a_s", ylab="meff", main=title)
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
        listresults <- list(effmass, x, TRUE, uwerrresults)

        potentialnew <- data.frame(R = x, space = TRUE, median = effmass$effmassfit$t0, errtotmean = effmass$effmassfit$se,
                m16 = effmass$effmassfit$m16, m84 = effmass$effmassfit$m84,
                bootmedian = effmass$boot$mean, booterrtot = effmass$boot$errtot, bootstat = effmass$boot$errstat,
                bootsys = sqrt(effmass$boot$errtot^2 - effmass$boot$errstat^2),
                bias = effmass$effmassfit$t0 - effmass$boot$mean)
        potential <- rbind(potential, potentialnew)
    }
}
    listfits[[x]] <- listresults
    if (opt$drawbootstrap) listtauint[[x]] <- uwerrresults


    #calculate W(x, t + 1)/W(x, t)
    if (TRUE) {
# temporal, loops

    filename <- sprintf(
        "result2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%ffinedistance",
        Nt, Ns, beta, xi, nape, alpha)

    title <- sprintf("beta = %.2f fine x = %d skipped %d",
            beta, x, skip * opt$nsave)
    if (opt$drawbootstrap) {
        WL <- calcplotWloopnormaltemporal(file = filename, skip = skip, Ns = Ns,
                Nt = Nt, x = x, bootsamples = bootsamples, title = title,
                path = opt$respath, zerooffset = opt$zerooffset, every = opt$every,
                nsave = opt$nsave, l = opt$bootl, fraction = opt$fraction,
                maxrows = opt$maxrows)
        uwerrresults <- uwerr.cf(WL)
        negatives[x] <- negatives[x] + sum(WL$cf0 < 0)

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

        plot(WL, log = "y", xlab = "t/a_t", ylab = "C(t)",
                main = sprintf("%s, logscale", title))
        WL <- bootstrap.effectivemass(WL, type = opt$effmasstype)
    } else {
        WL <- readinres[[Ns / 2 + x]][[1]]
        if (opt$aic) WL <- readinres[[Ns / 2 + x]][[1]]$effmass
    }


# m_eff
    if(!opt$aic) {
    t1 <- listbounds$lower[listbounds$spacial == FALSE & listbounds$yt == x]
    t2 <- listbounds$upper[listbounds$spacial == FALSE & listbounds$yt == x]
    WL.effmasslist <- deteffmass(WL.effmass = WL, t1 = t1, t2 = t2, yt = x,
        isspatial = 0, usecov = opt$usecov)
    #plot results, save results off meff in
    #short and long form (with and without bootstrapsamples)
    if (WL.effmasslist[[2]][[2]] != 0) {
        listresults <- list(WL.effmasslist[[1]], x, TRUE)
        potential <- rbind(potential, WL.effmasslist[[3]])
    }
    names(listresults) <- c("effmass", "x", "spacial")

    try(plot(WL.effmasslist[[1]], xlab = "t/a_t", ylab = "Meff",
            main = sprintf("%s, t1 = %d, t2 = %d, p = %f",
            title, t1, t2, WL.effmasslist[[1]]$effmassfit$Qval)))
    try(plot(WL.effmasslist[[1]], xlab = "t/a_t", ylab = "Meff",
            main = sprintf("%s, t1 = %d, t2 = %d, p = %f",
            title, t1, t2, WL.effmasslist[[1]]$effmassfit$Qval),
            ylim = WL.effmasslist[[2]]))
    }
    if (opt$aic) {
        effmass <- deteffmassaic(WL)
        plot(effmass$effmass, xlab="t/a_t", ylab="meff", main=title)
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
        listresults <- list(effmass, x, TRUE, uwerrresults)

        potentialnew <- data.frame(R = x, space = FALSE, median = effmass$effmassfit$t0, errtotmean = effmass$effmassfit$se,
                m16 = effmass$effmassfit$m16, m84 = effmass$effmassfit$m84,
                bootmedian = effmass$boot$mean, booterrtot = effmass$boot$errtot, bootstat = effmass$boot$errstat,
                bootsys = sqrt(effmass$boot$errtot^2 - effmass$boot$errstat^2),
                bias = effmass$effmassfit$t0 - effmass$boot$mean)
        potential <- rbind(potential, potentialnew)
    }
}
    listfits[[Ns / 2 + x]] <- listresults
    if (opt$drawbootstrap) listtauint[[Ns / 2 + x]] <- uwerrresults
    if (opt$drawbootstrap) message(negatives[x], " correlators are negative for x = ", x)
}
listfits[[Ns + 1]] <- githash

#write out results
t1 <- -1 #opt$lowerboundmeff
filenamepotential <- sprintf("%spotmeffnormal%s.csv",  opt$plotpath, endinganalysis)
filenamelist <- sprintf("%slistmeffnormal%s.RData", opt$plotpath, endinganalysis)
filenameuwerr <- sprintf("%slistuwerrtauintnormal%s.RData", opt$plotpath, endinganalysis)
filenamenegatives <- sprintf("%snegativesnormal%s.csv", opt$plotpath, endinganalysis)

write.table(potential, filenamepotential, row.names = FALSE)
saveRDS(listfits, file = filenamelist)
if (opt$drawbootstrap) {
    saveRDS(listtauint, file = filenameuwerr)
    message(sum(negatives), " correlators are negative overall")
    write.table(data.frame(x = seq(1, Ns / 2), neg = negatives),
        file = filenamenegatives, row.names = FALSE)
}

}

if (dofit) {
# use data from the potential to fit the expected form to the fine=temporal
# and coarse=spatial potential,
# then determine xi_ren by interpolation and subtraction,
# like in https://journals.aps.org/prd/pdf/10.1103/PhysRevD.70.014504
# then determine r_0 / a_s, determine the force, plot everything
# save results and bootstrapsamples of fits

#define functions for potential, force, force  =  -r^2 * derivation of potential
fnpot <- function (par, x, boot.r, ...) par[1] + par[2] * x + par[3] * log(x)
fnpotscal <- function (par, x, boot.r, finemask, ...) {
        return (par[1] + par[2] * x + par[3] * log(x) + par[4] * x[finemask])
}
fnforce <- function (par, x, boot.r, ...) {
        return ((- 1.0) * par[2] * x^2 - par[3] * x)
}
fnforceerr <- function (parerr, x, boot.r, correlation, ...) {
    return (sqrt(parerr[2]^2 * x^4 + parerr[3]^2 * x^2 + 2 * x^3 * correlation[2, 3] * parerr[2] * parerr[3]))
}

t1 <- opt$lowerboundmeff
filenamelist <- sprintf("%slistmeffnormal%s.RData", opt$plotpath, endinganalysis)


listmeff <- readRDS(file = filenamelist)
if (!opt$aic) bootsamples <- length(listmeff[[1]][[1]]$massfit.tsboot[, 1])
if (opt$aic) bootsamples <- listmeff[[1]][[1]]$effmass$boot.R

filenameforplots <- sprintf("%spotentialnormal%s.pdf", opt$plotpath, endingdofit)

pdf(file = filenameforplots, title = "")

# calculate plaquette as W(x = 1, y = 1, t = 0)
# with uwerr and with cf, to get bootstrapsamples
filename <- sprintf(
            "%sresult2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fcoarsedistance",
            opt$respath, Nt, Ns, beta, xi, nape, alpha)

alldata <- read.table(filename)
plaquettecolumn <- alldata[, (opt$zerooffset + Ns / 2) * opt$zerooffset + 1 + opt$zerooffset]
plaquettedata <- uwerrprimary(plaquettecolumn[seq(skip + 1, length(plaquettecolumn), opt$every)])
nom <- floor(length(plaquettecolumn))

column <- (opt$zerooffset + Ns / 2) * opt$zerooffset + 1 + opt$zerooffset

plaquettecf <- readloopfilecfonecolumn(file = filename, skip = skip,
        column = column, every = opt$every)
bootl <- max(opt$bootl, 4 * plaquettedata$tauint^2)
plaquettecf <- bootstrap.cf(plaquettecf,
                boot.R = bootsamples, boot.l = bootl)



#coarse potential
#read in data points, bootstrap samples: have to initialise empty vectors,
#because not every list element has a valid entry
xc <- rep(NA, Ns / 2)
yc <- rep(NA, Ns / 2)

# filter out points without valid entry in bootstrap.nlsfit,
# mask for combining vectors
maskc <- rep(FALSE, Ns / 2)
bsamplesc <- array(c(rep(NA, Ns / 2 * bootsamples)), dim = c(bootsamples, Ns / 2))

# also have vectors for combined potential
x <- rep(NA, Ns)
dx <- rep(NA, Ns)
y <- rep(NA, Ns)
mask <- rep(FALSE, Ns)
finemask <- rep(FALSE, Ns)
bsamples <- array(c(rep(NA, Ns * bootsamples)), dim = c(bootsamples, Ns))

for (i in seq(1, Ns / 2 - opt$omit, 1)) {
  if (listmeff[[i]][[2]]) {
      xc[i] <- listmeff[[i]][[2]]
      yc[i] <- listmeff[[i]][[1]]$effmassfit$t0[1]
      if(!opt$aic) {
        bsamplesc[, i] <- listmeff[[i]][[1]]$massfit.tsboot[, 1]
        bsamples[, i] <- listmeff[[i]][[1]]$massfit.tsboot[, 1]
      }
      else if (opt$aic) {
        bsamplesc[, i] <- listmeff[[i]][[1]]$boot$m50 - (listmeff[[i]][[1]]$effmassfit$t0 - listmeff[[i]][[1]]$boot$m50) * (listmeff[[i]][[1]]$effmassfit$se / listmeff[[i]][[1]]$boot$errstat - 1)
        bsamples[, i] <- bsamplesc[, i]
      }
      maskc[i] <- TRUE

      x[i] <- listmeff[[i]][[2]]
      y[i] <- listmeff[[i]][[1]]$effmassfit$t0[1]
      mask[i] <- TRUE
      finemask[i] <- FALSE
  }
}


#fine potential
#read in data points, bootstrap samples: have to initialise empty vectors,
#because not every list element has a valid entry
xf <- rep(NA, Ns / 2)
yf <- rep(NA, Ns / 2)

#filter out points without valid entry in bootstrap.nlsfit,
# finemask: needed later for combining potentials
maskf <- rep(FALSE, Ns / 2)
bsamplesf <- array(c(rep(NA, Ns / 2 * bootsamples)), dim = c(bootsamples, Ns / 2))

for (i in seq(1, Ns / 2 - opt$omit, 1)) {
  if (listmeff[[Ns / 2 + i]][[2]]) {
      xf[i] <- listmeff[[Ns / 2 + i]][[2]]
#~       message("xc =  ",xc[i])
      yf[i] <- listmeff[[Ns / 2 + i]][[1]]$effmassfit$t0[1]
      maskf[i] <- TRUE
      if (!opt$aic) {
        bsamplesf[, i] <- listmeff[[Ns / 2 + i]][[1]]$massfit.tsboot[, 1]
        bsamples[, Ns / 2 + i] <- listmeff[[Ns / 2 + i]][[1]]$massfit.tsboot[, 1]
      }
      else if (opt$aic) {
        bsamplesf[, i] <- listmeff[[Ns / 2 + i]][[1]]$boot$m50 - (listmeff[[Ns / 2 + i]][[1]]$effmassfit$t0 - listmeff[[Ns / 2 + i]][[1]]$boot$m50) * (listmeff[[Ns / 2 + i]][[1]]$effmassfit$se / listmeff[[Ns / 2 + i]][[1]]$boot$errstat - 1)
        bsamples[, Ns / 2 + i] <- bsamplesf[, i]
      }

      x[Ns / 2 + i] <- listmeff[[Ns / 2 + i]][[2]]
      y[Ns / 2 + i] <- listmeff[[Ns / 2 + i]][[1]]$effmassfit$t0[1]
      mask[Ns / 2 + i] <- TRUE
      finemask[Ns / 2 + i] <- TRUE
  }
}

#determine parameters of potential by bootstrap, save results
fit.resultcoarse <- bootstrap.nlsfit(fnpot, c(0.2, 0.2, 0.2),
                    yc, xc, bsamplesc, mask = maskc)
filenamecoarse <- sprintf("%sfitresultcoarsenormal%s.RData", opt$plotpath, endingdofit)
saveRDS(fit.resultcoarse, file = filenamecoarse)

fit.resultfine <- bootstrap.nlsfit(fnpot, c(0.2, 0.2, 0.2),
                    yf, xf, bsamplesf, mask = maskf)
filenamefine <- sprintf("%sfitresultfinenormal%s.RData", opt$plotpath, endingdofit)
saveRDS(fit.resultfine, file = filenamefine)

# plot potentials
plot(fit.resultfine, main  =  sprintf(
            "beta = %f, xi = %f, (2 + 1)D, Ns = %d, fine\n
            %d measurements of which skipped %d",
            beta, xi, Ns, nom * opt$nsave, skip * opt$nsave),
            ylab = "a_t V_{xt}(x)", xlab = "x / a_s")
plot(fit.resultcoarse, main  =  sprintf(
            "beta = %f, xi = %f, (2 + 1)D, Ns = %d, coarse\n
            %d measurements of which skipped %d",
            beta, xi, Ns, nom * opt$nsave, skip * opt$nsave),
            ylab = "a_s V_{xy}(x)", xlab = "x / a_s")



alldata <- data.frame(x = x, V = y, spatial = finemask)


xilist <- data.frame(xi = NA, rdiff = NA, r1 = NA)

# determine xi by averaging over all possible distances of the potentials
# between r_min and r_max
# like in https://journals.aps.org/prd/pdf/10.1103/PhysRevD.70.014504
for (bs in seq(1, bootsamples)) {
    for (r1 in seq(opt$rmin, opt$rmax - 1)) {
        for (r2 in seq(r1 + 1, opt$rmax)) {
            num <- (bsamples[bs, Ns / 2 + r1] - bsamples[bs, Ns / 2 + r2])
            den <- (bsamples[bs, r1]  - bsamples[bs, r2])
            xicalc <- (num  /  den)
            xilist <- rbind(xilist,
                    data.frame(xi = xicalc, rdiff = abs(r1 - r2), r1 = r1))
        }
}
}
xilist <- na.omit(xilist)

xicalcsub <- mean(xilist$xi)
dxicalcsub <- sd(xilist$xi)

xisquaredsub <- mean(xilist$xi^2)
dxisquaredsub <- sd(xilist$xi^2)

# determine xi by fitting a linear function to V_spatial (V_temporal)
matchpot <- function(par, x, boot.r, ...) par[1] + x / par[2]

bsamplesmatch <- array(c(rep(NA, Ns * bootsamples)), dim = c(bootsamples, Ns))

for (bs in seq(1, bootsamples, 1)) {
    bsamplesmatch[bs, 1:(Ns / 2)] <- bsamples[bs, (Ns / 2 + 1):Ns]
    bsamplesmatch[bs, (Ns / 2 + 1):Ns] <- bsamples[bs, 1:(Ns / 2)]
}

maskc <- maskc & c(rep(FALSE, opt$lowlim), rep(TRUE, Ns / 2 - opt$lowlim))

fit.match <- bootstrap.nlsfit(matchpot, c(0.1, xi),
            y = yc, x = yf, bsamples[, 1:Ns], mask = maskc)
plot(fit.match, xlab = "a_tV_s(x)", ylab = "a_sV_s(x)",
            main = "matching potentials to determine xi")
# print(fit.match)
filenamematch <- sprintf("%sfitresultmatchnormal%s.RData", opt$plotpath, endingdofit)
saveRDS(fit.match, file = filenamematch)

# rescale temporal potential
xicalc <- fit.match$t0[2]
dxicalc <- fit.match$se[2]
xisquared <- mean(fit.match$t^2)
dxisquared <- sd(fit.match$t^2)

y[finemask] <- y[finemask] / fit.match$t0[2] + fit.match$t0[1]
for (i in seq(1, bootsamples, 1)) {
    bsamples[i, finemask] <- bsamples[i, finemask] / fit.match$t[i, 2] + fit.match$t[i, 1]
}

# fit to overall potential
fit.resultscaled <- bootstrap.nlsfit(fnpot, c(0.1, 0.1, 0.1),
                    y, x, bsamples, mask = mask)
filenamescaled <- sprintf("%sfitresultscalednormal%s.RData", opt$plotpath, endingdofit)
saveRDS(fit.resultscaled, file = filenamescaled)

#if wanted determine correlation matrix seperately, as determined by the summary.bootstrapfit function
if (!is.null(fit.resultscaled$cov)) {
    cov.to.cor <- diag(1 / fit.resultscaled$se)
    correlation <- cov.to.cor % * % fit.resultscaled$cov % * % cov.to.cor
}else if (!is.null(fit.resultscaled$t)) {
    correlation <- cor(fit.resultscaled$t, fit.resultscaled$t, use = "na.or.complete")
}


#graphically represent force: for each bootstrap sample of paramters, determine value of force on a sequence,
#determine error as standard deviation of values on each point of the sequence
xx <- seq(-0.2, Ns / 2 + 0.2, 0.05)
bootsforce <- data.frame(matrix(rep(NA, length(xx)), ncol = length(xx), nrow = 1))

for (i in seq(1, bootsamples, 1)) {
    bootsforce <- rbind(bootsforce, fnforce(fit.resultscaled$t[i, ], xx, 0))
}
bootsforce <- na.omit(bootsforce)

forceerrs <- c()
for (i in seq(1, length(xx), 1)) {
    forceerrs[i] <- sd(bootsforce[, i])
}


#Determine r0 as solution of equation -r^2 d/dr V(r) = c, solve for each bootstrapsample, r0 = mean pm sd
# V(r) = a + sigma * r + b * ln(r)
# -d/dr V(r)  =  -sigma -b/r
# c  =  -sigma * r^2 - b * r, use p-q-formula for p = b/sigma, q = c/sigma
rzeroofc <- data.frame(r0 = NA, dr0 = NA, c = NA)

for (c in list(opt$crzero)) {
rzerolist <- data.frame(r1 = NA, r2 = NA)
for (bs in seq(1, bootsamples)) {
#~     a <- fit.resultscaled$t[bs, 1]
    sigma <- fit.resultscaled$t[bs, 2]
    b <- fit.resultscaled$t[bs, 3]
    r1 <- -b / (2 * sigma) + sqrt((b / (2 * sigma))^2 - c / sigma)
    r2 <- -b / (2 * sigma) - sqrt((b / (2 * sigma))^2 - c / sigma)
    rzerolist <- rbind(rzerolist, data.frame(r1 = r1, r2 = r2))
}
rzerolist <- na.omit(rzerolist)

rzero <- mean(rzerolist$r1)
drzero <- sd(rzerolist$r1)
rzeroofc <- rbind(rzeroofc, data.frame(r0 = rzero, dr0 = drzero, c = c))
bsrzero <- rzerolist[, 1]
}
rzeroofc <- na.omit(rzeroofc)

title <- sprintf(
        "beta = %.3f, xi = %.3f + /-%.3f, (2 + 1)D, Ns = %d\n
        %d measurements of which skipped %d",
        beta, xicalc, dxicalc, Ns, nom * opt$nsave, skip * opt$nsave)


#plot visible combination of temporal and spatial potential
plot(fit.resultscaled, main  =  title, ylab = "V(r)", xlab = "r/a_s",
            pch = 3, cex = 0.1)
pointsyerr(x = fit.resultscaled$x[finemask], y = fit.resultscaled$y[finemask],
            dy = fit.resultscaled$dy[finemask], col = 1, pch = 1)
pointsyerr(x = fit.resultscaled$x[!finemask], y = fit.resultscaled$y[!finemask],
            dy = fit.resultscaled$dy[!finemask], col = 2, pch = 2)
legend(legend = c("spatial", "temporal"), col = c(2, 1),
            pch = c(2, 1), x = "topleft")

#force with error and r0, zoomed in

polyval <- c(fnforce(fit.resultscaled$t0, xx, 0) + forceerrs,
            rev(fnforce(fit.resultscaled$t0, xx, 0) - forceerrs))
pcol <- col2rgb("gray", alpha = TRUE) / 255
pcol[4] <- 0.65
pcol <- rgb(red = pcol[1], green = pcol[2], blue = pcol[3], alpha = pcol[4])
xlim <- c(0.6 * min(rzeroofc$r0), 1.4 * max(rzeroofc$r0))
ylim <- c(0.6 * min(abs(rzeroofc$c)), 1.4 * max(abs(rzeroofc$c)))
ylim <- rev(-ylim)
lowylim <- -2 * max(abs(rzeroofc$c))
plot(xx, fnforce(fit.resultscaled$t0, xx, 0), type = "l", xlab = "r/a_s",
            ylab = "-r^2 d/dr V(r)", main  =  title, xlim = xlim, ylim = ylim)
polygon(x = c(xx, rev(xx)), y = polyval, col = pcol, lty = 0,
            lwd = 0.001, border = pcol)

for (i in seq(1, length(rzeroofc$c))) {
    c <- rzeroofc$c[i]
    rzero <- rzeroofc$r0[i]
    drzero <- rzeroofc$dr0[i]
#also make lines for position of r0, polygon for error of r0
#arrows(xhi, yhi, xlo, ylo)
    arrows(rzero, fnforce(fit.resultscaled$t0, rzero, 0), rzero, lowylim,
                    angle = 90, length = 0.1, code = 0)
    arrows(rzero, c, 0, c, angle = 90, length = 0.1, code = 0)
    r0polyval1 <- c(fnforce(fit.resultscaled$t0, rzero - drzero, 0) - fnforceerr(fit.resultscaled$se,
                    rzero - drzero, 0, correlation),
                    fnforce(fit.resultscaled$t0, rzero + drzero, 0) + fnforceerr(fit.resultscaled$se,
                    rzero + drzero, 0, correlation))
    r0polyval1 <- c(r0polyval1, lowylim, lowylim)
    r0polyval2 <- c(rzero - drzero, rzero + drzero)
    r0polyval2 <- c(r0polyval2, rev(r0polyval2))
    polygon(x = r0polyval2, y = r0polyval1, col = pcol,
            lty = 0, lwd = 0.001, border = pcol)
}

#plot thermalisation
plot(plaquettecolumn, main = "Thermalisation",
        xlab = sprintf("MCMC-steps/%d", opt$nsave), ylab = "P")

}

if (dofit) {
    # save bootstrapsamples and summary of results
resultssummary <- list(st = fit.resultscaled$t0[2], dst = fit.resultscaled$se[2], bsst = fit.resultscaled$t[, 2],
        rzeros = rzero, drzeros = drzero, crz = opt$crzero, bsrzeros = bsrzero,
        p = plaquettecf$cf.tsboot$t0[1], dp = plaquettecf$tsboot.se[1], bsp = plaquettecf$cf.tsboot$t[, 1],
        beta = beta, xiin = xi, Nt = Nt, Ns = Ns, bootsamples = bootsamples,
        nom = nom, skip = opt$skip, lowlim = opt$lowlim,
        xicalc  =  fit.match$t0[2], dxicalc  =  fit.match$se[2], bsxicalc = fit.match$t[, 2])
class(resultssummary) <- "resultssummary"


resultlist <- data.frame(xi = NA, beta = NA, xicalc = NA, dxicalc = NA,
    xi2 = NA, dxi2 = NA, r0 = NA, dr0 = NA, st = NA, dst = NA, p = NA,
    dp = NA, chi = NA, diff = NA, ddiff = NA, bs = NA,
    xicalcsub = NA, dxicalcsub = NA, xi2sub = NA, dxi2 = NA,
    c = NA, rmin = NA, rmax = NA, Ns = NA, Nt = NA, nape = NA, alpha = NA,
    omit = NA, nom = NA, skip = NA, t1 = NA, job = NA, hash = NA,
    every = NA, tauint = NA, dtauint = NA, bootl = NA, lowlim = NA, aic=NA, scaletauint=NA)

for (i in seq(1, max(1, length(rzeroofc$c)))) {
newline <- data.frame(xi = xi, beta = beta, xicalc = xicalc, dxicalc = dxicalc,
                      xi2 = xisquared, dxi2 = dxisquared,
                      r0 = rzeroofc$r0[i], dr0 = rzeroofc$dr0[i],
                      st = fit.resultscaled$t0[2], dst = fit.resultscaled$se[2],
                      p = plaquettecf$cf.tsboot$t0[1], dp = plaquettecf$tsboot.se[1],
                      chi = fit.resultscaled$chisqr / fit.resultscaled$dof,
                      diff = fit.match$t0[1], ddiff = fit.match$se[1],
                      bs = bootsamples, xicalcsub = xicalcsub,
                      dxicalcsub = dxicalcsub, xi2sub = xisquaredsub,
                      dxi2 = dxisquaredsub, c = rzeroofc$c[i], rmin = opt$rmin,
                      rmax = opt$rmax, Ns = Ns, Nt = Nt, nape = nape, alpha = alpha,
                      omit = opt$omit, nom = nom, skip = skip,
                      t1 = t1, job = opt$job, hash = githash, every = opt$every,
                      tauint = plaquettedata$tauint, dtauint = plaquettedata$dtauint,
                      bootl = opt$bootl, lowlim = opt$lowlim, aic=opt$aic, scaletauint=opt$scaletauint)
resultlist <- rbind(resultlist, newline)
}

resultlist <- resultlist[-1, ]
print(resultlist)
filename <- sprintf("%sresultsummary2p1dnormalb%.3fNs%d.csv",
        opt$plotpath, opt$betaone, opt$Ns)

columnnames <- FALSE
if (!file.exists(filename)) {
    columnnames <- TRUE
}

write.table(resultlist, filename, append = TRUE,
        row.names = FALSE, col.names = columnnames)
nameresults <- sprintf("%sresultsnormal%s.RData", opt$plotpath, endingdofit)
saveRDS(resultssummary, file = nameresults)
}

message(warnings())