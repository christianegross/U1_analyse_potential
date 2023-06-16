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
    make_option(c("--smearing"), action = "store_true", default = FALSE,
    help = "are the calculations done based on
            smeared lattices? (changes filenames) [default %default]"),

    make_option(c("--xidiff"), action = "store_true", default = FALSE,
    help = "Is xi different to L/T? [default %default]"),

    make_option(c("--respath"), type = "character", default = "",
    help = "path to where the resultfiles are stored [default %default]"),
    make_option(c("--plotpath"), type = "character", default = "",
    help = "path to where the plots are stored [default %default]"),
    make_option(c("--effmasstype"), type = "character", default = "log",
    help = "type of effective mass [default %default]"),

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

# boundaries for effective masses if no table can be found
t2 <- Ns / 2 - 2
t1 <- opt$lowerboundmeff
}


# set names for plot, tables, open, lists for saving results
filenameforplots <- sprintf(
            "%stauintplotsNt%dNs%dbeta%fxi%f.pdf",
            opt$plotpath, Nt, Ns, beta, xi)
if (opt$smearing) {
    filenameforplots <- sprintf(
            "%stauintplotsNt%dNs%dbeta%fxi%fnape%dalpha%f.pdf",
            opt$plotpath, Nt, Ns, beta, xi, nape, alpha)
}
pdf(file = filenameforplots, title = "")
listtauint <- list()
# determine Wilson Loop correlators and effective masses
# for all possible combinations Meff(x0/a_s) (use W(x0, t=0, y0+1)/W(x0, t=0, y0))
# and Meff(x0/a_t) (use W(x0, t0+1, y=0)/W(x0, t0, y=0))
# try method as in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.63.074501
# x here as z there, y here as x there
# measures a_s V(r / a_s)
# also determine number of correlators smaller than zero per distance
tauintdf <- data.frame()
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

    WL <- readloopfilecfrotated(file = filename, path = opt$respath, skip = skip,
            Nsmax = Ns * opt$fraction, yt = y, zerooffset = opt$zerooffset, every = opt$every, maxrows = opt$maxrows)

    autotime <- c()
    errors <- c()
    derrors <- c()
    for (boot.l in c(2, 4, 8, 16, 32, 64)) {
        WL <- bootstrap.cf(WL, boot.R = bootsamples, boot.l = boot.l)
            for (j in seq(1, opt$Ns*opt$fraction)){
                uwerrresults <- uwerr(data=WL$cf.tsboot$t[, j])
                autotime[j] <- uwerrresults$tauint
                errors[j] <- uwerrresults$dvalue
                derrors[j] <- uwerrresults$ddvalue
            }
        tauintdf <- rbind(tauintdf, data.frame(yt=y, spatial=T, l=boot.l, times=as.data.frame(t(autotime)),
                errors=as.data.frame(t(errors)), derrors=as.data.frame(t(derrors))))
    }
    # print(names(tauintdf))
    for (j in seq(1, opt$Ns*opt$fraction)){
        columnname <- paste("errors.V", j, sep="")
        dcolumnname <- paste("derrors.V", j, sep="")
        # print(columnname)
        # print(tauintdf[1:2, columnname])
        mask <- tauintdf$yt==y & tauintdf$spatial==T
        plotwitherror(x=1./tauintdf$l[mask], y=tauintdf[mask, columnname], dy=tauintdf[mask, dcolumnname], xlab="1/l", ylab="error", main=paste("x=", j, "y=", y, "Errors"))

        columnname <- paste("times.V", j, sep="")
        plotwitherror(x=1./tauintdf$l[mask], y=tauintdf[mask, columnname], xlab="1/l", ylab="tauint", main=paste("x=", j, "y=", y, "Autocorrelation time"))
        lines(x=seq(0, 1, by=0.1), y=rep(1, 11), col=2)
        }
    # listtauint[[y]] <- uwerrresults
}}
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
    WL <- readloopfilecfrotated(file = filename, path = opt$respath, skip = skip,
            Nsmax = Ns * opt$fraction, yt = t, zerooffset = opt$zerooffset, every = opt$every, maxrows = opt$maxrows)

    autotime <- c()
    errors <- c()
    derrors <- c()
    for (boot.l in c(2, 4, 8, 16, 32, 64)) {
        WL <- bootstrap.cf(WL, boot.R = bootsamples, boot.l = boot.l)
            for (j in seq(1, opt$Ns*opt$fraction)){
                uwerrresults <- uwerr(data=WL$cf.tsboot$t[, j])
                autotime[j] <- uwerrresults$tauint
                errors[j] <- uwerrresults$dvalue
                derrors[j] <- uwerrresults$ddvalue
            }
        tauintdf <- rbind(tauintdf, data.frame(yt=t, spatial=F, l=boot.l, times=as.data.frame(t(autotime)),
                errors=as.data.frame(t(errors)), derrors=as.data.frame(t(derrors))))
    }
    # print(names(tauintdf))
    for (j in seq(1, opt$Ns*opt$fraction)){
        columnname <- paste("errors.V", j, sep="")
        dcolumnname <- paste("derrors.V", j, sep="")
        # print(columnname)
        # print(tauintdf[1:2, columnname])
        mask <- tauintdf$yt==t & tauintdf$spatial==F
        plotwitherror(x=1./tauintdf$l[mask], y=tauintdf[mask, columnname], dy=tauintdf[mask, dcolumnname], xlab="1/l", ylab="error", main=paste("x=", j, "y=", y, "Errors"))

        columnname <- paste("times.V", j, sep="")
        plotwitherror(x=1./tauintdf$l[mask], y=tauintdf[mask, columnname], xlab="1/l", ylab="tauint", main=paste("x=", j, "y=", y, "Autocorrelation time"))
        lines(x=seq(0, 1, by=0.1), y=rep(1, 11), col=2)
        }
    # listtauint[[Ns / 2 + t]] <- uwerrresults
}}
# listtauint[[Ns/2 + Nt/2 + 1]] <- githash

# write out results
filenameuwerr <- sprintf(
            "%stauintNt%dNs%dbeta%fxi%fbsamples%d.csv",
            opt$plotpath, Nt, Ns, beta, xi, bootsamples)

if (opt$smearing) {
    filenameuwerr <- sprintf(
            "%stauintNt%dNs%dbeta%fxi%fnape%dalpha%fbsamples%d.csv",
            opt$plotpath, Nt, Ns, beta, xi, nape, alpha, bootsamples)
}
write.table(tauintdf, file = filenameuwerr, row.names=F, col.names=T)


mask <- tauintdf$spatial==T
plotwitherror(x=1./tauintdf$l[mask], y=tauintdf[mask, columnname], xlab="1/l", ylab="tauint", main=paste("Autocorrelation time spatial"))
mask <- tauintdf$spatial==F
plotwitherror(x=1./tauintdf$l[mask], y=tauintdf[mask, columnname], xlab="1/l", ylab="tauint", main=paste("Autocorrelation time temporal"))
