library(hadron)
library(optparse)

if (TRUE) {
    # set option list
option_list <- list(
    make_option(c("-b", "--beta"), type  =  "double", default  =  1.7,
    help  =  "beta-Parameter of simulation [default %default]"),
    make_option(c("--skip"), type  =  "integer", default  =  1000,
    help  =  "how many lines are skipped when reading in [default %default]"),

    make_option(c("-s", "--bootsamples"), type  =  "integer", default  =  500,
    help  =  "how many bootstrapsamples should be drawn [default %default]"),
    make_option(c("--maxrows"), type  =  "integer", default  =  -1,
    help  =  "Maximum of configurations that are read in, -1 = all [default %default]"),
    make_option(c("-r", "--Ns"), type  =  "integer", default  =  16,
    help  =  "L of lattice (spatial extent) [default %default]"),

    make_option(c("-t", "--Nt"), type  =  "integer", default  =  16,
    help  =  "T of lattice (temporal extent) [default %default]"),
    make_option(c("-n", "--nape"), type  =  "integer", default  =  0,
    help  =  "Number of APE-smears that were done [default %default]"),

    make_option(c("-a", "--alpha"), type  =  "double", default  =  1.0,
    help  =  "alpha used in APE-smearing [default %default]"),
    make_option(c("-S", "--sparam"), type  =  "double", default  =  1.5,
    help  =  "S parameter for uwerr [default %default]"),

    make_option(c("--betaone"), type  =  "double", default  =  0,
    help  =  "input beta at corresponding xi  =  1 [default %default]"),
    make_option(c("-x", "--xi"), type  =  "double", default  =  0,
    help  =  "xi used in lattice, only used if xidiff  =  TRUE, else
            xi is assumed to be L / T [default %default]"),

    make_option(c("--fraction"), type  =  "double", default  =  0.5,
    help  =  "Maximal extent of the Wilson-Loops as a fraction of the lattice size [default %default]"),
    make_option(c("--zerooffset"), type  =  "integer", default  =  0,
    help  =  "Offset for selecting the correct loops
            if W(x  =  0, y  =  0) was measured, set to 1 [default %default]"),

    make_option(c("-e", "--every"), type  =  "integer", default  =  1,
    help  =  "only reads every eth line [default %default]"),
    make_option(c("--bootl"), type  =  "integer", default  =  2,
    help  =  "block length in bootstrapping configurations [default %default]"),

    make_option(c("--nsave"), type  =  "integer", default  =  0,
    help  =  "steps between saved configs [default %default]"),
    make_option(c("--smearing"), action  =  "store_true", default  =  FALSE,
    help  =  "are the calculations done based on
            smeared lattices? (changes filenames) [default %default]"),

    make_option(c("--xidiff"), action  =  "store_true", default  =  FALSE,
    help  =  "Is xi different to L / T? [default %default]"),
    make_option(c("--errorall"), action  =  "store_true", default  =  FALSE,
    help  =  "Plot errors for different l for all y,t or only y,t = 1? [default %default]"),

    make_option(c("--respath"), type  =  "character", default  =  "",
    help  =  "path to where the resultfiles are stored [default %default]"),
    make_option(c("--plotpath"), type  =  "character", default  =  "",
    help  =  "path to where the plots are stored [default %default]"),
    make_option(c("--effmasstype"), type  =  "character", default  =  "log",
    help  =  "type of effective mass [default %default]"),

    make_option(c("--myfunctions"), type  =  "character",
    default  =  " / hiskp4 / gross / masterthesis / analyse / code / U1_analyse_potential / ",
#~     make_option(c("--myfunctions"), type  =  "character", default  =  "myfunctions.R",
    help  =  "path to where additional functions are stored,
            relative to folder where script is executed [default %default]")
)
parser <- OptionParser(usage  =  "%prog [options]", option_list  =  option_list)
args <- parse_args(parser, positional_arguments  =  0)
opt <- args$options
}

if (TRUE) {
# set some constants
source(paste(opt$myfunctions, "myfunctions.R", sep  =  ""))
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
options(digits = 7)
}

## set larger S
# set names for plot, for saving results
filenameforplots <- sprintf(
            "%stauintplotsuwerrNt%dNs%dbeta%fxi%f.pdf",
            opt$plotpath, Nt, Ns, beta, xi)
if (opt$smearing) {
    filenameforplots <- sprintf(
            "%stauintplotsuwerrNt%dNs%dbeta%fxi%fnape%dalpha%f.pdf",
            opt$plotpath, Nt, Ns, beta, xi, nape, alpha)
}
pdf(file  =  filenameforplots, title  =  "")
listtauint <- list()
tauintdf <- data.frame()
maxtimes <- c()
therm <- c()
## Determine autocorrelation for each point.
## Determine effect of blocksize on error for y=1, t=1, ans all y, t if errorall is set
## Determine unblocked for l=1, then with increasing blocksize.
## Determine error naively with standard deviation over sqrt(N), where N=number of blocks
## as well as with uwerr method
## To be able to compare to a method with autocorrelation
## First, spatial loops
## For y==1, save x=1, spatial-spatial plaquette, to be able to look at thermalisation
for (y in seq(1, Ns * opt$fraction, 1)) {
    if (y > 0) {
    filename <- sprintf(
        "result2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fcoarsedistance",
        Nt, Ns, beta, xi, nape, alpha)

    WL <- readloopfilecfrotated(file  =  filename, path  =  opt$respath, skip  =  skip,
            Nsmax  =  Ns * opt$fraction, yt  =  y, zerooffset  =  opt$zerooffset, every  =  opt$every, maxrows  =  opt$maxrows)
    if (y == 1) therm <- WL$cf[, 1]
    autotime <- c()
    dautotime <- c()
    errors <- c()
    derrors <- c()
    errorsnaive <- c()
    for (j in seq(1, opt$Ns * opt$fraction)){
            uwerrresults <- uwerr(data = WL$cf[, j], S = opt$sparam)
            autotime[j] <- uwerrresults$tauint
            dautotime[j] <- uwerrresults$dtauint
            maxtimes <- append(maxtimes, autotime[j])
            errors[j] <- uwerrresults$dvalue
            derrors[j] <- uwerrresults$ddvalue
            errorsnaive[j] <- sd(WL$cf[, j]) / sqrt(length(WL$cf[, j]))
            try(plot(uwerrresults, main = paste("spatial unblocked, x  = ", j, "y = ", y)))
       }
    tauintdf <- rbind(tauintdf, data.frame(yt = y, spatial = T, l = 1, times = as.data.frame(t(autotime)), dtimes = as.data.frame(t(dautotime)),
            errors = as.data.frame(t(errors)), derrors = as.data.frame(t(derrors)), errorsnaive  =  as.data.frame(t(errorsnaive))))

    if (opt$errorall || y == 1) {
    for (boot.l in c(2, 4, 8, 16, 32, 64, 128, 256)) {
        ## do not do bootstrapping yet, analyse the result of blocking first
        ## WL <- bootstrap.cf(WL, boot.R  =  bootsamples, boot.l  =  boot.l)
            for (j in seq(1, opt$Ns * opt$fraction)){
                myblocks <- blockts(data = WL$cf[, j], l = boot.l)
                autotime[j] <- NA
                dautotime[j] <- NA
                errors[j] <- NA
                derrors[j] <- NA
                errorsnaive[j] <- sd(myblocks) / sqrt(length(myblocks))
            }
        tauintdf <- rbind(tauintdf, data.frame(yt = y, spatial = T, l = boot.l, times = as.data.frame(t(autotime)), dtimes = as.data.frame(t(dautotime)),
                errors = as.data.frame(t(errors)), derrors = as.data.frame(t(derrors)), errorsnaive  =  as.data.frame(t(errorsnaive))))
    }}
}}
# same procedure as for y loop, here for temporal potential
# measured a_s V(r / a_t)
for (t in seq(1, Nt * opt$fraction, 1)) {
    if (t < Nt / 2 + 1) {
    filename <- sprintf(
            "result2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%ffinedistance",
            Nt, Ns, beta, xi, nape, alpha)

    # loops
    WL <- readloopfilecfrotated(file  =  filename, path  =  opt$respath, skip  =  skip,
            Nsmax  =  Ns * opt$fraction, yt  =  t, zerooffset  =  opt$zerooffset, every  =  opt$every, maxrows  =  opt$maxrows)

    autotime <- c()
    dautotime <- c()
    errors <- c()
    derrors <- c()
    errorsnaive <- c()
    for (j in seq(1, opt$Ns * opt$fraction)){
        uwerrresults <- uwerr(data = WL$cf[, j], S = opt$sparam)
        autotime[j] <- uwerrresults$tauint
        dautotime[j] <- uwerrresults$dtauint
        maxtimes <- append(maxtimes, autotime[j])
        errors[j] <- uwerrresults$dvalue
        derrors[j] <- uwerrresults$ddvalue
        errorsnaive[j] <- sd(WL$cf[, j]) / sqrt(length(WL$cf[, j]))
        try(plot(uwerrresults, main = paste("temporal unblocked, x  = ", j, "t = ", t)))
    }
    tauintdf <- rbind(tauintdf, data.frame(yt = t, spatial = F, l = 1, times = as.data.frame(t(autotime)), dtimes = as.data.frame(t(dautotime)),
            errors = as.data.frame(t(errors)), derrors = as.data.frame(t(derrors)), errorsnaive  =  as.data.frame(t(errorsnaive))))
    if (opt$errorall || t == 1) {
    for (boot.l in c(2, 4, 8, 16, 32, 64, 128, 256)) {
        ## do not do bootstrapping yet, analyse the result of blocking first
        ## WL <- bootstrap.cf(WL, boot.R  =  bootsamples, boot.l  =  boot.l)
            for (j in seq(1, opt$Ns * opt$fraction)){
                myblocks <- blockts(data = WL$cf[, j], l = boot.l)
                autotime[j] <- NA
                dautotime[j] <- NA
                errors[j] <- NA
                derrors[j] <- NA
                errorsnaive[j] <- sd(myblocks) / sqrt(length(myblocks))
            }
        tauintdf <- rbind(tauintdf, data.frame(yt = t, spatial = F, l = boot.l, times = as.data.frame(t(autotime)), dtimes = as.data.frame(t(dautotime)),
                errors = as.data.frame(t(errors)), derrors = as.data.frame(t(derrors)), errorsnaive  =  as.data.frame(t(errorsnaive))))
    }}
}}

# write out results
filenameuwerr <- sprintf(
            "%stauintNt%dNs%dbeta%fxi%fbsamples%d.csv",
            opt$plotpath, Nt, Ns, beta, xi, bootsamples)

if (opt$smearing) {
    filenameuwerr <- sprintf(
            "%stauintNt%dNs%dbeta%fxi%fnape%dalpha%fbsamples%d.csv",
            opt$plotpath, Nt, Ns, beta, xi, nape, alpha, bootsamples)
}
write.table(tauintdf, file  =  filenameuwerr, row.names = F, col.names = T)


## Plot errors and autocorrelation time
## into separate file, so one is not overwhelmed by the output of uwerr

filenameforplots <- sprintf(
            "%stauintplotsNt%dNs%dbeta%fxi%f.pdf",
            opt$plotpath, Nt, Ns, beta, xi)
if (opt$smearing) {
    filenameforplots <- sprintf(
            "%stauintplotsNt%dNs%dbeta%fxi%fnape%dalpha%f.pdf",
            opt$plotpath, Nt, Ns, beta, xi, nape, alpha)
}
pdf(file  =  filenameforplots, title  =  "")


for (y in seq(1, Ns * opt$fraction, 1)) {
    if (opt$errorall || y == 1) {
    ## error of 1/l
    for (j in seq(1, opt$Ns * opt$fraction)){
        columnname <- paste("errors.V", j, sep = "")
        dcolumnname <- paste("derrors.V", j, sep = "")
        columnname2 <- paste("errorsnaive.V", j, sep = "")
        mask <- tauintdf$yt == y & tauintdf$spatial == T
        try(plotwitherror(x = 1. / tauintdf$l[mask], y = tauintdf[mask, columnname], dy = tauintdf[mask, dcolumnname],
        xlab = "1 / l", ylab = "error", main = paste("x = ", j, "y = ", y, "Errors"),
        ylim = c(min(tauintdf[mask, columnname2]), max(tauintdf[mask, columnname2])),
        xlim = c(0, 1)))
        try(plotwitherror(x = 1. / tauintdf$l[mask], y = tauintdf[mask, columnname2], rep = T, col = 2, pch = 2))
        reldev <- max(tauintdf[mask, columnname2]) / tauintdf[mask & tauintdf$l == 1, columnname2]
        reldev <- sprintf("%.3f", (reldev - 1) * 100)

        ## line for result of uwerr, l=1
        mask <- tauintdf$yt == y & tauintdf$spatial == T & tauintdf$l == 1
        try(lines(x = seq(-1, 2, by = 1), y = rep(tauintdf[mask, columnname], 4), col = 1))
        try(lines(x = seq(-1, 2, by = 1), y = rep(tauintdf[mask, columnname] + rep(tauintdf[mask, dcolumnname], 4)), col = 1, lty = 2))
        try(lines(x = seq(-1, 2, by = 1), y = rep(tauintdf[mask, columnname] - rep(tauintdf[mask, dcolumnname], 4)), col = 1, lty = 2))

        try(legend(x = "topright", title = "max rel dev to l = 1:", legend = c(paste(reldev, "%")), pch = NA, col = NA))
        try(legend(x = "bottomright", legend = c("uwerr", "sd / sqrt(N)"), col = c(1, 2), pch = c(1, 2)))
        }}
        ## autocorrelation time
        try(plot(NA, xlim = c(1, Ns / 2), ylim = c(0, max(maxtimes)), xlab = "j", ylab = "tauint",
        main = paste("y = ", y, "Autocorrelation time")))
        try(lines(x = seq(0, opt$Nt, by = 1), y = rep(0.5, opt$Nt+1), col = 2))
        for (j in seq(1, opt$Ns * opt$fraction)){
            mask <- tauintdf$spatial == T & tauintdf$l == 1 & tauintdf$yt == y
            columnname <- paste("times.V", j, sep = "")
            dcolumnname <- paste("dtimes.V", j, sep = "")
            try(plotwitherror(x = j, y = tauintdf[mask, columnname], dy = tauintdf[mask, dcolumnname], rep = TRUE))
        }
}
for (t in seq(1, Nt * opt$fraction, 1)) {
    if (opt$errorall || t == 1) {
    ## error as a function of 1/l
    for (j in seq(1, opt$Ns * opt$fraction)){
        columnname <- paste("errors.V", j, sep = "")
        dcolumnname <- paste("derrors.V", j, sep = "")
        columnname2 <- paste("errorsnaive.V", j, sep = "")
        mask <- tauintdf$yt == t & tauintdf$spatial == F
        try(plotwitherror(x = 1. / tauintdf$l[mask], y = tauintdf[mask, columnname], dy = tauintdf[mask, dcolumnname],
        xlab = "1 / l", ylab = "error", main = paste("x = ", j, "t = ", t, "Errors"),
        ylim = c(min(tauintdf[mask, columnname2]), max(tauintdf[mask, columnname2])),
        xlim = c(0, 1)))

        try(plotwitherror(x = 1. / tauintdf$l[mask], y = tauintdf[mask, columnname2], rep = T, col = 2, pch = 2))
        reldev <- max(tauintdf[mask, columnname2]) / tauintdf[mask & tauintdf$l == 1, columnname2]
        reldev <- sprintf("%.3f", (reldev - 1) * 100)

        ## line for result of uwerr, l=1
        mask <- tauintdf$yt == t & tauintdf$spatial == F & tauintdf$l == 1
        try(lines(x = seq(-1, 2, by = 1), y = rep(tauintdf[mask, columnname], 4), col = 1))
        try(lines(x = seq(-1, 2, by = 1), y = rep(tauintdf[mask, columnname] + rep(tauintdf[mask, dcolumnname], 4)), col = 1, lty = 2))
        try(lines(x = seq(-1, 2, by = 1), y = rep(tauintdf[mask, columnname] - rep(tauintdf[mask, dcolumnname], 4)), col = 1, lty = 2))

        try(legend(x = "topright", title = "max rel dev to l = 1:", legend = c(paste(reldev, "%")), pch = NA, col = NA))
        try(legend(x = "bottomright", legend = c("uwerr", "sd / sqrt(N)"), col = c(1, 2), pch = c(1, 2)))
        }}
        ## autocorrelation time
        try(plot(NA, xlim = c(1, Ns / 2), ylim = c(0, max(maxtimes)), xlab = "x", ylab = "tauint",
        main = paste("t = ", t, "Autocorrelation time")))
        try(lines(x = seq(0, opt$Nt, by = 1), y = rep(0.5, opt$Nt+1), col = 2))
        for (j in seq(1, opt$Ns * opt$fraction)){
            mask <- tauintdf$spatial == F & tauintdf$l == 1 & tauintdf$yt == t
            columnname <- paste("times.V", j, sep = "")
            dcolumnname <- paste("dtimes.V", j, sep = "")
            try(plotwitherror(x = j, y = tauintdf[mask, columnname], dy = tauintdf[mask, dcolumnname], rep = TRUE))
        }
}

try(plot(therm, main=paste("thermalised values, skip=", opt$skip), xlab="MCMC-steps", ylab="P"))

print(warnings())
