## Determine effective masses for every single plateau choice
## probably not needed

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

    make_option(c("--useCov"), action = "store_true", default = FALSE,
    help = "Use covariamnt matrix for determining effective masses? [default %default]"),
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
}


if (TRUE) {
# set names for plot, tables, open, lists for saving results
filenameforplots <- sprintf(
            "%splotseffmass2p1dNt%dNs%dbeta%fxi%fbootl%duseCov%d.pdf",
            opt$plotpath, Nt, Ns, beta, xi, opt$bootl, opt$useCov)
if (opt$smearing) {
    filenameforplots <- sprintf(
            "%splotseff2p1dNt%dNs%dbeta%fxi%fnape%dalpha%fbootl%duseCov%d.pdf",
            opt$plotpath, Nt, Ns, beta, xi, opt$bootl, nape, alpha, opt$useCov)
}
pdf(file = filenameforplots, title = "")
potential <- data.frame(R = NA, space = NA, t1 = NA, t2 = NA, m = NA, dm = NA, p = NA, chi = NA)

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
    if (y > 0) {

    # loops
    print(y)
    filename <- sprintf(
        "result2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fcoarsedistance",
        Nt, Ns, beta, xi, nape, alpha)

    title <- sprintf("beta = %.2f coarse y = %d skipped %d",
                beta, y, skip * opt$nsave)
    WL <- calcplotWloopsideways(file = filename,
            path = opt$respath, skip = skip, Ns = Ns, yt = y,
            title = title, bootsamples = bootsamples,
            zerooffset = opt$zerooffset, nsave = opt$nsave,
            every = opt$every, l = opt$bootl, fraction = opt$fraction,
            maxrows = opt$maxrows, sim = "fixed")

    for (t1 in seq(0, WL$Time - 3, by = 1)) {
        for (t2 in seq(t1+1, WL$Time - 2, by = 1)) {
            WL.effmasslist <- deteffmass(WL = WL, yt = y,
            t1 = t1, t2 = t2, isspatial = 0, type = opt$effmasstype, useCov = opt$useCov)
            if (WL.effmasslist[[2]][[2]] != 0) {
                newline <- cbind(data.frame(R=y, space=TRUE, t1=t1, t2=t2), WL.effmasslist[[3]][, 2:3], WL.effmasslist[[3]][, 5:6])
                potential <- rbind(potential, newline)
            }
        }
    }
    try(plot(WL.effmasslist[[1]], xlab = "x / a_s", ylab = "Meff",
            main = title))
    try(plot(WL.effmasslist[[1]], xlab = "x / a_s", ylab = "Meff",
            main = title,
            ylim = WL.effmasslist[[2]]))
}
}
# same procedure as for y loop, here for temporal potential
# measured a_s V(r/a_t)
for (t in seq(1, Nt / 2, 1)) {
    if (t < Nt / 2 + 1) {

#~     message("\nt = ", t)
    print(t)
    filename <- sprintf(
            "result2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%ffinedistance",
            Nt, Ns, beta, xi, nape, alpha)

    # loops
    title <- sprintf("beta = %.2f fine t = %d skipped %d",
                beta, t, skip * opt$nsave)
    WL <- calcplotWloopsideways(file = filename,
            path = opt$respath, skip = skip, Ns = Ns, yt = t, title = title,
            bootsamples = bootsamples, zerooffset = opt$zerooffset,
            nsave = opt$nsave, every = opt$every, l = opt$bootl, fraction = opt$fraction,
            maxrows = opt$maxrows, sim = "fixed")

    for (t1 in seq(0, WL$Time - 3, by=1)) {
        for (t2 in seq(t1+1, WL$Time - 2, by=1)) {
            WL.effmasslist <- deteffmass(WL = WL, yt = t,
            t1 = t1, t2 = t2, isspatial = 0, type = opt$effmasstype, useCov = opt$useCov)
            if (WL.effmasslist[[2]][[2]] != 0) {
                newline <- cbind(data.frame(R=t, space=FALSE, t1=t1, t2=t2), WL.effmasslist[[3]][, 2:3], WL.effmasslist[[3]][, 5:6])
                potential <- rbind(potential, newline)
            }
        }
    }
    try(plot(WL.effmasslist[[1]], xlab = "x / a_s", ylab = "Meff",
            main = title))
    try(plot(WL.effmasslist[[1]], xlab = "x / a_s", ylab = "Meff",
            main = title,
            ylim = WL.effmasslist[[2]]))

}
}
#write out results
potential <- na.omit(potential)
filenamepotential <- sprintf(
            "%spotentialmeff2p1dmeffchosenNt%dNs%dbeta%fxi%fbsamples%dbootl%duseCov%d.csv",
            opt$plotpath, Nt, Ns, beta, xi, bootsamples, opt$bootl, opt$useCov)

if (opt$smearing) {
    filenamepotential <- sprintf(
            "%spotentialmeff2p1dmeffchosenNt%dNs%dbeta%fxi%fnape%dalpha%fbsamples%dbootl%duseCov%d.csv",
            opt$plotpath, Nt, Ns, beta, xi, nape, alpha, bootsamples, opt$bootl, opt$useCov)
}

write.table(potential, filenamepotential, row.names = FALSE)
}
