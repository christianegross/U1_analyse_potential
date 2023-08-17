library(optparse)

# TODO. implement possibility for differnt lower/upper for temporal/spatial

# makes tables in the right format for the boundaries of the effective masses.
# options are rotated for the sideways potential,
# subtracted for the normal potential and small for the non-integer potential.

if (TRUE) {
option_list <- list(
    make_option(c("-b", "--beta"), type = "double", default = 1.7,
    help = "beta-Parameter of simulation [default %default]"),
     make_option(c("-r", "--Ns"), type = "integer", default = 3,
    help = "Ns of lattice [default %default]"),

    make_option(c("-t", "--Nt"), type = "integer", default = 3,
    help = "Nt of lattice [default %default]"),
    make_option(c("-n", "--nape"), type = "integer", default = 0,
    help = "Number of APE-smears that were done [default %default]"),
     make_option(c("--lower"), type = "integer", default = 2,
    help = "lower boundary for determining meff [default %default]"),
     make_option(c("--upper"), type = "integer", default = 0,
    help = "upper boundary for determining meff,
            if 0, choose maximum extent [default %default]"),
     make_option(c("--lowertemporal"), type = "integer", default = 2,
    help = "lower boundary for determining meff [default %default]"),
     make_option(c("--uppertemporal"), type = "integer", default = 0,
    help = "upper boundary for determining meff,
            if 0, choose maximum extent [default %default]"),

    make_option(c("-a", "--alpha"), type = "double", default = 1.0,
    help = "alpha used in APE-smearing [default %default]"),
     make_option(c("-x", "--xi"), type = "double", default = 1,
    help = "xi used in lattice [default %default]"),

    make_option(c("--rotated"), action = "store_true", default = FALSE,
    help = "make table for analysis of sideways potential [default %default]"),
    make_option(c("--subtracted"), action = "store_true", default = FALSE,
    help = "make table for analysis of normal potential [default %default]"),
    make_option(c("--small"), action = "store_true", default = FALSE,
    help = "make table for analysis of small,
            nonplanar potential [default %default]"),

    make_option(c("--smearing"), action = "store_true", default = FALSE,
    help = "are the calculations done based on smeared lattices?
            (changes filenames) [default %default]"),
    make_option(c("--xidiff"), action = "store_true", default = FALSE,
    help = "Is xi different to Ns / Nt? [default %default]"),
   make_option(c("--path"), type = "character", default = ".",
    help = "path to where tables should be stored [default %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

Ns <- opt$Ns
Nt <- opt$Nt
xi <- Ns / (Nt * 1.0)
if (opt$xidiff) {
    xi <- opt$xi
}
nape <- opt$nape
alpha <- opt$alpha
upper <- opt$upper
if(opt$upper==0) upper <- opt$Ns / 2 - 2
uppertemporal <- opt$temporal
if(opt$uppertemporal==0) uppertemporal <- opt$Nt / 2 - 2
}

if (opt$rotated) {
    filename <- sprintf(
            "%s/tableanalysisrotatedNs%dNt%dbeta%fxi%f.csv",
            opt$path, opt$Ns, opt$Nt, opt$beta, xi)
    if (opt$smearing) {
        filename <- sprintf(
            "%s/tableanalysisrotatedNs%dNt%dbeta%fxi%fnap%dalpha%f.csv",
            opt$path, opt$Ns, opt$Nt, opt$beta, xi, nape, alpha)
    }
    print(filename)
    if (file.exists(filename)) {
        print("file already exists, doing nothing")
    } else {
        spacial <- data.frame(spacial =
                c(rep(TRUE, opt$Ns / 2), rep(FALSE, opt$Nt / 2)))
        yt <- data.frame(yt = c(seq(1, opt$Ns / 2), seq(1, opt$Nt / 2)))
        lowerdf <- data.frame(lower = rep(opt$lower, length(spacial)))
        upperdf <- data.frame(upper = rep(upper, length(spacial)))
        table <- cbind(spacial, yt, lowerdf, upperdf)
        write.table(table, filename, append = TRUE,
                row.names = FALSE, col.names = TRUE, sep=",")
    }
}

if (opt$subtracted) {
    filename <- sprintf(
                "%s/tableanalysissubtractedNs%dNt%dbeta%fxi%f.csv",
                opt$path, opt$Ns, opt$Nt, opt$beta, xi)
    if (opt$smearing) {
        filename <- sprintf(
                "%s/tableanalysissubtractedNs%dNt%dbeta%fxi%fnap%dalpha%f.csv",
                opt$path, opt$Ns, opt$Nt, opt$beta, xi, nape, alpha)
    }
    print(filename)
    if (file.exists(filename)) {
        print("file already exists, doing nothing")
    } else {
        spacial <- data.frame(spacial = rep(c(TRUE, FALSE), opt$Ns / 2))
        yt <- data.frame(yt = rep(NA, opt$Ns))
        yt$yt[seq(1, opt$Ns - 1, by = 2)] <- seq(1, opt$Ns / 2)
        yt$yt[seq(2, opt$Ns, by = 2)] <- seq(1, opt$Ns / 2)
        lowerdf <- data.frame(lower = rep(c(opt$lower, opt$lowertemporal), opt$Ns / 2))
        upperdf <- data.frame(upper = rep(c(upper, uppertemporal), opt$Ns / 2))
        table <- cbind(spacial, yt, lowerdf, upperdf)
        write.table(table, filename,
        append = TRUE, row.names = FALSE, col.names = TRUE, sep=",")
    }
}

if (opt$small) {
    filename <- sprintf(
                "%s/tableanalysissmallNs%dNt%dbeta%fxi%f.csv",
                opt$path, opt$Ns, opt$Nt, opt$beta, xi)
    if (opt$smearing) {
        filename <- sprintf(
                "%s/tableanalysissmallNs%dNt%dbeta%fxi%fnap%dalpha%f.csv",
                opt$path, opt$Ns, opt$Nt, opt$beta, xi, nape, alpha)
    }
    print(filename)
    if (file.exists(filename)) {
        print("file already exists, doing nothing")
    } else {
        x <- data.frame(x = c(rep(0, 4), rep(1, 4), rep(2, 4), rep(3, 4)))
        y <- data.frame(y = rep(c(0, 1, 2, 3), 4))
        lower <- data.frame(lower = rep(opt$lower, 16))
        upper <- data.frame(upper = rep(upper, 16))
        table <- cbind(x, y, lower, upper)

        write.table(table, filename,
            append = TRUE, row.names = FALSE, col.names = TRUE, sep=",")
    }
}
