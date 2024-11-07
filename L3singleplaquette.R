library(optparse)

if (TRUE) {
    # set option list
    option_list <- list(
        make_option(c("-b", "--beta"),
            type = "double", default = 1.7,
            help = "beta-Parameter of simulation [default %default]"
        ),
        make_option(c("--skip"),
            type = "integer", default = 1000,
            help = "how many lines are skipped when reading in [default %default]"
        ),
        make_option(c("-s", "--bootsamples"),
            type = "integer", default = 500,
            help = "how many bootstrapsamples should be drawn [default %default]"
        ),
        make_option(c("--maxrows"),
            type = "integer", default = -1,
            help = "Maximum of configurations that are read in, -1=all [default %default]"
        ),
        make_option(c("-r", "--Ns"),
            type = "integer", default = 3,
            help = "L of lattice (spatial extent) [default %default]"
        ),
        make_option(c("-t", "--Nt"),
            type = "integer", default = 16,
            help = "T of lattice (temporal extent) [default %default]"
        ),
        make_option(c("-n", "--nape"),
            type = "integer", default = 0,
            help = "Number of APE-smears that were done [default %default]"
        ),
        make_option(c("-a", "--alpha"),
            type = "double", default = 1.0,
            help = "alpha used in APE-smearing [default %default]"
        ),
        make_option(c("-x", "--xi"),
            type = "double", default = 0,
            help = "xi used in lattice, only used if xidiff = TRUE, else
            xi is assumed to be L/T [default %default]"
        ),
        make_option(c("-e", "--every"),
            type = "integer", default = 1,
            help = "only reads every eth line [default %default]"
        ),
        make_option(c("--uwerrs"),
            type = "double", default = 6,
            help = "parameter S used for the uwerr analysis [default %default]"
        ),
        make_option(c("--smearing"),
            action = "store_true", default = FALSE,
            help = "are the calculations done based on
            smeared lattices? (changes filenames) [default %default]"
        ),
        make_option(c("--xidiff"),
            action = "store_true", default = FALSE,
            help = "Is xi different to L/T? [default %default]"
        ),
        make_option(c("--respath"),
            type = "character", default = "",
            help = "path to where the resultfiles are stored [default %default]"
        ),
        make_option(c("--plotpath"),
            type = "character", default = "",
            help = "path to where the plots are stored [default %default]"
        ),
        make_option(c("--myfunctions"),
            type = "character",
            default = "/hiskp4/gross/masterthesis/analyse/code/U1_analyse_potential/",
            # ~     make_option(c("--myfunctions"), type = "character", default = "myfunctions.R",
            help = "path to where additional functions are stored,
            relative to folder where script is executed [default %default]"
        )
    )
    parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
    args <- parse_args(parser, positional_arguments = 0)
    opt <- args$options
}

# we have to use a custom installation of R to be able to change the parameter S in uwerr
require("hadron", lib.loc = "/hiskp4/gross/masterthesis/analyse/code/")
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

endinganalysis <- sprintf("Nt%dNs%dbeta%fxi%f", Nt, Ns, beta, xi)
if (opt$smearing) {
    endinganalysis <- sprintf(
        "Nt%dNs%dbeta%fxi%fnape%dalpha%f",
        Nt, Ns, beta, xi, nape, alpha
    )
}

pdf(sprintf("%s/plaqsingle%s.pdf", opt$plotpath, endinganalysis), title = "")


filename3 <- sprintf("%s/result2p1d.u1potential.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fnonplanar", opt$respath, Nt, Ns, beta, xi, opt$nape, opt$alpha)

## spatial-spatial-plaquette
data3 <- readloopfilecfonecolumn(file = filename3, path = "", skip = opt$skip, column = 6, memsafe = TRUE, every = opt$every, maxrows = opt$maxrows)
uwerr3 <- uwerrprimary(data3$cf[, 1], S = opt$uwerrs)
summary(uwerr3)
plot(x = seq(1, length(data3$cf[, 1])), y = data3$cf[, 1], main = paste("L = 3, T = ", Nt, "spatial-spatial"))
plot(uwerr3, main = paste("L = 3, T = ", Nt, "spatial-spatial"))
bootl <- ceiling(4 * uwerr3$tauint^2)
plaquettecf <- bootstrap.cf(data3,
    boot.R = bootsamples, boot.l = bootl
)
res <- list(uwerr = uwerr3, cf = plaquettecf, opt = opt, githash = githash)
saveRDS(res, sprintf("%s/plaqsingle%s.RData", opt$plotpath, endinganalysis))
