library("hadron")
library(optparse)

if (TRUE) {
    # set option list
option_list <- list(
    make_option(c("-b", "--beta"), type = "double", default = 1.7,
    help = "beta-Parameter of simulation [default %default]"),
    make_option(c("--skip"), type = "integer", default = 1000,
    help = "how many configurations are skipped when reading in [default %default]"),

    make_option(c("-r", "--Ns"), type = "integer", default = 16,
    help = "L of lattice (spatial extent) [default %default]"),

    make_option(c("-t", "--Nt"), type = "integer", default = 16,
    help = "T of lattice (temporal extent) [default %default]"),
    make_option(c("--betaone"), type = "double", default = 0,
    help = "input beta at corresponding xi = 1 [default %default]"),
    make_option(c("-x", "--xi"), type = "double", default = 0,
    help = "xi used in lattice, only used if xidiff = TRUE, else
            xi is assumed to be L/T [default %default]"),

    make_option(c("--nsave"), type = "integer", default = 0,
    help = "steps between saved configs [default %default]"),
    make_option(c("--xidiff"), action = "store_true", default = FALSE,
    help = "Is xi different to L/T? [default %default]"),
    make_option(c("--datapath"), type = "character", default = "./",
    help = "path to where the datafiles are stored [default %default]"),
    make_option(c("--plotpath"), type = "character", default = "",
    help = "path to where the plots are stored [default %default]"),
    make_option(c("--basename"), type = "character", default = "gradient_flow",
    help = "basename of resultfiles [default %default]"),

    make_option(c("--myfunctions"), type = "character",
        default = "/hiskp4/gross/masterthesis/analyse/code/U1_analyse_potential",
#~     make_option(c("--myfunctions"), type = "character", default = "myfunctions.R",
    help = "path to where additional functions are stored,
            relative to folder where script is executed [default %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
}


t0 <- 1.628e-3 # as determied in https://arxiv.org/pdf/2212.09627.pdf
dt0 <- 9.1e-5

tsqE <- c()
xi <- c()

filelist <- getorderedfilelist(path=opt$datapath, basename=opt$basename, last.digits=6, ending="")
# print(head(filelist))
for (index in seq(opt$skip, length(filelist))) {
    data <- read.table(file=filelist[index], header=F, skip=1,
    colClasses=c("numeric", "numeric", "NULL", "NULL", "numeric", rep("NULL", 4)),
    col.names=c("t", "xi", NA, NA, "E", rep(NA, 4)))
    data$tsqE <- data$t^2 * data$E
    print(head(data))
}
