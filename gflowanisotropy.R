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

    make_option(c("--xidiff"), action = "store_true", default = FALSE,
    help = "Is xi different to L/T? [default %default]"),
    make_option(c("--datapath"), type = "character", default = "./",
    help = "path to where the datafiles are stored [default %default]"),
    make_option(c("--plotpath"), type = "character", default = "",
    help = "path to where the plots are stored [default %default]"),
    make_option(c("--basename"), type = "character", default = "gradient_flow",
    help = "basename of resultfiles [default %default]"),

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

xiin <- opt$Ns / opt$Nt
if (opt$xidiff) {
    xiin <- opt$xi
}

source(paste(opt$myfunctions, "myfunctions.R", sep = ""))
githash <- printgitcommit(opt$myfunctions)

## as determied in https://arxiv.org/pdf/2212.09627.pdf
## take c0 to reproduce the result from the potential measurements
c0 <- 1.628e-3
dc0 <- 9.1e-5


## get all filenames
filelist <- getorderedfilelist(path = opt$datapath, basename = opt$basename, last.digits = 6, ending = "")

## prepare emtpy containers
xilist <- c()
times <- c()

data <- read.table(file = filelist[1], header = F, skip = 1,
    colClasses = c("numeric", "numeric", "NULL", "NULL", "numeric", rep("NULL", 4)),
    col.names = c("t", "xi", NA, NA, "E", rep(NA, 4)))

timesteps <- length(data$t)
timelist <- data$t

resultxi <- array(data=rep(NA, timesteps * length(filelist) - opt$skip),
dim=c(timesteps, length(filelist) - opt$skip))

resulttsqE <- resultxi

for (index in seq(opt$skip, length(filelist))) {
    ## for each file, read in necessary columns: t, xi, E
    data <- read.table(file = filelist[index], header = F, skip = 1,
    colClasses = c("numeric", "numeric", "NULL", "NULL", "numeric", rep("NULL", 4)),
    col.names = c("t", "xi", NA, NA, "E", rep(NA, 4)))
    ## determine t^2E
    data$tsqE <- data$t^2 * data$E
    ## determine indices for which t^2E in c0 +/- dc0
    match <- which(abs(data$tsqE - c0) < dc0)
    ## select anisotropies corresponding to these indices
    xilist <- append(xilist, data$xi[match])
    times <- append(times, data$t[match])
    ## save all values for plot
    resultxi[, index-opt$skip] <- data$xi
    resulttsqE[, index-opt$skip] <- data$t^2 * data$E
}
# print(xi)

## get mean and error
xi <- mean(xilist)
dxi <- sd(xilist)


## save result
result <- data.frame(beta = opt$beta, L = opt$Ns, T = opt$Nt, xiin = xiin,
xi = xi, dxi = dxi, time = mean(times), dtime = sd(times),
c0 = c0, dc0 = dc0, githash = githash,
noc = length(filelist) - opt$skip, nom = length(xilist))
print(result)

filename <- sprintf("%ssummaryanisotropygflowb%.3f.csv",
                    opt$plotpath, opt$betaone)
columnnames <- FALSE
if (!file.exists(filename)) {
    columnnames <- TRUE
}
write.table(result, filename,
        append = TRUE, row.names = FALSE, col.names = columnnames)

## prepare list of xi(t), E(t), plot

xitime <- apply(resultxi, MARGIN=1, FUN=mean)
dxitime <- apply(resultxi, MARGIN=1, FUN=sd)
Etime <- apply(resulttsqE, MARGIN=1, FUN=mean)
dEtime <- apply(resulttsqE, MARGIN=1, FUN=sd)

timeresult <- data.frame(t=timelist, xi=xitime, dxi = dxitime, tsqE = Etime, dtsqE = dEtime)

filename <- sprintf("%s/gflowb%fx%fL%dT%d", opt$plotpath, opt$beta, xiin, opt$Ns, opt$Nt)
write.table(x=timeresult, file=paste(filename, ".csv", sep=""), row.names=F, col.names=T)

pdf(paste(filename, ".pdf", sep=""), title="")

plotwitherror(x=timeresult$t, y=timeresult$tsqE, dy=timeresult$dtsqE,
main="Energy", xlab="t", ylab="t^2E")

plotwitherror(x=timeresult$t, y=timeresult$xi, dy=timeresult$dxi,
main="Anisotropy", xlab="t", ylab="xi")
