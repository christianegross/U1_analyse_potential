library(optparse)


if (TRUE) {
option_list <- list(
    make_option(c("-s", "--bootsamples"), type = "integer", default = 500,
    help = "how many bootstrapsamples should be drawn [default %default]"),
    make_option(c("-b", "--beta"), type = "double", default = 1.7,
    help = "beta at xi = 1 [default %default]"),
    make_option(c("--sparam"), type = "double", default = 6,
    help = "S parameter for the uwerr-function [default %default]"),
    make_option(c("--myfunctions"), type = "character",
        default = "/hiskp4/gross/masterthesis/analyse/code/U1_analyse_potential/",
#~     make_option(c("--myfunctions"), type = "character", default = "myfunctions.R",
    help = "path to where additional functions are stored [default %default]"),

    make_option(c("--respath"), type = "character", default = "plotstikz/",
    help = "path to where the resulting plots and data are stored [default %default]"),
    make_option(c("--datapath"), type = "character", default = "./",
    help = "path to where the data for the analyzed configs are stored [default %default]"),
    make_option(c("--type"), type = "character", default = "normal",
    help = "type of ensembles that shoule be analysed, one of normal, sideways [default %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
source(paste(opt$myfunctions, "myfunctions.R", sep = ""))
}

githash <- printgitcommit(opt$myfunctions)
bootsamples <-  500

beta <-  opt$beta
S <-  opt$sparam
readindata <- opt$readindata
if (opt$type ==  "normal") {
if (beta ==  1.65) {

betas <-  c(1.65, 1.615, 1.5614, 1.515, 1.48, 1.45, 1.42) # beta = 1.65
}

if (beta ==  1.7) {

betas <-  c(1.7, 1.6521, 1.6075, 1.5514, 1.525, 1.5, 1.49) # beta = 1.7
}
} else if (opt$type == "sideways") { #sideways potentials
if (beta ==  1.65) {

betas <-  c(1.65, 1.615, 1.565, 1.49, 1.4638, 1.4406, 1.42) # beta = 1.65
}

if (beta ==  1.7) {
betas <-  c(1.7, 1.6521, 1.6, 1.5381, 1.525, 1.4814, 1.49) # beta = 1.7
}
} else {
    stop(paste("type", opt$type, "is not supported"))
}

Nt <-  c(16, 20, 24, 32, 40, 48, 64)
therm16 <-  c(500, 500, 500, 500, 500, 1000, 2000)
therm3 <-  c(1333, 1333, 1333, 1333, 2666, 5333, 10666)

# xiconststr <- ""
# if (opt$xiconst) xiconststr <- "xiconst"
# print(opt$singlemulti)
# if (opt$singlemulti = =  "single") endname <- sprintf("%sbeta%fomit%d%slowlim%d", opt$type, opt$beta, opt$omit, xiconststr, opt$lowlimfitpot)
# if (opt$singlemulti = =  "multi") endname <- sprintf("multi%sbeta%fomit%d%slowlim%d", opt$type, opt$beta, opt$omit, xiconststr, opt$lowlimfitpot)
# if (opt$aic) endname <- sprintf("%saic", endname)
# if (opt$scaletauint) endname <- sprintf("%sscaletauintetp%d", endname, opt$errortotpot)


endname <- sprintf("%sbeta%f", opt$type, opt$beta)

## Read in result for L = 3


list3 <- list()
P3 <- array(NA, dim = c(opt$bootsamples, length(Nt)))


pdf(sprintf("L3%s.pdf", endname), title = "")

for (i in seq(1, length(betas))) {
    print(i)
    filengthame3 <-  sprintf("%s/result2p1d.u1potential.Nt%d.Ns3.b%f.xi%f.nape0.alpha1.000000nonplanar", opt$datapath, Nt[i], betas[i], 16 / Nt[i])
    data3 <-  readloopfilecfonecolumn(file = filengthame3, path = "", skip = therm3[i], column = 6, memsafe = TRUE)
    uwerr3 <-  uwerrprimary(data3$cf[, 1], S = S)
    summary(uwerr3)
    plot(x = seq(1, length(data3$cf[, 1])), y = data3$cf[, 1], main = paste("L = 3, T = ", Nt[i]))
    plot(uwerr3, main = paste("L = 3, T = ", Nt[i]))
    bootl <- ceiling(4 * uwerr3$tauint^2)
    plaquettecf <- bootstrap.cf(data3,
                boot.R = opt$bootsamples, boot.l = bootl)
    list3[[i]] <-  list(uwerr = uwerr3, cf = plaquettecf)
    P3[, i] <- plaquettecf$cf.tsboot$t[, 1]

}

names(list3) <- c("xi1", "xi2", "xi3", "xi4", "xi5", "xi6", "xi7")
list3[["plaquette"]] <- P3

saveRDS(object = list3, file = sprintf("plaquetteL3%s.RData", endname))

## check for normal distribution

for (i in seq(1, length(betas))) {
    shapiro3 <- shapiro.test(x = P3[, i])
    print(shapiro3)
    qqnorm(y = P3[, i], main = paste("L = 3, xi", i, shapiro3$p.value))
    qqline(y = P3[, i])
}



## completely rewrite programm:
## bootstrapsamples are determined such that the blocklength should cancel out the effects of autocorrelation,
## therefore for L = 16 just read in the results from chosepredict for P and xi.
## For L = 3, read in plaquette data, do the same as was done for L = 16:
## determine autocorrelation time with uwerr, determine P and bootstrapsamples with cf.
## save: beta, xi, P16, P3, ratio, with errors as well as bootstrapsamples
## Test if plaquette is normal distributed, save coefficients?
## beta should still be given as a fixed list, at some point this could be changed with an input file, not worth the effort now.
