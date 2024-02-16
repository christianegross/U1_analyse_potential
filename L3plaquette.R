library(optparse)


if (TRUE) {
option_list <- list(
    make_option(c("-s", "--bootsamples"), type = "integer", default = 500,
    help = "how many bootstrapsamples should be drawn [default %default]"),
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
    make_option(c("-i", "--inputfile"), type = "character", default = "./",
    help = "inputfile for further parameters (list of betas, types, ...) [default %default]"),
    make_option(c("-l", "--line"), type = "integer", default = "1",
    help = "line of inputfile to read [default %default]")
)

## the inputfile should be formatted the same as for chosepredict
## datapath is given separately, the datapath from the inputfile is ignored

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
source(paste(opt$myfunctions, "myfunctions.R", sep = ""))
}

input <- read.table(opt$input, header = TRUE, sep=",")
input <- input[opt$line, ]
print(input)

githash <- printgitcommit(opt$myfunctions)
bootsamples <-  500


xiconststr <- ""
if (input$xiconst) xiconststr <- "xiconst"
print(input$singlemulti)
if (input$singlemulti == "single") endname <- sprintf("%sbeta%fomit%d%sllxi%dllr0%d", input$type, input$beta, input$omit, xiconststr, input$lowlimxi, input$lowlimpot)
if (input$singlemulti == "multi") endname <- sprintf("multi%sbeta%fomit%d%slowlim%d", input$type, input$beta, input$omit, xiconststr, input$lowlimxi)
if (input$crzero != -1.65) endname <- sprintf("%sc%.2f", endname, input$crzero)
if (input$aic) endname <- sprintf("%saic", endname)
if (input$scaletauint) endname <- sprintf("%sscaletauintetp%d", endname, input$errortotpot)


S <-  opt$sparam
Nt <-  c(16, 20, 24, 32, 40, 48, 64, 80)
therm3 <-  c(1333, 1333, 1333, 1333, 2666, 5333, 10666, 15000)

# xiconststr <- ""
# if (input$xiconst) xiconststr <- "xiconst"
# print(input$singlemulti)
# if (input$singlemulti = =  "single") endname <- sprintf("%sbeta%fomit%d%slowlim%d", input$type, input$beta, input$omit, xiconststr, input$lowlimfitpot)
# if (input$singlemulti = =  "multi") endname <- sprintf("multi%sbeta%fomit%d%slowlim%d", input$type, input$beta, input$omit, xiconststr, input$lowlimfitpot)
# if (input$aic) endname <- sprintf("%saic", endname)
# if (input$scaletauint) endname <- sprintf("%sscaletauintetp%d", endname, input$errortotpot)


## Read in result for L = 3


list3 <- list()
P3 <- array(NA, dim = c(input$bootsamples, length(Nt)))
Punbiased <- c()


pdf(sprintf("L3%s.pdf", endname), title = "")

for (i in seq(1, 8)) {
    print(i)
    betaname <- paste0("xi", i, "beta")
    if (i==1) betaname <- "beta"
    beta <- input[betaname]
    filengthame3 <-  sprintf("%s/result2p1d.u1potential.Nt%d.Ns3.b%f.xi%f.nape0.alpha1.000000nonplanar", opt$datapath, Nt[i], beta, 16 / Nt[i])
    data3 <-  readloopfilecfonecolumn(file = filengthame3, path = "", skip = therm3[i], column = 6, memsafe = TRUE)
    uwerr3 <-  uwerrprimary(data3$cf[, 1], S = S)
    summary(uwerr3)
    plot(x = seq(1, length(data3$cf[, 1])), y = data3$cf[, 1], main = paste("L = 3, T = ", Nt[i]))
    plot(uwerr3, main = paste("L = 3, T = ", Nt[i]))
    bootl <- ceiling(4 * uwerr3$tauint^2)
    plaquettecf <- bootstrap.cf(data3,
                boot.R = input$bootsamples, boot.l = bootl)
    list3[[i]] <-  list(uwerr = uwerr3, cf = plaquettecf, beta = beta)
    P3[, i] <- plaquettecf$cf.tsboot$t[, 1]
    Punbiased[i] <- uwerr3$value

}

names(list3) <- c("xi1", "xi2", "xi3", "xi4", "xi5", "xi6", "xi7", "xi8")
list3[["plaquette"]] <- P3
list3$Punbiased <- Punbiased
list3$githash <- githash

saveRDS(object = list3, file = sprintf("plaquetteL3%s.RData", endname))

## check for normal distribution

for (i in seq(1, 8)) {
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
