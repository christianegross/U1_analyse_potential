library("hadron")
library(optparse)

# reads in data for the slope, reads in data for the


if (TRUE) {
option_list <- list(
    make_option(c("-s", "--bootsamples"), type = "integer", default = 500,
    help = "how many bootstrapsamples should be drawn [default %default]"),
    make_option(c("-b", "--beta"), type = "double", default = 1.7,
    help = "beta of ensemble [default %default]"),
    make_option(c("-x", "--xi"), type = "double", default = 1.0,
    help = "input xi of ensemble [default %default]"),
    make_option(c("-L", "--length"), type = "integer", default = 16,
    help = "spatial extent of ensemble [default %default]"),
    make_option(c("--lengthlarge"), type = "integer", default = 16,
    help = "spatial extent of ensemble from which the anisotropies were determined [default %default]"),

    make_option(c("-T", "--timeextent"), type = "integer", default = 16,
    help = "time extent of ensemble [default %default]"),
    make_option(c("--myfunctions"), type = "character",
        default = "/hiskp4/gross/masterthesis/su2/build/debug/analysisscripts/myfunctions.R",
#~     make_option(c("--myfunctions"), type = "character", default = "myfunctions.R",
    help = "path to where additional functions are stored [default %default]"),
    make_option(c("--respath"), type = "character", default = "",
    help = "path to where the resulting plots and data are stored [default %default]"),
    make_option(c("--datapath"), type = "character", default = "plotstikz",
    help = "path to where data for renormalized anisotropy [default %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
githash <- printgitcommit(opt$myfunctions)
# print(opt)
}


print(sprintf("resultspotNs%dNt%dbeta%fxi%fbs%d.RData", opt$length, opt$timeextent, opt$beta, opt$xi, opt$bootsamples))
result <- readRDS(sprintf("%sresultspotNs%dNt%dbeta%fxi%fbs%d.RData", opt$respath, opt$length, opt$timeextent, opt$beta, opt$xi, opt$bootsamples))

if(opt$xi!=1){

xifits <- readRDS(sprintf("%slistresultsrenormalization.RData", opt$datapath))

xifits <- xifits[["fitsxi"]]

# print(xifits[[2]])
xis <- c(1, 0.8, 2/3, 0.5, 0.4, 1/3, 0.25)
distances <- abs(xis-opt$xi)<0.01 ##roundabout way to take 0.33, 0.333333333333333333 equally into account

index <- match(TRUE, distances)
# print(index)

xiscale <- predict(xifits[[index]], opt$beta)
# print(xiscale$val)


newlist <- list(slopescaled=result[["slope"]] / xiscale$val,
                dslopescaled=sqrt((result[["dslope"]] / xiscale$val)^2+(result[["slope"]] * xiscale$err / xiscale$val^2)^2),
                bsslopescaled=result[["bsslope"]] / xiscale$boot,
                xiren=xiscale$val, dxiren=xiscale$err, bsxiren=xiscale$boot)

# print(newlist)

# resultscaled <- c(result, newlist)
} else{
    respotential <- readRDS(sprintf("%sresultssubtractedNs%dNt%dbeta%fxi%fbs%d.RData", opt$datapath, opt$lengthlarge, opt$timeextent, opt$beta, opt$xi, opt$bootsamples))
    newlist <- list(slopescaled=result[["slope"]] / respotential[["xicalc"]],
                    dslopescaled=sqrt((result[["dslope"]] / respotential[["xicalc"]])^2+(result[["slope"]] * respotential[["dxicalc"]] / respotential[["xicalc"]]^2)^2),
                    bsslopescaled=result[["bsslope"]],
                    xiren=respotential[["xicalc"]], dxiren=respotential[["dxicalc"]], bsxiren=respotential[["bsxicalc"]])
}
resultsmall <- c(newlist, list(p=as.double(result[["p"]]), dp=result[["dp"]], bsp=result[["bsp"]]))

saveRDS(resultsmall, sprintf("%sresultsmallscaledNs%dNt%dbeta%fxi%fbs%d.RData", opt$respath, opt$length, opt$timeextent, opt$beta, opt$xi, opt$bootsamples))

# print(xifits)