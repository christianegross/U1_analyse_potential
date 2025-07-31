library("hadron")

source("~/Documents/masterthesis/analyse_potential/myfunctions.R")
source("~/Documents/masterthesis/analyse_potential/matchwithellipse.R")

hamiltoniandata <- read.table("/home/gross/Documents/masterthesis/more_measurements/hamiltonian/interpolate_ham.csv", header = T)
modes <- c("xi0.20v2", "xi0.18v2")
types <- c("plaq", "beta", "plaqinter", "betainter", "plaqsmall", "ratio")

averagedrestable <- data.frame(
    beta = c(), etp = c(),
    contlimAIC = c(), dcontlimAIC = c(), dcontlimspreadAIC = c(),
    contlimerror = c(), dcontlimerror = c(), dcontlimspreaderror = c(),
    contlimunweighted = c(), dcontlimunweighted = c(), dcontlimspreadunweighted = c(),
    naratio = c(), nainbs = c(),
    mode = c(), type = c(), xiinter = c()
)

for (mode in modes) {
    for (type in types) {
        for (xiinter in c(T, F)) {
            tmp <- try(read.table(sprintf("%s/contlimittype%smode%s%saveraged.csv", "plotstikz", type, mode, ifelse(xiinter, "xiinter", ""))))
            if (!inherits(tmp, "try-error")) {
                tmp$mode <- mode
                tmp$type <- type
                tmp$xiinter <- xiinter
                averagedrestable <- rbind(averagedrestable, tmp)
            } else {
                message("problems with mode", mode, "type", type, "xiinter", xiinter)
            }
        }
    }
}
averagedrestable$interpolated_xi <- averagedrestable$xiinter | averagedrestable$type == "plaqinter" | averagedrestable$type == "betainter"
averagedrestable$interpolated_beta <- averagedrestable$type == "plaqinter" | averagedrestable$type == "betainter"
# averagedrestable[averagedrestable$mode=="xi0.18" & averagedrestable$type=="beta", ]
# averagedrestable[averagedrestable$mode=="xi0.18wo0.19" & averagedrestable$type=="beta", ]
averagedrestable
# stop()

pdf(sprintf("rescombined.pdf"), title = "")
modes <- c("all", "xi0.20", "xi0.19", "xi0.18", "xi0.18wo0.19", "xi0.25")
modes <- c("xi0.20v2", "xi0.18v2")
# modes <- c("xi0.19", "xi0.18", "xi0.18wo0.19")
cols <- c("black", "red", "blue", "darkgreen", "firebrick", "green")
pchs <- c(21, 22, 23, 24, 21, 22)

errtypes <- list(
    AIC = list(mean = "contlimAIC", sd = "dcontlimAIC", spread = "dcontlimspreadAIC", bs = "bsAIC", bsspread = "bsAICspread"),
    error = list(mean = "contlimerror", sd = "dcontlimerror", spread = "dcontlimspreaderror", bs = "bserror", bsspread = "bserrorspread"),
    unweighted = list(mean = "contlimunweighted", sd = "dcontlimunweighted", spread = "dcontlimspreadunweighted", bs = "bsunweighted", bsspread = "bsunweightedspread")
)
for (xiinter in c(T, F)) {
    for (err in c("AIC", "error", "unweighted")) {
        plot(
            NA,
            xlab = "beta", ylab = "P",
            xlim = range(
                averagedrestable[, errtypes[[err]]$mean][averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter] + averagedrestable[, errtypes[[err]]$spread][averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter],
                averagedrestable[, errtypes[[err]]$mean][averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter] - averagedrestable[, errtypes[[err]]$spread][averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter]
            ),
            ylim = range(
                averagedrestable[, errtypes[[err]]$mean][averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter] - averagedrestable[, errtypes[[err]]$spread][averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter],
                averagedrestable[, errtypes[[err]]$mean][averagedrestable$type == "plaqsmall" & averagedrestable$xiinter == xiinter] + averagedrestable[, errtypes[[err]]$spread][averagedrestable$type == "plaqsmall" & averagedrestable$xiinter == xiinter]
            ),
            main = paste("errtpye", err, "xiinter", xiinter)
        )
        for (mode_index in seq_along(modes)) {
            mode <- modes[mode_index]
            # print(errtypes[[err]]$mean)
            # print(names(averagedrestable))
            plaqmask <- averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter & averagedrestable$mode == mode
            plaqsmallmask <- averagedrestable$type == "plaqsmall" & averagedrestable$xiinter == xiinter & averagedrestable$mode == mode
            betamask <- averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter & averagedrestable$mode == mode
            plotwitherror(
                y = averagedrestable[, errtypes[[err]]$mean][plaqmask & averagedrestable$etp == 1],
                x = averagedrestable[, errtypes[[err]]$mean][betamask & averagedrestable$etp == 1],
                dy = array(
                    c(
                        averagedrestable[, errtypes[[err]]$sd][plaqmask & averagedrestable$etp == 0],
                        sqrt(averagedrestable[, errtypes[[err]]$sd][plaqmask & averagedrestable$etp == 1]^2 - averagedrestable[, errtypes[[err]]$sd][plaqmask & averagedrestable$etp == 0]^2),
                        sqrt(averagedrestable[, errtypes[[err]]$spread][plaqmask & averagedrestable$etp == 1]^2 - averagedrestable[, errtypes[[err]]$sd][plaqmask & averagedrestable$etp == 1]^2)
                    ),
                    dim = c(2, 3)
                ),
                dx = array(
                    c(
                        averagedrestable[, errtypes[[err]]$sd][betamask & averagedrestable$etp == 0],
                        sqrt(averagedrestable[, errtypes[[err]]$sd][betamask & averagedrestable$etp == 1]^2 - averagedrestable[, errtypes[[err]]$sd][betamask & averagedrestable$etp == 0]^2),
                        sqrt(averagedrestable[, errtypes[[err]]$spread][betamask & averagedrestable$etp == 1]^2 - averagedrestable[, errtypes[[err]]$sd][betamask & averagedrestable$etp == 1]^2)
                    ),
                    dim = c(2, 3)
                ),
                col = cols[mode_index], pch = pchs[mode_index], bg = cols[mode_index],
                rep = T, errsum.method = "quadrature"
            )
            plotwitherror(
                y = averagedrestable[, errtypes[[err]]$mean][plaqsmallmask & averagedrestable$etp == 1],
                x = averagedrestable[, errtypes[[err]]$mean][betamask & averagedrestable$etp == 1],
                dy = array(
                    c(
                        averagedrestable[, errtypes[[err]]$sd][plaqsmallmask & averagedrestable$etp == 0],
                        sqrt(averagedrestable[, errtypes[[err]]$sd][plaqsmallmask & averagedrestable$etp == 1]^2 - averagedrestable[, errtypes[[err]]$sd][plaqsmallmask & averagedrestable$etp == 0]^2),
                        sqrt(averagedrestable[, errtypes[[err]]$spread][plaqsmallmask & averagedrestable$etp == 1]^2 - averagedrestable[, errtypes[[err]]$sd][plaqsmallmask & averagedrestable$etp == 1]^2)
                    ),
                    dim = c(2, 3)
                ),
                dx = array(
                    c(
                        averagedrestable[, errtypes[[err]]$sd][betamask & averagedrestable$etp == 0],
                        sqrt(averagedrestable[, errtypes[[err]]$sd][betamask & averagedrestable$etp == 1]^2 - averagedrestable[, errtypes[[err]]$sd][betamask & averagedrestable$etp == 0]^2),
                        sqrt(averagedrestable[, errtypes[[err]]$spread][betamask & averagedrestable$etp == 1]^2 - averagedrestable[, errtypes[[err]]$sd][betamask & averagedrestable$etp == 1]^2)
                    ),
                    dim = c(2, 3)
                ),
                col = cols[mode_index], pch = pchs[mode_index], bg = cols[mode_index],
                rep = T, errsum.method = "quadrature"
            )
        }
        points(x = hamiltoniandata$lowerx, y = hamiltoniandata$lower, type = "o")
        legend("topleft", legend = c(modes, "H"), col = c(cols, "black"), pch = c(pchs, 1))
    }
}


warnings()


# resbs <- list(
#     beta = array(NA, dim = c(500, length(modes) * 16 * 3)), plaq3 = array(NA, dim = c(500, length(modes) * 16 * 3)),
#     plaq16 = array(NA, dim = c(500, length(modes) * 16 * 3)), plaq3ratio = array(NA, dim = c(500, length(modes) * 16 * 3))
# )
resbs <- list(
    beta = array(NA, dim = c(500, length(modes) * 16 * 1)), plaq3 = array(NA, dim = c(500, length(modes) * 16 * 1)),
    plaq16 = array(NA, dim = c(500, length(modes) * 16 * 1)), plaq3ratio = array(NA, dim = c(500, length(modes) * 16 * 1))
)
data <- data.frame(
    betacontlim = c(), dbetacontlim = c(),
    plaq3contlim = c(), dplaq3contlim = c(),
    plaq16contlim = c(), dplaq16contlim = c(),
    plaq3ratiocontlim = c(), dplaq3ratiocontlim = c(),
    betaiso = c(), mode = c(), xiinter = c(), spread = c(), errtype = c()
)

# error ratio result has to be collected from the bootstrap samples first
dplaq3ratio <- c()
index <- 1
bsnames <- c()
# for (err in c("AIC", "error", "unweighted")) {
for (err in c("unweighted")) {
    for (mode in modes) {
        for (xiinter in c(T, F)) {
        # for (xiinter in c(F)) {
            for (etp in c(1, 0)) {
                data <- rbind(data, data.frame(
                    betacontlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dbetacontlim = averagedrestable[, errtypes[[err]]$sd][averagedrestable$mode == mode & averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    plaq3contlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "plaqsmall" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dplaq3contlim = averagedrestable[, errtypes[[err]]$sd][averagedrestable$mode == mode & averagedrestable$type == "plaqsmall" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    plaq16contlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dplaq16contlim = averagedrestable[, errtypes[[err]]$sd][averagedrestable$mode == mode & averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    plaq3ratiocontlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp] * averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "ratio" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dplaq3ratiocontlim = NA,
                    betaiso = averagedrestable$beta[averagedrestable$mode == mode & averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    mode = mode, xiinter = xiinter, spread = F, errtype = err, etp = etp
                ))
                data <- rbind(data, data.frame(
                    betacontlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dbetacontlim = averagedrestable[, errtypes[[err]]$spread][averagedrestable$mode == mode & averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    plaq3contlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "plaqsmall" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dplaq3contlim = averagedrestable[, errtypes[[err]]$spread][averagedrestable$mode == mode & averagedrestable$type == "plaqsmall" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    plaq16contlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dplaq16contlim = averagedrestable[, errtypes[[err]]$spread][averagedrestable$mode == mode & averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    plaq3ratiocontlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp] * averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "ratio" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dplaq3ratiocontlim = NA,
                    betaiso = averagedrestable$beta[averagedrestable$mode == mode & averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    mode = mode, xiinter = xiinter, spread = T, errtype = err, etp = etp
                ))



                tmp <- try(readRDS(sprintf("%s/contlimittype%smode%s%saveraged.RData", "plotstikz", "beta", mode, ifelse(xiinter, "xiinter", ""))))
                resbs$beta[, 4 * index - 3] <- tmp[[paste0("etp", etp, "b1.65")]][[errtypes[[err]]$bs]]
                resbs$beta[, 4 * index - 2] <- tmp[[paste0("etp", etp, "b1.7")]][[errtypes[[err]]$bs]]
                resbs$beta[, 4 * index - 1] <- tmp[[paste0("etp", etp, "b1.65")]][[errtypes[[err]]$bsspread]]
                resbs$beta[, 4 * index + 0] <- tmp[[paste0("etp", etp, "b1.7")]][[errtypes[[err]]$bsspread]]
                tmp <- try(readRDS(sprintf("%s/contlimittype%smode%s%saveraged.RData", "plotstikz", "plaqsmall", mode, ifelse(xiinter, "xiinter", ""))))
                resbs$plaq3[, 4 * index - 3] <- tmp[[paste0("etp", etp, "b1.65")]][[errtypes[[err]]$bs]]
                resbs$plaq3[, 4 * index - 2] <- tmp[[paste0("etp", etp, "b1.7")]][[errtypes[[err]]$bs]]
                resbs$plaq3[, 4 * index - 1] <- tmp[[paste0("etp", etp, "b1.65")]][[errtypes[[err]]$bsspread]]
                resbs$plaq3[, 4 * index + 0] <- tmp[[paste0("etp", etp, "b1.7")]][[errtypes[[err]]$bsspread]]
                tmp <- try(readRDS(sprintf("%s/contlimittype%smode%s%saveraged.RData", "plotstikz", "plaq", mode, ifelse(xiinter, "xiinter", ""))))
                resbs$plaq16[, 4 * index - 3] <- tmp[[paste0("etp", etp, "b1.65")]][[errtypes[[err]]$bs]]
                resbs$plaq16[, 4 * index - 2] <- tmp[[paste0("etp", etp, "b1.7")]][[errtypes[[err]]$bs]]
                resbs$plaq16[, 4 * index - 1] <- tmp[[paste0("etp", etp, "b1.65")]][[errtypes[[err]]$bsspread]]
                resbs$plaq16[, 4 * index + 0] <- tmp[[paste0("etp", etp, "b1.7")]][[errtypes[[err]]$bsspread]]
                tmp2 <- try(readRDS(sprintf("%s/contlimittype%smode%s%saveraged.RData", "plotstikz", "ratio", mode, ifelse(xiinter, "xiinter", ""))))
                resbs$plaq3ratio[, 4 * index - 3] <- tmp[[paste0("etp", etp, "b1.65")]][[errtypes[[err]]$bs]] * tmp2[[paste0("etp", etp, "b1.65")]][[errtypes[[err]]$bs]]
                resbs$plaq3ratio[, 4 * index - 2] <- tmp[[paste0("etp", etp, "b1.7")]][[errtypes[[err]]$bs]] * tmp2[[paste0("etp", etp, "b1.7")]][[errtypes[[err]]$bs]]
                resbs$plaq3ratio[, 4 * index - 1] <- tmp[[paste0("etp", etp, "b1.65")]][[errtypes[[err]]$bsspread]] * tmp2[[paste0("etp", etp, "b1.65")]][[errtypes[[err]]$bsspread]]
                resbs$plaq3ratio[, 4 * index + 0] <- tmp[[paste0("etp", etp, "b1.7")]][[errtypes[[err]]$bsspread]] * tmp2[[paste0("etp", etp, "b1.7")]][[errtypes[[err]]$bsspread]]
                bsnames <- append(bsnames, paste0("errmode", err, "mode", mode, "xiinter", xiinter, "etp", etp, "beta", rep(c(1.65, 1.7), 2), "spread", c(F, F, T, T)))
                dplaq3ratio <- append(dplaq3ratio, apply(resbs$plaq3ratio[, 4 * index + (-3:0)], MARGIN = 2, FUN = sd, na.rm = T))
                index <- index + 1
            }
        }
    }
}
data$dplaq3ratiocontlim <- dplaq3ratio
indices <- data.frame(uprange = 45, lowrange = 35, bsindex = seq_along(data$betaiso), resindex = seq_along(data$betaiso))
# indices
mask <- seq_along(data$betaiso)
# mask <- 100:110
indices$resindex[mask] <- seq_along(mask)
indices$bsindex[mask] <- seq_along(mask)
resellipsemean <- getmatchingellipse(data = data, bsdata = resbs, hamres = hamiltoniandata, indices = indices[mask, ], verbose = F)
resellipsemean$mode <- data$mode[mask]
resellipsemean$xiinter <- data$xiinter[mask]
resellipsemean$spread <- data$spread[mask]
resellipsemean$errtype <- data$errtype[mask]
resellipsemean$etp <- data$etp[mask]
resellipsemean

data$bootstrapmeanplaq3 <- apply(resbs$plaq3, MARGIN = 2, FUN = mean, na.rm = T)
data$medianplaq3 <- apply(resbs$plaq3, MARGIN = 2, FUN = median, na.rm = T)
data$q16plaq3 <- data$medianplaq3 - apply(resbs$plaq3, MARGIN = 2, FUN = quantile, probs = 0.16, na.rm = T)
data$q84plaq3 <- apply(resbs$plaq3, MARGIN = 2, FUN = quantile, probs = 0.84, na.rm = T) - data$medianplaq3
data$bootstrapmeanplaq16 <- apply(resbs$plaq16, MARGIN = 2, FUN = mean, na.rm = T)
data$medianplaq16 <- apply(resbs$plaq16, MARGIN = 2, FUN = median, na.rm = T)
data$q16plaq16 <- data$medianplaq16 - apply(resbs$plaq16, MARGIN = 2, FUN = quantile, probs = 0.16, na.rm = T)
data$q84plaq16 <- apply(resbs$plaq16, MARGIN = 2, FUN = quantile, probs = 0.84, na.rm = T) - data$medianplaq16
data$bootstrapmeanbeta <- apply(resbs$beta, MARGIN = 2, FUN = mean, na.rm = T)
data$medianbeta <- apply(resbs$beta, MARGIN = 2, FUN = median, na.rm = T)
data$q16beta <- data$medianbeta - apply(resbs$beta, MARGIN = 2, FUN = quantile, probs = 0.16, na.rm = T)
data$q84beta <- apply(resbs$beta, MARGIN = 2, FUN = quantile, probs = 0.84, na.rm = T) - data$medianbeta
data$bootstrapmeanplaq3ratio <- apply(resbs$plaq3ratio, MARGIN = 2, FUN = mean, na.rm = T)
data$medianplaq3ratio <- apply(resbs$plaq3ratio, MARGIN = 2, FUN = median, na.rm = T)
data$q16plaq3ratio <- data$medianplaq3ratio - apply(resbs$plaq3ratio, MARGIN = 2, FUN = quantile, probs = 0.16, na.rm = T)
data$q84plaq3ratio <- apply(resbs$plaq3ratio, MARGIN = 2, FUN = quantile, probs = 0.84, na.rm = T) - data$medianplaq3ratio
data$cor <- onlydiagonalcorelements(resbs$plaq3, resbs$beta)
# bsnames
head(data)
# data
# stop()

## repeat the ellipse parameters with median and quantiles and median and sd


resellipsemedian <- getmatchingellipse(
    data = data.frame(betacontlim = data$medianbeta, dbetacontlim = data$q16beta, plaq3contlim = data$medianplaq3, dplaq3contlim = data$q84plaq3, betaiso = data$betaiso),
    bsdata = resbs, hamres = hamiltoniandata, indices = indices[mask, ], verbose = F
)
dim(resellipsemedian)
resellipsemedian$mode <- data$mode[mask]
resellipsemedian$xiinter <- data$xiinter[mask]
resellipsemedian$spread <- data$spread[mask]
resellipsemedian$errtype <- data$errtype[mask]
resellipsemedian$etp <- data$etp[mask]
resellipsemedian

resellipsemedianwithsd <- getmatchingellipse(
    data = data.frame(
        betacontlim = data$medianbeta, dbetacontlim = data$dbetacontlim,
        plaq3contlim = data$medianplaq3, dplaq3contlim = data$dplaq3contlim, betaiso = data$betaiso
    ),
    bsdata = resbs, hamres = hamiltoniandata, indices = indices[mask, ], verbose = F
)
dim(resellipsemedianwithsd)
resellipsemedianwithsd$mode <- data$mode[mask]
resellipsemedianwithsd$xiinter <- data$xiinter[mask]
resellipsemedianwithsd$spread <- data$spread[mask]
resellipsemedianwithsd$errtype <- data$errtype[mask]
resellipsemedianwithsd$etp <- data$etp[mask]
resellipsemedianwithsd

resellipseratio <- getmatchingellipse(
    data = data.frame(betacontlim = data$betacontlim, dbetacontlim = data$dbetacontlim, plaq3contlim = data$plaq3ratiocontlim, dplaq3contlim = data$dplaq3ratiocontlim, betaiso = data$betaiso),
    bsdata = list(beta = resbs$beta, plaq3 = resbs$plaq3ratio), hamres = hamiltoniandata, indices = indices[mask, ], verbose = F
)
dim(resellipseratio)
resellipseratio$mode <- data$mode[mask]
resellipseratio$xiinter <- data$xiinter[mask]
resellipseratio$spread <- data$spread[mask]
resellipseratio$errtype <- data$errtype[mask]
resellipseratio$etp <- data$etp[mask]
resellipseratio

resellipseratiomedian <- getmatchingellipse(
    data = data.frame(betacontlim = data$medianbeta, dbetacontlim = data$q16beta, plaq3contlim = data$medianplaq3ratio, dplaq3contlim = data$q84plaq3ratio, betaiso = data$betaiso),
    bsdata = list(beta = resbs$beta, plaq3 = resbs$plaq3ratio), hamres = hamiltoniandata, indices = indices[mask, ], verbose = F
)
dim(resellipseratiomedian)
resellipseratiomedian$mode <- data$mode[mask]
resellipseratiomedian$xiinter <- data$xiinter[mask]
resellipseratiomedian$spread <- data$spread[mask]
resellipseratiomedian$errtype <- data$errtype[mask]
resellipseratiomedian$etp <- data$etp[mask]
resellipseratiomedian

resellipseratiomedianwithsd <- getmatchingellipse(
    data = data.frame(
        betacontlim = data$medianbeta, dbetacontlim = data$dbetacontlim,
        plaq3contlim = data$medianplaq3ratio, dplaq3contlim = data$dplaq3ratiocontlim, betaiso = data$betaiso
    ),
    bsdata = list(beta = resbs$beta, plaq3 = resbs$plaq3ratio), hamres = hamiltoniandata, indices = indices[mask, ], verbose = F
)
dim(resellipseratiomedianwithsd)
resellipseratiomedianwithsd$mode <- data$mode[mask]
resellipseratiomedianwithsd$xiinter <- data$xiinter[mask]
resellipseratiomedianwithsd$spread <- data$spread[mask]
resellipseratiomedianwithsd$errtype <- data$errtype[mask]
resellipseratiomedianwithsd$etp <- data$etp[mask]
resellipseratiomedianwithsd

saveRDS(
    list(
        data = data, bs = resbs, resellipsemean = resellipsemean, resellipseratio = resellipseratio,
        resellipsemedian = resellipsemedian, resellipseratiomedian = resellipseratiomedian,
        resellipsemedianwithsd = resellipsemedianwithsd, resellipseratiomedianwithsd = resellipseratiomedianwithsd,
        bsnames = bsnames
    ),
    file = "plotstikz/ellipseparameters.RData"
)

data <- data[mask, ]

listellipse <- split(cbind(resellipsemean, data), f = seq(nrow(resellipsemean)))
pdf("plotellipse_prettypicture.pdf")
plot(NA, xlim = c(1.35, 1.52), ylim = c(0.6, 0.67))
plotlims <- lapply(listellipse, FUN = function(x) draw_ellipse_general(meanx = x$betacontlim, meany = x$plaq3contlim, radx = x$dbetacontlim, rady = x$dplaq3contlim, phi = x$theta, nstd = x$devstd, rep = T, cex = 0.01, points = 5000))
# plotwitherror(x = data$betacontlim, y = data$plaq3contlim, dy = data$dplaq3contlim, dx = data$dbetacontlim, col = "blue", rep = T)
points(x = hamiltoniandata$lowerx, y = hamiltoniandata$lower, type = "o", col = "red")

warnings()
