library("hadron")

source("~/Documents/masterthesis/analyse_potential/myfunctions.R")
source("~/Documents/masterthesis/analyse_potential/matchwithellipse.R")

hamiltoniandata <- read.table("/home/gross/Documents/masterthesis/more_measurements/hamiltonian/interpolate_ham.csv", header = T)
modes <- c("all", "xi0.20", "xi0.19", "xi0.18")
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
            }
        }
    }
}
averagedrestable$interpolated_xi <- averagedrestable$xiinter | averagedrestable$type == "plaqinter" | averagedrestable$type == "betainter"
averagedrestable$interpolated_beta <- averagedrestable$type == "plaqinter" | averagedrestable$type == "betainter"


pdf(sprintf("rescombined.pdf"), title = "")
modes <- c("all", "xi0.20", "xi0.19", "xi0.18")
cols <- c("black", "red", "blue", "darkgreen")
pchs <- c(21, 22, 23, 24)

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


resbs <- list(beta = array(NA, dim = c(500, length(modes) * 16 * 3)), plaq3 = array(NA, dim = c(500, length(modes) * 16 * 3)))
data <- data.frame(
    betacontlim = c(), dbetacontlim = c(),
    plaq3contlim = c(), dplaq3contlim = c(),
    plaq16contlim = c(), dplaq16contlim = c(),
    betaiso = c(), mode = c(), xiinter = c(), spread = c(), errtype = c()
)
index <- 1
bsnames <- c()
for (err in c("AIC", "error", "unweighted")) {
    for (mode in modes) {
        for (xiinter in c(T, F)) {
            for (etp in c(1, 0)) {
                data <- rbind(data, data.frame(
                    betacontlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dbetacontlim = averagedrestable[, errtypes[[err]]$sd][averagedrestable$mode == mode & averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    plaq3contlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "plaqsmall" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dplaq3contlim = averagedrestable[, errtypes[[err]]$sd][averagedrestable$mode == mode & averagedrestable$type == "plaqsmall" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    plaq16contlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dplaq16contlim = averagedrestable[, errtypes[[err]]$sd][averagedrestable$mode == mode & averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    betaiso = c(1.65, 1.7), mode = mode, xiinter = xiinter, spread = F, errtype = err, etp = etp
                ))
                data <- rbind(data, data.frame(
                    betacontlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dbetacontlim = averagedrestable[, errtypes[[err]]$spread][averagedrestable$mode == mode & averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    plaq3contlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "plaqsmall" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dplaq3contlim = averagedrestable[, errtypes[[err]]$spread][averagedrestable$mode == mode & averagedrestable$type == "plaqsmall" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    plaq16contlim = averagedrestable[, errtypes[[err]]$mean][averagedrestable$mode == mode & averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    dplaq16contlim = averagedrestable[, errtypes[[err]]$spread][averagedrestable$mode == mode & averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter & averagedrestable$etp == etp],
                    betaiso = c(1.65, 1.7), mode = mode, xiinter = xiinter, spread = T, errtype = err, etp = etp
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
                bsnames <- append(bsnames, paste0("errmode", err, "mode", mode, "xiinter", xiinter, "etp", etp, "beta", rep(c(1.65, 1.7), 2), "spread", c(F, F, T, T)))
                index <- index + 1
            }
        }
    }
}

indices <- data.frame(uprange = 45, lowrange = 35, bsindex = seq_along(data$betaiso), resindex = seq_along(data$betaiso))
# data
# indices
mask <- seq_along(data$betaiso)
# mask <- 60:96
indices$resindex[mask] <- seq_along(mask)
indices$bsindex[mask] <- seq_along(mask)
resellipse <- getmatchingellipse(data = data[mask, ], bsdata = resbs[, mask], hamres = hamiltoniandata, indices = indices[mask, ], verbose = F)
resellipse$mode <- data$mode[mask]
resellipse$xiinter <- data$xiinter[mask]
resellipse$spread <- data$spread[mask]
resellipse$errtype <- data$errtype[mask]
resellipse$etp <- data$etp[mask]
resellipse
# bsnames
saveRDS(list(data = data, bs = resbs, resellipse = resellipse, bsnames = bsnames), file = "plotstikz/ellipseparameters.RData")

data <- data[mask, ]

listellipse <- split(cbind(resellipse, data), f = seq(nrow(resellipse)))
pdf("plotellipse_prettypicture.pdf")
plot(NA, xlim = c(1.35, 1.52), ylim = c(0.6, 0.67))
plotlims <- lapply(listellipse, FUN = function(x) draw_ellipse_general(meanx = x$betacontlim, meany = x$plaq3contlim, radx = x$dbetacontlim, rady = x$dplaq3contlim, phi = x$theta, nstd = x$devstd, rep = T, cex = 0.01, points = 5000))
# plotwitherror(x = data$betacontlim, y = data$plaq3contlim, dy = data$dplaq3contlim, dx = data$dbetacontlim, col = "blue", rep = T)
points(x = hamiltoniandata$lowerx, y = hamiltoniandata$lower, type = "o", col = "red")

warnings()
