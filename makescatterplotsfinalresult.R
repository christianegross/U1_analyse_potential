library("hadron")

source("~/Documents/masterthesis/analyse_potential/myfunctions.R")

if (FALSE) {
    for (mode in c("all", "xi0.20", "xi0.19", "xi0.18")) {
        pdf(sprintf("res%s.pdf", mode), title = "")
        print(mode)
        plaq <- read.table(sprintf("%s/contlimittypeplaqmode%snotaveraged.csv", "plotstikz", mode))
        plaqsmall <- read.table(sprintf("%s/contlimittypeplaqsmallmode%snotaveraged.csv", "plotstikz", mode))
        beta <- read.table(sprintf("%s/contlimittypebetamode%snotaveraged.csv", "plotstikz", mode))
        plaqaveraged <- read.table(sprintf("%s/contlimittypeplaqmode%saveraged.csv", "plotstikz", mode))
        plaqsmallaveraged <- read.table(sprintf("%s/contlimittypeplaqsmallmode%saveraged.csv", "plotstikz", mode))
        betaaveraged <- read.table(sprintf("%s/contlimittypebetamode%saveraged.csv", "plotstikz", mode))
        number <- length(plaq$beta) / 2
        # print(dim(plaq))
        # print(dim(beta))
        # print(dim(plaqsmall))
        # print(data.frame(x = beta$contlim[beta$etp==1], y = plaq$contlim[beta$etp==1],
        #     dx = array(c(beta$dcontlim[beta$etp==0], sqrt(beta$dcontlimspread[beta$etp==1]^2 - beta$dcontlim[beta$etp==0]^2)), dim=c(number, 2)),
        #     dy = array(c(plaq$dcontlim[plaq$etp==0], sqrt(plaq$dcontlimspread[plaq$etp==1]^2 - plaq$dcontlim[plaq$etp==0]^2)), dim=c(number, 2))))
        # print(names(beta))
        # print(names(plaq))
        # print(range(beta$contlim + beta$dcontlim, beta$contlim - beta$dcontlim))
        # print(range(plaqsmall$contlim + plaqsmall$dcontlim, plaq$contlim - plaq$dcontlim))
        # print(beta)
        # print(plaq)
        plotwitherror(
            x = beta$contlim[beta$etp == 1], y = plaq$contlim[beta$etp == 1],
            dx = array(c(beta$dcontlim[beta$etp == 0], sqrt(beta$dcontlim[beta$etp == 1]^2 - beta$dcontlim[beta$etp == 0]^2)), dim = c(number, 2)),
            dy = array(c(plaq$dcontlim[plaq$etp == 0], sqrt(plaq$dcontlim[plaq$etp == 1]^2 - plaq$dcontlim[plaq$etp == 0]^2)), dim = c(number, 2)),
            errsum.method = "quadrature",
            xlim = range(beta$contlim + beta$dcontlim, beta$contlim - beta$dcontlim),
            ylim = range(plaqsmall$contlim + plaqsmall$dcontlim, plaq$contlim - plaq$dcontlim),
            col = ifelse(plaq$beta[beta$etp == 1] == 1.65, 1, 2),
            xlab = "beta", ylab = "P"
        )
        # print(data.frame(x = beta$contlim[beta$etp==1], y = plaqsmall$contlim[beta$etp==1],
        #     dx = array(c(beta$dcontlim[beta$etp==0], sqrt(beta$dcontlim[beta$etp==1]^2 - beta$dcontlim[beta$etp==0]^2)), dim=c(number, 2)),
        #     dy = array(c(plaqsmall$dcontlim[plaqsmall$etp==0], sqrt(plaqsmall$dcontlim[plaqsmall$etp==1]^2 - plaqsmall$dcontlim[plaqsmall$etp==0]^2)), dim=c(number, 2))))
        plotwitherror(
            x = beta$contlim[beta$etp == 1], y = plaqsmall$contlim[beta$etp == 1],
            dx = array(c(beta$dcontlim[beta$etp == 0], sqrt(beta$dcontlim[beta$etp == 1]^2 - beta$dcontlim[beta$etp == 0]^2)), dim = c(number, 2)),
            dy = array(c(plaqsmall$dcontlim[plaqsmall$etp == 0], sqrt(plaqsmall$dcontlim[plaqsmall$etp == 1]^2 - plaqsmall$dcontlim[plaqsmall$etp == 0]^2)), dim = c(number, 2)),
            errsum.method = "quadrature",
            col = ifelse(plaqsmall$beta[beta$etp == 1] == 1.65, 1, 2),
            pch = 2, rep = T
        )
        par(lwd = 3)
        # print(data.frame(
        #     x = betaaveraged$contlim[betaaveraged$etp==1], y = plaqaveraged$contlim[betaaveraged$etp==1],
        #     dx = array(c(betaaveraged$dcontlim[betaaveraged$etp==0], sqrt(betaaveraged$dcontlimspread[betaaveraged$etp==1]^2 - betaaveraged$dcontlim[betaaveraged$etp==0]^2)), dim=c(2, 2)),
        #     dy = array(c(plaqaveraged$dcontlim[plaqaveraged$etp==0], sqrt(plaqaveraged$dcontlimspread[plaqaveraged$etp==1]^2 - plaqaveraged$dcontlim[plaqaveraged$etp==0]^2)), dim=c(2, 2))))
        plotwitherror(
            x = betaaveraged$contlim[betaaveraged$etp == 1], y = plaqaveraged$contlim[betaaveraged$etp == 1],
            dx = array(c(betaaveraged$dcontlim[betaaveraged$etp == 0], sqrt(betaaveraged$dcontlimspread[betaaveraged$etp == 1]^2 - betaaveraged$dcontlim[betaaveraged$etp == 0]^2)), dim = c(2, 2)),
            dy = array(c(plaqaveraged$dcontlim[plaqaveraged$etp == 0], sqrt(plaqaveraged$dcontlimspread[plaqaveraged$etp == 1]^2 - plaqaveraged$dcontlim[plaqaveraged$etp == 0]^2)), dim = c(2, 2)),
            errsum.method = "quadrature",
            col = ifelse(plaqaveraged$beta[beta$etp == 1] == 1.65, 1, 2),
            rep = T, lwd = 3, cex = 3
        )
        # print(data.frame( x = betaaveraged$contlim[betaaveraged$etp==1], y = plaqsmallaveraged$contlim[betaaveraged$etp==1],
        #     dx = array(c(betaaveraged$dcontlim[betaaveraged$etp==0], sqrt(betaaveraged$dcontlimspread[betaaveraged$etp==1]^2 - betaaveraged$dcontlim[betaaveraged$etp==0]^2)), dim=c(2, 2)),
        #     dy = array(c(plaqsmallaveraged$dcontlim[plaqsmallaveraged$etp==0], sqrt(plaqsmallaveraged$dcontlimspread[plaqsmallaveraged$etp==1]^2 - plaqsmallaveraged$dcontlim[plaqsmallaveraged$etp==0]^2)), dim=c(2, 2))))
        plotwitherror(
            x = betaaveraged$contlim[betaaveraged$etp == 1], y = plaqsmallaveraged$contlim[betaaveraged$etp == 1],
            dx = array(c(betaaveraged$dcontlim[betaaveraged$etp == 0], sqrt(betaaveraged$dcontlimspread[betaaveraged$etp == 1]^2 - betaaveraged$dcontlim[betaaveraged$etp == 0]^2)), dim = c(2, 2)),
            dy = array(c(plaqsmallaveraged$dcontlim[plaqsmallaveraged$etp == 0], sqrt(plaqsmallaveraged$dcontlimspread[plaqsmallaveraged$etp == 1]^2 - plaqsmallaveraged$dcontlim[plaqsmallaveraged$etp == 0]^2)), dim = c(2, 2)),
            errsum.method = "quadrature",
            col = ifelse(plaqsmallaveraged$beta[beta$etp == 1] == 1.65, 1, 2),
            rep = T, lwd = 3, cex = 3, pch = 2
        )
        par(lwd = 1)
    }
}



hamiltoniandata <- read.table("/home/gross/Documents/masterthesis/more_measurements/hamiltonian/interpolate_ham.csv", header=T)
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
    AIC = list(mean = "contlimAIC", sd = "dcontlimAIC", spread = "dcontlimspreadAIC"),
    error = list(mean = "contlimerror", sd = "dcontlimerror", spread = "dcontlimspreaderror"),
    unweighted = list(mean = "contlimunweighted", sd = "dcontlimunweighted", spread = "dcontlimspreadunweighted")
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
            print(errtypes[[err]]$mean)
            # print(names(averagedrestable))
            plaqmask <- averagedrestable$type == "plaq" & averagedrestable$xiinter == xiinter & averagedrestable$mode == mode
            plaqsmallmask <- averagedrestable$type == "plaqsmall" & averagedrestable$xiinter == xiinter & averagedrestable$mode == mode
            betamask <- averagedrestable$type == "beta" & averagedrestable$xiinter == xiinter & averagedrestable$mode == mode
            # print(averagedrestable[plaqmask, ])
            # print(array(c(averagedrestable[, errtypes[[err]]$sd][plaqmask & averagedrestable$etp == 0],
            #     sqrt(averagedrestable[, errtypes[[err]]$sd][plaqmask & averagedrestable$etp == 1]^2 - averagedrestable[, errtypes[[err]]$sd][plaqmask & averagedrestable$etp == 0]^2)),
            #     dim=c(2, 2)))
            # print(data.frame(
            #     x = averagedrestable[, errtypes[[err]]$mean][plaqmask & averagedrestable$etp == 1],
            #     y = averagedrestable[, errtypes[[err]]$mean][betamask & averagedrestable$etp == 1],
            #     dx=array(c(averagedrestable[, errtypes[[err]]$sd][plaqmask & averagedrestable$etp == 0],
            #     sqrt(averagedrestable[, errtypes[[err]]$sd][plaqmask & averagedrestable$etp == 1]^2 - averagedrestable[, errtypes[[err]]$sd][plaqmask & averagedrestable$etp == 0]^2)),
            #     dim=c(2, 2))))
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
        points(x=hamiltoniandata$lowerx, y=hamiltoniandata$lower, type="o")
        legend("topleft", legend = c(modes, "H"), col = c(cols, "black"), pch = c(pchs, 1))
    }
}


# pdf(sprintf("rescombined_old.pdf"), title = "")
# modes <- c("all", "xi0.20", "xi0.19", "xi0.18")
# cols <- c("black", "red", "blue", "darkgreen")
# pchs <- c(21, 22, 23, 24)
# for(xiinter in c(T, F)) {
# plot(NA, xlim=c(1.35, 1.5), ylim=c(0.58, 0.68), xlab="beta", ylab="P", main=paste("xiinter", xiinter))
# for (index in seq_along(modes)) {
#     mode <- modes[index]
#     print(mode)
#     plaqaveraged <- read.table(sprintf("%s/contlimittypeplaqmode%s%saveraged.csv", "plotstikz", mode, ifelse(xiinter, "xiinter", "")))
#     plaqsmallaveraged <- read.table(sprintf("%s/contlimittypeplaqsmallmode%s%saveraged.csv", "plotstikz", mode,ifelse(xiinter, "xiinter", "")))
#     betaaveraged <- read.table(sprintf("%s/contlimittypebetamode%s%saveraged.csv", "plotstikz", mode,ifelse(xiinter, "xiinter", "")))
#     # print(data.frame(
#     #     x = betaaveraged$contlim[betaaveraged$etp==1], y = plaqaveraged$contlim[betaaveraged$etp==1],
#     #     dx = array(c(betaaveraged$dcontlim[betaaveraged$etp==0], sqrt(betaaveraged$dcontlimspread[betaaveraged$etp==1]^2 - betaaveraged$dcontlim[betaaveraged$etp==0]^2)), dim=c(2, 2)),
#     #     dy = array(c(plaqaveraged$dcontlim[plaqaveraged$etp==0], sqrt(plaqaveraged$dcontlimspread[plaqaveraged$etp==1]^2 - plaqaveraged$dcontlim[plaqaveraged$etp==0]^2)), dim=c(2, 2))))
#     plotwitherror(
#         x = betaaveraged$contlim[betaaveraged$etp==1], y = plaqaveraged$contlim[betaaveraged$etp==1],
#         dx = array(c(betaaveraged$dcontlim[betaaveraged$etp==0], sqrt(betaaveraged$dcontlimspread[betaaveraged$etp==1]^2 - betaaveraged$dcontlim[betaaveraged$etp==0]^2)), dim=c(2, 2)),
#         dy = array(c(plaqaveraged$dcontlim[plaqaveraged$etp==0], sqrt(plaqaveraged$dcontlimspread[plaqaveraged$etp==1]^2 - plaqaveraged$dcontlim[plaqaveraged$etp==0]^2)), dim=c(2, 2)),
#         errsum.method = "quadrature",
#         col = cols[index], pch=pchs[index], bg=cols[index],
#         rep = T
#     )
#     # print(data.frame( x = betaaveraged$contlim[betaaveraged$etp==1], y = plaqsmallaveraged$contlim[betaaveraged$etp==1],
#     #     dx = array(c(betaaveraged$dcontlim[betaaveraged$etp==0], sqrt(betaaveraged$dcontlimspread[betaaveraged$etp==1]^2 - betaaveraged$dcontlim[betaaveraged$etp==0]^2)), dim=c(2, 2)),
#     #     dy = array(c(plaqsmallaveraged$dcontlim[plaqsmallaveraged$etp==0], sqrt(plaqsmallaveraged$dcontlimspread[plaqsmallaveraged$etp==1]^2 - plaqsmallaveraged$dcontlim[plaqsmallaveraged$etp==0]^2)), dim=c(2, 2))))
#     plotwitherror(
#         x = betaaveraged$contlim[betaaveraged$etp==1], y = plaqsmallaveraged$contlim[betaaveraged$etp==1],
#         dx = array(c(betaaveraged$dcontlim[betaaveraged$etp==0], sqrt(betaaveraged$dcontlimspread[betaaveraged$etp==1]^2 - betaaveraged$dcontlim[betaaveraged$etp==0]^2)), dim=c(2, 2)),
#         dy = array(c(plaqsmallaveraged$dcontlim[plaqsmallaveraged$etp==0], sqrt(plaqsmallaveraged$dcontlimspread[plaqsmallaveraged$etp==1]^2 - plaqsmallaveraged$dcontlim[plaqsmallaveraged$etp==0]^2)), dim=c(2, 2)),
#         errsum.method = "quadrature",
#         col = cols[index], pch=pchs[index], bg=cols[index],
#         rep = T
#     )
# }
# legend("topleft", legend=modes, col=cols, pch=pchs)
# }


warnings()
