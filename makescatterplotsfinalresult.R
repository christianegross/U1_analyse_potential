library("hadron")

source("~/Documents/masterthesis/analyse_potential/myfunctions.R")

for (mode in c("all", "xi0.20", "xi0.19", "xi0.18")) {
    pdf(sprintf("res%s.pdf", mode), title = "")
    plaq <- read.table(sprintf("%s/contlimittypeplaqmode%snotaveraged.csv", "plotstikz", mode))
    plaqsmall <- read.table(sprintf("%s/contlimittypeplaqsmallmode%snotaveraged.csv", "plotstikz", mode))
    beta <- read.table(sprintf("%s/contlimittypebetamode%snotaveraged.csv", "plotstikz", mode))
    plaqAIC <- read.table(sprintf("%s/contlimittypeplaqmode%saveraged.csv", "plotstikz", mode))
    plaqsmallAIC <- read.table(sprintf("%s/contlimittypeplaqsmallmode%saveraged.csv", "plotstikz", mode))
    betaAIC <- read.table(sprintf("%s/contlimittypebetamode%saveraged.csv", "plotstikz", mode))
    plotwitherror(
        x = beta$contlim, dx = beta$dcontlim, y = plaq$contlim, dy = plaq$dcontlim,
        xlim = range(beta$contlim + beta$dcontlim, beta$contlim - beta$dcontlim),
        ylim = range(plaqsmall$contlim + plaqsmall$dcontlim, plaq$contlim - plaq$dcontlim),
        col = ifelse(plaq$beta == 1.65, 1, 2),
        xlab = "beta", ylab = "P"
    )
    plotwitherror(
        x = beta$contlim, dx = beta$dcontlim, y = plaqsmall$contlim, dy = plaqsmall$dcontlim,
        xlim = range(beta$contlim + beta$dcontlim, beta$contlim - beta$dcontlim),
        col = ifelse(plaqsmall$beta == 1.65, 1, 2),
        pch = 2, rep = T
    )
    par(lwd = 3)
    plotwitherror(
        x = betaAIC$contlim, dx = betaAIC$dcontlim, y = plaqAIC$contlim, dy = plaqAIC$dcontlim,
        xlim = range(betaAIC$contlim + betaAIC$dcontlim, betaAIC$contlim - betaAIC$dcontlim),
        ylim = range(plaqAIC$contlim + plaqAIC$dcontlim, plaqAIC$contlim - plaqAIC$dcontlim),
        col = ifelse(plaqAIC$beta == 1.65, 1, 2),
        rep = T, lwd = 3, cex = 3
    )
    plotwitherror(
        x = betaAIC$contlim, dx = betaAIC$dcontlim, y = plaqsmallAIC$contlim, dy = plaqsmallAIC$dcontlim,
        xlim = range(betaAIC$contlim + betaAIC$dcontlim, betaAIC$contlim - betaAIC$dcontlim),
        ylim = range(plaqsmallAIC$contlim + plaqsmallAIC$dcontlim, plaqsmallAIC$contlim - plaqsmallAIC$dcontlim),
        col = ifelse(plaqsmallAIC$beta == 1.65, 1, 2),
        rep = T, lwd = 3, cex = 3, pch = 2
    )
    par(lwd = 1)
}
