path <- "~/Documents/masterthesis/more_measurements/heatbath/contlim/plotstikz"
pdf(sprintf("%s/determineallcorrelationscontlimit.pdf", path), title = "")
res <- data.frame(beta = c(), etp = c(), omit = c(), pot = c(), xiinter = c(), uplim = c(), lowlim = c(), degree = c(), covpb = c(), covp3b = c(), covrb = c())
for (betaiso in c(1.65, 1.70)) {
    for (etp in c(0, 1)) {
        for (omit in c(0, 1)) {
            for (pot in c("normal", "sideways")) {
                for (xiinter in c(T, F)) {
                    filename <- sprintf("beta%fomit%dxiconstllxi1llr00fl0.30aicscaletauintetp%d%scont1maxiter300", betaiso, omit, etp, ifelse(xiinter, "xiinter", ""))
                    dataplaqbeta <- readRDS(sprintf("%s/%s%s%s.RData", path, "listpolynomialrenormalizationbetachosen", pot, filename))
                    dataplaq3 <- readRDS(sprintf("%s/%s%s%s.RData", path, "listpolynomialrenormalizationbetachosenL3", pot, filename))
                    stopifnot(length(dataplaqbeta$fits) == length(dataplaq3$fits))
                    for (i in seq_along(dataplaq3$uplimfit)) {
                        for (degree in seq(1, length(dataplaqbeta$fits[[i]]) / 3)) {
                            plaq <- dataplaqbeta$fits[[i]][[paste0("plaq", degree)]]
                            beta <- dataplaqbeta$fits[[i]][[paste0("plaq", degree)]]
                            plaq3 <- dataplaq3$fits[[i]][[paste0("plaqsmall", degree)]]
                            ratio <- dataplaq3$fits[[i]][[paste0("ratio", degree)]]
                            # print(dataplaqbeta$lowlim[[i]])
                            # print(plaq3)
                            newline <- try(data.frame(
                                beta = betaiso, etp = etp, omit = omit, pot = pot, xiinter = xiinter, uplim = dataplaqbeta$uplim[i], lowlim = dataplaqbeta$lowlim[i], degree = degree,
                                covpb = cor(plaq$t[, 1], beta$t[, 1], use = "na.or.complete"), covp3b = cor(plaq3$t[, 1], beta$t[, 1], use = "na.or.complete"),
                                covrb = cor(ratio$t[, 1], beta$t[, 1], use = "na.or.complete")
                            ))
                            # print(lapply(X=c(
                            #         betaiso, etp, omit, pot, xiinter, dataplaqbeta$uplim[[i]], dataplaqbeta$lowlim[[i]], degree), FUN=class))
                            if (inherits(newline, "try-error")) {
                                print(sprintf(
                                    "problem with beta %.2f etp %d omit %d pot %s xiinter %d uplim %d lowlim %d degree %d",
                                    betaiso, etp, omit, pot, xiinter, dataplaqbeta$uplim[[i]], dataplaqbeta$lowlim[[i]], degree
                                ))
                            } else {
                                res <- rbind(res, newline)
                                for (dat in list(plaq$t[, 1], beta$t[, 1], plaq3$t[, 1], ratio$t[, 1])) {
                                    qqnorm(y = dat, main = sprintf(
                                        "beta %.2f etp %d omit %d pot %s xiinter %d uplim %d lowlim %d degree %d",
                                        betaiso, etp, omit, pot, xiinter, dataplaqbeta$uplim[[i]], dataplaqbeta$lowlim[[i]], degree
                                    ))
                                    qqline(y = dat)
                                }
                            }
                        }
                    }
                    # print(res)
                }
            }
        }
    }
}
res
write.table(res, file = sprintf("%s/determineallcorrelationscontlimit.csv", path), row.names = F)
