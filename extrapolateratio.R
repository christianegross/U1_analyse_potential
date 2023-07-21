library("hadron")
source("/hiskp4/gross/masterthesis/analyse/code/U1_analyse_potential/myfunctions.R")
args <-  commandArgs(trailingOnly = TRUE)

bootsamples <-  500

beta <-  as.numeric(args[1])
S <-  as.numeric(args[2])

if (beta == 1.65) {

betas <-  c(1.65, 1.615, 1.58, 1.515, 1.48, 1.46, 1.4378) # beta = 1.65
columns <-  c(1, 1, 11, 1, 1, 1, 1) # beta = 1.65
# betaextrapolate <-  c(1.45, 1.47) # beta = 1.65
# Ntextrapolate <-  48 # beta = 1.65
betaextrapolate <- c()
Ntextrapolate <- 0
}

if (beta == 1.7) {

betas <-  c(1.7, 1.65, 1.6075, 1.55, 1.525, 1.5, 1.49) # beta = 1.7
columns <-  c(11, 11, 1, 11, 1, 1, 1) # beta = 1.7
betaextrapolate <-  c() # beta = 1.7
Ntextrapolate <-  0 # beta = 1.7
}

Nt <-  c(16, 20, 24, 32, 40, 48, 64)
therm16 <-  c(500, 500, 500, 500, 500, 1000, 2000)
therm3 <-  c(1333, 1333, 1333, 1333, 1333, 2667, 2667)

P16 <-  c()
dP16 <-  c()
P3 <-  c()
dP3 <-  c()
resxi <-  read.table(sprintf("plotstikz/resultsrenormalizationbeta%f.csv", betas[1]), header = TRUE)
# TODO: read in xi from file
# print(resxi)
# print(resxi$xiphys)
# print(resxi$dxiphys)
# tableresname <-  "resultsummary2p1dnormalb1.700Ns16_170723.csv"
tableresname <-  "resultsummary2p1dnormalb1.700Ns16.csv"
tableres <-  read.table(tableresname, header = T)

resultxi <-  data.frame()


for (i in seq(1, length(betas))) {
    # print(Nt[i])
    if (Nt[i] !=  Ntextrapolate) {
    mask <-  tableres$beta == betas[i] & abs(tableres$xi - 16 / Nt[i]) < 1e-5
    resultxi <-  rbind(resultxi, data.frame(xiphys = tableres$xicalc[mask], dxiphys = tableres$dxicalc[mask]))
    } else {
    mask1 <-  tableres$beta == betaextrapolate[1] & abs(tableres$xi - 16 / Nt[i]) < 1e-5
    mask2 <-  tableres$beta == betaextrapolate[2] & abs(tableres$xi - 16 / Nt[i]) < 1e-5
    resultxi <-  rbind(resultxi, data.frame(xiphys = (tableres$xicalc[mask1] + tableres$xicalc[mask2]) * 0.5,
        dxiphys = sqrt(tableres$dxicalc[mask1]^2 + tableres$dxicalc[mask2]^2) * 0.5))
    # print("do extrapolation")
    }
}

if(FALSE){
pdf(sprintf("extrapolatecontroluwerrbeta%f.pdf", betas[1]), title = "")
list3 <-  list()
list16 <-  list()
for (i in seq(1, length(betas))) {
    print(Nt[i])
    if (Nt[i] !=   Ntextrapolate) {
    filengthame16 <-  sprintf("result2p1d.u1potential.rotated.Nt%d.Ns16.b%f.xi%f.nape0.alpha1.000000coarsedistance", Nt[i], betas[i], 16 / Nt[i])
    data16 <-  readloopfilecfonecolumn(file = filengthame16, path = "", skip = therm16[i], column = columns[i])
    data16 <-  data16$cf
    mask <-  tableres$beta == betas[i] & abs(tableres$xi - 16 / Nt[i]) < 1e-5
    } else {
    filengthame161 <-  sprintf("result2p1d.u1potential.rotated.Nt%d.Ns16.b%f.xi%f.nape0.alpha1.000000coarsedistance", Nt[i], betaextrapolate[1], 16 / Nt[i])
    data161 <-  readloopfilecfonecolumn(file = filengthame161, path = "", skip = therm16[i], column = columns[i])
    filengthame162 <-  sprintf("result2p1d.u1potential.rotated.Nt%d.Ns16.b%f.xi%f.nape0.alpha1.000000coarsedistance", Nt[i], betaextrapolate[2], 16 / Nt[i])
    data162 <-  readloopfilecfonecolumn(file = filengthame162, path = "", skip = therm16[i], column = columns[i])
    # print(length(data161$cf))
    # print(length(data162$cf))
    data16 <-  (data161$cf[1:2000, 1] + data162$cf[1:2000, 1]) * 0.5
    mask1 <-  tableres$beta == betaextrapolate[1] & abs(tableres$xi - 16 / Nt[i]) < 1e-5
    mask2 <-  tableres$beta == betaextrapolate[2] & abs(tableres$xi - 16 / Nt[i]) < 1e-5
    # print("do extrapolation")
    }
    uwerr16 <-  uwerrprimary(data16, S = S)
    list16[[i]] <-  uwerr16
    P16[i] <-  uwerr16$value
    dP16[i] <-  uwerr16$dvalue
    plot(x = seq(1, length(data16)), y = data16, main = paste("L = 16, T = ", Nt[i]))
    plot(uwerr16, main = paste("L = 16, T = ", Nt[i]))
}


print(P16)
print(resultxi)


for (i in seq(1, length(betas))) {
    print(i)
    filengthame3 <-  sprintf("result2p1d.u1potential.Nt%d.Ns3.b%f.xi%f.nape0.alpha1.000000nonplanar", Nt[i], betas[i], 16 / Nt[i])
    data3 <-  readloopfilecfonecolumn(file = filengthame3, path = "", skip = therm3[i], column = 6)
    data3 <-  data3$cf
    uwerr3 <-  uwerrprimary(data3, S = S)
    list3[[i]] <-  uwerr3
    P3[i] <-  uwerr3$value
    dP3[i] <-  uwerr3$dvalue
    plot(x = seq(1, length(data3)), y = data3, main = paste("L = 3, T = ", Nt[i]))
    plot(uwerr3, main = paste("L = 3, T = ", Nt[i]))
}

result <-  data.frame(P3 = P3, dP3 = dP3, P16 = P16, dP16 = dP16, xiphyssq = resultxi$xiphys^2, dxiphyssq = resultxi$dxiphys^2, betas = betas)
filename <-  sprintf("extrapolatebeta%fbs%dS%f.csv", betas[1], bootsamples, S)
write.table(result, filename, row.names = F, col.names = T)
write.table(resultxi, sprintf("resultxisingleextrapolatebeta%f.csv", betas[i]), row.names = F, col.names = T)
}

filename <-  sprintf("extrapolatebeta%fbs%dS%f.csv", betas[1], bootsamples, S)
result <-  read.table(filename, header = T)
P3 <-  result$P3
dP3 <-  result$dP3

P16 <-  result$P16
dP16 <-  result$dP16

bootstraps3 <-  parametric.bootstrap(bootsamples, P3, dP3, seed = 12345 + 1000 * betas[1])
bootstraps16 <-  parametric.bootstrap(bootsamples, P16, dP16, seed = 24680 + 1000 * betas[1])
bootstrapsxi <-  parametric.bootstrap(bootsamples, resultxi$xiphys, resultxi$dxiphys, seed = 67890 + 1000 * betas[1])
bootstrapsxi <-  bootstrapsxi^2
bootstrapsbeta <-  parametric.bootstrap(bootsamples, betas, rep(1e-2, length(betas)), seed = 54321 + 1000 * betas[1])
# bootstrapsbeta <- t(array(rep(betas, bootsamples), dim = c(length(betas), bootsamples)))
# print(head(bootstrapsbeta))


pdf(sprintf("extrapolateratiobeta%f.pdf", betas[1]), title = "")

fit.contlim16 <-  bootstrap.nlsfit(fnpar, par.guess = c(1, 1, 1), y = P16, x = resultxi$xiphys^2, bsamples = cbind(bootstraps16, bootstrapsxi))
# summary(fit.contlim16)
plot(fit.contlim16, plot.range = c( - 0.1, 1.3), xlab = "xi_ren^2", ylab = "P(16)", main = "anisotropy from data")

fit.contlim3 <-  bootstrap.nlsfit(fnpar, par.guess = c(1, 1, 1), y = P3, x = resultxi$xiphys^2, bsamples = cbind(bootstraps3, bootstrapsxi))
# summary(fit.contlim3)
plot(fit.contlim3, plot.range = c( - 0.1, 1.3), xlab = "xi_ren^2", ylab = "P(3)", main = "anisotropy from data at L = 16")

fit.xibeta <-  bootstrap.nlsfit(fncub, par.guess = c(1.4, 1, -0.1, 1), y = betas, x = resultxi$xiphys^2, bsamples = cbind(bootstrapsbeta, bootstrapsxi))
# summary(fit.xibeta)
plot(fit.xibeta, plot.range = c( - 0.1, 1.3), xlab = "xi_ren^2", ylab = "beta", main = "anisotropy from data")



# print(resxi$xiphys)
# print(resxi$dxiphys)


# print(P3)
# print(P16)
# print(P3 / P16)
# print(sqrt((dP3 / P16)^2 + (P3 * dP16 / P16^2)^2))
# plotwitherror(x = resxi$xiphys^2, dx = 2 * resxi$dxiphys, y = P3 / P16,
#   dy = sqrt((dP3 / P16)^2 + (P3 * dP16 / P16^2)^2), xlab = "xi_ren^2", ylab = "P(3) / P(16)", main = "anisotropy from fit")
# plotwitherror(x = resultxi$xiphys^2, dx = 2 * resultxi$dxiphys, y = P3 / P16,
#   dy = sqrt((dP3 / P16)^2 + (P3 * dP16 / P16^2)^2), xlab = "xi_ren^2", ylab = "P(3) / P(16)", main = "anisotropy from data")



fit.ratio <-  bootstrap.nlsfit(fnpar, par.guess = c(1, 1, 1), y = P3 / P16, x = resultxi$xiphys^2, bsamples = cbind(bootstraps3 / bootstraps16, bootstrapsxi))
# summary(fit.ratio)
plot(fit.ratio, plot.range = c( - 0.1, 1.3), xlab = "xi_ren^2", ylab = "P(3) / P(16)", main = "anisotropy from data")


string <-  sprintf("continuum limit at L = 3 for beta = %f:\nbeta = %s, P = %s\ncalculated from P(16) = %s and ratio = %s", betas[1],
tex.catwitherror(x = fit.xibeta$t0[1], dx = fit.xibeta$se[1], digits = 2, with.dollar = F),
tex.catwitherror(x = mean(fit.contlim16$t[, 1] * fit.ratio$t[, 1]), dx = sd(fit.contlim16$t[, 1] * fit.ratio$t[, 1]), digits = 2, with.dollar = F),
tex.catwitherror(x = fit.contlim16$t0[1], dx = fit.contlim16$se[1], digits = 2, with.dollar = F),
tex.catwitherror(x = fit.ratio$t0[1], dx = fit.ratio$se[1], digits = 2, with.dollar = F)
)
print(string)

inversepar <- function(par, x, boot.r, ...){
    return( (-par[2]/2./par[3]) - sqrt((par[2]/2./par[3])^2-(par[1]-x)/par[3]))
}

fitcubicroot <- function(par, x, boot.r, ...){
    return( par[1] + par[2]*x + par[3]*sqrt(x) + par[4]*x^(1./3.))
}



fitinverse <- bootstrap.nlsfit(fitcubicroot, par.guess = c(1, 1, 1, 1), y = resultxi$xiphys^2, x = betas, bsamples = bootstrapsxi)
print(fitinverse)
plot(fitinverse, main = "xi(beta) for L = 16", xlab = "beta", ylab = "xi_ren^2", plot.range = c(1, 2), ylim = c(-0.1, 1.1))
lines(x = seq(0, 3), y = rep(0, 4), col = "red")

# print(names(fitinverse))
# print(fitinverse$fn)

rootpar <- function(x, par, bootr = 1) fnpar(par, x, bootr)
rootcub <- function(x, par, bootr = 1) fncub(par, x, bootr)
rootinv <- function(x, par, bootr = 1) fitcubicroot(par, x, bootr)

find_root <- function(func, fit.result, interval = c(1, 2)){
    # mean value
    root <- uniroot(func, interval, par = fit.result$t0, tol = 1e-12)
    res <- list()
    res$t0 <- root$root
    # print(root)
    for(bs in seq(1, fit.result$boot.R)){
        root <- uniroot(func, interval, par = fit.result$t[bs, ], tol = 1e-12, extendInt = "yes")
        res$t[bs] <- root$root
    }
    res$se <- sd(res$t)
    return(res)
}

betalim <- find_root(rootinv, fitinverse, interval = c(1.3, 1.6))
string <- sprintf("root-finding: %f  + /- %f", betalim$t0, betalim$se)
print(string)
string <- sprintf("fit with xy: %f  + /- %f", fit.xibeta$t0[1], fit.xibeta$se[1])
print(string)


resultlist <- list(fit.contlim16 = fit.contlim16, fit.contlim3 = fit.contlim3,
        fit.xibeta = fit.xibeta, fit.ratio = fit.ratio,
        P3 = list(P3val = P3, P3err = dP3, P3boot = bootstraps3),
        betainv = betalim, fitinv = fitinverse)
saveRDS(resultlist, file = sprintf("listextrapolatebeta%f.RData", betas[1]))
write.table(resultxi, sprintf("resultxisingleextrapolatebeta%f.csv", betas[1]), row.names = F, col.names = T)


xilist <- resultxi$xiphys^2
dxilist <- 2*resultxi$xiphys*resultxi$dxiphys
xilistboot <- bootstrapsxi

