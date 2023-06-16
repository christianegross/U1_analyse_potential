library("hadron")
source("/hiskp4/gross/masterthesis/analyse/code/U1_analyse_potential/myfunctions.R")

bootsamples <- 100

betas <- c(1.65, 1.615, 1.58, 1.515, 1.48, 1.46, 1.43)
columns <- c(1, 1, 11, 1, 1, 1, 1)
betaextrapolate <- c(1.45, 1.47)
Nt <- c(16, 20, 24, 32, 40, 48, 64)
Ntextrapolate <- 48
therm16 <- c(500, 500, 500, 500, 500, 1000, 2000)
therm3 <- c(1333, 1333, 1333, 1333, 1333, 2667, 2667)

P16 <- c()
dP16 <- c()
P3 <- c()
dP3 <- c()
resxi <- read.table("plotstikz/resultsrenormalizationbeta1.650000.csv", header=TRUE)
print(resxi)
print(resxi$xiphys)
print(resxi$dxiphys)


pdf("extrapolateratio.pdf", title="")
list3 <- list()
list16 <- list()
for (i in seq(1:length(betas))) {
    if (Nt[i] != Ntextrapolate) {
    filengthame16 <- sprintf("result2p1d.u1potential.rotated.Nt%d.Ns16.b%f.xi%f.nape0.alpha1.000000coarsedistance", Nt[i], betas[i], 16/Nt[i])
    data16 <- readloopfilecfonecolumn(file = filengthame16, path = "", skip = therm16[i], column = columns[i])
    data16 <- data16$cf
    } else {
    filengthame161 <- sprintf("result2p1d.u1potential.rotated.Nt%d.Ns16.b%f.xi%f.nape0.alpha1.000000coarsedistance", Nt[i], betaextrapolate[1], 16/Nt[i])
    data161 <- readloopfilecfonecolumn(file = filengthame161, path = "", skip = therm16[i], column = columns[i])
    filengthame162 <- sprintf("result2p1d.u1potential.rotated.Nt%d.Ns16.b%f.xi%f.nape0.alpha1.000000coarsedistance", Nt[i], betaextrapolate[2], 16/Nt[i])
    data162 <- readloopfilecfonecolumn(file = filengthame162, path = "", skip = therm16[i], column = columns[i])
    print(length(data161$cf))
    print(length(data162$cf))
    data16 <- (data161$cf[1:2000, 1] + data162$cf[1:2000, 1]) * 0.5
    print("do extrapolation")
    }
    uwerr16 <- uwerrprimary(data16)
    list16[[i]] <- uwerr16
    P16[i] <- uwerr16$value
    dP16[i] <- uwerr16$dvalue
    plot(x=seq(1:length(data16)), y=data16)
    plot(uwerr16)
}

for (i in seq(1:length(betas))) {
    filengthame3 <- sprintf("result2p1d.u1potential.Nt%d.Ns3.b%f.xi%f.nape0.alpha1.000000nonplanar", Nt[i], betas[i], 16/Nt[i])
    data3 <- readloopfilecfonecolumn(file = filengthame3, path = "", skip = therm3[i], column = 6)
    data3 <- data3$cf
    uwerr3 <- uwerrprimary(data3)
    list3[[i]] <- uwerr3
    P3[i] <- uwerr3$value
    dP3[i] <- uwerr3$dvalue
    plot(x=seq(1:length(data3)), y=data3)
    plot(uwerr3)
}

bootstraps3 <- parametric.bootstrap(bootsamples, P3, dP3)
bootstraps16 <- parametric.bootstrap(bootsamples, P16, dP16)
# print(head(bootstraps3))


# pdf("extrapolateratio.pdf", title="")
# print(resxi$xiphys)
# print(resxi$dxiphys)
print(P3)
print(P16)
print(P3 / P16)
print(sqrt((dP3 / P16)^2 + (P3 * dP16 / P16^2)^2))
plotwitherror(x = resxi$xiphys, dx = resxi$dxiphys, y = P3 / P16,
  dy = sqrt((dP3 / P16)^2 + (P3 * dP16 / P16^2)^2))
