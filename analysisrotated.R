library(hadron)
library(optparse)

if(TRUE){
option_list <- list(
    make_option(c("-b", "--beta"), type="double", default=1.7, 
    help="beta-Parameter of simulation [default %default]"),
    make_option(c("--skip"), type="integer", default=1000, 
    help="how many lines are skipped when reading in [default %default]"),
     
    make_option(c("-s", "--bootsamples"), type="integer", default=500, 
    help="how many bootstrapsamples should be drawn [default %default]"), 
    make_option(c("-r", "--Ns"), type="integer", default=16, 
    help="Ns of lattice [default %default]"),
    
    make_option(c("-t", "--Nt"), type="integer", default=16, 
    help="Nt of lattice [default %default]"),
    make_option(c("-n", "--nape"), type="integer", default=0, 
    help="Number of APE-smears that were done [default %default]"), 
    
    make_option(c("-l", "--lowerboundmeff"), type="integer", default=0, 
    help="meff is determined in region [lowerboundmeff, Ns/2-2] [default %default]"), 
    make_option(c("-j", "--job"), type="integer", default=0, 
    help="qbig job number [default %default]"), 
    
    make_option(c("-o", "--omit"), type="integer", default=0, 
    help="how many points should be omitted in calculating the potentials? [default %default]"), 
    make_option(c("-a", "--alpha"), type="double", default=1.0, 
    help="alpha used in APE-smearing [default %default]"), 
    
    make_option(c("--betaone"), type="double", default=0, 
    help="input beta at corresponding xi=1 [default %default]"), 
    make_option(c("-x", "--xi"), type="double", default=0, 
    help="xi used in lattice [default %default]"),   
    make_option(c("--zerooffset"), type="integer", default=0, 
    help="Offset for selecting the correct loops if W(x=0, y=0) was measured [default %default]"),
    
    make_option(c("-e", "--every"), type="integer", default=1, 
    help="only reads every eth line [default %default]"), 
    make_option(c("--bootl"), type="integer", default=2, 
    help="block length in bootstrapping configurations [default %default]"), 
    
    make_option(c("--nsave"), type="integer", default=0, 
    help="steps between saved configs [default %default]"),  
    
    make_option(c("--analyse"), action="store_true", default=FALSE, 
    help="if true, correlators and effective masses are determined from provided table [default %default]"),
    make_option(c("--determinemeff"), action="store_true", default=FALSE, 
    help="if true, correlators and effective masses are determined for given t1 [default %default]"),
    
    make_option(c("--dofit"), action="store_true", default=FALSE, 
    help="if true, potentials and xi are calculated [default %default]"),
    make_option(c("--smearing"), action="store_true", default=FALSE, 
    help="are the calculations done based on smeared lattices? (changes filenames) [default %default]"),
    
    make_option(c("--xidiff"), action="store_true", default=FALSE, 
    help="Is xi different to Ns/Nt? [default %default]"),
    make_option(c("--plotonlymeff"), action="store_true", default=FALSE, 
    help="The plots of the fitted meff are plotted seperately [default %default]"),
    
    make_option(c("--respath"), type="character", default="", 
    help="path to where the resultfiles are stored [default %default]"),
    make_option(c("--plotpath"), type="character", default="", 
    help="path to where the plots are stored [default %default]"),
    
    make_option(c("--myfunctions"), type="character", default="/hiskp4/gross/masterthesis/su2/build/debug/analysisscripts/myfunctions.R", 
#~     make_option(c("--myfunctions"), type="character", default="myfunctions.R", 
    help="path to where additional functions are stored [default %default]")
)
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments=0)
opt <- args$options
}

#~ print(opt)

if(TRUE){
source(opt$myfunctions)
beta <- opt$beta
skip <- opt$skip
bootsamples <- opt$bootsamples
Ns <- opt$Ns
Nt <- opt$Nt
xi <- Ns/Nt
if(opt$xidiff){
    xi <- opt$xi
}

nape <- opt$nape
alpha <- opt$alpha

potential <- data.frame(R=NA, m=NA, dm=NA, space=NA, p=NA, chi=NA)

#determining the effective masses takes the most times, so do not do this if it has already been done
#fit also does not have to be done every time
analyse <- opt$analyse
dofit <- opt$dofit

t2 <- Ns/2-2
t1 <- opt$lowerboundmeff
}

if(analyse){
filenameforplots <- sprintf("%splotseffectivemass2p1dmeffchosenNs%dbeta%fxi%f.pdf", opt$plotpath, Ns, beta, xi)
if(opt$smearing){
    filenameforplots <- sprintf("%splotseffectivemass2p1dmeffchosenNs%dbeta%fxi%fnape%dalpha%f.pdf", opt$plotpath, Ns, beta, xi, nape, alpha)
}
pdf(file=filenameforplots, title="")
listfits <- list()
listtauint <- list()
namelistbounds <- sprintf("%stableanalysisrotatedNs%dNt%dbeta%fxi%f.csv", opt$respath, Ns, Nt, beta, xi)
if(opt$smearing){
    namelistbounds <- sprintf("%stableanalysisrotatedNs%dNt%dbeta%fxi%fnape%dalpha%f.csv", opt$respath, Ns, Nt, beta, xi, nape, alpha)
}
if (file.exists(namelistbounds)){    
    listbounds <- read.table(namelistbounds, sep=",", header=TRUE)
} else{
    listbounds <- cbind(data.frame(yt=c(seq(1, Ns/2), seq(1, Nt/2))), data.frame(spacial=c(rep(TRUE, Ns/2), rep(FALSE, Nt/2))), data.frame(lower=rep(opt$lowerboundmeff, Ns/2+Nt/2)), data.frame(upper=rep(Ns/2-2, Ns/2+Nt/2)))
}
print(listbounds)
negatives <- c()
#determine Wilson Loop correlators and effective masses for all possible combinations Meff(x0, y0) (use W(x0, t, y0)) 
#and Meff(x0, t0) (use W(x0, t0, y)) 
#try method as in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.63.074501
#x here as z there, y here as x there
for (y in seq(1,Ns/2,1)){
    #save results, list has to have an entry everywhere
    listresults <- list(FALSE, FALSE, FALSE, FALSE)
    uwerrresults <- FALSE
    if(y>0){
    
#~     message("\ny=", y)
    print(y)
#~     filename <- sprintf("resultsqbig/meff2p1d/result2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fcoarsedistance",
#~                         Nt, Ns, beta, xi, nape, alpha)
#~     if(onqbig){
    filename <- sprintf("result2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fcoarsedistance",Nt, Ns, beta, xi, nape, alpha)
#~     }
    
    title<-sprintf("beta=%.2f coarse y=%d skipped %d", beta, y, skip*opt$nsave)
    WL <- calcplotWloopsideways(file=filename, path=opt$respath, skip=skip, Ns=Ns, yt=y, title=title, bootsamples=bootsamples,
    zerooffset=opt$zerooffset, nsave=opt$nsave, every=opt$every, l=opt$bootl)
    
    uwerrresults <- uwerr.cfboot(WL)
#~     print(uwerrresults)
    
    #determine effective masses, fit one effective mass, therefore set bordes. borders may have to be adjusted for anisotropic lattices
    t1 <- listbounds$lower[listbounds$spacial==TRUE&listbounds$yt==y]
    t2 <- listbounds$upper[listbounds$spacial==TRUE&listbounds$yt==y]
    WL.effmasslist <- deteffmass(WL=WL, yt=y, potential=potential, t1=t1, t2=t2, isspatial=1)
    
    #plot results, save results off meff in short and long form (with and without bootstrapsamples
#~     message(WL.effmasslist[[2]][[1]])
    if(WL.effmasslist[[2]][[2]]!=0){
        listresults <- list(WL.effmasslist[[1]], y, FALSE, uwerrresults)
    }
    #fine=!spacial
    names(listresults) <- c("effmass", "yt", "fine", "uwerr")
    
    title <- sprintf("%s, %d configs\n", title, (length(WL$cf[,1])+skip)*opt$nsave)
    try(plot(WL.effmasslist[[1]], xlab="x/a_s", ylab="Meff", main=sprintf("%s, t1=%d, t2=%d, p=%f", title, t1, t2, WL.effmasslist[[1]]$effmassfit$Qval)))
    try(plot(WL.effmasslist[[1]], xlab="x/a_s", ylab="Meff", main=sprintf("%s, t1=%d, t2=%d, p=%f, zoomed", title, t1, t2, WL.effmasslist[[1]]$effmassfit$Qval), ylim=WL.effmasslist[[2]]))
}
    listfits[[y]] <- listresults
    listtauint[[y]] <- uwerrresults
    negatives[y] <- sum(WL$cf0 < 0)
}
#same procedure as for y loop
for (t in seq(1,Nt/2,1)){
    listresults <- list(FALSE, FALSE, FALSE, FALSE)
    uwerrresults <- FALSE
    if(t<Nt/2+1){
        
#~     message("\nt=", t)
    print(t)
    filename <- sprintf("result2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%ffinedistance",Nt, Ns, beta, xi, nape, alpha)
    
    
    title<-sprintf("beta=%.2f fine t=%d skipped %d", beta, t, skip*opt$nsave)
    WL <- calcplotWloopsideways(file=filename, path=opt$respath, skip=skip, Ns=Ns, yt=t, title=title, bootsamples=bootsamples,
    zerooffset=opt$zerooffset, nsave=opt$nsave, every=opt$every, l=opt$bootl)
    
    uwerrresults <- uwerr.cfboot(WL)
#~     print(uwerrresults)
    
    t1 <- listbounds$lower[listbounds$spacial==FALSE&listbounds$yt==t]
    t2 <- listbounds$upper[listbounds$spacial==FALSE&listbounds$yt==t]
    WL.effmasslist <- deteffmass(WL=WL, yt=t, potential=potential, t1=t1, t2=t2, isspatial=0)
#~     message(WL.effmasslist[[2]][[1]])
    if(WL.effmasslist[[2]][[2]]!=0){
        listresults <- list(WL.effmasslist[[1]], t, TRUE, uwerrresults)
    }
    names(listresults) <- c("effmass", "yt", "fine", "uwerr")

    title <- sprintf("%s, %d configs\n", title, (length(WL$cf[,1])+skip)*opt$nsave)
    try(plot(WL.effmasslist[[1]], xlab="x/a_s", ylab="Meff", main=sprintf("%s, t1=%d, t2=%d, p=%f", title, t1, t2, WL.effmasslist[[1]]$effmassfit$Qval)))
    try(plot(WL.effmasslist[[1]], xlab="x/a_s", ylab="Meff", main=sprintf("%s, t1=%d, t2=%d, p=%f", title, t1, t2, WL.effmasslist[[1]]$effmassfit$Qval), ylim=WL.effmasslist[[2]]))
}
    listfits[[Ns/2+t]] <- listresults
    listtauint[[Ns/2+t]] <- uwerrresults
    negatives[Ns/2+t] <- sum(WL$cf0 < 0)
}
#write out results
potential <- na.omit(potential)
filenamepotential <- sprintf("%spotentialmeff2p1dmeffchosenNs%dbeta%fxi%fbsamples%d.csv", opt$plotpath, Ns, beta, xi, bootsamples)
filenamelist <- sprintf("%slistmeff2p1dmeffchosenNs%dbeta%fxi%fbsamples%d.RData", opt$plotpath, Ns, beta, xi, bootsamples)
filenameuwerr <- sprintf("%slistuwerrtauintmeffchosenNs%dbeta%fxi%fbsamples%d.RData", opt$plotpath, Ns, beta, xi, bootsamples)
    filenamenegatives <- sprintf("%snegativessidewaysNs%dbeta%fxi%fbsamples%d.csv", opt$plotpath, Ns, beta, xi, bootsamples)
#~ filenamelistwl <- sprintf("listwl2p1drotatedNs%dbeta%fxi%fbsamples%d.RData", Ns, beta, xi, bootsamples)

if(opt$smearing){
    filenamepotential <- sprintf("%spotentialmeff2p1dmeffchosenNs%dbeta%fxi%fnape%dalpha%fbsamples%d.csv", opt$plotpath, Ns, beta, xi, nape, alpha, bootsamples)
    filenamelist <- sprintf("%slistmeff2p1dmeffchosenNs%dbeta%fxi%fnape%dalpha%fbsamples%d.RData", opt$plotpath, Ns, beta, xi, nape, alpha, bootsamples)
    filenameuwerr <- sprintf("%slistuwerrtauintmeffchosenNs%dbeta%fxi%fnape%dalpha%fbsamples%d.RData", opt$plotpath, Ns, beta, xi, nape, alpha, bootsamples)
    filenamenegatives <- sprintf("%snegativessisewaysNs%dbeta%fxi%fnape%dalpha%fbsamples%d.csv", opt$plotpath, Ns, beta, xi, nape, alpha, bootsamples)
#~     filenamelistwf <- sprintf("listwl2p1drotatedNs%dbeta%fxi%fnape%dalpha%fbsamples%d.RData", Ns, beta, xi, nape, alpha, bootsamples)
}
write.table(potential, filenamepotential, row.names=FALSE)
#~ print(listfits)
saveRDS(listfits, file=filenamelist)
saveRDS(listtauint, file=filenameuwerr)
write.table(data.frame(x=c(seq(1, Ns/2), seq(1, Nt/2)), neg=negatives), file=filenamenegatives, row.names=FALSE)
#~ saveRDS(listallWL, file=filenamelistwl)
#~ print(listfits[[1]][[1]])
}

if(dofit){
filenamepotential <- sprintf("%spotentialmeff2p1dmeffchosenNs%dbeta%fxi%fbsamples%d.csv", opt$plotpath, Ns, beta, xi, bootsamples)
filenamelist <- sprintf("%slistmeff2p1dmeffchosenNs%dbeta%fxi%fbsamples%d.RData", opt$plotpath, Ns, beta, xi, bootsamples)

if(opt$smearing){
    filenamepotential <- sprintf("%spotentialmeff2p1dmeffchosenNs%dbeta%fxi%fnape%dalpha%fbsamples%d.csv", opt$plotpath, Ns, beta, xi, nape, alpha, bootsamples)
    filenamelist <- sprintf("%slistmeff2p1dmeffchosenNs%dbeta%fxi%fnape%dalpha%fbsamples%d.RData", opt$plotpath, Ns, beta, xi, nape, alpha, bootsamples)
}
potential <- read.table(filenamepotential, header=TRUE)
if(opt$determinemeff){
    t1 <- opt$lowerboundmeff
} else{
    t1 <- -1
}


listmeff <- readRDS(file=filenamelist)
#~ bootsamples <- length(listmeff[[1]][[1]]$massfit.tsboot[,1])
message("bootsamples= ", bootsamples)
#~ message("read in files")



filenameforplots <- sprintf("%spotentialmeff2p1dmeffchosenNs%dbeta%fxi%fomit%d.pdf", opt$plotpath, Ns, beta, xi, opt$omit)
if(opt$smearing){
filenameforplots <- sprintf("%spotentialmeff2p1dmeffchosenNs%dbeta%fxi%fomit%dnape%dalpha%f.pdf", opt$plotpath, Ns, beta, xi, opt$omit, nape, alpha)
}
pdf(file=filenameforplots, title="")

#define functions for potential, force, force = -r^2*derivation of potential
fnpot <- function (par, x, boot.r, ...) par[1] + par[2] * x + par[3] * log(x)
fnforce <- function (par, x, boot.r, ...) (- 1.0) * par[2] * x^2 - par[3] * x
fnforceerr <- function (parerr, x, boot.r, correlation, ...) sqrt( parerr[2]^2 * x^4 + parerr[3]^2 * x^2 + 2 * x^3 * correlation[2,3] * parerr[2] * parerr[3] )


potall <- function(par, x, boot.r, maskpot, ...){
    x[mask] <- par[4]*x[mask]
    ret <- par[1] + par[2] * x + par[3] * log(x)
    return (invisible(ret))
}

#~ message("defined functions")


#coarse potential
#read in data points, bootstrap samples: have to initialise empty vectors, because not every list element has a valid entry
xc <- rep(NA, Ns/2)
yc <- rep(NA, Ns/2)

#filter out points without valid entry in bootstrap.nlsfit, mask for combining vectors
maskc <- rep(FALSE, Ns/2)
bsamplesc <- array(c(rep(NA, Ns/2*bootsamples)), dim=c(bootsamples, Ns/2))

#also have vectors for combined potential
x <- rep(NA, Ns/2+Nt/2)
dx <- rep(NA, Ns/2+Nt/2)
y <- rep(NA, Ns/2+Nt/2)
mask <- rep(FALSE, Ns/2+Nt/2)
finemask <- rep(FALSE, Ns/2+Nt/2)
bsamples <- array(c(rep(NA, 2*(Ns/2+Nt/2)*bootsamples)), dim=c(bootsamples, 2*(Ns/2+Nt/2)))
#~ message("defined empty arrays")

for (i in seq(1, Ns/2-opt$omit, 1)){
  if(listmeff[[i]][[2]]){
      xc[i] <- listmeff[[i]][[2]]
#~       message("xc= ",xc[i])
      yc[i] <- listmeff[[i]][[1]]$effmassfit$t0[1]
      bsamplesc[, i] <- listmeff[[i]][[1]]$massfit.tsboot[,1]
      maskc[i] <- TRUE
      
      x[i] <- listmeff[[i]][[2]]
      y[i] <- listmeff[[i]][[1]]$effmassfit$t0[1]
      bsamples[, i] <- listmeff[[i]][[1]]$massfit.tsboot[,1]
      mask[i] <- TRUE
      finemask[i] <- FALSE
  }
}
xall <- x
yall <- y
bsamplesall <- bsamples[, 1:Ns/2+Nt/2]



#determine parameters by bootstrap
fit.resultcoarse <- bootstrap.nlsfit(fnpot, c(0.2, 0.2, 0.2), yc, xc, bsamplesc, mask=maskc)
filenamecoarse <- sprintf("%sfitresultcoarseb%fxi%fNs%dbsamples%domit%dl%d.RData", opt$plotpath, beta, xi, Ns, bootsamples, opt$omit, t1)
if(opt$smearing){
filenamecoarse <- sprintf("%sfitresultcoarseb%fxi%fNs%dnape%dalpha%fbsamples%domit%dl%d.RData", opt$plotpath, beta, xi, Ns, nape, alpha, bootsamples, opt$omit, t1)
}
saveRDS(fit.resultcoarse, file=filenamecoarse)

#fine potential
#read in data points, bootstrap samples: have to initialise empty vectors, because not every list element has a valid entry
xf <- rep(NA, Nt/2)
yf <- rep(NA, Nt/2)

#filter out point without valid entry in bootstrap
maskf <- rep(FALSE, Nt/2)
bsamplesf <- array(c(rep(NA, Nt/2*bootsamples)), dim=c(bootsamples, Nt/2))
for (i in seq(1, Nt/2-opt$omit/xi, 1)){
  if(listmeff[[Ns/2+i]][[2]]){
      xf[i]     <- listmeff[[Ns/2+i]][[2]]
#~       message("xf= ",xf[i])
      yf[i]     <- listmeff[[Ns/2+i]][[1]]$effmassfit$t0[1]
      bsamplesf[, i]        <- listmeff[[Ns/2+i]][[1]]$massfit.tsboot[,1]
      maskf[i]      <- TRUE
      
      x[Ns/2+i] <- listmeff[[Ns/2+i]][[2]]
      y[Ns/2+i] <- listmeff[[Ns/2+i]][[1]]$effmassfit$t0[1]
      bsamples[, Ns/2+i]    <- listmeff[[Ns/2+i]][[1]]$massfit.tsboot[,1]
      mask[Ns/2+i]  <- TRUE
      finemask[Ns/2+i] <- TRUE
    
  }
}

#head(bsamples)
#~ warnings() 

#determine parameters by bootstrap
fit.resultfine <- bootstrap.nlsfit(fnpot, c(0.2, 0.2, 0.2), yf, xf, bsamplesf, mask=maskf)
filenamefine <- sprintf("%sfitresultfineb%fxi%fNs%dbsamples%domit%dl%d.RData", opt$plotpath, beta, xi, Ns, bootsamples, opt$omit, t1)
if(opt$smearing){
filenamefine <- sprintf("%sfitresultfineb%fxi%fNs%dnape%dalpha%fbsamples%domit%dl%d.RData", opt$plotpath, beta, xi, Ns, nape, alpha, bootsamples, opt$omit, t1)
}
saveRDS(fit.resultfine, file=filenamefine)

#~ message("starting to calculate xi")
singlexis <- data.frame(xi=NA)
for(bs in seq(1, bootsamples)){
    val <- fit.resultfine$t[bs, 2]/fit.resultcoarse$t[bs, 2]
#~     message(val)
    singlexis <- rbind(singlexis, data.frame(xi=val))
}
singlexis <- na.omit(singlexis)
xisingle <- mean(singlexis$xi)
dxisingle <- sd(singlexis$xi)
strings <- sprintf("xi from single potentials: %f +/- %f", xisingle, dxisingle)
#~ print(strings)
#~ message("single xi calculated")

#calculate plaquette as W(x=1, y=1, t=0)

filename <- sprintf("%sresult2p1d.u1potential.rotated.Nt%d.Ns%d.b%f.xi%f.nape%d.alpha%fcoarsedistance", opt$respath, Nt, Ns, beta, xi, nape, alpha)

alldata <- read.table(filename)
plaquettecolumn <- alldata[, (opt$zerooffset+Ns/2)*opt$zerooffset+1+opt$zerooffset]
plaquettedata <- uwerrprimary(plaquettecolumn[seq(skip+1, length(plaquettecolumn), opt$every)])
newline <- data.frame(value=plaquettedata$value, dvalue=plaquettedata$dvalue, ddvalue=plaquettedata$ddvalue, 
                      tauint=plaquettedata$tauint, dtauint=plaquettedata$dtauint)
nom <- floor(length(plaquettecolumn)/opt$every)
#print("Value for plaquette")
#print(newline)
column <- (opt$zerooffset+Ns/2)*opt$zerooffset+1+opt$zerooffset

plaquettecf <- readloopfilecfonecolumn(file=filename, skip=skip, column=column, every=opt$every)
#~ summary(plaquettecf)
plaquettecf <- bootstrap.cf(plaquettecf, boot.R=bootsamples, boot.l=opt$bootl)
#~ summary(plaquettecf)


uwerrresults <- uwerr.cf(plaquettecf)
#~ print(uwerrresults)
#~ print(uwerrresults$uwcf$tauint[1])



try(plot(fit.resultfine, main = sprintf("beta=%.3f, xi=%.3f, (2+1)D, Ns=%d, fine\n%d measurements of which skipped %d", beta, xi, Ns, nom*100, skip*100), ylab="a_s V(t)", xlab="t/a_t"))
try(plot(fit.resultcoarse, main = sprintf("beta=%.3f, xi=%.3f, (2+1)D, Ns=%d, coarse\n%d measurements of which skipped %d", beta, xi, Ns, nom*100, skip*100), ylab="a_s V(y)", xlab="y/a_s"))


# determine xi by interpolating between different V(t) to get V(x)=V(t xi), 
#like in https://journals.aps.org/prd/abstract/10.1103/PhysRevD.63.074501

xis <- data.frame(xi=NA, t=NA, y=NA, j=NA)
xibootsamples <- c()
#Determine parameters of interpolation for V_t, V_t(t)=a*t+b, defined by two points, lower and upper, 
#with distance 1 a_t, V(i)=lower, V(i+1)=upper
#calculate the xi for all Vs for all bootstrapsamples, overall xi = mean pm sd
for (bs in seq(1, bootsamples, 1)){
#~     #message("bs=", bs)
interpolation <- data.frame(lower=NA, upper=NA, a=NA, b=NA, bootsample=NA, lowerx=NA)
for (i in seq(1, Nt/2-1-opt$omit/xi, 1)){
    if(listmeff[[Ns/2+i+1]][[2]]){
    lower     <- listmeff[[Ns/2+i]][[1]]$massfit.tsboot[bs,1]
    upper     <- listmeff[[Ns/2+i+1]][[1]]$massfit.tsboot[bs,1]
    a <- upper-lower
#~     #message("i=", i, " a=", a)
    b <- lower - i * a
    newline <- data.frame(lower=lower, upper=upper, a=a, b=b, bootsample=bs, lowerx=i)
    interpolation <- rbind(interpolation, newline)
}
}

interpolation <- interpolation[-1,]
#~ #print(interpolation)

#solve equation V_s(x)=V_t(xi t)=a*t+b -> t=(V_s-b)/a -> xir=x/t=xiresult for all possible V_s with x >= 2 (variable name xi is already taken)
#potential not dominated by linear part for y>2
xiboots <- c()
for ( i in seq(2, Ns/2-opt$omit, 1)){
  if(listmeff[[i]][[2]]){
    Vs <- listmeff[[i]][[1]]$massfit.tsboot[bs,1]
    index <- NA
#~     #print(Vs)
#determine in which interval of the V_t V_s is sitting by selecting lowest possible interval for which the upper limit is bigger than V_s
    for (j in seq(1, Nt/2-1-opt$omit/xi)){
#~         print(interpolation$lowerx==j&interpolation$bootsample==bs)
        if (length(interpolation$upper[interpolation$lowerx==j&interpolation$bootsample==bs])>0){
        if (Vs < interpolation$upper[interpolation$lowerx==j&interpolation$bootsample==bs]){
            index <- interpolation$lowerx==j&interpolation$bootsample==bs
#~             #message(index)
            break
        }
        }
    }
    #highest V_s could be bigger than highest V_t, so there is not an index for every y 
    if( !is.na(index) ){
        t <- (Vs-interpolation$b[index])/interpolation$a[index]
        xir <- i/t
        newline <- data.frame(xi=xir, t=t, y=i, j=interpolation$lowerx[index])
        xis <- rbind(xis, newline)
#~         xiboots[i-1] <- xir
        xiboots <- append(xiboots, xir)
    }
  }
}
#~ print(xiboots)
xibootsamples[bs] <- mean(xiboots)
}
#~ print(length(xibootsamples))
#~ write.table(xis, sprintf("weirdxisnape%dalpha%f.csv", nape, alpha), row.names=FALSE)
#~ write.table(interpolation, sprintf("weirdxisnape%dalpha%f.csv", nape, alpha), row.names=FALSE, append=TRUE)

#do matching directly, without interpolation
#~ findmatch <- function(x, par, V, boot.r, ...) fnpot(par, x, boot.r, ...) - V

#~ for(bs in seq(1, bootsamples)){
#~     for(i in seq(1, Nt/2-opt$omit/xi)){
#~         V <- listmeff[[Ns/2+i]][[1]]$massfit.tsboot[bs,1]
#~         root <- uniroot(findmatch, interval=c(0.01, Nt), par=fit.resultcoarse$t[bs,], V=V, boot.r=0)
#~         t <- root$root
#~         xir <- t/i
#~         newline <- data.frame(xi=xir, t=t, y=i)
#~         xis <- rbind(xis, newline)
#~     }
#~ }
    
    

#calculate xi as mean
xis <- xis[-1, ]
xicalc <- mean(xis$xi) 
dxicalc <-  sd(xis$xi)
xisquared <- mean(xis$xi^2)
dxisquared <- sd(xis$xi^2)
#~ xibootsamples <- na.omit(xibootsamples)
#~ print(xibootsamples)
strings <- sprintf("%e +/- %e %f", xicalc-mean(xibootsamples), dxicalc-sd(xibootsamples), mean(xibootsamples))
#~ message("deviation different means: ", strings)
#~ message(xicalc)
message(warnings())

}   

if(dofit){
#~ xicalc <- fit.resultfine$t0[2] / fit.resultcoarse$t0[2]
#~ dxicalc <- sqrt( (fit.resultfine$se[2] / fit.resultcoarse$t0[2])^2 + 
#~ (fit.resultfine$t0[2] * fit.resultcoarse$se[2] / fit.resultcoarse$t0[2]^2)^2)

#rescale t so potential, r0 can be determined, all r are given in units of a_s


x[finemask] <- x[finemask] * xicalc
dx[finemask] <- x[finemask] * dxicalc
#For the coarse potential, no rescaling is necessary, but an error has to be given to the fit-function. could error be zero?
dx[!finemask] <- 1e-8


bsamplesx <- parametric.bootstrap(bootsamples, c(x[!finemask]), c(dx[!finemask]))
#~ bsamplesx <- parametric.bootstrap(bootsamples, c(x), c(dx))
#~ head(bsamplesx)

#join bootstrapsamples for x, y together
#~ for (i in seq(1, length(x))){
for (i in seq(1, Ns/2-opt$omit)){
    bsamples[, i+Ns/2+Nt/2] <- bsamplesx[, i]
}
for (i in seq(1, Nt/2-opt$omit/xi)){
    bsamples[, i+Ns+Nt/2] <- i * array(xibootsamples, dim=c(bootsamples, 1))
}
#~ print(x[finemask])
#~ print(head(bsamples))


title <- sprintf("beta=%f, xi=%f+/-%f, (2+1)D, Ns=%d\n%d measurements of which skipped %d", beta, xicalc, dxicalc, Ns, nom*100, skip*100)
#~ cat(title) 


fit.resultscaled <- bootstrap.nlsfit(fnpot, c(0.1, 0.13, 0.05), y, x, bsamples, mask=mask)
#~ summary(fit.resultscaled)
#~ print(mask)
#~ print(names(fit.resultscaled))
#~ print(mean(fit.resultscaled$t[,1]))
#~ print(mean(fit.resultscaled$t[,2]))
#~ print(mean(fit.resultscaled$t[,3]))
#~ print(fit.resultscaled$t0)

#~ plot(fit.resultscaled)
filenamescaled <- sprintf("%sfitresultscaledb%fxi%fNs%dbsamples%domit%dl%d.RData", opt$plotpath, beta, xi, Ns, bootsamples, opt$omit, t1)
if(opt$smearing){
filenamescaled <- sprintf("%sfitresultscaledb%fxi%fNs%dnape%dalpha%fbsamples%domit%dl%d.RData", opt$plotpath, beta, xi, Ns, nape, alpha, bootsamples, opt$omit, t1)
}
saveRDS(fit.resultscaled, file=filenamescaled)

strings <- sprintf("\nxi = %f +/- %f\neta= %f +/- %f\n\n", xicalc, dxicalc, xicalc/xi, dxicalc/xi)
#~ cat(strings)

#~ filenamescaled <- sprintf("fitresultscaledb%fxi%fNs%dbsamples%d.RData", beta, xi, Ns, bootsamples)
#~ if(opt$smearing){
#~ filenamescaled <- sprintf("fitresultscaledb%fxi%fNs%dnape%dalpha%fbsamples%d.RData", beta, xi, Ns, nape, alpha, bootsamples)
#~ }

#if wanted determine correlation matrix seperately, as determined by the summary.bootstrapfit function
if(!is.null(fit.resultscaled$cov)) {
    cov.to.cor <- diag(1 / fit.resultscaled$se)
    correlation <- cov.to.cor %*% fit.resultscaled$cov %*% cov.to.cor
}else if(!is.null(fit.resultscaled$t)) {
    correlation <- cor(fit.resultscaled$t, fit.resultscaled$t, use="na.or.complete")
}

#~ print(correlation)

#graphically represent force: for each bootstrap sample of paramters, determine value of force on a sequence, 
#determine error as standard deviation of values on each point of the sequence
xx <- seq(-0.2, Ns/2+0.2, 0.05)
bootsforce <- data.frame(matrix(rep(NA,length(xx)), ncol=length(xx), nrow=1))

for (i in seq(1, bootsamples, 1)){
    bootsforce <- rbind(bootsforce, fnforce(fit.resultscaled$t[i,], xx, 0))
}
bootsforce <- na.omit(bootsforce)
#~ print(head(bootsforce))

forceerrs <- c()
for (i in seq(1, length(xx), 1)){
    forceerrs[i] <- sd(bootsforce[,i])
}

#~ print(forceerrs)

#~ c <- opt$c

#Determine r0 as solution of equation -r^2 d/dr V(r)=c, solve for each bootstrapsample, r0=mean pm sd
# V(r)=a + sigma*r + b*ln(r)
# -d/dr V(r) = -sigma -b/r
# c = -sigma*r^2 - b*r, use p-q-formula for p=b/sigma, q=c/sigma
rzeroofc <- data.frame(r0=NA, dr0=NA, c=NA)

for (c in list(-1.65)){
#~ for (c in seq(-1.5, 0, 0.15)){
rzerolist <- data.frame(r1=NA, r2=NA)
for (bs in seq(1, bootsamples)){
#~     a <- fit.resultscaled$t[bs, 1]
    sigma <- fit.resultscaled$t[bs, 2]
    b <- fit.resultscaled$t[bs, 3]
    r1 <- -b/(2.0*sigma) + sqrt( (b/(2.0*sigma))^2 - c/sigma)
    r2 <- -b/(2.0*sigma) - sqrt( (b/(2.0*sigma))^2 - c/sigma)
    rzerolist <- rbind(rzerolist , data.frame(r1=r1, r2=r2))
}
rzerolist <- rzerolist[-1,]

#~ print(head(cbind(rzerolist, data.frame(sigma=fit.resultscaled$t[,2], b=fit.resultscaled$t[,3]))))
rzero <- mean(rzerolist$r1)
drzero <- sd(rzerolist$r1)
rzeroofc <- rbind(rzeroofc, data.frame(r0=rzero, dr0=drzero, c=c))
bsrzero <- rzerolist[, 1]
}

#~ print(bsrzero)
#~ print(length(bsrzero))

#~ print(rzeroofc)
rzeroofc <- rzeroofc[-1,]
#~ rzeroofc <- na.omit(rzeroofc)
#~ print(rzeroofc)

#plot results
title <- sprintf("beta=%.3f, xi=%.3f+/-%.3f, (2+1)D, Ns=%d\n%d measurements of which skipped %d", beta, xicalc, dxicalc, Ns, nom*opt$nsave, skip*opt$nsave)

try(plot(fit.resultscaled, main = title, ylab="V(r)", xlab="r/a_s"))

#force with error
try(plot(xx, fnforce(fit.resultscaled$t0, xx, 0), type="l", xlab="r/a_s", ylab="-r^2 d/dr V(r)", main = title, xlim=c(0,8)))
#make outline of errors, code from bootstrap.nlsfit
polyval <- c(fnforce(fit.resultscaled$t0, xx, 0)+forceerrs, rev(fnforce(fit.resultscaled$t0, xx, 0)-forceerrs))
    pcol <- col2rgb("gray", alpha=TRUE)/255 
    pcol[4] <- 0.65
    pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])
try(polygon(x=c(xx, rev(xx)), y=polyval, col=pcol, lty=0, lwd=0.001, border=pcol))
#~ lines(x=xx, y=fnforce(fit.resultscaled$t0, xx, 0)+fnforceerr(fit.resultscaled$se, xx, 0, correlation), col=2)
#~ lines(x=xx, y=fnforce(fit.resultscaled$t0, xx, 0)-fnforceerr(fit.resultscaled$se, xx, 0, correlation), col=2)

#force with error and r0, zoomed in

printarrows <- function(x) return (x)

xlim <- c(0.6*min(rzeroofc$r0), 1.4*max(rzeroofc$r0))
ylim <- c(0.6*min(abs(rzeroofc$c)), 1.4*max(abs(rzeroofc$c)))
ylim <- rev(-ylim)
lowylim <- -2*max(abs(rzeroofc$c))
#~ print(xlim)
#~ print(ylim)
try(plot(xx, fnforce(fit.resultscaled$t0, xx, 0), type="l", xlab="r/a_s", ylab="-r^2 d/dr V(r)", main = title, xlim=xlim, ylim=ylim))
try(polygon(x=c(xx, rev(xx)), y=polyval, col=pcol, lty=0, lwd=0.001, border=pcol))


for (i in seq(1,length(rzeroofc$c))){
    c <- rzeroofc$c[i]
    rzero <- rzeroofc$r0[i]
    drzero <- rzeroofc$dr0[i]
#~ lines(x=xx, y=fnforce(fit.resultscaled$t0, xx, 0)+fnforceerr(fit.resultscaled$se, xx, 0, correlation), col=2)
#~ lines(x=xx, y=fnforce(fit.resultscaled$t0, xx, 0)-fnforceerr(fit.resultscaled$se, xx, 0, correlation), col=2)

#also make lines for position of r0, polygon for error of r0
#arrows(xhi, yhi, xlo, ylo)
try(arrows(rzero, fnforce(fit.resultscaled$t0, rzero, 0), rzero, lowylim, angle=90, length=0.1, code=0))
try(arrows(rzero, c, 0, c , angle=90, length=0.1, code=0))
r0polyval1 <- c(fnforce(fit.resultscaled$t0, rzero-drzero, 0)-fnforceerr(fit.resultscaled$se ,rzero-drzero, 0, correlation ),
                fnforce(fit.resultscaled$t0, rzero+drzero, 0)+fnforceerr(fit.resultscaled$se ,rzero+drzero, 0, correlation ))
r0polyval1 <- c(r0polyval1, lowylim, lowylim)
r0polyval2 <- c(rzero-drzero, rzero+drzero)
r0polyval2 <- c(r0polyval2, rev(r0polyval2))
try(polygon(x=r0polyval2, y=r0polyval1 , col=pcol, lty=0, lwd=0.001, border=pcol))
#~ axis(3)
#~ lines(xx, fnforce(fit.resultscaled$t0, xx, 0))
}

#plot visible combination of temporal and spatial potential
try(plot(fit.resultscaled, main = title, ylab="V(r)", xlab="r/a_s", pch=3, cex=0.1))
try(pointsyerr(x=fit.resultscaled$x[finemask], y=fit.resultscaled$y[finemask], dy=fit.resultscaled$dy[finemask], col=1, pch=1))
try(pointsyerr(x=fit.resultscaled$x[!finemask], y=fit.resultscaled$y[!finemask], dy=fit.resultscaled$dy[!finemask], col=2, pch=2))
legend(legend=c("spatial", "temporal"), col=c(2,1), pch=c(2,1), x="topleft")

#plot thermalisation
plot(plaquettecolumn, main="Thermalisation", xlab=sprintf("MCMC-steps/%d", opt$nsave), ylab="P")

#save values of force
resultforce <- data.frame(R=NA, force=NA, forceerr=NA)
for(i in seq (1, length(xx), 1)){
    resultforce <- rbind(resultforce, data.frame(R=xx[i], force=fnpot(fit.resultscaled$t0, xx[i], 0), forceerr=forceerrs[i]))
}
resultforce <- na.omit(resultforce)

rzerocoarse <- determinerzero(fit.resultcoarse, bootsamples=bootsamples)
rzerofine <- determinerzero(fit.resultfine, bootsamples=bootsamples)

differentxis <- data.frame(xirzero=NA, xist=NA)

for(bs in seq(1, bootsamples)){
xirzero <- rzerocoarse[[3]][bs] / rzerofine[[3]][bs]
xist <- fit.resultfine$t[bs, 2] / fit.resultcoarse$t[bs, 2]
differentxis <- rbind(differentxis, data.frame(xirzero=xirzero, xist=xist))
}
#~ print(differentxis)
differentxis <- differentxis[-1,]
differentxiresults <- c(mean(differentxis$xirzero), sd(differentxis$xirzero), mean(differentxis$xist), sd(differentxis$xist))
#~ print(differentxiresults)


}
if(dofit){
write.table(resultforce, file=sprintf("%sresultsforcebeta%fxi%fNs%dnape%dalpha%fomit%dmeffchosen.csv", opt$plotpath, beta, xi, Ns, nape, alpha, opt$omit), row.names=FALSE)


if(opt$determinemeff){
    t1 <- opt$lowerboundmeff
} else{
    t1 <- -1
    }

#~ print(length(bsrzero))
resultssummary <- list(st=fit.resultscaled$t0[2], dst=fit.resultscaled$se[2], bsst=fit.resultscaled$t[,2], 
        rzeros=rzero, drzeros=drzero, crz=-1.65, bsrzeros=bsrzero,
        p=plaquettedata$value, dp=plaquettedata$dvalue, bsp=plaquettecf$cf.tsboot$t[,1],
        beta=beta, xiin=xi, Nt=Nt, Ns=Ns, bootsamples=bootsamples,
        nom=nom, skip=opt$skip, 
        xicalc=xicalc, dxicalc=dxicalc, bsxicalc=xibootsamples)
class(resultssummary) <- "resultssummary"


resultlist <- data.frame(xi=NA, beta=NA, xicalc=NA, dxicalc=NA, r0=NA, dr0=NA, st=NA, dst=NA, p=NA, dp=NA, 
                         chi=NA, chifine=NA, chicoarse=NA, c=NA, bs=NA, Ns=NA, Nt=NA, nape=NA, alpha=NA, omit=NA, nom=NA, skip=NA, xi2=NA,
                         dxi2=NA, xisingle=NA, dxisingle=NA, t1=NA,  job=NA, every=NA, tauint=NA, dtauint=NA, bootl=NA, 
                         xirzero=NA, dxirzero=NA, xist=NA, dxist=NA)
#~ print(length(resultlist))

for(i in seq(1, max(1, length(rzeroofc$c)))){
#~ message(xi, " ",  beta, " ",  xicalc, " ",  dxicalc, " ", 
#~                       rzeroofc$r0[i], " ",  rzeroofc$dr0[i], " ",  
#~                       fit.resultscaled$t0[2], " ",  fit.resultscaled$se[2], " ",  plaquettedata$value, " ",  plaquettedata$dvalue, " ", 
#~                       fit.resultscaled$chisqr/fit.resultscaled$dof, " ",  rzeroofc$c[i], " ",  bootsamples, " ",  
#~                       Ns, " ",  Nt, " ",  nape, " ",  alpha, " ",  opt$omit, " ",  nom, " ",  skip, " ", 
#~                       xisquared, " ",  dxisquared, " ",  xisingle, " ",  dxisingle, " ",  t1, " ",  opt$job)
newline <- data.frame(xi=xi, beta=beta, xicalc=xicalc, dxicalc=dxicalc,
                      r0=rzeroofc$r0[i], dr0=rzeroofc$dr0[i], 
                      st=fit.resultscaled$t0[2], dst=fit.resultscaled$se[2], p=uwerrresults$uwcf$value[1], dp=uwerrresults$uwcf$dvalue[1],
                      chi=fit.resultscaled$chisqr/fit.resultscaled$dof, chifine=fit.resultfine$chisqr/fit.resultfine$dof,
                      chicoarse=fit.resultcoarse$chisqr/fit.resultcoarse$dof, c=rzeroofc$c[i], bs=bootsamples, 
                      Ns=Ns, Nt=Nt, nape=nape, alpha=alpha, omit=opt$omit, nom=nom, skip=skip,
                      xi2=xisquared, dxi2=dxisquared, xisingle=xisingle, dxisingle=dxisingle, t1=t1, job=opt$job,
                      every=opt$every, tauint=uwerrresults$uwcf$tauint[1], dtauint=uwerrresults$uwcf$dtauint[1], bootl=opt$bootl,
                      xirzero=differentxiresults[1], dxirzero=differentxiresults[2], xist=differentxiresults[3], dxist=differentxiresults[4])
                      #-1: chosen individually

resultlist <- rbind(resultlist, newline)
}

resultlist <- resultlist[-1, ]
print(resultlist)
filename <- sprintf("%sresultsummary2p1dsidewaysb%.3fNs%d.csv", opt$plotpath, opt$betaone, opt$Ns)
columnnames <- FALSE
if(!file.exists(filename)){
    columnnames <- TRUE
#~     newline <- data.frame(xi="#input xi", beta="#input beta", xicalc="#xi calculated from limear interpolation", dxicalc="standard deviation",
#~                       r0="#calculated r0", dr0="sd", 
#~                       st="#string tension", dst="#sd", p="#plaquette", dp="#sd",
#~                       chi="#chi2red of fit of combined potential", c="#c used in definition of r0", bs="#bootstrapsamples", 
#~                       Ns="#spatial extent of lattice", Nt="#temporal extent", nape="#number of APE-smears", alpha="#alpha used for APE-smears",
#~                       omit="#last omit points in potential are ignored for fitting", nom="#number of measurements",
#~                       xi2="#xisquared", dxi2="#sd", xisingle="#ratio between linear coefficients of single fits", dxisingle="#sd",
#~                       t1="#lower bound for meff", job="#job number on qbig")
#~     resultlist <- rbind(newline, resultlist)
}
write.table(resultlist, filename, append=TRUE, row.names=FALSE, col.names=columnnames)
nameresults <- sprintf("%sresultsrotatedNs%dNt%dbeta%fxi%fbs%d.RData", opt$plotpath, Ns, Nt, beta, xi, bootsamples)
saveRDS(resultssummary, file=nameresults)


#~ fit.all <- bootstrap.nlsfit(potall, c(0.2, 0.2, 0.2, xi), yall, xall, bsamplesall, maskpot=finemask)
#~ summary(fit.all)
}
message(warnings())

if(opt$plotonlymeff){

filenamelist <- sprintf("%slistmeff2p1dmeffchosenNs%dbeta%fxi%fbsamples%d.RData", opt$plotpath, Ns, beta, xi, bootsamples)

if(opt$smearing){
    filenamelist <- sprintf("%slistmeff2p1dmeffchosenNs%dbeta%fxi%fnape%dalpha%fbsamples%d.RData", opt$plotpath, Ns, beta, xi, nape, alpha, bootsamples)
}
if(file.exists(filenamelist)){
listmeff <- readRDS(file=filenamelist)
} else{
    filenamelist <- sprintf("%slistmeff2p1drotatedNs%dbeta%fxi%fbsamples%dl%d.RData", opt$plotpath, Ns, beta, xi, bootsamples, t1)
if(opt$smearing){
    filenamelist <- sprintf("%slistmeff2p1drotatedNs%dbeta%fxi%fnape%dalpha%fbsamples%dl%d.RData", opt$plotpath, Ns, beta, xi, nape, alpha, bootsamples, t1)
}
listmeff <- readRDS(file=filenamelist)
}


filenameforplots <- sprintf("%splotsonlymeff2p1dmeffchosenNs%dbeta%fxi%f.pdf", opt$plotpath, Ns, beta, xi)
if(opt$smearing){
    filenameforplots <- sprintf("%splotsonlymeff2p1dmeffchosenNs%dbeta%fxi%fnape%dalpha%f.pdf", opt$plotpath, Ns, beta, xi, nape, alpha)
}
pdf(file=filenameforplots, title="")

for(i in seq(1,Ns/2)){
    if(listmeff[[i]][[2]]){
        titleall <- sprintf("beta=%.2f y=%d", beta, i)
        if(!is.null(listmeff[[i]][[1]]$effmassfit)){
            title <- sprintf("%s coarse p=%f", titleall, listmeff[[i]][[1]]$effmassfit$Qval)
        } else{
            title <- sprintf("%s coarse", titleall)
        }
        message(title)
        max <- min(max(na.omit(listmeff[[i]][[1]]$effMass[2:(Ns/2-1)])),4)
        if(i>7){
            max <- min(max(na.omit(listmeff[[i]][[1]]$effMass[2:(Ns/2-1)])),5)
        }
        min <- max(min(na.omit(listmeff[[i]][[1]]$effMass[2:(Ns/2-1)])),0)
        ylim <- c(min, max)
        message(ylim)
        try(plot(listmeff[[i]][[1]], xlab="x/a_s", ylab="Meff", main=title))
        try(plot(listmeff[[i]][[1]], xlab="x/a_s", ylab="Meff", ylim=ylim, main=title))
}
}

for(i in seq(1,Nt/2)){
#~     print(listmeff[[Ns/2+i]][[2]])
    if(listmeff[[Ns/2+i]][[2]]){
        titleall <- sprintf("beta=%.2f t=%d", beta, i)
        if(!is.null(listmeff[[i]][[1]]$effmassfit)){
            title <- sprintf("%s fine p=%f", titleall, listmeff[[i+Ns/2]][[1]]$effmassfit$Qval)
        } else{
            title <- sprintf("%s fine", titleall)
        }
        message(title)
        max <- min(max(na.omit(listmeff[[i+Ns/2]][[1]]$effMass[2:(Ns/2-1)])),4)
        if(i>7){
            max <- min(max(na.omit(listmeff[[i+Ns/2]][[1]]$effMass[2:(Ns/2-1)])),5)
        }
        min <- max(min(na.omit(listmeff[[i+Ns/2]][[1]]$effMass[2:(Ns/2-1)])),0)
        ylim <- c(min, max)
        message(ylim)
        try(plot(listmeff[[i+Ns/2]][[1]], xlab="x/a_s", ylab="Meff", main=title))
        try(plot(listmeff[[i+Ns/2]][[1]], xlab="x/a_s", ylab="Meff", ylim=ylim, main=title))
}
}
}
