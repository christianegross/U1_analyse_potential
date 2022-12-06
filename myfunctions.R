library(hadron)


readloopfilecfrotated <- function(file, path="", skip=0, Nsmax, yt, zerooffset=0, every=1) {
    # Nsmax: maximum extent measured in the spacial x-direction, typically half of spatial extent of lattice
    # yt: read columns with R=sqrt(x^2+yt^2) for x starting at zero into a cf-container
    # set zerooffset to 1 if W(x=0, y=0) and so on has also been measured
  tmp <- as.matrix(read.table(paste(path, file, sep=""), skip=skip))
  Dm <- dim(tmp)
  confno <- tmp[[Dm[2]]]
  nom <- length(tmp[, 1])
#~   print(Dm)
#~   print(confno)
#~   print(nom)
  Time <- Nsmax
  ii <- seq((yt+zerooffset-1)*(Nsmax+zerooffset)+1+zerooffset, (yt+zerooffset)*(Nsmax+zerooffset), 1)
#~   print(ii)
  ret <- cf_meta(nrObs = 1, Time=Time, nrStypes = 1)
  ret <- cf_orig(ret, cf = tmp[seq(1, nom, every) ,ii])
  ret$conf.index <- confno
  
  return (invisible(ret))
}

calcplotWloopsideways <- function(file, skip, Ns, yt, bootsamples, title, path="", nsave=100, zerooffset=0, every=1, l=2){
    WL <- readloopfilecfrotated(file=filename, path=path, skip=skip, Nsmax=Ns/2, yt=yt, zerooffset=zerooffset, every=every)
    WL <- bootstrap.cf(WL, boot.R=bootsamples, boot.l=l)
    title <- sprintf("%s, %d configs\n used every %d, boot.l=%d\n", title, (length(WL$cf[,1])+skip)*nsave, nsave*every, l)
    
#~     title<-sprintf("beta=%.2f, %d configs, skipped %d\ncoarse  y=%d", beta, (length(WL$cf[,1])+skip)*100, skip*100, y)
    plot(WL, log="y", xlab="x/a_s", ylab="C(x)", main=sprintf("%s, logscale", title))
    if (-1 %in% sign(WL$cf0)){
      plot(WL, xlab="x/a_s", ylab="C(x)", main=sprintf("%s, linear scale", title), ylim=c(-WL$cf0[3], WL$cf0[3]))
      points(x=seq(-1,Ns/2+1,1), y=rep(0, times=length(seq(-1,Ns/2+1,1))), col=2, type="l")
      message(WL$cf0, " ")
    }
    return (WL)
}

readloopfilecfsub <- function(file, path="", skip=0, Nytmax, x=x, Nsmax=8, zerooffset=0, every=1) {
    # want to calculate ratios W(x, t)/W(x, t+1)
    # Nytmax: maximum extent measured in secondary direction
    # Nsmax: maximum extent measured in primary direction
    # x: read columns with x=x 
    # set zerooffset to 1 if W(x=0, y=0) and so on has also been measured
  tmp <- as.matrix(read.table(paste(path, file, sep=""), skip=skip))
  Dm <- dim(tmp)
  confno <- tmp[[Dm[2]]]
  nom <- length(tmp[, 1])
  Time <- Nytmax
#~   print(Time)
  ii <- seq((Nsmax+zerooffset)*zerooffset+zerooffset+x, (Nsmax+zerooffset)*(Nytmax+zerooffset), Nsmax+zerooffset)
#~     print(ii)
  ret <- cf_meta(nrObs = 1, Time=Time, nrStypes = 1)
  ret <- cf_orig(ret, cf = tmp[seq(1, nom, every) ,ii])
  ret$conf.index <- confno
  
  return (invisible(ret))
}

readloopfilecfsmall <- function(file, path="", skip=0, Nsmax=4, x, y, Ntmax=8, start=0) {
    # Nsmax: maximum extent measured in the spacial x-direction, typically min(4, lattice size)
    # x, y: read columns with W(t, x, y) for all possible t 
    # Ntmax: temporal extent of the lattice
  tmp <- as.matrix(read.table(paste(path, file, sep=""), skip=skip))
#~   print(tmp[1,])
#~   print(tmp[5,])
  Dm <- dim(tmp)
  confno <- tmp[[Dm[2]]]
  Time <- Ntmax+1-start
#~   message("Time=", Time)
  ii <- seq((x)*(Nsmax+1)+y+1+start*(Nsmax+1)*(Nsmax+1), (Nsmax+1)*(Nsmax+1)*(Ntmax+1), (Nsmax+1)*(Nsmax+1))
#~   print(ii)
  ret <- cf_meta(nrObs = 1, Time=Time, nrStypes = 1)
  ret <- cf_orig(ret, cf = tmp[,ii])
  ret$conf.index <- confno
  
  return (invisible(ret))
}

readloopfilecfplaquette <- function(file, path="", skip=0, Nsmax=4, Ntmax=8, every=1) {
    # Nsmax: maximum extent measured in the spacial x-direction, typically min(4, lattice size)
    # Ntmax: temporal extent of the lattice
    # reads W(x=1, t=1, y=0) and W(x=1, y=1, t=0) into one cf-container 
  tmp <- as.matrix(read.table(paste(path, file, sep=""), skip=skip))
  Dm <- dim(tmp)
  confno <- tmp[[Dm[2]]]
  Time <- 2
  nom <-length(tmp[, 1])
  ii <- c((Nsmax+1)+2, 2*(Nsmax+1)*(Nsmax+1), (Nsmax+1)*(Nsmax+1))
#~   print(ii)
  ret <- cf_meta(nrObs = 1, Time=Time, nrStypes = 1)
  ret <- cf_orig(ret, cf = tmp[seq(1, nom, every) ,ii])
  ret$conf.index <- confno
  
  return (invisible(ret))
}

readloopfilecfonecolumn <- function(file, path="", skip=0, column=1, every=1) {
    #reads column into a cf container
  tmp <- as.matrix(read.table(paste(path, file, sep=""), skip=skip))
  Dm <- dim(tmp)
  confno <- tmp[[Dm[2]]]
  Time <- 1
  nom <- length(tmp[, 1])
  ii <- c(column, column)
#~   print(ii)
#~   print(dim(tmp[,ii]))
  ret <- cf_meta(nrObs = 1, Time=Time, nrStypes = 1)
  ret <- cf_orig(ret, cf = tmp[,ii])
  ret$conf.index <- confno
  
  return (invisible(ret))
}

calcplotWloopnormalspatial <- function(file, skip, Ns, x, bootsamples, title, path="", zerooffset=0, every=1, nsave=100, l=2){
#~     print(Ns)
    WL <- readloopfilecfsub(file=filename, skip=skip, Nytmax=Ns/2, x=x, Nsmax=Ns/2, zerooffset=zerooffset, every=every)
    WL <- bootstrap.cf(WL, boot.R=bootsamples, boot.l=l)
    title <- sprintf("%s, %d configs\n used every %d, boot.l=%d\n", title, (length(WL$cf[,1])+skip)*nsave, nsave*every, l)
    
#~     title<-sprintf("beta=%.2f, %d configs, skipped %d\ncoarse  y=%d", beta, (length(WL$cf[,1])+skip)*100, skip*100, y)
    plot(WL, log="y", xlab="y/a_s", ylab="C(y)", main=sprintf("%s, logscale", title))
    if (-1 %in% sign(WL$cf0)){
      plot(WL, xlab="y/a_s", ylab="C(y)", main=sprintf("%s, linear scale", title), ylim=c(-WL$cf0[3], WL$cf0[3]))
      points(x=seq(-1,Ns/2+1,1), y=rep(0, times=length(seq(-1,Ns/2+1,1))), col=2, type="l")
#~       message(WL$cf0, " ")
    }
    return (WL)
}

calcplotWloopnormaltemporal <- function(file, skip, Ns, Nt, x, bootsamples, title, path="", zerooffset=0, every=1, nsave=100, l=2){
    WL <- readloopfilecfsub(file=filename, skip=skip, Nytmax=Nt/2, x=x, Nsmax=Ns/2, zerooffset=zerooffset, every=every)
    WL <- bootstrap.cf(WL, boot.R=bootsamples, boot.l=l)
    title <- sprintf("%s, %d configs\n used every %d, boot.l=%d\n", title, (length(WL$cf[,1])+skip)*nsave, nsave*every, l)
    
#~     title<-sprintf("beta=%.2f, %d configs, skipped %d\ncoarse  y=%d", beta, (length(WL$cf[,1])+skip)*100, skip*100, y)
    plot(WL, log="y", xlab="t/a_s", ylab="C(t)", main=sprintf("%s, logscale", title))
    if (-1 %in% sign(WL$cf0)){
      plot(WL, xlab="t/a_s", ylab="C(t)", main=sprintf("%s, linear scale", title), ylim=c(-WL$cf0[3], WL$cf0[3]))
      points(x=seq(-1,Nt/2+1,1), y=rep(0, times=length(seq(-1,Nt/2+1,1))), col=2, type="l")
#~       message(WL$cf0, " ")
    }
    return (WL)
}

calcplotWloopsmall <- function(filename, skip, Nsmax, Ntmax, Nt, x, y, start, bootsamples, title, nsave=0, l=2){
    WL <- readloopfilecfsmall(file=filename, skip=skip, Nsmax=Nsmax, x=x, y=y, Ntmax=Ntmax, start=start)
    WL <- bootstrap.cf(WL, boot.R=bootsamples, boot.l=l)
#~     message("WL calculated")
    title <- sprintf("%s%d configs, skipped %d",title, (length(WL$cf[,1])+skip-1)*nsave, skip*nsave)
#~         print(names(WL))
    try(plot(WL, log="y", xlab="t/a_t", ylab="C(x)", main=sprintf("%s, logscale", title)))
    arrows((Nt-start)/2, 10, (Nt-start)/2, 1e-08, angle=90, length=0.1, code=0, col=2)
#~     message("logscale plotted")
    if (-1 %in% sign(WL$cf0)){
      try(plot(WL, xlab="t/a_t", ylab="C(x)", main=sprintf("%s, linear scale", title), ylim=c(-WL$cf0[3], WL$cf0[3])))
      points(x=seq(-1,Ntmax+1,1), y=rep(0, times=length(seq(-1,Ntmax+1,1))), col=2, type="l")
      arrows((Nt-start)/2, 10, (Nt-start)/2, -10, angle=90, length=0.1, code=0, col=2)
#~       message("linear plotted")
      message(WL$cf0, " ")
    }
    return (WL)
}

deteffmass <- function(WL, yt, potential, t1, t2, isspatial){
    ylim <- c(0,0)
    WL.effmass <- bootstrap.effectivemass(WL, type="log")
#~     message("meff determined")
    if(length(na.omit(as.vector(WL.effmass$effMass))) > 1) {
        WL.effmass <- try(fit.effectivemass(WL.effmass, t1=t1, t2=t2, useCov=TRUE), silent=FALSE)
        if (!inherits(WL.effmass, "try-error")) {
#~             message("fit not successful")
            newline <- data.frame(R=yt, m=WL.effmass$effmassfit$t0[1], dm=WL.effmass$effmassfit$se[1], space=isspatial, p=WL.effmass$effmassfit$Qval, 
                                chi=WL.effmass$effmassfit$chisqr/WL.effmass$effmassfit$dof)
#~             print(newline)
            potential <- rbind(potential, newline)
            regionplotzoom <- seq(ceiling(0.2*length(WL.effmass$effMass))+1, length(WL.effmass$effMass)-1)
#~             print(WL.effmass$effMass)
#~             print(regionplotzoom)
            if(length(regionplotzoom > 0)){
                max <- min(max(na.omit(WL.effmass$effMass[regionplotzoom])),2)
                min <- max(min(na.omit(WL.effmass$effMass[regionplotzoom])),0)
                ylim <- c(min, max)
            }
#~             print(ylim)
        } else {
#~             message("fit done and successful")
            #have to recalculate so the correct effective masses are stored instead of error
            WL.effmass <- bootstrap.effectivemass(WL, type="log")
        }
    }
    return (list(WL.effmass, ylim))
}

deteffmasssmall <- function(WL, x, y, t1, t2){
    ylim <- c(0,0)
    listresults <- list(FALSE, FALSE, FALSE)
    newline <- data.frame(R=NA, m=NA, dm=NA)
    WL.effmass <- bootstrap.effectivemass(WL, type="log")
#~     message("meff determined")
    if(length(na.omit(as.vector(WL.effmass$effMass))) > 1) {
        WL.effmass <- try(fit.effectivemass(WL.effmass, t1=t1, t2=t2, useCov=TRUE), silent=FALSE)
        listresults <- list(WL.effmass, x, y)
        if (inherits(WL.effmass, "try-error")) {
#~             message("fit not successful")
            #have to recalculate so the correct effective masses are stored instead of error
            WL.effmass <- bootstrap.effectivemass(WL, type="log")
            listresults <- list(FALSE, FALSE, FALSE)
            newline <- data.frame(R=sqrt(x^2+y^2), m=NA, dm=NA)
        } 
        else{
#~             message("fit done and successful")
            regionplotzoom <- seq(floor(length(WL.effmass$effMass))*0.25, length(WL.effmass$effMass)-1)
            if(length(regionplotzoom > 0)){
                max <- min(max(na.omit(WL.effmass$effMass[regionplotzoom])),2)
                min <- max(min(na.omit(WL.effmass$effMass[regionplotzoom])),0)
                ylim <- c(min, max)
            }
            newline <- data.frame(R=sqrt(x^2+y^2), m=WL.effmass$effmassfit$t0[1], dm=WL.effmass$effmassfit$se[1])
        }
    }
    return (list(newline, ylim, listresults))
}

plotyerr <- function (x, y, dy, arlength=0.02, log="", ...){
    xlo <- x
    xhi <- x
    ylo <- y-dy
    yhi <- y+dy
    dots <- list(...)
    if (!'ylim' %in% names(dots)){
        ylim <- c(min(na.omit(ylo)), max(na.omit(yhi)))
        plot(x=x, y=y, ylim=ylim, log=log, ...)
    } else{
        plot(x=x, y=y, log=log, ...)
    }
    arrows(xhi, yhi, xlo, ylo, angle=90, length=arlength, code=3, ...)
}

pointsyerr <- function (x, y, dy, arlength=0.02, ...){
    xlo <- x
    xhi <- x
    ylo <- y-dy
    yhi <- y+dy
    points(x=x, y=y, ...)
    arrows(xhi, yhi, xlo, ylo, angle=90, length=arlength, code=3, ...)
}

plotxyerr <- function (x, y, dx, dy, arlength=0.02, ...){
    xlo <- x-dx
    xhi <- x+dx
    ylo <- y
    yhi <- y
    plot(x=x, y=y, ...)
    arrows(xhi, yhi, xlo, ylo, angle=90, length=arlength, code=0, ...)
    xlo <- x
    xhi <- x
    ylo <- y-dy
    yhi <- y+dy
    arrows(xhi, yhi, xlo, ylo, angle=90, length=arlength, code=0, ...)
}

pointsxyerr <- function (x, y, dx, dy, arlength=0.02, ...){
    xlo <- x-dx
    xhi <- x+dx
    ylo <- y
    yhi <- y
    points(x=x, y=y, ...)
    arrows(xhi, yhi, xlo, ylo, angle=90, length=arlength, code=0, ...)
    xlo <- x
    xhi <- x
    ylo <- y-dy
    yhi <- y+dy
    arrows(xhi, yhi, xlo, ylo, angle=90, length=arlength, code=0, ...)
}

errorpolygon <- function(X, fitresult, col.p, col.band="gray", polygon=TRUE, arlength=0.1, pch=1,...){
    prediction <- predict(fitresult, X)
    if(missing(pch)){ 
        p.pch<-col.p
    } else{
        p.pch <- pch
    }
    if(fitresult$errormodel == "yerrors") {
      limits <- pointsyerr(x=fitresult$x, y=fitresult$y, dy=fitresult$dy, col=col.p, pch=p.pch, arlength=arlength, ...)
    }
    else {
      limits <- pointsxyerr(x=fitresult$x, y=fitresult$y, dy=fitresult$dy, dx=fitresult$dx, col=col.p, pch=p.pch, arlength=arlength, ...)
    }
#~     limits <- plotwitherror(x=fitresult$x, y=fitresult$y, dy=fitresult$dy, col=col.p, pch=p.pch, rep=TRUE, ...)
    
    xlim <- limits$xlim
    ylim <- limits$ylim
    
    if(polygon){
        polyval <- c(prediction$val + prediction$err, rev(prediction$val - prediction$err))
        if(any(polyval < ylim[1]) || any(polyval > ylim[2])) {
          polyval[polyval < ylim[1]] <- ylim[1]
          polyval[polyval > ylim[2]] <- ylim[2]
        }
        col.band="gray"
        pcol <- col2rgb(col.band, alpha=TRUE)/255 
        pcol[4] <- 0.65
        pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])
        polygon(x=c(X, rev(X)), y=polyval, col=pcol, lty=0, lwd=0.001, border=pcol)
    }
    lines(x=X, y=prediction$val, col=col.p, ...)
}

points.effectivemass <- function (x, ..., ref.value, col, col.fitline) {
  effMass <- x
  if(missing(col)) {
    col <- c("black", rainbow(n=(effMass$nrObs-1)))
  }
  if(missing(col.fitline)) {
	col.fitline <- col[1]
  }
  # BaKo: is this also valid for acosh type effective masses?
  t <- effMass$t.idx
  suppressWarnings(pointsyerr(x=t-1, y=effMass$effMass[t], dy=effMass$deffMass[t], col=col[1], ...))
  if(effMass$nrObs > 1) {
    for(i in 1:(effMass$nrObs-1)) {
      suppressWarnings(pointsyerr(x=t-1, y=effMass$t0[t+i*length(t)], dy=effMass$se[t+i*length(t)], rep=TRUE, col=col[i+1], ...))
    }
  }
  if(!missing(ref.value)) {
    abline(h=ref.value, col=c("darkgreen"), lwd=c(3))
  }
  if(!is.null(effMass$effmassfit)){
    lines(x=c(effMass$t1,effMass$t2),
          y=c(effMass$effmassfit$t0[1],effMass$effmassfit$t0[1]),
          col=col.fitline,
          lwd=1.3)
      pcol <- col2rgb(col.fitline,alpha=TRUE)/255                                                                                                   
      pcol[4] <- 0.65
      pcol <- rgb(red=pcol[1],green=pcol[2],blue=pcol[3],alpha=pcol[4])
      rect(xleft=effMass$t1, ybottom=effMass$effmassfit$t0[1]-effMass$effmassfit$se[1],
           xright=effMass$t2, ytop=effMass$effmassfit$t0[1]+effMass$effmassfit$se[1],
           col=pcol,
           border=NA)
  }
}

determinerzero <- function(fit.result, bootsamples=2000, c=-1.65, xi=1){
    #Determine r0 as solution of equation -r^2 d/dr V(r)=c, solve for each bootstrapsample, r0=mean pm sd
    # V(r)=a + sigma*r + b*ln(r)
    # -d/dr V(r) = -sigma -b/r
    # c = -sigma*r^2 - b*r, use p-q-formula for p=b/sigma, q=c/sigma
#~     rzerolist <- data.frame(r1=NA, r2=NA)
    # if a_t V is determined, c has to be resacled so effectively a_s V is used
    c <- c*xi
    rzerolist <- c()
    for (bs in seq(1, bootsamples)){
    #~     a <- fit.resultscaled$t[bs, 1]
        sigma <- fit.result$t[bs, 2]
        b <- fit.result$t[bs, 3]
        r1 <- -b/(2*sigma) + sqrt( (b/(2*sigma))^2 - c/sigma)
#~         r2 <- -b/(2*sigma) - sqrt( (b/(2*sigma))^2 - c/sigma)
#~         rzerolist <- rbind(rzerolist , data.frame(r1=r1, r2=r2))
        rzerolist[bs] <- r1
    }
    rzerolist <- na.omit(rzerolist)
    
#~     rzero <- mean(rzerolist$r1)
#~     drzero <- sd(rzerolist$r1)
    rzero <- mean(rzerolist)
    drzero <- sd(rzerolist)
    return (list(rzero, drzero, rzerolist))
}

drawallticks <- function(all=FALSE, inward=FALSE){
    # want to draw tick labels on every side, but do not always want tick labels
    if(!inward){ 
        tick=-0.02
    } else{
        tick=0.02
    }
    axis(1, labels=TRUE, tck=tick)
    axis(2, labels=TRUE, tck=tick)
    axis(3, labels=all, tck=tick)
    axis(4, labels=all, tck=tick)
}

find_bootsamples <- function(Ns, Nt, beta, xi, plotpath="", smearing=FALSE, nape=0, alpha=1){
    basestring <- sprintf("%slistmeffsmallNs%dNt%dbeta%fxi%fbsamples*", plotpath, Ns, Nt, beta, xi)
    if(smearing){
        basestring <- sprintf("%slistmeffsmallNs%dNt%dbeta%fxi%fnape%dalpha%fbsamples*", plotpath, Ns, Nt, beta, xi, nape, alpha)
    }
#~     print(basestring)
    len <- nchar(basestring)
    files <- list.files(pattern=basestring)
    if(length(files)==0){
        stop("There are no data for this set of parameters")
    }
#~     print(files)
    numbers <- c()
    for(i in seq(1, length(files))){
        numbers[i] <- substring(files[i], len, len+4)
    }
#~     print(numbers)
    numbers <- suppressWarnings(as.integer(numbers))
    for(i in seq(1, length(numbers))){
#~         print(i)
        for(j in seq(3, 1)){
#~             print(j)
            if(!is.na(numbers[i])){
#~                 print("fixed")
                break
            } else{
                numbers[i] <- suppressWarnings(as.integer(substring(files[i], len, len+j)))
            }
        }
    }
    print("simulations with the following amount of bootstrapsamples present:")
    print(numbers)
    return(max(numbers))
}

fitplotfunctions <- function(fun1, fun2, x1, y1, dy1, x2, y2, dy2, par1, par2, title="", boot.R=500, bs=FALSE, bsamples1=array(), bsamples2=array(), verbose=FALSE, polygon=TRUE, ...){
    #generate bootstrapsamples, do the fit
    if (!bs){
        bsamples1 <- parametric.bootstrap(boot.R, c(y1), c(dy1))
        bsamples2 <- parametric.bootstrap(boot.R, c(y2), c(dy2))
    }
    fit.result1 <- try(bootstrap.nlsfit(fun1, par1, y1, x1, bsamples1))
    fit.result2 <- try(bootstrap.nlsfit(fun2, par2, y2, x2, bsamples2))

    
    #determine proper limits
    dots <- list(...)
    my.ylim <- FALSE
    my.xlim <- FALSE
    if (!'ylim' %in% names(dots)){
        ylimnew <- c(min(na.omit(y1-dy1), na.omit(y2-dy2)), max(na.omit(y1+dy1), na.omit(y2+dy2)))
        my.ylim <- TRUE
    }
    if (!'xlim' %in% names(dots)){
        xlimnew <- c(min(na.omit(x1), na.omit(x2)), max(na.omit(x1), na.omit(x2)))
        my.xlim <- TRUE
        xmin <- min(xlimnew)
        xmax <- max(xlimnew)
    } else{
        xmin <- min(dots$xlim)
        xmax <- max(dots$xlim)
    }
    if(verbose) print(c(xmin, xmax))
    if(verbose) print("dy1")
    if(verbose) print(dy1)
    if(verbose) print("dy2")
    if(verbose) print(dy2)
    
    #points for plotting the result
    xvalues <- seq(from=xmin-1, to=xmax+1, length.out=200)
    #plot everything together
    if(my.ylim && my.xlim){
        plotyerr(x=x1, y=y1, dy=dy1, main=title, col=1, pch=1, xlim=xlimnew, ylim=ylimnew, ...)
        if(verbose) print("no limits given")
    }
    if(!my.ylim && my.xlim){
        plotyerr(x=x1, y=y1, dy=dy1, main=title, col=1, pch=1, xlim=xlimnew, ...)
        if(verbose) print("y-limits given")
    }
    if(my.ylim && !my.xlim){
        plotyerr(x=x1, y=y1, dy=dy1, main=title, col=1, pch=1, ylim=ylimnew, ...)
        if(verbose) print("x-limits given")
    }
    if(!my.ylim && !my.xlim){
        plotyerr(x=x1, y=y1, dy=dy1, main=title, col=1, pch=1, ...)
        if(verbose) print("both limits given")
    }
    pointsyerr(x=x2, y=y2, dy=dy2, col=2, pch=2, ...)
    
    if(!inherits(fit.result1, "try-error")){
        errorpolygon(xvalues, fit.result1, col.p=1, col.band="gray", polygon=polygon)
        if(verbose) print("res1")
        if(verbose) print(names(fit.result1))
        if(verbose) print("bsamples")
        if(verbose) print(apply(bsamples1, 2, mean))
        if(verbose) print(apply(bsamples1, 2, sd))
        if(verbose) print("fit")
        if(verbose) print(fit.result1$x)
        if(verbose) print(fit.result1$y)
        if(verbose) print(fit.result1$dy)
    }
    if(!inherits(fit.result2, "try-error")){
        errorpolygon(xvalues, fit.result2, col.p=2, col.band="red", polygon=polygon)
        if(verbose) print("res2")
        if(verbose) print("bsamples")
        if(verbose) print(apply(bsamples2, 2, mean))
        if(verbose) print(apply(bsamples2, 2, sd))
        if(verbose) print("fit")
        if(verbose) print(fit.result2$x)
        if(verbose) print(fit.result2$y)
        if(verbose) print(fit.result2$dy)
    }
        
#~     print("dy1")
#~     print(dy1)
#~     print("dy2")
#~     print(dy2)
    return (list(fit.result1, fit.result2))
}

readinbootstrapsamples <- function(beta, Ns, Nt, xi, bootsamples=500, path="", names=c('bsslope', 'bsst', 'bsratio', 'bsp', 'bsrzeros'), columns=c(1,1,1,1,5), filename="resultspot"){
    if(length(columns)!=length(names)){
        stop("The columns per name and names do not have the same length! check input")
    }
    filename <- sprintf("%s%sNs%dNt%dbeta%fxi%fbs%d.RData", path, filename, Ns, Nt, beta, xi, bootsamples)
    listres <- readRDS(filename)
#~     print(head(listres))
#~     print(listres)
#~     print(names(listres))
    bsamples <- array(c(rep(NA, sum(columns)*bootsamples)), dim=c(bootsamples, sum(columns)))
    index <- 1
    for(i in seq(1, length(names))){
        name <- names[i]
#~         print(names[i])
        for(j in seq(1, columns[i])){
#~             print(j)
#~             print(length(bsamples[,index]))
            if (columns[i] > 1){
#~                 print(length(listres[[name]][[j]]))
                bsamples[, index] <- listres[[name]][[j]]
            }else{
#~                 print(length(listres[[name]]))
#~                 print(listres[[name]])
                bsamples[, index] <- listres[[name]]
            }
            index <- index+1
        }
    }
    return (bsamples)
}


mmtoinches <- function(size) return (size/25.4)


uwerr.cfboot <- function(cf){
  stopifnot(inherits(cf, 'cf_boot'))

  uw_wrapper <- function(x){
    uw_tmp <- try(uwerrprimary(data=x), silent=TRUE)
    if( any(class(uw_tmp) == "try-error") ){
      c(value=NA, dvalue=NA, ddvalue=NA, tauint=NA, dtauint=NA)
    } else {
      c(value=uw_tmp$value, dvalue=uw_tmp$dvalue, ddvalue=uw_tmp$ddvalue,
        tauint=uw_tmp$tauint, dtauint=uw_tmp$dtauint)
    }
  }

  res <- list()
  res[["uwcf"]] <- cbind(as.data.frame(t(apply(X=cf$cf.tsboot$t, MARGIN=2L, FUN=uw_wrapper))),
                         t=(1:ncol(cf$cf.tsboot$t)))
  if( has_icf(cf) ){
    res[["uwicf"]] <- cbind(as.data.frame(t(apply(X=cf$icf.tsboot$t, MARGIN=2L, FUN=uw_wrapper))),
                            t=(1:ncol(cf$icf.tsboot$t)))
  }
  return(res)
}
