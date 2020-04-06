
colF <- colorRampPalette( c('#8c510a','#d8b365','#c7eae5','#5ab4ac','#01665e','#2166ac') )

colT <- colorRampPalette( c('#8c510a','#d8b365','#c7eae5','#5ab4ac','#01665e','#2166ac') )

mergeList <- function( list1, list2 ){
  
  # update elements of list1 if contained list2
  
  stopifnot(is.list(list1), is.list(list2))
  
  n2 <- names(list2)
  
  for(k in n2){
    if(!k %in% names(list1)){
      list1 <- append( list1, list2[k] )
    }else{
      list1[k] <- list2[k]
    }
  }
  list1
}

.pasteCols <- function(mm){
  tmp <- apply(mm,1,paste0,collapse='-')
  names(tmp) <- NULL
  tmp
}


gjamPredTime <- function(output, nsim = 10, quant = .05, minw = .1){
  
  # evaluates carrying capacity for TIME
  # must supply either 'groups' or 'ngrid'
  # groups - xdata$groups column to predict for groups
  # ngrid  - no. values per covariate if !BYGROUP
  # nsim   - simulation steps per covariate combination
  # NOTE: x is centered and standardized
  
  #  breaks <- optim <- rhoStandXmu <- sensAlpha <- sensBeta <- sensRho <- 
  #    wA <- wL <- xStand <- xUnstand <- NULL
  
  x          <- output$inputs$xStand
  xl         <- output$inputs$xRho
  #xnames     <- output$inputs$xnames
  xlnames    <- output$inputs$xlnames
  factorBeta <- output$inputs$factorBeta
  factorRho  <- output$inputs$factorRho
  ydata      <- output$inputs$y
  effMat     <- output$inputs$effMat
  uindex     <- output$parameters$uindex
  gindex     <- output$parameters$gindex
  other      <- output$inputs$other
  notOther   <- output$inputs$notOther
  bg         <- output$parameters$betaMu
  xnames     <- rownames(bg)
  Rmat       <- output$parameters$RmatStandXmu
  Amat       <- output$parameters$Amu
  sigMu      <- output$parameters$sigMu
  wB         <- output$parameters$wB
  wL         <- output$parameters$wL
  wA         <- output$parameters$wA
  bgibbs     <- output$chains$bgibbs
  lgibbs     <- output$chains$lgibbs
  alphaGibbs <- output$chains$alphaGibbs
  groups     <- output$inputs$xdata$groups
  ng         <- output$modelList$ng
  burnin     <- output$modelList$burnin
  timeList   <- output$inputs$timeList
  groupCol   <- output$inputs$timeList$groups
  timeCol    <- output$inputs$timeList$times
  group      <- output$inputs$xdata[,groupCol]
 # time       <- output$inputs$xdata[,timeCol]
  assign(timeCol, time)
  time       <- output$inputs$xdata$times
  timeZero   <- output$inputs$timeList$timeZero
 # year       <- output$inputs$xdata[,groupCol]
  
  cx <- colnames(xl)[!colnames(xl) %in% colnames(x)]
  if(length(cx) > 0)x <- cbind(x, xl[,cx])
  
  S <- ncol(ydata)
  snames <- colnames(ydata)
  termB  <- termR <- termA <- FALSE
  
  if(!is.null(bgibbs))termB <- TRUE
  if(!is.null(lgibbs))termR <- TRUE
  if(!is.null(alphaGibbs))termA <- TRUE
  
  w <- ydata/effMat
  
 # wstart <- apply( ydata/effMat, 2, mean, na.rm=T)
  
  fn <- function(wt, termA, termB, termR, xmu, xlmu,
                 gindex, uindex, notOther,
                 amat, bmat, rmat, epsilon){
    wz <- wt
    wz[ wz < 0 ] <- 0
    ff <- wz*0
    
    if(termB)ff[,notOther] <- t(xmu[,rownames(bmat)]%*%bmat[,notOther])
    if(termR){
      Vmat <- wz[ drop=F, ,gindex[,'colW']]*xlmu[ drop=F, ,gindex[,'colX']]
      colnames(Vmat) <- rownames(gindex)
      ff[,notOther] <- ff[,notOther] + Vmat%*%rmat[colnames(Vmat),notOther]
    }
    if(termA){
      Umat <- wz[ drop=F, ,uindex[,1]]*wz[ drop=F, ,uindex[,2]]
      ff[,notOther] <- ff[,notOther] + Umat%*%amat[,notOther]
    }
    ff[,notOther] <-  ff[,notOther] + wz[,notOther] + epsilon
    ff
  }
  
  
  ii <- burnin:ng
  
  times  <- sort(unique(time))
  groups <- sort(unique(group))
  nt     <- length(times)
  ng     <- length(groups)
  
  pbar <- txtProgressBar(min = 1,max = nt, style=1)
  
  ylo <- yhi <- ymu <- w*0
  
  yj <- array(NA, dim = c( ng, ncol(w), nsim) )
  dimnames(yj)[[1]] <- groups
  dimnames(yj)[[2]] <- colnames(w)
  
  wj <- w[timeZero+1,] + .01
  for(k in 1:nsim)yj[,,k] <- wj
  
  for(j in 1:nt){
    
    gj <- which(time == j)  # x[gj,] is one step ahead of w[timeZero,]
    gi <- group[gj]
    
    zj <- yj <- yj[gi,,] 
    zj[ zj < 0 ] <- 0
    zj <- zj + minw
    
    if(j> 1){
      
      for(k in 1:nsim){
        
        ig <- sample(ii, 1)   
        
        muw <- w*0
        
        if(termB){
          bmat <- bg*0
          bmat[ 1:length(bg) ] <- bgibbs[ig,]
        }
        
        if(termA){
          amat <- Amat*0
          amat[ wA ] <- alphaGibbs[ig,]
        }
        
        if(termR){
          rmat <- Rmat*0
          rmat[ wL ] <- lgibbs[ig,]
        }
        
        epsilon <- .rMVN(1, 0, sigma = sigMu )[notOther]
        
        yj[,,k] <- fn(wt = zj[ , ,k], termA, termB, termR, 
                      xmu = x[drop=F, gj, ], xlmu = xl[drop=F, gj, ],
                      gindex, uindex, notOther,
                      amat, bmat, rmat, epsilon)
      }
    }
    
    yj[ yj < 0 ] <- 0
    ymu[gj,] <- apply( yj, c(1,2), mean )
    ylo[gj,] <- apply( yj, c(1,2), quantile, quant )
    yhi[gj,] <- apply( yj, c(1,2), quantile, 1 - quant )
    
    setTxtProgressBar(pbar,j)
    
  }
  
  ymu <- ymu*effMat[,notOther]
  ylo <- ylo*effMat[,notOther]
  yhi <- yhi*effMat[,notOther]
  
  list(ymu = ymu, ylo = ylo, yhi = yhi)
}
##########################


##########################

.wrapperEquilAbund <- function(output, covars = NULL, nsim = 10, ngrid = NULL, 
                               BYFACTOR = FALSE, BYGROUP = TRUE, verbose = FALSE){
  
  # evaluates equil abundances for TIME
  # covars - names of covariates in x to use as gradients
  # must supply 'ngrid' or BYGROUP
  # groups - xdata$groups column to predict for groups
  # ngrid  - no. values per covariate if !BYGROUP
  # nsim   - simulation steps per covariate combination
  # NOTE: x is centered and standardized
  
  #  breaks <- optim <- rhoStandXmu <- sensAlpha <- sensBeta <- sensRho <- 
  #    wA <- wL <- xStand <- xUnstand <- NULL
  
  x          <- output$inputs$xStand
  xl         <- output$inputs$xRho
  xnames     <- output$inputs$xnames
  xlnames    <- output$inputs$xlnames
  factorBeta <- output$inputs$factorBeta
  factorRho  <- output$inputs$factorRho
  ydata      <- output$inputs$y
  effMat     <- output$inputs$effMat
  uindex     <- output$parameters$uindex
  gindex     <- output$parameters$gindex
  other      <- output$inputs$other
  notOther   <- output$inputs$notOther
  bg         <- output$parameters$betaMu
  Rmat       <- output$parameters$RmatStandXmu
  Amat       <- output$parameters$Amu
  sigMu      <- output$parameters$sigMu
  wB         <- output$parameters$wB
  wL         <- output$parameters$wL
  wA         <- output$parameters$wA
  bgibbs     <- output$chains$bgibbs
  lgibbs     <- output$chains$lgibbs
  alphaGibbs <- output$chains$alphaGibbs
  groups     <- output$inputs$xdata$groups
  ng         <- output$modelList$ng
  burnin     <- output$modelList$burnin
  
  if( !is.null(ngrid) )BYGROUP <- F
  
  cx <- colnames(xl)[!colnames(xl) %in% colnames(x)]
  if(length(cx) > 0)x <- cbind(x, xl[,cx])
  
  S <- ncol(ydata)
  snames <- colnames(ydata)
  termB  <- termR <- termA <- FALSE
  
  if(!is.null(bgibbs))termB <- TRUE
  if(!is.null(lgibbs))termR <- TRUE
  if(!is.null(alphaGibbs))termA <- TRUE
  
  
  xnames  <- colnames(x)
  vxnames <- xnames[ !xnames %in% factorBeta$isFactor ][-1] # exclude intercept
  vxnames <- vxnames[ !vxnames %in% factorRho$isFactor ]
  
  gnames <- vxnames
  
  if( !is.null(covars) )gnames <- vxnames[ vxnames %in% covars ]
  
  ix <- grep(':', xnames)
  intb <- xnames[ix]
  if (length(ix) > 0)xnames <- xnames[ -ix ]
  
  ix <- grep(':', xlnames)
  intr <- xlnames[ix]
  if (length(ix) > 0)xlnames <- xlnames[ -ix ]
  
  ix <- grep(':', vxnames)
  intv <- vxnames[ix]
  if (length(ix) > 0)vxnames <- vxnames[ -ix ]
  
  if( !BYGROUP ){   # a prediction grid for covariates
    
    sgrid <- 0
    if(ngrid > 1)sgrid <- seq(-2.2, 2.2, length = ngrid)  # std devs for standardized x
    
    grid <- vector( 'list', length(vxnames) )
    for(k in 1:length(vxnames)){
      if(vxnames[k] %in% gnames){
        grid[[k]] <- sgrid
      }else{
        grid[[k]] <- 0
      }
    }
    names(grid) <- vxnames
    xgrid <- expand.grid(grid)
    
  }else{                # predict for groups
    ii <- rep( groups, ncol(x) )
    jj <- rep( colnames(x), each = nrow(x) )
    xgrid <- tapply( as.vector(x), list(groups = ii, x = jj), mean, na.rm=T)
    xgrid <- xgrid[,vxnames, drop=F]
  }
  
  fnames <- unique( c(factorBeta$facNames,  factorRho$facNames) )
  
  nfact <- length(fnames)
  factorColumns <- character(0)

  if(nfact > 0){
    
    flist  <- vector( 'list', length(fnames))
    
    for(k in 1:length(fnames)){ # find factor in either factorBeta or factorRho
      
      
      fk <- fl <- numeric(0)
      knames <- lnames <- character(0)
      
      if ( !is.null(factorBeta$factorList) ) {
        wf <- grep(fnames[k], names(factorBeta$factorList))
        if(length(wf) > 0){
          fk <- factorBeta$contrast[fnames[k]][[wf]]
          knames <- factorBeta$factorList[fnames[k]][[1]]
          fcc   <- factorBeta$factorList[[wf]]
        }
      }
      if (!is.null(factorRho$factorList) ) {
        wf <- grep(fnames[k], names(factorRho$factorList))
        if(length(wf) > 0){
          fk <- factorRho$contrast[fnames[k]][[1]]
          knames <- factorRho$factorList[fnames[k]][[1]]
          fcc   <- factorRho$factorList[[wf]]
        }
      }
      if(!BYFACTOR)fk <- fk[drop = F, 1, ]
      
      xindex <- rep( 1:nrow(xgrid), each = nrow(fk) )
      kindex <- rep( 1:nrow(fk), nrow(xgrid) )
      
      colnames(fk) <- fcc
      
      xgrid <- cbind(xgrid[drop=F, xindex,], fk[drop=F, kindex,])
      
      factorColumns <- c( factorColumns, fcc )
    }
  }
      
  
  intercept <- 1
  xgrid <- as.matrix( cbind(intercept, xgrid) )
  xgrid <- xgrid[drop=F, ,xnames]
  
  attr(xgrid,'factors') <- factorColumns
  
  wstart <- apply( ydata/effMat, 2, mean, na.rm=T)
  
  fn <- function(w, termA, termB, termR, xmu, xlmu,
                 gindex, uindex, notOther,
                 amat, bmat, rmat, epsilon){
    wz <- w
    wz[ wz < 0 ] <- 0
    ff <- wz*0
    
    if(termB)ff[notOther] <- t(xmu[,rownames(bmat)]%*%bmat[,notOther])
    if(termR){
      Vmat <- wz[ gindex[,'colW']]*xlmu[ gindex[,'colX']]
      names(Vmat) <- rownames(gindex)
      ff[notOther] <- ff[notOther] + t(Vmat%*%rmat[names(Vmat),notOther])
    }
    if(termA){
      Umat <- matrix( wz[uindex[,1]]*wz[uindex[,2]], 1)
      ff[notOther] <- ff[notOther] + t(Umat%*%amat[,notOther])
    }
    ff[notOther] <-  ff[notOther] - epsilon
    sum( ff^2 )
  }
  
  lo <- rep(0, S)
  hi <- 2*max(ydata/effMat)
  
  if( verbose ) print(dim(xgrid))
  
  if(nrow(xgrid) > 10)pbar <- txtProgressBar(min = 1,max = nrow(xgrid), style=1)
  
  ccMu <- matrix(0, nrow(xgrid), S)
  colnames(ccMu) <- snames
  rownames(ccMu) <- rownames(xgrid)
  ccSd <- ccMu
  
  xnames <- rownames(bg)
  
  ii <- burnin:ng
  
  for(j in 1:nrow(xgrid)){
    
    xmu   <- xgrid[drop=F, j, xnames]
    xlmu  <- xgrid[drop=F, j, xlnames]
    wstar <- matrix(NA, nsim, S)
    
    for(k in 1:nsim){
      
      ig <- sample(ii, 1)
      
      if(termB){
        bmat <- bg*0
        bmat[ 1:length(bg) ] <- bgibbs[ig,]
        
        cnn <- colnames(xmu)
        if(length(intb) > 0){
          for(m in 1:length(intb)){
            cii <- unlist( strsplit(intb[m], ':') )
            xmu <- cbind(xmu, xmu[1,cii[1]]*xmu[1,cii[2]])
          }
          colnames(xmu) <- c(cnn,intb)
        }
      }
      
      if(termA){
        amat <- Amat*0
        amat[ wA ] <- alphaGibbs[ig,]
      }
      
      if(termR){
        rmat <- Rmat*0
        rmat[ wL ] <- lgibbs[ig,]
        cnn <- colnames(xlmu)
        if(length(intr) > 0){
          for(m in 1:length(intr)){
            cii <- unlist( strsplit(intr[m], ':') )
            xlmu <- cbind(xlmu, xlmu[1,cii[1]]*xlmu[1,cii[2]])
          }
          colnames(xlmu) <- c(cnn,intr)
        }
      }
      
      epsilon <- .rMVN(1, 0, sigma = sigMu )[notOther]
      
      wstar[k,] <- optim(wstart, fn, method = "L-BFGS-B", 
                         termA = termA, termB = termB, termR = termR, 
                         xmu = xmu, xlmu = xlmu,
                         gindex = gindex, uindex = uindex, 
                         notOther = notOther,
                         amat = amat, bmat = bmat, rmat = rmat, 
                         epsilon = epsilon, lower = lo, upper  = hi )$par
    }
    ccMu[j,] <- colMeans(wstar)
    ccSd[j,] <- apply(wstar, 2, sd)
    
    if(nrow(xgrid) > 10)setTxtProgressBar(pbar,j)
    
  }
  ccMu[,other] <- ccSd[,other] <- NA
  
  list(x = xgrid, ccMu = ccMu, ccSd = ccSd)
}

plotEquilAbund <- function(output, nsim = 20, ngrid = NULL, BYFACTOR = FALSE, 
                           verbose = T){
  
  wstar <- .wrapperEquilAbund(output, nsim, ngrid, BYFACTOR, 
                              verbose = verbose)
  ccMu <- wstar$ccMu[,notOther]
  ccSd <- wstar$ccSd[,notOther]
  ccx  <- wstar$x
  ccx  <- ccx[, !colnames(ccx) %in% attributes(ccx)$factors, drop=F]
  ccx  <- ccx[,-1, drop=F]
  
  np <- ncol(ccMu)
  
  npage <- 1
  o   <- 1:np
  if(np > 16){
    npage <- ceiling(np/16)
    np    <- 16
  }
  
  mfrow <- .getPlotLayout(np)
  
  nbin <- 12
  ylimit <- c(0, max(ccMu) )
  
  for(m in 1:ncol(ccx)){   # loop over predictors
    
    xm  <- colnames(ccx)[m]
    xx  <- ccx[,m]
    atx <- quantile(xx,seq(0, 1, length=nbin))
    
    xlimit <- c( sum( c(.7, .3)*atx[1:2] ), sum( c(.3, .7)*atx[(nbin-1):nbin] ) )
    
    k   <- 0
    add <- F
    o   <- 1:np
    o   <- o[o <= 16]
    
    for(p in 1:npage){
      
      file <- paste('equilAbund_', xm, '_', p,'.pdf',sep='')
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
      
      npp <- ncol(ccMu) - k
      if(npp > np)npp <- np
      mfrow <- .getPlotLayout(np)
      par(mfrow=mfrow, bty='n', omi=c(.3,.3,0,0), mar=c(3,3,2,1), 
          tcl= tcl, mgp=mgp)
      
      for(j in o){
        
        yy  <- ccMu[,j]
        
        minbin <- 5
        xmids  <- 0
        while( length(xmids) == 1){
          minbin <- minbin - 1
          tt  <- .getBin(xx, yy, minbin = minbin)
          xbin <- tt$xbin
          xmids <- tt$xmids
        }
        
        c95 <- tapply( yy, list(bin = xbin), quantile, pnorm( c(-1.96, -1, 1, 1.96) ))
        ci <- matrix( unlist( c95 ), ncol = 4, byrow = T )
        rownames(ci) <- names(c95)
        
        xmids <- xmids[ as.numeric(names(c95)) ]
        
        .shadeInterval( xmids, loHi = ci[,c(1, 4)],col=specColor[j],PLOT=T,add=F,
                        xlab=' ',ylab=' ', xlim = xlimit,  
                        LOG=F, trans = .3)
        .shadeInterval( xmids, loHi = ci[,c(2, 3)],col=specColor[j], PLOT=T, add=T, 
                        trans = .3)
        mu <- tapply( yy, xbin, mean )
        lines(xmids,  mu, lty=2, lwd=2, col = specColor[j])
        
        
        k <- k + 1
        if(k > 26)k <- 1
        
        lab <- colnames(ccMu)[j]
        
        .plotLabel( lab,above=T )
      }
      mtext(xm, 1, outer=T)
      mtext('Equilibrium abundance', 2, outer=T)
      
      if(!SAVEPLOTS){
        readline('equilibrium abundance -- return to continue ')
      } else {
        dev.off()
      }
      o <- o + 16
      o <- o[o <= SO]
    }
  }
}

.rMVN <- function (nn, mu, sigma = NULL, sinv = NULL){
  
  # nn - no. samples from one mu vector or nrow(mu) for matrix
  
  if(!is.null(sigma)){
    m <- ncol(sigma)
  }else if(!is.null(sinv)){
    m <- ncol(sinv)
  }else{
    stop( '.rMNV requires either sigma or sinv' )
  }
  
  if(length(mu) > 1){
    if( !is.matrix(mu) ) mu <- matrix( mu, nn, length(mu) )  # mu is a vector of length m
    if( ncol(mu) == 1 & nn == 1 )  mu <- t(mu)
    if( length(mu) == m & nn > 1) mu <- matrix( mu, nn, length(mu), byrow=T )
  }
  
  if(is.null(sinv)){          # from sigma
    
    vv <- try(svd(sigma),T)
    
    if( inherits(vv,'try-error') ){
      ev <- eigen(sigma, symmetric = TRUE)
      rr <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
    } else {
      rr <- vv$v %*% (t(vv$u) * sqrt(vv$d)) 
    }
    
  }else{     # from sinv
    
    L  <- chol(sinv)  
    rr <- backsolve(t(L), diag(m), upper.tri = F) 
  }
  ps <- matrix(rnorm(nn * m), nn) %*% rr
  ps + mu 
}

gjamSimTime <- function(S, Q = 0, nsite, ntime = 50, termB, termR, termA, obsEffort = 100,
                        predPrey = NULL, zeroAlpha = NULL, PLOT = FALSE){
  
  # observations are counts ('DA' in gjam)
  # S           - no. species in Y
  # Q           - no. predictors in X
  # nsite       - no. groups/plots/sites
  # ntime       - mean length of time series
  # termB, termR, termA are logical for inclusion of immigration/emigration X%*%B, 
  #               DI growth VL, DD growth UA
  # obsEffort   - effort used in gjam
  # predPrey    - interactions assumed to be competition, unless predPrey = c(pred, prey), where
  #               prey and pred are the rows/columns in alpha 
  # zeroAlpha   - second column does not affect first column
  # PLOT        - logical to plot series
  
  cols <- colF( S )
  
  wdata <- NULL
  if(Q < 1)Q <- 1
  
  if(Q <= 1 | !termB){
    form <- '~ 1'
  }else{
    form <- paste0( 'x', c(2:Q), collapse = '+')
    form <- as.formula( paste(' ~', form ) )
  }
  
  if(!is.null(predPrey) & length(predPrey) == 2)   predPrey <- matrix(predPrey, nrow=1 )
  if(!is.null(zeroAlpha) & length(zeroAlpha) == 2)zeroAlpha <- matrix(zeroAlpha, nrow=1 )
  
  beta <- rhoTrue <- alphaTrue <- wstar <- NULL
  
  nt     <- 2 + rpois(nsite, ntime - 2)
  ntot   <- sum(nt)
  groups <- rep(1:nsite, nt)
  
  times  <- 1:nt[1]
  if(nsite > 1){
    for(k in 2:nsite)times <- c(times, c(1:nt[k]))
  }
  
  w <- matrix(1000,ntot,S)
  colnames(w) <- paste('s', 1:S, sep='')
  
  if(termB){ # environmental immigration/emigration
    
    bsd <- 2
    if(termR)bsd <- 2
    if(termA)bsd <- 10
    
    x     <- matrix( 0, ntot, Q)
    beta  <- matrix( rnorm(Q*S, 0, bsd), Q, S )
    beta[1,] <- rnorm(S, 0, 1)
    colnames(x) <- rownames(beta) <- paste('x',1:Q,sep='')
    rownames(beta)[1] <- 'intercept'
    colnames(beta) <- paste('s', 1:S, sep='')
    
  }else{
    
    x <- matrix(1, ntot, 1)
    colnames(x) <- 'intercept'
    
  }
  
  if(termR){ # growth rate rho
    
    gam <- runif(S, -.03, .03)       
    if(termB){
      gam <- runif(S, -.03, .05)
    }
    if(termA){
      gam <- runif(S, .01, .1)
    }
    rhoTrue <- gam
  }
  
  if(termA){ # alpha matrix
    
    intrWt <- 1.2
    if(S > 10)intrWt <- 10
    
    allPos <- F
    np     <- 0
    
    while(!allPos){
      daa <- runif(S,-1/80000,-1/150000)  # competition
      ir  <- runif(S^2, min(daa), 0)
      aa  <- matrix(ir, S, S)
      diag(aa) <- intrWt*daa
      
      if(!is.null(predPrey))aa[ predPrey ] <- -aa[ predPrey ] # prey has positive effect on pred
      if(!is.null(zeroAlpha))aa[ zeroAlpha ] <- 0             # no effect
      alphaTrue <- aa
      
      wstar <- -solve(crossprod(aa))%*%t(aa)%*%gam # carrying capacity
      if( all(wstar > 0) )allPos <- T
      np <- np + 1
    }
    print( paste(np, 'iterations for termA') )
  }
  
  # residual covariance
  sigma <- diag(1, S)
  XB    <- 0
  
  for(k in 1:nsite){
    
    ww     <- matrix(10, nt[k], S)
    ww[1,] <- runif(S, 200, 400)
    if(!termR & !termA)ww[1,] <- beta[1,]
    if(termR & termA)  ww[1,] <- runif(S, wstar*.1, wstar*2)
    
    xx <- matrix( rnorm(nt[k]*Q), nt[k], Q)
    xx[,1] <- 1
    if(termB)XB <- xx%*%beta
    
    for(t in 2:nt[k]){
      
      ep    <- .rMVN(1, 0, sigma)
      ww[t,] <- ww[t-1,] + t( ep )
      if(termB) ww[t,] <- ww[t,] + XB[t,] 
      if(termR) ww[t,] <- ww[t,] + ww[t-1,]*gam
      if(termA) ww[t,] <- ww[t,] + diag(ww[t-1,])%*%aa%*%t(ww[drop=F,t-1,]) 
      
      if(sum(ww[t,]) > 1000000)stop('try again')
    }
    wk <- which(groups == k)
    w[wk,] <- ww
    if(termB)x[wk,] <- xx
  }
  
  wkeep <- which( is.finite(rowSums(w)) & apply(w, 1, min) > 0 )
  w <- w[wkeep,]
  if(termB){
    x <- x[wkeep,]
    x <- x[,!colnames(x) == 'x']
  }
  times <- times[wkeep]
  groups <- groups[wkeep]
  
  y <- round( w*obsEffort )
  colnames(y) <- colnames(w)
  
  if(PLOT){
    ylim <- c(0, 1.5*max(y))
    xlim <- c(0, max(nt))
    par(bty='n', cex=1.5)
    xlab <- expression( paste("Time ", italic(t), sep=''))
    ylab <- expression( paste("Count ", italic(y[s][t]), sep=''))
    
    wk <- which(groups == 1)
    plot(times[wk], y[wk,1],type='l', xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab)
    
    for(k in 1:nsite){
      wk <- which(groups == k)
      for(s in 1:S){
        lines(times[wk], y[wk,s],col = cols[s],lwd=2)
        if(termA)abline(h=wstar[s]*obsEffort, lty=2, col= cols[s], lwd=2)
      }
    }
  }
  rownames(w) <- paste(groups, times, sep='-')
  
  if(termA){
    print( 'eigenvalues and carrying capacities' )
    wdata <- data.frame( eigenvalues = eigen(aa)$values, carryingCapacity = wstar )
    print( wdata )
  }
  
  rho <- matrix(0, S, S)
  diag(rho) <- rhoTrue
  
  trueValues <- list(beta = beta, rho = rho, alpha = alphaTrue, 
                     sigma = sigma, w = w)
  
  wt <- which( sapply(trueValues, length) > 0 )
  trueValues <- trueValues[ wt ]
  
  xdata <- data.frame( groups = groups, times = times, x, stringsAsFactors = F)
  
  list(xdata = xdata, ydata = y, edata = y*0 + obsEffort, formula = form, 
       groups = groups, times = times, trueValues = trueValues,
       wdata = wdata)
}

.cleanNames <- function(xx){
  
  xx <- .replaceString(xx,'-','')
  xx <- .replaceString(xx,'_','')
  xx <- .replaceString(xx,' ','')
  xx <- .replaceString(xx,"'",'')
  xx
}


.replaceString <- function(xx,now='_',new=' '){  #replace now string in vector with new
  
  ww <- grep(now,xx,fixed=T)
  if(length(ww) == 0)return(xx)
  
  for(k in ww){
    s  <- unlist( strsplit(xx[k],now,fixed=T) )
    ss <- s[1]
    if(length(s) == 1)ss <- paste( ss,new,sep='')
    if(length(s) > 1)for(kk in 2:length(s)) ss <- paste( ss,s[kk],sep=new)
    xx[k] <- ss
  }
  xx
}

gjamTimePrior <- function( xdata, ydata, edata, priorList, minSign = 5, 
                           betaMax = 30, rhoMax = 1 ){
  
  # minSign - minimum number of co-occurrences to estimate alpha
  
  bp <- lp <- ap <- NULL
  formulaBeta <- formulaRho <- alphaSign <- NULL
  termB <- termR <- termA <- FALSE
  
  betaPrior <- rhoPrior <- alphaPrior <- NULL
  
  for(k in 1:length(priorList))assign( names(priorList)[k], priorList[[k]] )
  
  S <- ncol(ydata)
  w <- ydata/edata
  n <- nrow(w)
  
  colnames(w) <- .cleanNames( colnames(w ) )
  
  
  if(!is.null(betaPrior))termB <- TRUE
  if(!is.null(rhoPrior)) termR <- TRUE
  if(!is.null(alphaSign))termA <- TRUE
  
  if( 'formulaBeta' %in% names(priorList) )termB <- TRUE
  
  if(termB){  # only XB
    
    timeZero <- grep('-0', rownames(ydata)) # rownames from gjamFillMissing
    timeLast <- c( timeZero - 1, nrow(ydata))[-1]
    
    tmp <- model.frame(formulaBeta, data=xdata, na.action=NULL) # standardized design
    
    x   <- model.matrix(formulaBeta, data=tmp)
    
    # do not center intercept or factors
    wf <- names( which( attr(attributes(tmp)$terms, 'dataClasses') == 'factor' ) )
    wf <- which( startsWith(colnames(x), wf ) )
    wf <- c(1, wf)
    wp <- c(1:ncol(x))[-wf]
    
    if(length(wp) > 0){
      xm  <- colMeans(x, na.rm=T)
      xs  <- apply(x, 2, sd, na.rm=T)
      x[,wp]   <- t( (t(x[,wp]) - xm[wp])/xs[wp]  )
    }
    
    kname <- as.vector( outer( colnames(x), colnames(w), FUN = paste, sep='-') )
    
    xnames <- attributes(x)$dimnames[[2]]
    blo <- matrix(Inf, ncol(x), ncol(ydata))
    rownames(blo) <- colnames(x)
    colnames(blo) <- colnames(ydata)
    bhi <- -blo
    
    # variance in dy
    gr  <- unique(xdata$groups)
    nr  <- length(gr)
    
    sumw <- sumw2 <- rep(0, ncol(w))
    sumx <- sumx2 <- rep(0, ncol(x))
    
    for(j in gr){
      wj <- which(xdata[,'groups'] == j & is.finite(rowSums(x)) &
                    is.finite(rowSums(w)))
      wj <- wj[ !wj %in% timeZero ]
      
      if(length(wj) <= ncol(x)) next
      dw <- w[ wj[ -1 ], ] - w[ wj[ -length(wj) ], ]
      xj <- x[wj,][-1,]
      
      wv <- which( apply(xj, 2, var) > 0 )
      if(length(wv) == 0)next
      bb <- solve( crossprod(xj[,wv]) )%*%crossprod(xj[,wv], dw)
      
      bx <- blo
      bx[rownames(bb),] <- bb
      blo[bx < blo] <- bx[bx < blo]
      
      bx <- bhi
      bx[rownames(bb),] <- bb
      bhi[bx > bhi] <- bx[bx > bhi]
      
      sumw  <- sumw + colSums(dw)
      sumw2 <- sumw2 + colSums(dw^2)
      sumx  <- sumx + colSums(xj)
      sumx2 <- sumx2 + colSums(xj^2)
    }
    
    # variables that are fixed within a group
    wsd <- sqrt( sumw2/nrow(x) - (sumw/nrow(x))^2 )
    xsd <- sqrt( sumx2/nrow(x) - (sumx/nrow(x))^2 )
    xsd[1] <- 1
    
    wx  <- matrix(wsd, length(xsd), length(wsd), byrow=T)/matrix(xsd, length(xsd), length(wsd) )
    
    blo[ blo == 0 | !is.finite(blo) ] <- -wx[ blo == 0 | !is.finite(blo) ]
    bhi[ bhi == 0 | !is.finite(bhi) ] <-  wx[ bhi == 0 | !is.finite(bhi) ]
    
    bl <- blo - 2*abs(blo)
    bh <- bhi + 2*abs(bhi)
    rownames(bl)[1] <- rownames(bh)[1] <- 'intercept'
    
    blo[ blo < -betaMax ] <- -betaMax
    bhi[ bhi > betaMax ]  <- betaMax
    
    if( is.list(betaPrior) ){
      gg  <- gjamPriorTemplate(formulaBeta, xdata, ydata = ydata, 
                               lo = betaPrior$lo, hi = betaPrior$hi)
      blo <- gg[[1]]
      bhi <- gg[[2]]
      
      blo[ !rownames(blo)%in% names(betaPrior$lo), colnames(blo) ] <- 
        bl[ !rownames(blo)%in% names(betaPrior$lo),  ]
      
      bhi[ !rownames(bhi)%in% names(betaPrior$hi), colnames(bhi) ] <- 
        bh[ !rownames(bhi)%in% names(betaPrior$hi),  ]
      
    }else{
      blo <- bl
      attr(blo,'formula') <- formulaBeta
      bhi <- bh
    }
    
    
    bp  <- list(lo = blo, hi = bhi )
  }
  
  if(termR){  # rho
    if( !is.list(rhoPrior) ){
      lp <- NULL
    }else{
      formulaRho <- priorList$formulaRho
      if(is.null(formulaRho)){
        cat('\nformulaRho omitted from priorList, used ~1\n')
      }
      gg <- gjamPriorTemplate(formulaRho, xdata, ydata = ydata, 
                              lo = rhoPrior$lo, hi = rhoPrior$hi)
      lp <- list(lo = gg[[1]], hi = gg[[2]])
    }
    
    lp$lo[ lp$lo < -betaMax ] <- -betaMax
    lp$hi[ lp$hi > betaMax ]  <- betaMax
  }
  
  if(termA){  # alpha
    
    if( !is.list(rhoPrior) ){
      if(is.null(lp))stop(' must have rhoPrior if there is alphaSign' )
      ap <- NULL
    }else{
      
      #must co-occur
      yind <- ydata
      rownames(alphaSign) <- colnames(alphaSign) <- 
        colnames(yind) <- .cleanNames( colnames(alphaSign ) )
      
      yind[yind > 1] <- 1
      yind <- crossprod(yind)
      
      yind[yind < minSign] <- 0  # minimum co-occurence
      yind[yind > 1] <- 1
      
      alphaSign <- alphaSign*yind[rownames(alphaSign),colnames(alphaSign)]
      
      rho <- (lp$lo['intercept',] + lp$hi['intercept',])/2 + .01
      
      timeZero <- which(xdata$times == 0)
      timeLast <- (timeZero - 1)[-1]
      wt <- unique( c(timeZero, timeLast) )

      wmu    <- colMeans( w, na.rm=T )
      wdelta <- apply( w, 2, diff )  # pop rate
      wdelta[!is.finite(wdelta) ] <- 0
      wdelta <- wdelta[-wt,]
      
      ##############
   #   rmat <- matrix(rho, n-1, S, byrow=T) 
   #   wdelta[ wdelta >  rmat ] <- rmat[ wdelta > rmat ]
   #   wi <- 1/w[-n,]
      
   #   wi[ !is.finite(wi) ] <- 0
   #   wd <- wdelta*wi
   #   wd[ !is.finite(wd) ] <- 0
   #   wd[ wd >  rmat ] <- rmat[ wd > rmat ]
      
      ww <- n/crossprod(w)
      wrange <- apply(wdelta, 2, quantile, c(.05,.95), na.rm=T)
      wrange[1, wrange[1,] >= 0 ] <- mean( wrange[1, wrange[1,] < 0 ] )
      wrange[2, wrange[2,] <= 0 ] <- mean( wrange[2, wrange[2,] > 0 ] )
      
      wlo <- matrix(wrange[1,], S, S, byrow=T)*ww   # E[dw]/E[w_s * w_s']
      whi <- matrix(wrange[2,], S, S, byrow=T)*ww   # E[dw]/E[w_s * w_s']
      
      rw  <- apply( w, 2, quantile, .8, na.rm = T)
      rlh <- matrix( rho/rw, S, S, byrow= T) # rho/w_s
      
      whi[ whi > rlh ]  <- rlh[ whi > rlh ] 
      wlo[ wlo < -rlh ] <- -rlh[ wlo < -rlh ] 
      
      alo <- ahi <- wlo*0
      ww  <- which(alphaSign < 0)
      
      scale <- 100
      if(max(whi) > .1)scale <- 1
      if(length(ww) > 0){
        alo[ww] <- scale*wlo[ww]
        ahi[ww] <- 0
      }
      ww <- which(alphaSign > 0)
      if(length(ww) > 0){
        alo[ww] <- 0
        ahi[ww] <- scale*whi[ww]
      }
      ww <- which(alphaSign == 0)
      if(length(ww) > 0)alo[ww] <- ahi[ww] <- NA
      
      ap <- list(lo = alo, hi = ahi)
    }
  }
  
  list(betaPrior = bp, rhoPrior = lp, alphaPrior = ap, formulaBeta = formulaBeta,
       formulaRho = formulaRho)
}

foodWebDiagram <- function(S, guildList = NULL, predPrey = NULL, zeroAlpha = NULL,
                           intraComp = 1:S, label = NULL, PLOT = TRUE, layout = 'rr'){
  
  # S - no. species
  # default interaction is negative, arrows only for negative interactions
  # guildList - overrides default, only members of the same guild compete
  # predPrey  - matrix with 2 columns, second column is prey of first column
  # zeroAlpha - matrix with 2 columns, second column does not affect first column
  # intraComp - needed for intraspecific comp if guildList is specified
  # layout can be 'tree', 'rr', ...
  
  require( DiagrammeR )
  
  pp <- numeric(0)
  
  fromTo <- as.matrix( expand.grid( c(1:S), c(1:S) ) )
  ft <- .pasteCols( fromTo )
  
  qq <- numeric(0)
  if( !is.null(guildList) ){
    
    fromTo <- numeric(0)
    for(k in 1:length(guildList)){
      ft <- as.matrix(expand.grid(guildList[[k]], guildList[[k]]))
      fromTo <- rbind(fromTo, ft)
    }
    if(!is.null(intraComp)){
      ft <- cbind(1:S, 1:S)
      fromTo <- rbind(fromTo, ft)
    }
    ft <- .pasteCols( fromTo )
    ww <- which(!duplicated(ft))
    ft <- ft[ww]
    fromTo <- fromTo[ww,]
    zeroAlpha <- NULL
  }
  
  if(!is.null(predPrey)){
    pp <- .pasteCols( predPrey[drop=F,,c(2,1)] )
    qq <- .pasteCols( predPrey)
    #   fromTo <- fromTo[ !ft %in% pp, ]
    fromTo <- rbind(fromTo, predPrey, predPrey[drop=F,,c(2,1)])
    ft <- .pasteCols( fromTo )
    ww <- which(!duplicated(ft))
    fromTo <- fromTo[ww,]
    ft <- ft[ww]
  }
  
  if(!is.null(zeroAlpha)){
    za <- .pasteCols( zeroAlpha[drop=F,,c(1,2)])
    
    fromTo <- fromTo[ !ft %in% za, ]
    ft <- ft[ !ft %in% za ]
  }
  
  if( length(qq) > 0 ){
    ft <- .pasteCols( fromTo )
    fromTo <- rbind( fromTo[ft %in% qq,], fromTo[!ft %in% qq, ] )
  }
  ft <- .pasteCols( fromTo )
  pp <- pp[ pp %in% ft ]
  
  if(!PLOT){
    return( fromTo )
  }else{
    ecol  <- rep( 'tan', length(ft) )
    #  ncol  <- rep( 'tan', S )
    ncol  <- colF(S)
    shape <- rep( 'rectangle', S)
    
    if( length(qq) > 0 ){
      
      wt <- which( ft %in% qq )
      ecol[ wt ] <- 'brown'
      
    }
    if( length(pp) > 0 ){
      
      wt <- which( ft %in% pp )
      ecol[ wt ] <- 'blue'
    }
    
    if( is.null(layout) )layout <- "tree"
    if(is.null(label))label <- paste('s', 1:S, sep='')
    nodes <- create_node_df(n = S, label = label, style = "filled", color = ncol, shape = shape)
    edges <- create_edge_df(from = fromTo[,1], to = fromTo[,2], color = ecol )
    graph <- create_graph(nodes_df = nodes, edges_df = edges)
    render_graph( graph,  layout = layout )
  }
}


