
#' Two-sample Q-Q plot
#'
#' Plot the quartiles of two distributions against each other: graphical test of similarity of fit
#' @param data1 First vector of angles to be tested.
#' @param data2 Second vector of angles to be tested.
#' @export
#' @examples
#' two.sample.QQ(q.a, q.b)
two.sample.QQ <- function(data1, data2) {
    
    n1 <- length(data1) 
    n2 <- length(data2) 
    nmin <- min(n1,n2) 
    nmax <- max(n1,n2)
    
    if (n2 < n1) {dataref <- data2 ; dataoth <- data1
    } else {dataref <- data1 ; dataoth <- data2}
    
    zref <- sin(0.5*(dataref - median.circular(dataref))) 
    szref <- sort(zref)
    
    zoth <- sin(0.5*(dataoth-median.circular(dataoth)))
    szoth <- sort(zoth)
    
    koth <- 0
    szothred <- 0
    szreffin <- 0
    
    for (k in 1:nmin) {
        koth[k] <- 1+nmax*(k-0.5)/nmin
        szothred[k] <- szoth[koth[k]] 
        szreffin[k] <- szref[k]
    }
    plot(szreffin, szothred, pch=16, xlim=c(-1,1), ylim=c(-1,1), xlab = "Smaller sample", ylab = "Larger sample")
    xlim <- c(-1,1) ; ylim <- c(-1,1) 
    lines(xlim, ylim, lwd=2, col = "lightseagreen")
}


#' Watson's two-sample test of common mean
#'
#' Watson's large-sample nonparametric test for the similarity of a common mean for two (or more) distributions. Does not assume that the distributions tested have a common concentration or shape.
#' @param samples A list containing the vectors of angles (in radians) to be tested.
#' @return List containing the test statistic, p-value and ratio of maximum to minimum sample dispersion.
#' @export
#' @examples
#' comp <- watson.common.mean.test(list(q.4.a, q.4.b))
watson.common.mean.test <- function(samples) {
    
    data <- unlist(samples)
    N <- length(data)
    g <- length(samples)
    sample.sizes <- 0
    for (i in 1:g) {sample.sizes[i] <- length(samples[[i]])}
    
    size.csum <- cumsum(sample.sizes) 
    
    delhat <- 0 
    tbar <- 0
    
    for (k in 1:g) {
        sample <- samples[[k]]
        
        tm1 <- trigonometric.moment(sample, p=1)
        tm2 <- trigonometric.moment(sample, p=2)
        
        Rbar1 <- tm1$rho
        Rbar2 <- tm2$rho 
        tbar[k] <- tm1$mu %% (2*pi)
        
        delhat[k] <- (1-Rbar2)/(2*Rbar1*Rbar1)
    }
    
    dhatmax <- max(delhat) 
    dhatmin <- min(delhat)
    
    if (dhatmax/dhatmin <= 4) { # use P procedure
        
        CP <- 0
        SP <- 0
        dhat0 <- 0
        
        for (k in 1:g) {
            CP <- CP + sample.sizes[k]*cos(tbar[k])
            SP <- SP + sample.sizes[k]*sin(tbar[k])
            dhat0 <- dhat0 + sample.sizes[k]*delhat[k] 
        }
        
        dhat0 <- dhat0/N
        RP <- sqrt(CP*CP+SP*SP)
        
        Yg <- 2*(N-RP)/dhat0
    } else {
        
        CM <- 0 
        SM <- 0
        Yg <- 0
        
        for (k in 1:g) {
            CM <- CM + (sample.sizes[k]*cos(tbar[k])/delhat[k])
            SM <- SM + (sample.sizes[k]*sin(tbar[k])/delhat[k])
            Yg <- Yg + (sample.sizes[k]/delhat[k]) 
        }
        RM <- sqrt(CM*CM+SM*SM)
        Yg <- 2*(Yg-RM)
    }
    
    pval = pchisq(Yg, g-1, lower.tail = F)
    list(Y.g = Yg, p.val = pval, disp.ratio = dhatmax/dhatmin)
}


#' Bootstrap version of Watson's two-sample test of common mean
#'
#' Bootstrap of Watson's nonparametric test for the similarity of a common mean for two (or more) distributions. Does not assume that the distributions tested have a common concentration or shape.
#' @param samples A list containing the vectors of angles (in radians) to be tested.
#' @param B Number of bootstrap samples to use to obtain the estimated p-value. Default is 9999.
#' @param show.progress Boolean indicating whether or not to display a progress bar as the bootstrap is run.
#' @return List containing the test statistic, p-value and ratio of maximum to minimum sample dispersion.
#' @export
#' @examples
#' comp <- watson.mean.test.boot(list(q.4.a, q.4.b), B = 99)
watson.mean.test.boot <- function(samples, B = 9999, show.progress = T) {
    
    data <- unlist(samples)
    N <- length(data)
    g <- length(samples)
    sample.sizes <- c()
    g.id <- c()
    for (i in 1:g) {sample.sizes[i] <- length(samples[[i]])
                    g.id <- c(g.id, rep(i, sample.sizes[i]))}
    
    
    delhat <- c() 
    tbar <- c() 
    centdat <- c()
    
    for (k in 1:g) {
        sample <- samples[[k]] 
        
        tm1 <- trigonometric.moment(sample, p=1)        
        tm2 <- trigonometric.moment(sample, p=2)
        
        Rbar1 <- tm1$rho
        Rbar2 <- tm2$rho 
        tbar[k] <- tm1$mu %% (2*pi)
        
        delhat[k] <- (1-Rbar2)/(2*Rbar1*Rbar1)
        
        centsamp <- circular(sample - tbar[k])    
        centdat <- circular(c(centdat,centsamp))
    }
    dhatmax <- max(delhat)
    dhatmin <- min(delhat)
    
    if (dhatmax/dhatmin <= 4) {
        PorM <- 1
        CP <- 0
        SP <- 0
        dhat0 <- 0
        
        for (k in 1:g) {
            CP <- CP + sample.sizes[k]*cos(tbar[k])
            SP <- SP + sample.sizes[k]*sin(tbar[k])
            dhat0 <- dhat0 + sample.sizes[k]*delhat[k] 
        }
        dhat0 <- dhat0/N
        RP <- sqrt(CP*CP+SP*SP)
        Yg <- 2*(N-RP)/dhat0
        
    } else {
        PorM <- 0 
        CM <- 0 
        SM <- 0 
        Yg <- 0
        
        for (k in 1:g) {
            CM <- CM + (sample.sizes[k]*cos(tbar[k])/delhat[k])
            SM <- SM + (sample.sizes[k]*sin(tbar[k])/delhat[k])
            Yg <- Yg + (sample.sizes[k]/delhat[k]) 
        }
        
        RM <- sqrt(CM*CM+SM*SM)
        Yg <- 2*(Yg-RM)
    }
    YgObs <- Yg
    nxtrm <- 1
    
    if (show.progress) {pb <- txtProgressBar(min = 0, max = B, style = 3)}
    for (b in 1:B) {
        for (k in 1:g) {
            centsamp <- centdat[g.id == k]         
            bootsamp <- sample(centsamp, size=sample.sizes[k], replace=TRUE)
            
            tm1 <- trigonometric.moment(bootsamp, p=1)
            tm2 <- trigonometric.moment(bootsamp, p=2)
            
            Rbar1 <- tm1$rho
            Rbar2 <- tm2$rho 
            tbar[k] <- tm1$mu
            
            delhat[k] <- (1-Rbar2)/(2*Rbar1*Rbar1)
        }
        
        if (PorM == 1) {
            CP <- 0 
            SP <- 0 
            dhat0 <- 0
            
            for (k in 1:g) {
                CP <- CP + sample.sizes[k]*cos(tbar[k])
                SP <- SP + sample.sizes[k]*sin(tbar[k])
                dhat0 <- dhat0 + sample.sizes[k]*delhat[k] 
            }
            
            dhat0 <- dhat0/N
            RP <- sqrt(CP*CP+SP*SP)
            Yg <- 2*(N-RP)/dhat0
        } else {
            CM <- 0 ; SM <- 0 ; Yg <- 0
            
            for (k in 1:g) {
                CM <- CM + (sample.sizes[k]*cos(tbar[k])/delhat[k])
                SM <- SM + (sample.sizes[k]*sin(tbar[k])/delhat[k])
                Yg <- Yg + (sample.sizes[k]/delhat[k]) 
            }
            RM <- sqrt(CM*CM+SM*SM)
            Yg <- 2*(Yg-RM)
        }
        YgBoot <- Yg
        
        if (YgBoot >= YgObs) {nxtrm <- nxtrm+1}
        if (show.progress) {setTxtProgressBar(pb, b)}
    }
    if (show.progress) {close(pb)}
    pval <- nxtrm/(B+1)
    list(statistic = YgObs, p.value = pval, disp.ratio = dhatmax/dhatmin)
}


#' Wallraff's two-sample test of common concentration
#'
#' Wallraff's nonparametric test of circular homoscedasticity (nonparametric test for similar concentration of multiple samples).
#' @param samples A list containing the vectors of angles (in radians) to be tested.
#' @return List containing the test statistic, p-value and ratio of maximum to minimum sample dispersion.
#' @export
#' @examples
#' comp <- watson.common.mean.test(list(q.4.a, q.4.b))
wallraff.concentration.test <- function(samples) {
    
    data <- unlist(samples)
    g <- length(samples)
    g.id <- c()
    for (i in 1:g) {
        g.id <- c(g.id, rep(i, length(samples[[i]])))
    }
    
    N <- length(data)
    #    sample.sizescsum <- cumsum(sample.sizes) 
    
    tbar <- circular(0) 
    distdat <- c()
    for (k in 1:g) {
        
        #        dist <- 0 
        sample <- samples[[k]]
        
        tm1 <- trigonometric.moment(sample, p=1) 
        tbar[k] <- tm1$mu %% (2*pi)
        
        dist <- pi-abs(pi-abs(sample-tbar[k]))
        distdat <- c(distdat, dist)
    }
    
    TestRes <- kruskal.test(distdat, g = g.id)
    
    list(p.val = TestRes$p.value, result = TestRes)
} 


#' Cosine & sine rank scores
#'
#' Support function: calculate cosine & sine rank scores for a sample of angles
#' @param samples A list containing the vectors of angles (in radians) to be evaluated.
#' @return List containing the cosine and sine rank scores.
#' @export
#' @examples
#' cs.scores <- cs.unif.scores(list(q.4.a, q.4.b))
cs.unif.scores <- function(samples) {
    
    data <- unlist(samples)
    N <- length(data)
    ranks <- rank(data, ties.method="random")
    cos.u.scores <- cos(ranks*2*pi/N)
    sin.u.scores <- sin(ranks*2*pi/N)
    
    list(cos.scores = cos.u.scores, sin.scores = sin.u.scores)    
}


#' Test for common distribution among multiple samples
#'
#' Large-sample Mardia-Wheeler-Watson test for common distribution.
#' @param cs.scores A list containing sine and cosine rank scores, output from the data of interest by \code{\link{cs.unif.scores}}.
#' @param sample.sizes A vector defining the sizes of the samples that were used to produce the \code{cs.scores} object.
#' @return List containing the test statistic and p-value.
#' @export
#' @examples
#' cd.res <- mww.common.dist.LS(cs.unif.scores(list(q.4.a, q.4.b)), c(length(q.4.a), length(q.4.b)))
mww.common.dist.LS <- function(cs.scores, sample.sizes) {
    
    N <- sum(sample.sizes)
    
    g <- length(sample.sizes)
    g.id <- c()
    for (i in 1:g) {
        g.id <- c(g.id, rep(i, sample.sizes[i]))
    }
    
    Wg <- 0
    
    for (k in 1:g) {
        cos.u.scores.k <- cs.scores$cos.scores[g.id == k] 
        sin.u.scores.k <- cs.scores$sin.scores[g.id == k] 
        
        sum.cos.k.sq <- (sum(cos.u.scores.k))^2
        sum.sin.k.sq <- (sum(sin.u.scores.k))^2
        
        Wg <- Wg+(sum.cos.k.sq+sum.sin.k.sq)/length(g.id[g.id == k])
    }
    Wg <- 2*Wg 
    p.val <- pchisq(Wg, 2*(g-1), lower.tail = F)
    
    list(statistic = Wg, p.val = p.val)    
}


#' Permutation test for common distribution among multiple samples
#'
#' Permutation-based Mardia-Wheeler-Watson test for common distribution, to be used when any of the samples has less than 10 angles.
#' @param cs.scores A list containing sine and cosine rank scores, output from the data of interest by \code{\link{cs.unif.scores}}.
#' @param sample.sizes A vector defining the sizes of the samples that were used to produce the \code{cs.scores} object.
#' @param NR Number of permutation samples to use to estimate the p-value. Default is 9999.
#' @param show.progress Boolean indicating whether or not to display a progress bar as the bootstrap is run.
#' @return p-value for the test.
#' @export
#' @examples
#' cd.res <- mww.common.dist.rand(list(q.4.a, q.4.b))
mww.common.dist.rand <- function(cs.scores, sample.sizes, NR = 9999, show.progress = T) {
    
    g <- length(sample.sizes)
    g.id <- c()
    for (i in 1:g) {
        g.id <- c(g.id, rep(i, sample.sizes[i]))
    }
    N <- sum(sample.sizes)
    
    cos.u.scores <- cs.scores$cos.scores
    sin.u.scores <- cs.scores$sin.scores
    
    WgObs <- mww.common.dist.LS(cs.scores, sample.sizes)$statistic
    nxtrm <- 1
    
    if (show.progress) {pb <- txtProgressBar(min = 0, max = NR, style = 3)}
    for (r in 1:NR) {
        cos.u.scores.rand <- c() 
        sin.u.scores.rand <- c()
        
        rand.ind <- sample(1:N, replace = F)
        
        for (k in 1:g) {            
            cos.u.scores.rand <- c(cos.u.scores.rand, cos.u.scores[rand.ind[g.id == k]])
            sin.u.scores.rand <- c(sin.u.scores.rand, sin.u.scores[rand.ind[g.id == k]])
        }
        
        cs.scores.rand <- list(cos.scores = cos.u.scores.rand, sin.scores = sin.u.scores.rand)
        
        WgRand <-  mww.common.dist.LS(cs.scores.rand, sample.sizes)$statistic
        
        if (WgRand >= WgObs) { nxtrm <- nxtrm+1 }
        if (show.progress) {setTxtProgressBar(pb, r)}
    }
    if (show.progress) {close(pb)}
    nxtrm/(NR+1)
}


#' Permutation test for common distribution between two samples
#'
#' Extension of \code{\link{watson.two.test}}, using permutation to obtain the p-value.
#' @param data1 First vector of angles to be tested.
#' @param data2 Second vector of angles to be tested.
#' @param NR Number of permutation samples to use to estimate the p-value. Default is 9999.
#' @param show.progress Boolean indicating whether or not to display a progress bar as the bootstrap is run.
#' @return List containing the observed value of the test statistic and the p-value.
#' @seealso \code{\link{watson.two.test}}
#' @export
#' @examples
#' cd.res <- watson.two.test.rand(list(q.4.a, q.4.b))
watson.two.test.rand <- function(data1, data2, NR = 9999, show.progress = T){
    
    wats.obs <- watson.two.test(data1, data2)$statistic 
    nxtrm <- 1
    
    n1 <- length(data1)
    n2 <- length(data2)
    N <- n1+n2
    
    combined <- c(data1, data2)
    
    if (show.progress) {pb <- txtProgressBar(min = 0, max = NR, style = 3)}
    
    for (r in 1:NR) {
        rs <- sample(combined)
        rs1 <- rs[1:n1]
        rs2 <- rs[(n1+1):N]
        
        wats.rand <- watson.two.test(rs1, rs2)$statistic
        
        if (wats.rand >= wats.obs) {nxtrm <- nxtrm+1}
        if (show.progress) {setTxtProgressBar(pb, r)}
    }
    
    pval <- nxtrm/(NR+1) 
    if (show.progress) {close(pb)}
    list(observed = wats.obs, p.val = pval)
}
