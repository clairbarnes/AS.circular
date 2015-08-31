#' Jones-Pewsey normalising constant
#'
#' Support function to calculate the normalising constant for the Jones-Pewsey distribution, based on specified values of \code{kappa} and \code{psi}.
#' @param kappa Concentration parameter.
#' @param psi Shape parameter.
#' @return A constant value to be passed to further functions.
#' @export
#' @examples
#' nc <- JP.NCon(kappa = 1, psi = -3)
JP.NCon <- function(kappa, psi){
    if (kappa < 0.001) {ncon <- 1/(2*pi) ; return(ncon)} 
    else {
        eps <- 10*.Machine$double.eps
        if (abs(psi) <= eps) {ncon <- 1/(2*pi*I.0(kappa)) ; return(ncon)}
        
        
        else {
            intgrnd <- function(x){ (cosh(kappa*psi)+sinh(kappa*psi)*cos(x))**(1/psi) }
            
            ncon <- 1/integrate(intgrnd, lower=-pi, upper=pi)$value
            return(ncon) } }
}


#' Jones-Pewsey probability density function
#'
#' Calculate the pdf of the Jones-Pewsey distribution at given points, based on specified values of \code{mu}, \code{kappa} and \code{psi}.
#' @param theta Angle, in radians, at which the density is to be calculated. Can be either numeric or vector.
#' @param mu Mean direction parameter.
#' @param kappa Concentration parameter.
#' @param psi Shape parameter.
#' @return The density function evaluated at \code{theta}: either a numeric or a vector, depending on the input.
#' @export
#' @examples
#' pdf <- JP.pdf(theta = c(k.2), mu = 3, kappa = 1, psi = -3, JP.NCon(kappa = 1, psi = -3))
JP.pdf <- function(theta, mu, kappa, psi, ncon){
    
    if (kappa < 0.001) {pdfval <- 1/(2*pi) ; return(pdfval)}
    else {
        eps <- 10*.Machine$double.eps
        if (abs(psi) <= eps) {
            pdfval <- ncon*exp(kappa*cos(theta-mu)) ; return(pdfval) }
        
        
        else { 
            pdfval <- (cosh(kappa*psi)+sinh(kappa*psi)*cos(theta-mu))**(1/psi)
            pdfval <- ncon*pdfval ; return(pdfval) } }
    
}


#' Jones-Pewsey distribution function
#'
#' Calculate the cdf of the Jones-Pewsey distribution at given points, based on specified values of \code{mu}, \code{kappa} and \code{psi}.
#' @param theta Angle, in radians, at which the distribution is to be calculated. Can be either numeric or vector.
#' @param mu Mean direction parameter.
#' @param kappa Concentration parameter.
#' @param psi Shape parameter.
#' @return The distribution function evaluated at \code{theta}: either a numeric or a vector, depending on the input.
#' @export
#' @examples
#' cdf <- JP.df(theta = c(k.2), mu = 3, kappa = 1, psi = -3, JP.NCon(kappa = 1, psi = -3))
JP.df <- function(theta, mu, kappa, psi, ncon) {
    
    eps <- 10*.Machine$double.eps
    if (theta <= eps) {dfval <- 0 ; return(dfval)}
    
    else 
        if (theta >= 2*pi-eps) {dfval <- 1 ; return(dfval)} else
            if (kappa < 0.001) {dfval <- theta/(2*pi) ; return(dfval)}
    else {
        if (abs(psi) <= eps) {
            vMPDF <- function(x){ ncon*exp(kappa*cos(x-mu)) }
            dfval <- integrate(vMPDF, lower=0, upper=theta)$value
            return(dfval) }
        
        
        else { 
            dfval <- integrate(JP.pdf, mu=mu, kappa=kappa, psi=psi, ncon = ncon, lower=0, upper=theta)$value
            return(dfval) }
    }
}


#' Jones-Pewsey quantile function
#'
#' Evaluate the quartile of the Jones-Pewsey distribution at given points, based on specified values of \code{mu}, \code{kappa} and \code{psi}.
#' @param u Quantile to be evaluated.
#' @param mu Mean direction parameter.
#' @param kappa Concentration parameter.
#' @param psi Shape parameter.
#' @return The quantile function evaluated at \code{u}: either a numeric or a vector, depending on the input.
#' @export
#' @examples
#' q <- JP.qf(u = 0.9, mu = 3, kappa = 1, psi = -3, JP.NCon(kappa = 1, psi = -3))
JP.qf <- function(u, mu, kappa, psi, ncon) {
    
    eps <- 10*.Machine$double.eps
    if (u <= eps) {theta <- 0 ; return(theta)}
    
    else 
        if (u >= 1-eps) {theta <- 2*pi-eps ; return(theta)} else
            if (kappa < 0.001) {theta <- u*2*pi ; return(theta)}
    else {
        roottol <- .Machine$double.eps**(0.6)
        qzero <- function(x) {
            y <- JP.df(x, mu, kappa, psi, ncon) - u ; return(y) }
        res <- uniroot(qzero, lower=0, upper=2*pi-eps, tol=roottol)
        theta <- res$root ; return(theta) }
}


#' Jones-Pewsey simulation function
#'
#' Simulate a number of angles from the Jones-Pewsey distribution specified by \code{mu}, \code{kappa} and \code{psi}.
#' @param n Number of angles to be simulated from the distribution.
#' @param mu Mean direction parameter.
#' @param kappa Concentration parameter.
#' @param psi Shape parameter.
#' @return A vector of n simulated angles from the specified distribution.
#' @export
#' @examples
#' sample <- JP.sim(n = 100, mu = 3, kappa = 1, psi = -3, JP.NCon(kappa = 1, psi = -3))
JP.sim <- function(n, mu, kappa, psi, ncon) {
    
    fmax <- JP.pdf(mu, mu, kappa, psi, ncon) ; theta <- 0
    for (j in 1:n) {
        stopgo <- 0
        while (stopgo == 0) {
            u1 <- runif(1, 0, 2*pi)
            pdfu1 <- JP.pdf(u1, mu, kappa, psi, ncon)
            u2 <- runif(1, 0, fmax)
            if (u2 <= pdfu1) { theta[j] <- u1 ; stopgo <- 1 }
        } }
    return(theta)
}


#================================================================================


#' Maximum likelihood estimator for Jones-Pewsey parameters
#'
#' Obtain maximum likelihood estimates of Jones-Pewsey parameters for a given data set.
#' @param data Vector of angles over which maximum likelihood estimation is to be performed.
#' @return A list containing the maximum likelihood achieved, estimates of \code{mu}, \code{kappa} and \code{psi}, and the Hessian matrix of the optimised parameters.
#' @export
#' @examples
#' mu <- JP.mle(k.2)$mu
JP.mle <- function(data) {
    n <- length(data)
    s <- sum(sin(data))
    c <- sum(cos(data))
    muvM <- atan2(s,c) %% (2*pi)
    kapvM <- A1inv(sqrt(s*s+c*c)/n)
    
    JPnll <- function(p){
        mu <- p[1] ; kappa <- p[2] ; psi <- p[3] ; parlim <- abs(kappa*psi)
        if (parlim > 10) {
            y <- 9999.0
            return(y)
        } else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(data, mu, kappa, psi, ncon)))
            return (y)
        }
    }
    
    out <- optim(par=c(muvM, kapvM, 0.001), fn=JPnll, gr = NULL, method = "L-BFGS-B", lower = c(muvM-pi, 0, -Inf), upper = c(muvM+pi, Inf, Inf), hessian = TRUE)
    
    list(maxll = -out$value,
         mu = out$par[1] %% (2*pi),
         kappa = out$par[2],
         psi = out$par[3],
         HessMat = out$hessian)
}


#' Normal-theory confidence interval for Jones-Pewsey parameter estimators
#'
#' Obtain nominal 100(1-\code{alpha})% confidence intervals for Jones-Pewsey estimators, by solving the Hessian matrix produced by \code{JP.mle}.
#' @param jp.ests List (as output by \code{JP.mle}) containing named estimates of \code{mu}, \code{kappa} and \code{psi}, and the Hessian matrix \code{HessMat} of the optimised parameters.
#' @param alpha Significance level of confidence interval to be obtained. Default is 0.05 (95% confidence interval).
#' @return A list containing estimates of \code{mu}, \code{kappa} and \code{psi}, and the upper and lower bounds of the confidence intervals calculated.
#' @export
#' @examples
#' ci <- JP.ci.nt(JP.mle(k.2), alpha = 0.05)
JP.ci.nt <- function(jp.ests, alpha = 0.05) {
    quant <- qnorm(1-alpha/2)
    infmat <- solve(jp.ests$HessMat)
    standerr <- sqrt(diag(infmat))
    
    list(alpha = alpha,
         mu = c(est = jp.ests$mu, lower = jp.ests$mu-(quant*standerr[1]), upper = jp.ests$mu+(quant*standerr[1])),
         kappa = c(est = jp.ests$kappa, lower = jp.ests$kappa-(quant*standerr[2]), upper = jp.ests$kappa+(quant*standerr[2])),
         psi = c(est = jp.ests$psi, lower = jp.ests$psi-(quant*standerr[3]), upper = jp.ests$psi+(quant*standerr[3]))    )
}


#' Bootstrap confidence interval for Jones-Pewsey parameter estimators
#'
#' Obtain nominal 100(1-\code{alpha})% confidence intervals for Jones-Pewsey estimators, using a bootstrap resampling method.
#' @param data Vector of angles over which maximum likelihood estimation is to be performed.
#' @param alpha Significance level of confidence interval to be obtained. Default is 0.05 (95% confidence interval).
#' @param B Number of bootstrap samples to use to obtain the confidence interval. Default is 9999.
#' @param show.progress Boolean indicating whether or not to display a progress bar as the bootstrap is run.
#' @return A list containing estimates of \code{mu}, \code{kappa} and \code{psi}, and the upper and lower bounds of the confidence intervals calculated.
#' @export
#' @examples
#' ci.boot <- JP.ci.boot(q)
JP.ci.boot <- function(data, alpha = 0.05, B = 9999, show.progress = T) {
    n <- length(data)
    
    JPmleRes <- JP.mle(data)
    mu.est <- JPmleRes$mu
    kappa.est <- JPmleRes$kappa
    psi.est <- JPmleRes$psi
    ncon.est <- JP.NCon(kappa.est, psi.est)
    
    # resample from distribution & get B parameter estimates
    if (show.progress) {pb <- txtProgressBar(min = 0, max = B, style = 3)}
    for (b in 2:(B+1)) {
        sample <- JP.sim(n, mu.est[1], kappa.est[1], psi.est[1], ncon.est[1])
        JPmleRes <- JP.mle(sample) 
        mu.est[b] <- JPmleRes$mu
        kappa.est[b] <- JPmleRes$kappa
        psi.est[b] <- JPmleRes$psi
        ncon.est[b] <- JP.NCon(kappa.est[b], psi.est[b])
        if (show.progress) {setTxtProgressBar(pb, b)}
    }
    dist <- pi-abs(pi-abs(mu.est-mu.est[1]))
    sdist <- sort(dist)
    mu.lower <- mu.est[1]-sdist[(B+1)*(1-alpha)] 
    mu.upper <- mu.est[1]+sdist[(B+1)*(1-alpha)]
    
    skappa.est <- sort(kappa.est) 
    kap.lower <- skappa.est[(B+1)*alpha/2] ; kap.upper <- skappa.est[(B+1)*(1-alpha/2)]
    
    spsi.est <- sort(psi.est)
    psi.lower <- spsi.est[(B+1)*alpha/2] ; psi.upper <- spsi.est[(B+1)*(1-alpha/2)]
    
    if (show.progress) {close(pb)}
    
    list(alpha = alpha,
         mu = c(est = mu.est, lower = jp.ests$mu-(quant*standerr[1]), upper = jp.ests$mu+(quant*standerr[1])),
         kappa = c(est = kappa.est, lower = jp.ests$kappa-(quant*standerr[2]), upper = jp.ests$kappa+(quant*standerr[2])),
         psi = c(est = psi.est, lower = jp.ests$psi-(quant*standerr[3]), upper = jp.ests$psi+(quant*standerr[3]))    )
}


#================================================================================


#' Jones-Pewsey P-P plot
#'
#' Produces a P-P plot of the data against a specified Jones-Pewsey distribution, to graphically assess the goodness of fit of the model.
#' @param data Vector of angles to be fitted against the Jones-Pewsey distribution.
#' @param mu Mean direction parameter.
#' @param kappa Concentration parameter.
#' @param psi Shape parameter.
#' @export
#' @examples
#' JP.PP(k.2, mu = 3, kappa = 1, psi = -3)
JP.PP <- function(data, mu, kappa, psi) {
    ncon <- JP.NCon(kappa, psi)
    n <- length(data)
    edf <- ecdf(data) 
    tdf <- 0 
    for (j in 1:n) {tdf[j] <- JP.df(data[j], mu, kappa, psi, ncon)}
    
    plot.default(tdf, edf(data), pch=20, xlim=c(0,1), ylim=c(0,1),
                 xlab = "Jones-Pewsey distribution function", ylab = "Empirical distribution function")
    lines(c(0,1), c(0,1), lwd=2, col = "lightseagreen")
    edf(data) - tdf
}


#' Jones-Pewsey Q-Q plot of data
#'
#' Produces a Q-Q plot of the data against a specified Jones-Pewsey distribution, to graphically assess the goodness of fit of the model.
#' @param data Vector of angles to be fitted against the Jones-Pewsey distribution.
#' @param mu Mean direction parameter.
#' @param kappa Concentration parameter.
#' @param psi Shape parameter.
#' @export
#' @examples
#' JP.QQ(k.2, mu = 3, kappa = 1, psi = -3)
JP.QQ <- function(data, mu, kappa, psi) {
    ncon <- JP.NCon(kappa, psi)
    n <- length(data)
    edf <- ecdf(data) 
    tqf <- 0 
    for (j in 1:n) {tqf[j] <- JP.qf(edf(data)[j], mu, kappa, psi, ncon)}
    
    plot.default(tqf, data, pch=20, xlim=c(0,2*pi), ylim=c(0,2*pi), xlab = "Jones-Pewsey quantile function", ylab = "Empirical quantile function") 
    lines(c(0,2*pi), c(0,2*pi), lwd=2, col = "lightseagreen")
    data - tqf
}


#================================================================================


#' Goodness-of-fit tests for Jones-Pewsey distribution
#'
#' Runs a bundle of goodness-of-fit tests on the data against the Jones-Pewsey distribution.
#' @details The tests rely on doubling the Jones-Pewsey distribution function of each angle in the data set, and testing the resulting distribution for uniformity. 
#' @param data Vector of angles to be tested against the Jones-Pewsey distribution.
#' @param mu Mean direction parameter of the Jones-Pewsey distribution to be fitted.
#' @param kappa Concentration parameter of the Jones-Pewsey distribution to be fitted.
#' @param psi Shape parameter of the Jones-Pewsey distribution to be fitted.
#' @param display Boolean specifying whether to print the results of the four tests to the console or not.
#' @return List of test statistics for the four tests, plus the p-value for the Rayleigh tests.
#' @seealso The individual tests included in the bundle: \code{\link{kuiper.test}}, \code{\link{watson.test}}, \code{\link{rao.spacing.test}}, and \code{\link{rayleigh.test}}. 
#' @export
#' @examples
#' GoF <- JP.GoF(q, mu, kappa, psi)
JP.GoF <- function(data, mu, kappa, psi, display = T) {
    n <- length(data)
    tdf <- 0
    ncon <- JP.NCon(kappa, psi)
    
    for (j in 1:n) {tdf[j] <- JP.df(data[j], mu, kappa, psi, ncon)}
    cunif <- circular(2*pi*tdf)
    
    if (display) {
        print(kuiper.test(cunif))
        print(watson.test(cunif))
        print(rao.spacing.test(cunif))
        print(rayleigh.test(cunif))
    }
    
    list(kuiper.statistic = kuiper.test(cunif)$statistic,
         watson.statistic = watson.test(cunif)$statistic,
         rao.statistic = rao.spacing.test(cunif)$statistic,
         rayleigh.statistic = rayleigh.test(cunif)$statistic,
         rayleigh.p = rayleigh.test(cunif)$p.value)
}


#' Bootstrap goodness-of-fit tests for Jones-Pewsey distribution
#'
#' Runs a bundle of bootstrap goodness-of-fit tests on the data against the Jones-Pewsey distribution.
#' @details The tests rely on doubling the Jones-Pewsey distribution function of each angle in the data set, and testing the resulting distribution for uniformity. 
#' @param data Vector of angles to be tested against the Jones-Pewsey distribution.
#' @param B Number of bootstrap samples to use to obtain the confidence interval. Default is 9999.
#' @param show.progress Boolean indicating whether or not to display a progress bar as the bootstrap is run.
#' @return List of p-values for the four tests.
#' @seealso The individual tests included in the bundle: \code{\link{kuiper.test}}, \code{\link{watson.test}}, \code{\link{rao.spacing.test}}, and \code{\link{rayleigh.test}}. 
#' @export
#' @examples
#' GoF.boot <- JP.GoF.boot(q, B = 99)
JP.GoF.boot <- function(data, B = 9999, show.progress = T) {
    n <- length(data)
    JPmleRes <- JP.mle(data) 
    muhat0 <- JPmleRes$mu
    kaphat0 <- JPmleRes$kappa
    psihat0 <- JPmleRes$psi
    ncon0 <- JP.NCon(kaphat0, psihat0)
    
    tdf <- 0 
    for (j in 1:n) {
        tdf[j] <- JP.df(data[j], muhat0, kaphat0, psihat0, ncon0)
    }
    
    cunif <- circular(2*pi*tdf)
    nxtrm <- rep(1,4)
    
    unif.test.0 <- c(kuiper.test(cunif)$statistic, 
                     watson.test(cunif)$statistic,
                     rao.spacing.test(cunif)$statistic,
                     rayleigh.test(cunif)$statistic)
    
    if (show.progress) {pb <- txtProgressBar(min = 0, max = B, style = 3)}
    for (b in 2:(B+1)) {
        bootstrap.sample <- JP.sim(n, muhat0, kaphat0, psihat0, ncon0)
        JPmleRes <- JP.mle(bootstrap.sample) 
        muhat1 <- JPmleRes$mu
        kaphat1 <- JPmleRes$kappa
        psihat1 <- JPmleRes$psi
        ncon1 <- JP.NCon(kaphat1, psihat1)
        tdf <- 0
        for (j in 1:n) {
            tdf[j] <- JP.df(bootstrap.sample[j], muhat1, kaphat1, psihat1, ncon1)
        }
        cunif <- circular(2*pi*tdf)
        
        nxtrm[1] <- nxtrm[1] + (kuiper.test(cunif)$statistic >= unif.test.0[1])
        
        nxtrm[2] <- nxtrm[2] + (watson.test(cunif)$statistic >= unif.test.0[2])
        nxtrm[3] <- nxtrm[3] + (rao.spacing.test(cunif)$statistic  >= unif.test.0[3])
        nxtrm[4] <- nxtrm[4] + (rayleigh.test(cunif)$statistic >= unif.test.0[4])
        if (show.progress) {setTxtProgressBar(pb, b)}
    }
    pval <- nxtrm/(B+1)
    names(pval) <- c("kuiper", "watson", "rao", "rayleigh")
    if (show.progress) {close(pb)}
    return(pval)
}


#' Likelihood ratio test for Jones-Pewsey against nested models
#'
#' Likelihood ratio test of Jones-Pewsey distribution with MLE parameters against a simpler model with a specified value of \code{psi}.
#' @param data Vector of angles to be tested.
#' @param psi.0 Value of \code{psi} specifying a simpler model, to test against the three-parameter Jones-Pewsey.
#' @param alpha Significance level of likelihood ratio test. Default is 0.05 (95% confidence interval).
#' @details von Mises has psi -> 0; cardioid has psi = 1; wrapped Cauchy has psi = -1.
#' @return List containing the deviance, p-value, matrix ofthe parameters tested, and a comment explaining the output.
#' @export
#' @examples
#' JP.LR.test <- JP.psi.LR.test(q, psi.0 = 1)
JP.psi.LR.test <- function(data, psi.0 = 0, alpha = 0.05) {
    
    null.model = "null"
    if (psi.0 == 0) {null.model <- "von Mises"}
    if (psi.0 == 1) {null.model <- "cardioid"}
    if (psi.0 ==-1) {null.model <- "wrapped Cauchy"}
    
    n <- length(data)
    s <- sum(sin(data))
    c <- sum(cos(data))
    mu.vM <- atan2(s,c) %% (2*pi)
    kappa.vM <- A1inv(sqrt(s*s+c*c)/n)
    
    JPnll.psi <- function(p){
        mu <- p[1]
        kappa <- p[2]
        psi <- p[3]
        parlim <- abs(kappa*psi)
        if (parlim > 10) {y <- 9999.0 ; return(y)}
        else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(data, mu, kappa, psi, ncon)))
            return (y)
        }
    }
    
    JPnll.psi.0 <- function(p){
        mu <- p[1]
        kappa <- p[2]
        psi <- psi.0 
        parlim <- abs(kappa*psi)
        if (parlim > 10) {y <- 9999.0 ; return(y)} 
        else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(data, mu, kappa, psi, ncon))) 
            return(y) 
        }
    }
    
    out <- optim(par=c(mu.vM, kappa.vM, 0), fn=JPnll.psi, gr = NULL, method = "L-BFGS-B", lower = c(mu.vM-pi, 0, -Inf), upper = c(mu.vM+pi, Inf, Inf))
    maxll1 <- -out$value
    muhat1 <- out$par[1]
    kaphat1 <- out$par[2]
    psihat1 <- out$par[3]
    
    out <- optim(par=c(muhat1, kaphat1), fn=JPnll.psi.0, gr = NULL, method = "L-BFGS-B", lower = c(muhat1-pi, 0), upper = c(muhat1+pi, Inf))
    maxll0 <- -out$value 
    muhat0 <- out$par[1] 
    kaphat0 <- out$par[2]
    
    D <- round(-2*(maxll0-maxll1),3)
    
    pval <- pchisq(D, df= 1, lower.tail=F)
    if (pval < alpha) {outcome = "gives"} else {outcome = "does not give"}
    pval <- round(pval, 3)
    
    comment = paste("Jones-Pewsey ", outcome, " an improvement on ", null.model,
                    " model at the ", alpha * 100, "% level.", sep = "")
    
    comparison <- round(rbind(max.ll = c(maxll0, maxll1),
                              mu = c(muhat0, muhat1),
                              kappa = c(kaphat0, kaphat1),
                              psi = c(psi.0, psihat1)),3)
    colnames(comparison) <- c(null.model, "Jones-Pewsey")
    
    list(D = D, p.val = pval, comment = comment, comparison = comparison)
}


#' Bootstrap likelihood ratio test for Jones-Pewsey against nested models
#'
#' Bootstrap likelihood ratio test of Jones-Pewsey distribution with MLE parameters against a simpler model with a specified value of \code{psi}.
#' @param data Vector of angles to be tested.
#' @param psi.0 Value of \code{psi} specifying a simpler model, to test against the three-parameter Jones-Pewsey.
#' @param alpha Significance level of likelihood ratio test. Default is 0.05 (95% confidence interval).
#' @param B Number of bootstrap samples to use to obtain the confidence interval. Default is 9999.
#' @details von Mises has psi -> 0; cardioid has psi = 1; wrapped Cauchy has psi = -1.
#' @return List containing the deviance, p-value, matrix ofthe parameters tested, and a comment explaining the output.
#' @export
#' @examples
#' JP.LR.test <- JP.psi.LR.boot(q, psi.0 = 1, B = 99)
JP.psi.LR.boot <- function(data, psi.0 = 0, B = 9999, alpha = 0.05, show.progress = T) {
    
    x <- data 
    n <- length(x)
    s <- sum(sin(x))
    c <- sum(cos(x))
    mu.vM <- atan2(s,c) %% (2*pi)
    kappa.vM <- A1inv(sqrt(s*s+c*c)/n)
    
    JPnll.psi <- function(p){
        mu <- p[1]
        kappa <- p[2]
        psi <- p[3] 
        parlim <- abs(kappa*psi)
        if (parlim > 10) {y <- 9999.0 ; return(y)}
        else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(x, mu, kappa, psi, ncon))) ; return (y) 
        }
    }
    
    JPnll.psi0 <- function(p){
        mu <- p[1] 
        kappa <- p[2] 
        psi <- psi.0
        parlim <- abs(kappa*psi)
        if (parlim > 10) {y <- 9999.0 ; return(y)} 
        else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(x, mu, kappa, psi, ncon))) ; return(y)
        }
    }
    
    out <- optim(par=c(mu.vM, kappa.vM, 0), fn=JPnll.psi, gr = NULL, method = "L-BFGS-B", lower = c(mu.vM-pi, 0, -Inf), upper = c(mu.vM+pi, Inf, Inf))
    maxll1 <- -out$value 
    muhat1 <- out$par[1]
    kaphat1 <- out$par[2]
    psihat1 <- out$par[3]
    
    out <- optim(par=c(muhat1, kaphat1), fn=JPnll.psi0, gr = NULL, method = "L-BFGS-B", lower = c(muhat1-pi, 0), upper = c(muhat1+pi, Inf))
    maxll0 <- -out$value 
    muhat0 <- out$par[1]
    kaphat0 <- out$par[2]
    
    D <- -2*(maxll0-maxll1) 
    
    if (show.progress) {pb <- txtProgressBar(min = 0, max = B, style = 3)}
    nxtrm <- 1
    for (j in 2:(B+1)) {
        stopgo <- 0
        while (stopgo == 0) {
            
            x <- JP.sim(n, muhat0, kaphat0, psi.0)
            
            out <- optim(par=c(muhat0, kaphat0, psi.0), fn=JPnll.psi, gr = NULL, method = "L-BFGS-B", lower = c(muhat0-pi, 0, -Inf), upper = c(muhat0+pi, Inf, Inf))
            maxll1 <- -out$value 
            
            out <- optim(par=c(muhat0, kaphat0), fn=JPnll.psi0, gr = NULL, method = "L-BFGS-B", lower = c(muhat0-pi, 0), upper = c(muhat0+pi, Inf))
            maxll0 <- -out$value
            
            D[j] <- -2*(maxll0-maxll1) 
            
            if (D[j] >= 0) {
                if (D[j] < 40) { 
                    if (D[j] >= D[1]) {nxtrm <- nxtrm+1}
                    stopgo <- 1 } }
        } 
        if (show.progress) {setTxtProgressBar(pb, j)}
    }
    if (show.progress) {close(pb)}
    
    pval <- nxtrm/(B+1) 
    if (pval < alpha) {outcome = "gives"} else {outcome = "does not give"}
    pval <- round(pval,3)
    
    null.model = "null"
    if (psi.0 == 0) {null.model <- "von Mises"}
    if (psi.0 == 1) {null.model <- "cardioid"}
    if (psi.0 ==-1) {null.model <- "wrapped Cauchy"}
    
    comment = paste("Jones-Pewsey ", outcome, " an improvement on ", null.model,
                    " model at the ", alpha * 100, "% level.", sep = "")
    
    comparison <- round(rbind(max.ll = c(maxll0, maxll1),
                              mu = c(muhat0, muhat1),
                              kappa = c(kaphat0, kaphat1),
                              psi = c(psi.0, psihat1)), 3)
    colnames(comparison) <- c(null.model, "Jones-Pewsey")
    
    list(p.val = pval, comment = comment, comparison = comparison)
}


#' AIC and BIC for Jones-Pewsey against nested models
#'
#' Obtain the AIC and BIC scores for Jones-Pewsey distribution with MLE parameters over a simpler model with a specified value of \code{psi}.
#' @param data Vector of angles to be tested.
#' @param psi.0 Value of \code{psi} specifying a simpler model, to test against the three-parameter Jones-Pewsey.
#' @details von Mises has psi -> 0; cardioid has psi = 1; wrapped Cauchy has psi = -1.
#' @return List containing a matrix ofthe parameters tested and their associated AIC and BIC scores, as well as a comment explaining the output.
#' @seealso \code{\link{AIC}}, \code{\link{BIC}}
#' @export
#' @examples
#' JP.info <- JP.psi.info(q)
JP.psi.info <- function(data, psi.0 = 0) {
    null.model = "null"
    if (psi.0 == 0) {null.model <- "von Mises"}
    if (psi.0 == 1) {null.model <- "cardioid"}
    if (psi.0 ==-1) {null.model <- "wrapped Cauchy"}
    
    x <- data
    n <- length(x)
    s <- sum(sin(x)) 
    c <- sum(cos(x))
    k <- c(2,3)
    mu.vM <- atan2(s,c) %% (2 * pi)
    kappa.vM <- A1inv(sqrt(s*s+c*c)/n)
    
    JPnll <- function(p){
        mu <- p[1] 
        kappa <- p[2] 
        psi <- p[3] 
        parlim <- abs(kappa*psi)
        if (parlim > 10) {y <- 9999.0 ; return(y)}
        else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(x, mu, kappa, psi, ncon))) 
            return (y)
        }
    }
    
    JPnllpsi0 <- function(p){
        mu <- p[1] 
        kappa <- p[2] 
        psi <- psi.0 
        parlim <- abs(kappa*psi)
        if (parlim > 10) {y <- 9999.0 ; return(y)} 
        else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(x, mu, kappa, psi, ncon))) 
            return(y)
        }
    }
    
    out <- optim(par=c(mu.vM, kappa.vM, 0), fn=JPnll, gr = NULL, method = "L-BFGS-B", lower = c(mu.vM-pi, 0, -Inf), upper = c(mu.vM+pi, Inf, Inf))
    maxll1 <- -out$value 
    muhat1 <- out$par[1]
    kaphat1 <- out$par[2] 
    psihat1 <- out$par[3]
    nu <- 3
    AIC1 <- (2*nu)-(2*maxll1)
    BIC1 <- (log(n)*nu)-(2*maxll1) 
    
    out <- optim(par=c(muhat1, kaphat1), fn=JPnllpsi0, gr = NULL, method = "L-BFGS-B", lower = c(muhat1-pi, 0), upper = c(muhat1+pi, Inf))
    maxll0 <- -out$value 
    nu <- 2
    AIC0 <- (2*nu)-(2*maxll0)
    BIC0 <- (log(n)*nu)-(2*maxll0) 
    
    if (AIC0 < AIC1) {
        comment.AIC <- paste("AIC favours ", null.model, " model over full Jones-Pewsey distributon.", sep = "")
    } else {
        comment.AIC <- paste("AIC favours full Jones-Pewsey distribution over ", null.model, " model.", sep = "")
    }
    if (BIC0 < BIC1) {
        comment.BIC <- paste("BIC favours ", null.model, " model over full Jones-Pewsey distributon.", sep = "")
    } else {
        comment.BIC <- paste("BIC favours full Jones-Pewsey distribution over ", null.model, " model.", sep = "")
    }
    comparison <- round(rbind(mu = c(mu.vM, muhat1),
                              kappa = c(kappa.vM, kaphat1),
                              psi = c(psi.0, psihat1),
                              llh = c(maxll0, maxll1),
                              AIC = c(AIC0, AIC1),
                              BIC = c(BIC0, BIC1),
                              AICc = c(AIC0, AIC1) + (2*k*(k+1)) / (n-k-1)),3)
    colnames(comparison) <- c(null.model, "Jones-Pewsey")
    
    list(comparison = comparison, AIC = comment.AIC, BIC = comment.BIC)
}