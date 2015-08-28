
#' von Mises P-P plot
#'
#' Produces a P-P plot of the data against a specified von Mises distribution, to graphically assess the goodness of fit of the model.
#' @param data Vector of angles (in radians) to be fitted against the von Mises distribution.
#' @param mu Mean direction parameter for von Mises distribution.
#' @param kappa Concentration parameter for von Mises distribution. Must be between 0 and 1.
#' @export
#' @examples
#' vM.PP(k.2, mu = 3, kappa = 1)
vM.PP <- function(data, mu, kappa)  {
    edf <- ecdf(data)
    tdf <- pvonmises(data, mu, kappa, from=circular(0), tol = 1e-06)
    
    plot.default(tdf, edf(data), pch=20, xlim=c(0,1), ylim=c(0,1), xlab = "von Mises distribution function", ylab = "Empirical distribution function")
    lines(c(0,1), c(0,1), lwd=2, col = "lightseagreen")
    edf(data) - tdf
}


#' von Mises Q-Q plot
#'
#' Produces a Q-Q plot of the data against a specified von Mises distribution, to graphically assess the goodness of fit of the model.
#' @param data Vector of angles (in radians) to be fitted against the von Mises distribution.
#' @param mu Mean direction parameter for von Mises distribution.
#' @param kappa Concentration parameter for von Mises distribution. Must be between 0 and 1.
#' @export
#' @examples
#' vM.QQ(k.2, mu = 3, kappa = 1)
vM.QQ <- function(data, mu, kappa)  {
    edf <- ecdf(data)
    tqf <- qvonmises(edf(data), mu, kappa, from=circular(0), tol = 1e-06)
    
    plot.default(tqf, data, pch=20, xlim=c(0,2*pi), ylim=c(0,2*pi), xlab = "von Mises quantile function", ylab = "Empirical quantile function")
    lines(c(0,2*pi), c(0,2*pi), lwd=2, col = "lightseagreen")
    data - tqf
}


#' Goodness-of-fit tests for von Mises distribution
#'
#' Runs a bundle of goodness-of-fit tests on the data against a specified von Mises distribution.
#' @details Four of the tests rely on doubling the von Mises distribution function of each angle in the data set, and testing the resulting distribution for uniformity. 
#' @param data Vector of angles to be tested against the von Mises distribution.
#' @param mu Mean direction parameter of the von Mises distribution to be fitted.
#' @param kappa Concentration parameter of the von Mises distribution to be fitted. Must be between 0 and 1.
#' @param display Boolean specifying whether to print the results of the four tests to the console or not.
#' @return List of test statistics for the four tests, plus the p-value for the Rayleigh tests.
#' @seealso The individual tests included in the bundle: \code{\link{kuiper.test}}, \code{\link{watson.test}}, \code{\link{rao.spacing.test}}, and \code{\link{rayleigh.test}}. 
#' @export
#' @examples
#' GoF <- vM.GoF(q, mu, kappa)
vM.GoF <- function(data, mu, kappa, display = T) {
    tdf <- pvonmises(data, circular(mu), kappa, from=circular(0), tol = 1e-06)
    cunif <- circular(2*pi*tdf)
    
    if (display) {
        print(watson.test(data, dist = "vonmises"))
        print(kuiper.test(cunif))
        print(watson.test(cunif))
        print(rao.spacing.test(cunif))
        print(rayleigh.test(cunif))
    }
    
    list(watson.vM.stat = watson.test(data, dist = "vonmises")$statistic,
         kuiper.unif.stat = kuiper.test(cunif)$statistic,
         watson.unif.stat = watson.test(cunif)$statistic,
         rao.unif.stat = rao.spacing.test(cunif)$statistic,
         rayleigh.unif.stat = rayleigh.test(cunif)$statistic,
         rayleigh.unif.pval = rayleigh.test(cunif)$p.val)
}


#' Bootstrap goodness-of-fit tests for von Mises distribution
#'
#' Runs a bundle of bbootstrap goodness-of-fit tests on the data against the von Mises distribution.
#' @details The tests rely on doubling the von Mises distribution function of each angle in the data set, and testing the resulting distribution for uniformity. 
#' @param data Vector of angles to be tested against the von Mises distribution.
#' @param B Number of bootstrap samples to use to obtain the confidence interval. Default is 9999.
#' @param show.progress Boolean indicating whether or not to display a progress bar as the bootstrap is run.
#' @return List of p-values for the four tests.
#' @seealso The individual tests included in the bundle: \code{\link{kuiper.test}}, \code{\link{watson.test}}, \code{\link{rao.spacing.test}}, and \code{\link{rayleigh.test}}. 
#' @export
#' @examples
#' GoF.boot <- vM.GoF.boot(q, B = 99)
vM.GoF.boot <- function(data, B = 9999, show.progress = T) {
    
    n <- length(data)
    vM.mle <- mle.vonmises(data, bias=TRUE)
    mu.0 <- vM.mle$mu
    kappa.0 <- vM.mle$kappa
    
    tdf <- pvonmises(data, mu.0, kappa.0, from=circular(0), tol = 1e-06)
    cunif <- circular(2*pi*tdf)
    
    unif.test.0 <- rep(0,4)
    nxtrm <- rep(1,4)
    
    unif.test.0[1] <- kuiper.test(cunif)$statistic
    unif.test.0[2] <- watson.test(cunif)$statistic
    unif.test.0[3] <- rao.spacing.test(cunif)$statistic
    unif.test.0[4] <- rayleigh.test(cunif)$statistic
    
    if (show.progress) {pb <- txtProgressBar(min = 0, max = B, style = 3)}
    for (b in 2:(B+1)) {
        bootstrap.sample <- rvonmises(n, mu.0, kappa.0)
        vM.mle <- mle.vonmises(bootstrap.sample, bias = T)
        mu.1 <- vM.mle$mu
        kappa.1 <- vM.mle$kappa
        
        tdf <- pvonmises(bootstrap.sample, mu.1, kappa.1, from=circular(0), tol = 1e-06)
        cunif <- circular(2*pi*tdf)
        
        nxtrm[1] <- nxtrm[1] + (kuiper.test(cunif)$statistic >= unif.test.0[1])
        nxtrm[2] <- nxtrm[2] + (watson.test(cunif)$statistic >= unif.test.0[2])
        nxtrm[3] <- nxtrm[3] + (rao.spacing.test(cunif)$statistic >= unif.test.0[3])
        nxtrm[4] <- nxtrm[4] + (rayleigh.test(cunif)$statistic >= unif.test.0[4])
        if (show.progress) {setTxtProgressBar(pb,b)}      
    }
    p.val <- nxtrm/(B+1)
    names(p.val) <- c("kuiper", "watson", "rao", "rayleigh")
    if (show.progress) {close(pb)}      
    p.val
}