
#' Circular sample moments
#'
#' Obtain the key sample trigonometric moments of the distribution.
#' @param data Vector of angles.
#' @return A list containing the sample mean, mean resultant length, and second, third and fourth sine and consine moments of the data.
#' @export
#' @examples
#' bar <- get.moments(q)
get.moments <- function(data) {
    
    t2 <- trigonometric.moment(data, p = 2, center = T)
    t3 <- trigonometric.moment(data, p = 3, center = T)
    t4 <- trigonometric.moment(data, p = 4, center = T)
    
    list(mu = as.numeric(mean(data)),
         r = rho.circular(data),
         b2 = t2$sin,
         b3 = t3$sin,
         b4 = t4$sin,
         a2 = t2$cos,
         a3 = t3$cos,
         a4 = t4$cos)
}


#' Bias-corrected circular sample statistics.
#'
#' Obtain the main circular sample statistics of the distribution, with large-sample normal-theory confidence intervals.
#' @param data Vector of angles.
#' @param symmetric Boolean indicating whether or not the data is assumed to be symmetric.
#' @return A list containing estimates of mu (location), rho (concentration), beta2 (skew) and alpha2 (kurtosis).
#' @export
#' @examples
#' bar <- bc.sample.statistics(q)
bc.sample.statistics <- function(data, symmetric = F) {
    
    n <- length(data)
    bar <- get.moments(data)    
    
    rho.est <- bar$r - ((1-bar$a2)/(4*n*bar$r))
    r2 <- bar$r^2
    r4 <- r2^2
    
    if (symmetric) {
        beta2.est <- 0
    } else {
        beta2.est <- bar$b2 + ((bar$b3/bar$r)+(bar$b2/r2)-(2*bar$a2*bar$b2/r4))/n
    }
    
    div <- 2*n*r2
    mu.est <- (bar$mu + (bar$b2/div)) %% (2*pi)
    
    alpha2.est <- bar$a2 - (1-(bar$a3/bar$r)-((bar$a2*(1-bar$a2)+bar$b2*bar$b2)/r2))/n
    
    list(mu = mu.est, rho = rho.est, beta2 = beta2.est, alpha2 = alpha2.est)   
}


#' Large-sample confidence interval for sample statistics.
#'
#' Obtain the main circular sample statistics of the distribution, with large-sample normal-theory confidence intervals.
#' @param data Vector of angles.
#' @param alpha Significance level of confidence interval to be obtained. Default is 0.05 (95% confidence interval).
#' @param symmetric Boolean indicating whether or not the data is assumed to be symmetric.
#' @return A list containing estimates of mu (location), rho (concentration), beta2 (skew) and alpha2 (kurtosis), and the upper and lower bounds of the confidence intervals calculated.
#' @export
#' @examples
#' bar <- bc.ci.LS(q)
bc.ci.LS <- function(data, alpha = 0.05, symmetric = F) {
    
    n <- length(data)
    ests <- bc.sample.statistics(data, symmetric = symmetric)
    
    bar <- get.moments(data)

    r2 <- bar$r^2 ; r4 <- r2^2
    
    qval <- qnorm(1-alpha/2)
    
    rho.bc <- ests$rho
    r.SE <- sqrt((1-2*r2+bar$a2)/(2*n))    
    rho.upper <- rho.bc + qval*r.SE
    rho.lower <- rho.bc - qval*r.SE
    
    rho <- c(estimate = rho.bc, lower = rho.lower, upper = rho.upper)
    
    if (symmetric) {
        beta2 <- c(estimate = 0, lower = 0, upper = 0)
    } else {    
        beta2.bc <- ests$beta2
        beta2.SE <- sqrt((((1-bar$a4)/2)-(2*bar$a2)-(bar$b2^2)+(2*bar$a2/bar$r)*(bar$a3+(bar$a2*(1-bar$a2)/bar$r)))/n)
        beta2.upper <- beta2.bc + qval*beta2.SE
        beta2.lower <- beta2.bc - qval*beta2.SE
        
        beta2 <- c(estimate = beta2.bc, lower = beta2.lower, upper = beta2.upper)
    }
    
    div <- 2*n*r2 
    
    mu.bc <- ests$mu
    mu.SE <- sqrt((1-bar$a2)/div)
    
    mu.upper <- mu.bc + qval * mu.SE
    mu.lower <- mu.bc - qval * mu.SE
    
    mu <- c(estimate = mu.bc, lower = mu.lower, upper = mu.upper)
    
    alpha2.bc <- ests$alpha2
    alpha2.SE <- sqrt((((1+bar$a4)/2)-(bar$a2*bar$a2)+(2*bar$b2/bar$r)*(bar$b3+(bar$b2*(1-bar$a2)/bar$r)))/n)
    alpha2.upper <- alpha2.bc + qval*alpha2.SE
    alpha2.lower <- alpha2.bc - qval*alpha2.SE
    
    alpha2 <- c(estimate = alpha2.bc, lower = alpha2.lower, upper = alpha2.upper)
    
    list(alpha = alpha, symmetric = symmetric, mu = mu, rho = rho, beta2 = beta2, alpha2 = alpha2)
}


#' Bootstrap confidence interval for sample statistics.
#'
#' Obtain the main circular sample statistics of the distribution, with large-sample normal-theory confidence intervals.
#' @param data Vector of angles.
#' @param alpha Significance level of confidence interval to be obtained. Default is 0.05 (95% confidence interval).
#' @param symmetric Boolean indicating whether or not the data is assumed to be symmetric.
#' @param B Number of bootstrap samples to use to estimate the confidence intervals. Default is 9999.
#' @param show.progress Boolean indicating whether or not to display a progress bar as the bootstrap is run.
#' @return A list containing estimates of mu (location), rho (concentration), beta2 (skew) and alpha2 (kurtosis), and the upper and lower bounds of the confidence intervals calculated.
#' @export
#' @examples
#' ests <- bc.ci.boot(q)
bc.ci.boot <- function(data, symmetric = F, alpha = 0.05, B = 9999, show.progress = T) {
    
    n <- length(data)
    
    ests <- bc.sample.statistics(data, symmetric = symmetric)
    
    mu.est <- ests$mu
    rho.est <- ests$rho
    beta2.est <- ests$beta2
    alpha2.est <- ests$alpha2
    
    if (symmetric) {
        reflection <- 2 * mu.est - data
        sample.data <- c(data, reflection)
    } else {
        sample.data <- data
    }
    
    if (show.progress) {pb <- txtProgressBar(min = 0, max = B, style = 3)}
    for (b in 2 : (B+1)) { 
        bootstrap.sample <- sample(sample.data, size = n, replace = TRUE)
        ests <- bc.sample.statistics(bootstrap.sample, symmetric)
        
        mu.est[b] <- ests$mu
        rho.est[b] <- ests$rho
        beta2.est[b] <- ests$beta2
        alpha2.est[b] <- ests$alpha2
        if (show.progress) {setTxtProgressBar(pb, b)}
    }
    
    dist <- 0
    
    if (symmetric) {
        dist <- pi - abs(pi - abs(mu.est - mu.est[1]))        
        sdist <- sort(dist)
        
        mu.lower <- mu.est[1] - sdist[(B+1)*(1-alpha)]
        mu.upper <- mu.est[1] + sdist[(B+1)*(1-alpha)]
    } else {
        if (mu.est[1] < pi) {
            ref <- mu.est[1] + pi
            for (b in 1:(B+1)) {
                dist[b] <- -(pi-abs(pi-abs(mu.est[b]-mu.est[1])))
                if (mu.est[b] > mu.est[1]) {
                    if (mu.est[b] < ref) {
                        dist[b] <- -dist[b]
                    }
                }
            }
        } else {    #  if (mu.est[1] >= pi)
            ref <- mu.est[1] - pi
            for (b in 1:(B+1)) { 
                dist[b] <- pi-abs(pi-abs(mu.est[b]-mu.est[1]))
                if (mu.est[b] > ref) {
                    if (mu.est[b] < mu.est[1]) {
                        dist[b] <- -dist[b]
                    }                    
                }  
            }
        }
        sdist <- sort(dist)
        mu.lower <- mu.est[1]+sdist[(B+1)*(alpha/2)]
        mu.upper <- mu.est[1]+sdist[(B+1)*(1-alpha/2)]
        
        sbeta2.est <- sort(beta2.est)
        
        beta2.lower <- sbeta2.est[(B+1)*(alpha/2)]
        beta2.upper <- sbeta2.est[(B+1)*(1-alpha/2)]
        
        beta2 <- c(beta2.est[1], beta2.lower, beta2.upper)
    }    
    
    mu <- c(mu.est[1], mu.lower, mu.upper) 
    
    srho.est <- sort(rho.est)
    
    rho.lower <- srho.est[(B+1)*(alpha/2)] ; rho.upper <- srho.est[(B+1)*(1-alpha/2)]
    
    salpha2.est <- sort(alpha2.est)
    
    alpha2.lower<- salpha2.est[(B+1)*(alpha/2)]
    alpha2.upper <- salpha2.est[(B+1)*(1-alpha/2)]
    
    rho <- c(rho.est[1], rho.lower, rho.upper)
    alpha2 <- c(alpha2.est[1], alpha2.lower, alpha2.upper)
    if (show.progress) {close(pb)}
    list(mu = mu, rho = rho, beta2 = beta2, alpha2 = alpha2)    
}


#' Reflective symmetry test
#'
#' Large-sample asymptotic-theory test of reflective symmetry about an unspecified mean.
#' @param data Vector of angles.
#' @param alpha Significance level of confidence interval to be obtained. Default is 0.05 (95% confidence interval).
#' @param symmetric Boolean indicating whether or not the data is assumed to be symmetric.
#' @return A list containing the test statistic and p-value of the test.
#' @seealso \code{\link{wilcox.test}} for a rank test of reflective symmetry about a specified direction.
#' @export
#' @examples
#' reject.symmetry <- r.symm.test.stat(q)$p.value < 0.05
r.symm.test.stat <- function(data) {
    n <- length(data)
    bar <- get.moments(data)
    
    var <- ((1-bar$a4)/2-(2*bar$a2)+(2*bar$a2/bar$r) * (bar$a3 + (bar$a2*(1-bar$a2)/bar$r)))/n
    ts <- abs(bar$b2/sqrt(var))
    pval <- 2*pnorm(ts, mean = 0, sd = 1, lower = F)
    list(test.statistic = ts, p.val = pval)
}


#' Bootstrap reflective symmetry test
#'
#' Bootstrap test of reflective symmetry about an unspecified mean.
#' @param data Vector of angles.
#' @param alpha Significance level of confidence interval to be obtained. Default is 0.05 (95% confidence interval).
#' @param B Number of bootstrap samples to use to estimate the standard error of the p-value. Default is 9999.
#' @param show.progress Boolean indicating whether or not to display a progress bar as the bootstrap is run.
#' @return A list containing the p-value of the test and its standard error.
#' @seealso \code{\link{wilcox.test}} for a rank test of reflective symmetry about a specified direction.
#' @export
#' @examples
#' reject.symmetry <- r.symm.test.boot(q)
r.symm.test.boot <- function(data, B = 9999, alpha = 0.05) {
    
    n <- length(data)
    absz <- r.symm.test.stat(data)$test.statistic
    
    tbar <- mean(data)
    reflection <- 2*tbar-data
    symm.data <- c(data, reflection)
    
    if (show.progress) {pb <- txtProgressBar(min = 0, max = B, style = 3)}
    for (b in 2:(B+1)) { 
        bootstrap.sample <- sample(symm.data, size = n, replace = TRUE)
        absz[b] <- r.symm.test.stat(bootstrap.sample)$test.statistic
        if (show.progress) {setTxtProgressBar(pb, b)}
    }
    p.val = length(absz[absz >= absz[1]]) / (B + 1)
    if (show.progress) {close(pb)}
    list(p.val = p.val,
         std.error = qnorm(1-(alpha/2)) * sqrt((p.val*(1-p.val))/(B + 1)))
}


#' Bundle of tests of uniformity
#'
#' Runs a bundle of goodness-of-fit tests on the data against the uniform distribution.
#' @param data Vector of angles.
#' @param display Boolean specifying whether to print the results of the four tests to the console or not.
#' @return List of test statistics for the four tests, plus the p-value for the Rayleigh tests.
#' @seealso The individual tests included in the bundle: \code{\link{kuiper.test}}, \code{\link{watson.test}}, \code{\link{rao.spacing.test}}, and \code{\link{rayleigh.test}}. 
#' @export
#' @examples
#' bar$mu <- bc.point.estimates(q)$mu
uniformity.tests <- function(data, display = T) {
    
    if (display) {
        print(kuiper.test(data))
        print(watson.test(data))
        print(rao.spacing.test(data))
        print(rayleigh.test(data))
    }
    list(kuiper.statistic = kuiper.test(data)$statistic,
         watson.statistic = watson.test(data)$statistic,
         rao.statistic = rao.spacing.test(data)$statistic,
         rayleigh.statistic = rayleigh.test(data)$statistic,
         rayleigh.p = rayleigh.test(data)$p.value)
}


#' Uniformity plot of data
#'
#' Produces a uniform P-P plot. Graphical tool to assess uniformity of data.
#' @param data Vector of angles, in radians.
#' @param extend Boolean specifying whether to extend the tails of the plot to make any pattern at the extremes of the data easier to see.
#' @export
#' @examples
#' uniformity.plot(q)
uniformity.plot <- function(data, extend = F) {

    data <- data %% (2 * pi)
    n <- length(data)
    
    x <- c(1:n)/(n+1)
    y <- sort(data) / (2 * pi)
    
    if (extend) {
            f <- floor(n/5)
            x.end <- x[1:f] + 1
            y.end <- y[1:f] + 1
            
            x.start <- x[(n - f):n] - 1
            y.start <- y[(n - f):n] - 1
            
            x <- c(x.start, x, x.end)
            y <- c(y.start, y, y.end)
        } 
    
    plot(x, y, pch = 20, asp = T, xlab = "Uniform quantiles", ylab = "Sample quantiles")
    abline(a = 0, b = 1, col = "lightseagreen")
}