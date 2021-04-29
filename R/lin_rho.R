#' Lin's Rho (or CCC)
#'
#' Lin's agreement/concordance correlation coefficient (using
#' epiR::epi.ccc as workhorse)
#'
#' @param x vector of first set of measurements
#' @param y vector of first set of measurements
#' @param ci method to be used for ci
#' @param conf.level magniture of the returned confidence interval
#' @param rep.measure if TRUE there are repeated observation across subject
#' @param subjectid factor providing id of the observer if rep.measure == TRUE
#'
#' @export
lin_rho <- function(x = NULL,
                    y = NULL,
                    ci = c('z-transform', 'asymptotic'),
                    conf.level = 0.95,
                    rep.measure = FALSE,
                    subjectid = NULL)
{
    if (is.null(x) || is.null(y))
        stop("x and y can't be NULL")
    ci <- match.arg(ci)
    if (conf.level >= 1 || conf.level <= 0)
        stop("conf.level must be between 0 and 1")
    if (rep.measure && is.null(subjectid))
        stop('if rep.measure == TRUE you must specify subjectid')

    rval <- epiR::epi.ccc(x = x,
                          y = y,
                          ci = ci,
                          conf.level = conf.level,
                          rep.measure = rep.measure,
                          subjectid = subjectid)
    setNames(rval$rho.c, c("estimate", "lower", "upper"))
}

#' Scatterplot with (0,0) to (1,1) line
#'
#' @param x generally the gold standard (on x axis)
#' @param y generally the testing measurement (on y axis)
#' @param stat which statistics to plot in a (topleft-placed) legend
#' 
#' @export
lin_plot <- function(x, y, stat = c("ci", 'coef', 'none'), ...){
    ## pch = 20
    stat <- match.arg(stat)
    mi <- min(c(x, y), na.rm = TRUE)
    ma <- max(c(x, y), na.rm = TRUE)
    xylim <- c(floor(mi), ceiling(ma))
    lin <- lin_rho(x = x, y = y)
    lin_string <- switch(stat,
                         ci = sprintf("CCC: %.2f (95%% CI: %.2f - %.2f)",
                                      lin$estimate,
                                      lin$lower,
                                      lin$upper),
                         coef = sprintf("CCC: %.2f", lin$estimate))

    plot(x = x, y = y,
         xlim = xylim, ylim = xylim,
         las = 1, pch = 20,
         ...)
    graphics::grid()
    if (stat %in% c("ci", "coef"))
        legend("topleft", legend = lin_string, bg = 'white')
    abline(a = 0, b = 1, col = "blue")
}




#' OCCC (overall concordance correlation coefficient)
#'
#' da citare barnhart2002occc
#' 
#' @param x data.frame of measurement one column for rater
#' @param ... further arguments passed to agRee::agree.ccc
#'
#' @export
occc <- function(x, ...){
    
    x <- NA_remove(x, quiet = TRUE)
    res <- agRee::agree.ccc(ratings = as.matrix(x), ...)
    setNames(unlist(res), c('estimate', 'lower', 'upper'))

}



#' computes difference between lin CCC (and confidence interval)
#' comparing x and y with a reference standard
#'
#' @param x first measurement (classical)
#' @param y second measurement (experimental)
#' @param refstd reference standard for measurement
#' @param conf confidence level
#' @param R bootstrap replicates
#' 
#' @export
lin_diff <- function(x, y, refstd, conf = 0.95, R = 10000){
    f <- function(a, b) lbagree::lin_rho(x = a, y = b)$estimate
    lin_x <- lbagree::lin_rho(x = x, y = refstd)
    lin_y <- lbagree::lin_rho(x = y, y = refstd)
    
    est <- f(y, refstd) - f(x, refstd) 
    ## percentile confidence interval
    db <- data.frame("x" = x, "y" = y, "refstd" = refstd)
    diff_f <- function(data, i) {
        f(data[i, 'y'], data[i, 'refstd']) - 
        f(data[i, 'x'], data[i, 'refstd']) 
    }
    res <-  boot::boot(data = db, statistic = diff_f, R = R)
    ci <- boot::boot.ci(res, type = 'bca', conf = conf)
    boot_ci <- ci[[4]][4:5]
    res <- c(est, boot_ci)
    nm <- c('CCC diff', 'Lower CI', "Upper CI")
    diff <- setNames(res, nm = nm)
    list('lin_x' = lin_x, 'lin_y' = lin_y, 'diff' = diff)
}
