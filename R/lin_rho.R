#' Lin's Rho
#'
#' Lin's agreement correlation coefficient (using epiR::epi.ccc as workhorse)
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
