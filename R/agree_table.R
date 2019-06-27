#' Calculate agreement between two nominal/ordinal variables
#' 
#' @param x first qualitative variable
#' @param y first second variable
#' @param cohen_k weight for cohen (uw = unweighted, 'wl = linearly weighted',
#' @param conf.level confidence level for estimate confidence interval
#' @param ... further arguments passed to biv_quali
#' 
#' @export
agree_table <- function(x, y, cohen_k = c('uw', 'wl', 'wq'), 
                        conf.level = 0.95,  
                        ...) 
{
    cohen_k <- match.arg(cohen_k) 
    cohen <- lbagree::cohen_k(x = x, y = y, conf.level = conf.level)    
    stat <- switch(cohen_k,
                   uw = cohen$unweighted[1],
                   wl = cohen$`weighted (linear)`[1],
                   wq = cohen$`weighted (quadratic)`[1])
    method <- switch(cohen_k,
                   uw = "Cohen's Kappa (unweighted)",
                   wl = "Cohen's Kappa (linearly weighted)",
                   wq = "Cohen's Kappa (quadratically weighted)")
    caption <- sprintf("Agreement table, %s = %.3f ", method, stat)
    biv_quali(x = x, y = y, 
              perc = 'none', 
              test = 'none', 
              useNA = 'no',
              caption = caption, 
              ...)
}
