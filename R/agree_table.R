#' Calculate agreement between two nominal/ordinal variables
#' 
#' @param x first qualitative variable
#' @param y first second variable
#' @param cohen_k weight for cohen (uw = unweighted, 'wl = linearly weighted',
#' @param caption_prefix prefix used in the latex caption
#' @param conf.level confidence level for estimate confidence interval
#' @param ... further arguments passed to biv_quali
#' 
#' @export
agree_table <- function(x, y,
                        cohen_k = c('uw', 'wl', 'wq'),
                        caption_prefix = "",
                        conf.level = 0.95,  
                        ...) 
{
    cohen_k <- match.arg(cohen_k) 
    cohen <- lbagree::cohen_k(x = x, y = y, conf.level = conf.level)    
    method <- switch(cohen_k,
                     uw = "Cohen's Kappa (unweighted)",
                     wl = "Cohen's Kappa (linearly weighted)",
                     wq = "Cohen's Kappa (quadratically weighted)")
    stat <- switch(cohen_k,
                   uw = cohen$unweighted,
                   wl = cohen$`weighted (linear)`,
                   wq = cohen$`weighted (quadratic)`)
    caption <- sprintf("%s%s table, %s = %.3f (%.2f CI: %.3f to %.3f)",
                       caption_prefix,
                       if (caption_prefix == '') "Agreement" else " agreement",
                       method,
                       stat[1], # coefficient
                       conf.level,
                       stat[2], # low.ci
                       stat[3])  # up.ci
    lbstat::biv_quali(x = x, y = y, 
                      perc = 'none', 
                      test = 'none', 
                      useNA = 'no',
                      caption = caption, 
                      ...)
}
