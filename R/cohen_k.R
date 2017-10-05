#' Cohen k
#'
#' Cohen k (psych version)
#' 
#' @param x a table of agreement, a 2 columns structure or a single
#'     vector (first rater) if y is a vector too
#' @param y a vector (second rater) or NULL
#' @param conf.level confidence level for confidence intervals
#' @examples
#'
#' rater_a <- c(rep(1,50), rep(2,30), rep(3,20))
#' rater_b <- c(rep(1, 44), rep(2,  5), rep(3, 1),
#'        rep(1,  7), rep(2, 20), rep(3, 3),
#'        rep(1,  9), rep(2,  5), rep(3, 6))
#' tab <- table(rater_a, rater_b)
#' rate <- data.frame("A" = rater_a, "B" = rater_b)
#'
#' cohen_k(rater_a, rater_b)
#' cohen_k(tab)
#' cohen_k(rate)
#'
#' ## check with
#' ## ----------
#' ## irr::kappa2(rate)$value
#' ## irr::kappa2(rate, weight = "equal")$value
#' ## irr::kappa2(rate, weight = "squared")$value
#' ## psych::cohen.kappa(rate)
#' 
#' @export
cohen_k <- function(x = NULL,
                    y = NULL,
                    ## weights = c('none', 'linear', 'squared'),
                    conf.level = 0.95) 
{

    if( !is.null(y) ){
        if (! (is.atomic(x) && is.atomic(y)))
            stop("x and y must be atomic")
        if (length(x) != length(y))
            stop("x and y must have the same length")
        ## se passa tutti i controlli uniformiamo al caso con un
        ## data.frame o matrice su cui si basa il codice a seguire
        x <- cbind(x, y)
    } else if(! (ncol(x) ==2 || (is.table(x))))
        stop('if y is NULL x must be a table or a 2 columns data structure')
    
    if(conf.level >= 1 || conf.level <= 0)
        stop('conf.level must be between 0 and 1')

    ## weights <- match.arg(weights)
    tab <- if (is.table(x)) x else table(x[, 1], x[, 2])

    ## la accetta un tipo di peso simile a quello dato nell'esempio
    ## su http://vassarstats.net/kappaexp.html
    ## (anzi è l'unico modo che ho velocemente trovato affinché combacino
    ## in weigthed le seguenti stime anche per quanto riguarda l'intervallo
    ## di confidenza
    ## 
    ## cohen_k(x = tab, weights = 'none')
    ## cohen_k(x = tab, weights = 'squared')
    ## 
    ## ossia le stime col peso generato coincidano con le stime del peso 
    ## di default (che appunto è quadratico
   
    w <- function(n = NULL, type = NULL){
        db <- data.frame(seq_len(n), seq_len(n))
        d <- as.matrix(dist(db, method = 'maximum', diag = TRUE, upper = TRUE))
        if (type == 'linear') 
            1 - (abs(d) / max(d))
        else if (type == 'squared') 
            1 - (d^2 / (max(d)^2))
    }

    ## stima normale (non pesata e quadratica)
    rval_std <- psych::cohen.kappa(x = tab,
                                   w = NULL,
                                   n.obs = NULL,
                                   alpha = 1 - conf.level,
                                   levels = NULL)
    ## stima con pesi lineari
    rval_lin <- psych::cohen.kappa(x = tab,
                                   w = w(n = ncol(tab), type = 'linear'),
                                   n.obs = NULL,
                                   alpha = 1 - conf.level,
                                   levels = NULL)

    ## return value
    rval <- list('unweighted' = (rval_std$confid)[1, ],
                 'weighted (linear)' = (rval_lin$confid)[2, ],
                 'weighted (quadratic)' = (rval_std$confid)[2, ])
                 
    rval
}
