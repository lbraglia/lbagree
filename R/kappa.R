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
        ## l'uso di data.frame affinché i factor non si perdano in una
        ## matrice nei casi particolari (in cui una modalità non sia
        ## rappresentata per un rater)
        x <- cbind(data.frame(x), y)
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
    var_order <- c('estimate', 'lower', 'upper')
    rval <- list(
        'table' = tab,
        'unweighted' = (rval_std$confid)[1, var_order],
        'weighted (linear)' = (rval_lin$confid)[2, var_order],
        'weighted (quadratic)' = (rval_std$confid)[2, var_order])
    
    rval
}

#' Fleiss' kappa
#'
#' Fleiss' kappa (Cohen's equivalent for multiple raters)
#' 
#' @param x a data.frame with one column per scale item/subject and
#'     one row per rater. Judgement can be a factors or integers
#' @param nlevels number of levels of possible rating/judgement
#' @param raters_par arguments passed to raters::concordance
#' @param irr_par arguments passed to irr::kappam.fleiss
#' @param gwet_par arguments passed to irrCAC::gwet.ac1.dist
#' @examples
#' ## reproducing the example from
#' ## https://en.wikipedia.org/wiki/Fleiss%27_kappa
#' item1  <- rep(5, 14)
#' item2  <- c(rep(1,0), rep(2,2), rep(3,6), rep(4,4), rep(5,2))
#' item3  <- c(rep(1,0), rep(2,0), rep(3,3), rep(4,5), rep(5,6))
#' item4  <- c(rep(1,0), rep(2,3), rep(3,9), rep(4,2), rep(5,0))
#' item5  <- c(rep(1,2), rep(2,2), rep(3,8), rep(4,1), rep(5,1))
#' item6  <- c(rep(1,7), rep(2,7), rep(3,0), rep(4,0), rep(5,0))
#' item7  <- c(rep(1,3), rep(2,2), rep(3,6), rep(4,3), rep(5,0))
#' item8  <- c(rep(1,2), rep(2,5), rep(3,3), rep(4,2), rep(5,2))
#' item9  <- c(rep(1,6), rep(2,5), rep(3,2), rep(4,1), rep(5,0))
#' item10 <- c(rep(1,0), rep(2,2), rep(3,2), rep(4,3), rep(5,7))
#' 
#' df <-  data.frame(item1, item2, item3, item4, item5,
#'                   item6, item7, item8, item9, item10)
#' df
#' fleiss_k(x = df, nlevels = 5)
#' 
#' @export
fleiss_k <- function(x,
                     nlevels = NULL,
                     raters_par = list(test = "MC", B = 1000, alpha = 0.05),
                     irr_par = list(exact = FALSE, detail = FALSE),
                     gwet_par = list(weights = "unweighted",
                                     categ = NULL,
                                     conflev = 0.95, N = Inf)
                     ){

    
    if (is.null(nlevels)) {
        stop("specify number of levels for the rating")
        ## tmp <- as.matrix(x)
        ## dim(tmp) <- NULL
        ## nlevels <- length(unique(tmp %without% NA))
    }

    ## todo. E' giusto eliminare un rater se ha dato missing ad un
    ## item o è più corretto eliminare un item? è più giusto un item

    ## eliminare rater
    not_na <- lbmisc::NA_remove(x, quiet = FALSE)

    type_norm <- function(y)
        factor(as.integer(y), levels = seq_len(nlevels))

    ## irr package results
    irr_params <-  c(list("ratings" = t(not_na)), irr_par)
    irr_res <- do.call(irr::kappam.fleiss, irr_params)
    
    ## raters package results
    data_norm <- as.data.frame(lapply(not_na, type_norm))
    freqs <- do.call(rbind, lapply(data_norm, table))
    raters_params <- c(list("db" = freqs), raters_par)
    capture.output(raters_res <- do.call(raters::concordance, raters_params))

    ## gwet AC1 (irrCAC package)
    gwet_params <- c(list("ratings" = freqs), gwet_par)
    gwet_res <- do.call(irrCAC::gwet.ac1.dist, gwet_params)
    
    ## reporting
    list('freqs' = freqs,
         'irr_res'    = irr_res,
         'raters_res' = raters_res,
         'gwet_ac1' = gwet_res)
    
}


