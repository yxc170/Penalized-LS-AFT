#' Function to fit...
#'
#' @param formula a formula object
#' @param data an optional data frame
#' @param subset an optional logical vector
#' @param id an optional
#' @param weights an optional
#' @param B a numeric value
#' @param lambda a numeric value
#' @param control a list
#'
#' @importFrom BB dfsane
#' @importFrom parallel makeCluster parSapply stopCluster clusterExport detectCores
#' @export
paftgee <- function(formula, data, subset, id = NULL, weights = NULL, lambda = NULL,
                    penalty = c("scad", "lasso", "alasso"), B = 0,
                    corstr = c("independence", "exchangeable", "AR1"),
                    control = list()) {
    penalty <- match.arg(penalty)
    corstr <- match.arg(corstr)
    Call <- match.call()
    if (missing(formula)) stop("Argument 'formula' is required.")
    if (missing(data)) data <- environment(formula)
    if (!missing(subset)) {
        sSubset <- substitute(subset)
        subIdx <- eval(sSubset, data, parent.frame())
        if (!is.logical(subIdx))
            stop("'subset' must be logical")
        subIdx <- subIdx & !is.na(subIdx)
        data <- data[subIdx, ]
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "id", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$data <- data
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    xname <- colnames(model.matrix(attr(mf, "terms"), mf))
    DF <- cbind(mf[,1], model.matrix(attr(mf, "terms"), mf))
    DF <- as.data.frame(DF)
    obj <- model.response(mf)
    if (!is.Surv(obj)) stop("Response must be a `Surv` object")
    if (attr(obj, "type") != "right") stop("Response must be a right censored object")
    ctrl <- paft.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl)))
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (ctrl$parallel & is.null(ctrl$parCl)) ctrl$parCl <- detectCores() / 2L
    if (any(names(mf) %in% "(id)")) id <- mf$"(id)"
    else id <- 1:nrow(DF)
    if (any(names(mf) %in% "(weights)")) w <- mf$"(weights)"
    else w <- rep(1, nrow(DF))
    if (is.null(ctrl$pindex)) {
        ctrl$pindex <- rep(1, ncol(DF) - 2)
        if (attr(attr(mf, "terms"), "intercept") > 0) ctrl$pindex[1] <- 0
    }
    DF[,1] <- log(DF[,1])
    y <- DF[,1]
    d <- DF[,2]
    x <- as.matrix(DF[,-(1:2)])
    if (is.character(ctrl$init)) b0 <- getInit(DF, ctrl$init, w)
    if (is.character(ctrl$init) && !any(c("gehan_s", "gehan_ns", "lm") == ctrl$init))
        stop("Undefined initial value type")
    if (is.numeric(ctrl$init)) {
        b0 <- ctrl$init
        if (length(b0) != ncol(x))
            stop("Initial value and the number of covaraites must be the same length")
    }
    output <- list(beta.init = b0)
    nt <- tabulate(id)
    if (missing(lambda)){#penalized model:
        foo <- paftgeeEst1(y = y, X = x, D = d, b0 = b0, nt = nt, pindex = ctrl$pindex,
                           w = w, penalty = penalty, nCV = ctrl$nfold, corstr = corstr,
                           lambda = ctrl$lambda.vec, CVprop = ctrl$CVprop, rule = ctrl$cv.rule, eps = 1e-6, maxit = ctrl$maxit, tol = ctrl$tol)
        output$beta <- foo$b
        output$lambdaBest <- foo$lambdaBest
        output$steps <- foo$iter
    }
    else if (lambda < 0) stop("Lambda must be a postive value")
    else if (lambda == 0) {#non-penalized model
        foo <- aftgeeEst(y = y, X = x, D = d, b0 = b0, nt = nt, w = w,
                            corstr = corstr, tol = ctrl$tol, maxit = ctrl$maxit)
        output$beta <- foo$b
        output$steps <- foo$iter
    }
    else if (lambda > 0) {#penalized model: lambda is already specified
        foo <- paftgeeEst(y = y, X = x, D = d, b0 = b0, nt = nt, pindex = ctrl$pindex,
                              w = w, penalty = penalty, corstr = corstr, lambda = lambda,
                              eps = 1e-6, maxit = ctrl$maxit, tol = ctrl$tol)
        output$beta <- foo$b
        output$steps <- foo$iter
        output$lambda <- lambda
    }
    ## Assumes id matches x and y, and is sorted
    return(output[order(names(output))])
}


#' Control list
#'
#' @export
paft.control <- function(tol = 1e-3, maxit = 500, init = NULL, nfold = 5,
                         cv.rule = c("cv", "1se", "bic"), lambda.vec = 1:100 / 50,
                         pindex = NULL, CVprop = FALSE, parallel = FALSE, parCl = NULL) {
    if (parallel & is.null(parCl)) parCl <- detectCores() / 2L
    cv.rule <- match.arg(cv.rule)
    if (is.null(init)) init <- "lm"
    list(tol = tol, maxit = maxit, init = init, nfold = nfold,
         cv.rule = cv.rule, lambda.vec = lambda.vec, pindex = pindex,
         CVprop = FALSE, parallel = parallel, parCl = parCl)
}


#' Function to calculate initial value
#' @noRd
getInit <- function(DF, init, w) {
    y <- DF[,1]
    d <- DF[,2]
    x <- as.matrix(DF[,-(1:2)])
    if (is.character(init)) {
        b0 <- unname(coef(lm(time ~ . - status - 1, subset = status > 0, data = DF, weights = w)))

        is.int <- colnames(x) == "(Intercept)"
        if (init == "gehan_s")
            b0 <- dfsane(b0[!is.int], fn = function(b) drop(gehan_s(b, x[,!is.int], y, d, w)),
                         quiet = TRUE, alertConvergence = FALSE)$par
        if (init == "gehan_ns")
            b0 <- dfsane(b0[!is.int], fn = function(b) drop(gehan_ns(b, x[,!is.int], y, d, w)),
                         quiet = TRUE, alertConvergence = FALSE)$par
        if (init != "lm" & sum(is.int) > 0)
            b0 <- c(mean(y - colSums(b0 * t(x[,!is.int]))), b0)
    }
    b0
}




