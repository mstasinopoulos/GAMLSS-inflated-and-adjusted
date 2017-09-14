require(gamlss)
gamlssInf0to1 <- function(y = NULL, mu.formula = ~1, sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1, 
    xi0.formula = ~1, xi1.formula = ~1, data = NULL, family = BE, weights = rep(1, length(Y_)), trace = FALSE, 
    ...) {
    call <- match.call()
    vars <- list(...)
    nopar <- eval(gamlss.family(family))$nopar
    Res <- if (!is.null(data)) 
        get(deparse(substitute(y)), envir = as.environment(data))
    else y
    WhatDist <- if (any(Res == 0) && !any(Res == 1)) 
        "Zero"
    else if (!any(Res == 0) && any(Res == 1)) 
        "One"
    else if (any(Res == 0) && any(Res == 1)) 
        "Zero&One"
    else stop("the response variable do not have zero or one")
    Y_ <- switch(WhatDist, Zero = (Res == 0) * 1 + (!Res == 0) * 0, One = (Res == 1) * 1 + (!Res == 1) * 
        0, `Zero&One` = (Res == 0) * 1 + (Res == 1) * 2 + (Res != 1) * (Res != 0) * 3)
    formula0 <- if (WhatDist == "Zero" || WhatDist == "Zero&One") {
        formula(paste("Y_", deparse(substitute(xi0.formula), width.cutoff = 500L), sep = ""))
    }
    else {
        formula(paste("Y_", deparse(substitute(xi1.formula), width.cutoff = 500L), sep = ""))
    }
    if (WhatDist == "Zero&One") {
        if (trace) 
            cat("******      The multinomial model        ******", "\n")
        m0 <- if (is.null(data)) {
            gamlss(formula0, sigma.formula = xi1.formula, weights = weights, family = MN3, trace = trace, 
                ...)
        }
        else {
            gamlss(formula0, sigma.formula = xi1.formula, weights = weights, data = data, family = MN3, trace = trace, 
                ...)
        }
    }
    if (WhatDist == "Zero" || WhatDist == "One") {
        if (trace) 
            cat("*****       The binomial model        *****", "\n")
        m0 <- if (is.null(data)) {
            gamlss(formula0, sigma.formula = xi1.formula, weights = weights, family = BI, trace = trace, 
                ...)
        }
        else {
            gamlss(formula0, sigma.formula = xi1.formula, weights = weights, data = data, family = BI, trace = trace, 
                ...)
        }
    }
    WW <- ifelse((Res != 1) & (Res != 0), 1, 0)
    WEIGHTS <- ifelse((Res != 1) & (Res != 0), 1, .Machine$double.eps^0.8)
    meanRes <- mean(Res)
    Res1 <- ifelse(WW, Res, meanRes)
    formulaM <- formula(paste("Res1", deparse(substitute(mu.formula), width.cutoff = 500L), sep = ""), envir = as.environment(data))
    if (trace) 
        cat("***** The continuous distribution model *****", "\n")
    mM <- if (is.null(data)) {
        gamlss(formulaM, sigma.formula = sigma.formula, nu.formula = nu.formula, tau.formula = tau.formula, 
            weights = WEIGHTS * weights, family = family, trace = trace, ...)
    }
    else {
        gamlss(formulaM, data = data, sigma.formula = sigma.formula, nu.formula = nu.formula, tau.formula = tau.formula, 
            weights = WEIGHTS * weights, family = family, trace = trace, ...)
    }
    mM$call$family <- substitute(family)
    m0$call$data <- substitute(data)
    mM$call$data <- substitute(data)
    m0$call$weights <- substitute(weights)
    mM$call$weights <- substitute(WEIGHTS)
    if ("mu" %in% mM$parameters) {
        formulaM[[2]] <- substitute(y)
        mM$call$formula <- formula(formulaM)
        mM$mu.terms[[2]] <- substitute(y)
        mM$mu.formula[[2]] <- substitute(y)
    }
    mM$call$sigma.formula <- formula(sigma.formula)
    mM$call$nu.formula <- formula(nu.formula)
    mM$call$tau.formula <- formula(tau.formula)
    m0$call$formula <- formula(formula0)
    m0$call$sigma.formula <- formula(xi1.formula)
    m0$mu.formula[[2]] <- substitute(y)
    fam <- paste("Inf", mM$family[1], sep = "")
    family.original <- mM$family[1]
    method <- mM$method
    parameters <- switch(WhatDist, Zero = c(mM$parameters, "xi0"), One = c(mM$parameters, "xi1"), `Zero&One` = c(mM$parameters, 
        "xi0", "xi1"))
    G.deviance <- deviance(m0) + deviance(mM)
    N <- length(Res)
    type <- "mixed"
    iter <- c(m0$iter, mM$iter)
    converged <- m0$converged & mM$converged
    noObs <- sum(weights)
    df.fit <- mM$df.fit + m0$df.fit
    df.residual <- noObs - df.fit
    aic <- G.deviance + 2 * df.fit
    sbc <- G.deviance + log(length(Res)) * df.fit
    out <- list(multinom = m0, dist = mM, family = fam, noObs = noObs, original.family = family.original, 
        method = method, parameters = parameters, call = call, y = Res, weights = weights, G.deviance = G.deviance, 
        N = N, type = type, typeInf = WhatDist, df.fit = df.fit, df.residual = df.residual, aic = aic, sbc = sbc)
    if ("mu" %in% parameters) {
        out$mu.fv <- mM$mu.fv
        out$mu.lp <- mM$mu.lp
        out$mu.wv <- mM$mu.wv
        out$mu.wt <- mM$mu.wt
        out$mu.link <- mM$mu.link
        out$mu.terms <- mM$mu.terms
        out$mu.x <- mM$mu.x
        out$mu.qr <- mM$mu.qr
        out$mu.coefficients <- mM$mu.coefficients
        out$mu.offset <- mM$mu.offset
        out$mu.xlevels <- mM$mu.xlevels
        out$mu.formula <- mM$mu.formula
        out$mu.df <- mM$mu.df
        out$mu.nl.df <- mM$mu.nl.df
        out$mu.pen <- mM$mu.pen
    }
    if ("sigma" %in% parameters) {
        out$sigma.fv <- mM$sigma.fv
        out$sigma.lp <- mM$sigma.lp
        out$sigma.wv <- mM$sigma.wv
        out$sigma.wt <- mM$sigma.wt
        out$sigma.link <- mM$sigma.link
        out$sigma.terms <- mM$sigma.terms
        out$sigma.x <- mM$sigma.x
        out$sigma.qr <- mM$sigma.qr
        out$sigma.coefficients <- mM$sigma.coefficients
        out$sigma.offset <- mM$sigma.offset
        out$sigma.xlevels <- mM$sigma.xlevels
        out$sigma.formula <- mM$sigma.formula
        out$sigma.df <- mM$sigma.df
        out$sigma.nl.df <- mM$sigma.nl.df
        out$sigma.pen <- mM$sigma.pen
    }
    if ("nu" %in% parameters) {
        out$nu.fv <- mM$nu.fv
        out$nu.lp <- mM$nu.lp
        out$nu.wv <- mM$nu.wv
        out$nu.wt <- mM$nu.wt
        out$nu.link <- mM$nu.link
        out$nu.terms <- mM$nu.terms
        out$nu.x <- mM$nu.x
        out$nu.qr <- mM$nu.qr
        out$nu.coefficients <- mM$nu.coefficients
        out$nu.offset <- mM$nu.offset
        out$nu.xlevels <- mM$nu.xlevels
        out$nu.formula <- mM$nu.formula
        out$nu.df <- mM$nu.df
        out$nu.nl.df <- mM$nu.nl.df
        out$nu.pen <- mM$nu.pen
    }
    if ("tau" %in% parameters) {
        out$tau.fv <- mM$tau.fv
        out$tau.lp <- mM$tau.lp
        out$tau.wv <- mM$tau.wv
        out$tau.wt <- mM$tau.wt
        out$tau.link <- mM$tau.link
        out$tau.terms <- mM$tau.terms
        out$tau.x <- mM$tau.x
        out$tau.qr <- mM$tau.qr
        out$tau.coefficients <- mM$tau.coefficients
        out$tau.offset <- mM$tau.offset
        out$tau.xlevels <- mM$tau.xlevels
        out$tau.formula <- mM$tau.formula
        out$tau.df <- mM$tau.df
        out$tau.nl.df <- mM$tau.nl.df
        out$tau.pen <- mM$tau.pen
    }
    if ("xi0" %in% parameters) {
        out$xi0.fv <- m0$mu.fv
        out$xi0.lp <- m0$mu.lp
        out$xi0.wv <- m0$mu.wv
        out$xi0.wt <- m0$mu.wt
        out$xi0.link <- m0$mu.link
        out$xi0.terms <- m0$mu.terms
        out$xi0.x <- m0$mu.x
        out$xi0.qr <- m0$mu.qr
        out$xi0.coefficients <- m0$mu.coefficients
        out$xi0.offset <- m0$mu.offset
        out$xi0.xlevels <- m0$mu.xlevels
        out$xi0.formula <- m0$mu.formula
        out$xi0.df <- m0$mu.df
        out$xi0.nl.df <- m0$mu.nl.df
        out$xi0.pen <- m0$mu.pen
    }
    if ("xi1" %in% parameters) {
        if (WhatDist == "Zero&One") {
            out$xi1.fv <- m0$sigma.fv
            out$xi1.lp <- m0$sigma.lp
            out$xi1.wv <- m0$sigma.wv
            out$xi1.wt <- m0$sigma.wt
            out$xi1.link <- m0$sigma.link
            out$xi1.terms <- m0$sigma.terms
            out$xi1.x <- m0$sigma.x
            out$xi1.qr <- m0$sigma.qr
            out$xi1.coefficients <- m0$sigma.coefficients
            out$xi1.offset <- m0$sigma.offset
            out$xi1.xlevels <- m0$sigma.xlevels
            out$xi1.formula <- m0$sigma.formula
            out$xi1.df <- m0$sigma.df
            out$xi1.nl.df <- m0$sigma.nl.df
            out$xi1.pen <- m0$sigma.pen
        }
        if (WhatDist == "One") {
            out$xi1.fv <- m0$mu.fv
            out$xi1.lp <- m0$mu.lp
            out$xi1.wv <- m0$mu.wv
            out$xi1.wt <- m0$mu.wt
            out$xi1.link <- m0$mu.link
            out$xi1.terms <- m0$mu.terms
            out$xi1.x <- m0$mu.x
            out$xi1.qr <- m0$mu.qr
            out$xi1.coefficients <- m0$mu.coefficients
            out$xi1.offset <- m0$mu.offset
            out$xi1.xlevels <- m0$mu.xlevels
            out$xi1.formula <- m0$mu.formula
            out$xi1.df <- m0$mu.df
            out$xi1.nl.df <- m0$mu.nl.df
            out$xi1.pen <- m0$mu.pen
        }
    }
    if (trace) 
        cat("            The Final Global Deviance  =", out$G.deviance, "\n")
    if (WhatDist == "Zero") {
        cdf <- Inf0to1.p(mM$family[1], type.of.Inflation = "Zero")
        prob <- m0$mu.fv
        uval <- ifelse(Res == 0, runif(length(Res), 0, prob), 0)
        uval <- switch(nopar, ifelse(Res > 0 & Res < 1, cdf(q = Res, mu = out$mu.fv, xi0 = out$xi0.fv), uval), 
            ifelse(Res > 0 & Res < 1, cdf(q = Res, mu = out$mu.fv, sigma = out$sigma.fv, xi0 = out$xi0.fv), 
                uval), ifelse(Res > 0 & Res < 1, cdf(q = Res, mu = out$mu.fv, sigma = out$sigma.fv, nu = out$nu.fv, 
                xi0 = out$xi0.fv), uval), ifelse(Res > 0 & Res < 1, cdf(q = Res, mu = out$mu.fv, sigma = out$sigma.fv, 
                nu = out$nu.fv, tau = out$tau.fv, xi0 = out$xi0.fv), uval))
    }
    if (WhatDist == "One") {
        cdf <- Inf0to1.p(mM$family[1], type.of.Inflation = "One")
        prob <- m0$mu.fv
        uval <- switch(nopar, ifelse(Res > 0 & Res < 1, cdf(q = Res, mu = out$mu.fv, xi1 = out$xi1.fv), 1), 
            ifelse(Res > 0 & Res < 1, cdf(q = Res, mu = out$mu.fv, sigma = out$sigma.fv, xi1 = out$xi1.fv), 
                1), ifelse(Res > 0 & Res < 1, cdf(q = Res, mu = out$mu.fv, sigma = out$sigma.fv, nu = out$nu.fv, 
                xi1 = out$xi1.fv), 1), ifelse(Res > 0 & Res < 1, cdf(q = Res, mu = out$mu.fv, sigma = out$sigma.fv, 
                nu = out$nu.fv, tau = out$tau.fv, xi1 = out$xi1.fv), 1))
        uval <- ifelse(Res == 1, runif(length(Res), 1 - prob, 1), uval)
    }
    if (WhatDist == "Zero&One") {
        cdf <- Inf0to1.p(mM$family[1], type.of.Inflation = "Zero&One")
        prob <- cbind(m0$mu.fv/(1 + m0$mu.fv + m0$sigma.fv), m0$sigma.fv/(1 + m0$mu.fv + m0$sigma.fv))
        uval <- ifelse(Res == 0, runif(length(Res), 0, prob[, 1]), 0)
        uval <- switch(nopar, ifelse(Res > 0 & Res < 1, cdf(q = Res, mu = out$mu.fv, xi0 = out$xi0.fv, xi1 = out$xi1.fv), 
            uval), ifelse(Res > 0 & Res < 1, cdf(q = Res, mu = out$mu.fv, sigma = out$sigma.fv, xi0 = out$xi0.fv, 
            xi1 = out$xi1.fv), uval), ifelse(Res > 0 & Res < 1, cdf(q = Res, mu = out$mu.fv, sigma = out$sigma.fv, 
            nu = out$nu.fv, xi0 = out$xi0.fv, xi1 = out$xi1.fv), uval), ifelse(Res > 0 & Res < 1, cdf(q = Res, 
            mu = out$mu.fv, sigma = out$sigma.fv, nu = out$nu.fv, tau = out$tau.fv, xi0 = out$xi0.fv, xi1 = out$xi1.fv), 
            uval))
        uval <- ifelse(Res == 1, runif(length(Res), 1 - prob[, 2], 1), uval)
    }
    out$residuals <- qnorm(uval)
    class(out) <- c("gamlssinf0to1", class(m0))
    out
}
fitted.gamlssinf0to1 <- function(object, parameter = c("mu", "sigma", "nu", "tau", "xi0", "xi1"), ...) {
    parameter <- match.arg(parameter)
    if (!parameter %in% object$par) 
        stop(paste(parameter, "is not a parameter in the gamlss object", "\n"))
    x <- if (is.null(object$na.action)) 
        object[[paste(parameter, "fv", sep = ".")]]
    else napredict(object$na.action, object[[paste(parameter, "fv", sep = ".")]])
    x
}
coef.gamlssinf0to1 <- function(object, parameter = c("mu", "sigma", "nu", "tau", "xi0", "xi1"), ...) {
    parameter <- match.arg(parameter)
    if (!parameter %in% object$par) 
        stop(paste(parameter, "is not a parameter in the object", "\n"))
    x <- object[[paste(parameter, "coefficients", sep = ".")]]
    x
}
print.gamlssinf0to1 <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nFamily: ", deparse(x$family), "\n")
    cat("Fitting method:", deparse(x$method), "\n")
    cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)
    cat("Mu Coefficients")
    if (is.character(co <- x$contrasts)) 
        cat("  [contrasts: ", apply(cbind(names(co), co), 1, paste, collapse = "="), "]")
    cat(":\n")
    if ("mu" %in% x$parameters) {
        print.default(format(coef(x, "mu"), digits = digits), print.gap = 2, quote = FALSE)
    }
    if ("sigma" %in% x$parameters) {
        cat("Sigma Coefficients:\n")
        print.default(format(coef(x, "sigma"), digits = digits), print.gap = 2, quote = FALSE)
    }
    if ("nu" %in% x$parameters) {
        cat("Nu Coefficients:\n")
        print.default(format(coef(x, "nu"), digits = digits), print.gap = 2, quote = FALSE)
    }
    if ("tau" %in% x$parameters) {
        cat("Tau Coefficients:\n")
        print.default(format(coef(x, "tau"), digits = digits), print.gap = 2, quote = FALSE)
    }
    if ("xi0" %in% x$parameters) {
        cat("xi0 Coefficients:\n")
        print.default(format(coef(x, "xi0"), digits = digits), print.gap = 2, quote = FALSE)
    }
    if ("xi1" %in% x$parameters) {
        cat("xi1 Coefficients:\n")
        print.default(format(coef(x, "xi1"), digits = digits), print.gap = 2, quote = FALSE)
    }
    cat("\n Degrees of Freedom for the fit:", x$df.fit, "Residual Deg. of Freedom  ", x$df.residual, "\n")
    cat("Global Deviance:    ", format(signif(x$G.deviance)), "\n            AIC:    ", format(signif(x$aic)), 
        "\n            SBC:    ", format(signif(x$sbc)), "\n")
    invisible(x)
}
deviance.gamlssinf0to1 <- function(object, ...) {
    x <- object$G.deviance
    x
}
vcov.gamlssinf0to1 <- function(object, type = c("vcov", "cor", "se", "coef", "all"), robust = FALSE, hessian.fun = c("R", 
    "PB"), ...) {
    Mdiag <- function(x) {
        if (is.null(x)) 
            M <- NULL
        else {
            x <- x[!sapply(x, is.null)]
            dimlist <- sapply(x, dim)
            a <- apply(dimlist, 1, cumsum)
            dimMat <- rowSums(dimlist)
            M <- array(0, dimMat)
            indexdim <- rbind(c(0, 0), a)
            for (i in 1:length(x)) M[(indexdim[i, 1] + 1):indexdim[i + 1, 1], (indexdim[i, 2] + 1):indexdim[i + 
                1, 2]] <- x[[i]]
        }
        M
    }
    type <- match.arg(type)
    hessian.fun <- match.arg(hessian.fun)
    M1 <- vcov(object$dist, robust = robust, hessian.fun = hessian.fun)
    M2 <- vcov(object$multin, robust = robust, hessian.fun = hessian.fun)
    VCOV <- Mdiag(list(M1, M2))
    colnames(VCOV) <- c(rownames(M1), rownames(M2))
    rownames(VCOV) <- c(rownames(M1), rownames(M2))
    coef <- unlist(lapply(object$parameter, function(x) coef(object, parameter = x)))
    se <- sqrt(diag(VCOV))
    corr <- cov2cor(VCOV)
    switch(type, vcov = VCOV, cor = corr, se = se, coef = coef, all = list(se = se, vcov = VCOV, coef = coef, 
        cor = corr))
}
summary.gamlssinf0to1 <- function(object, type = c("vcov", "qr"), robust = FALSE, save = FALSE, hessian.fun = c("R", 
    "PB"), digits = max(3, getOption("digits") - 3), ...) {
    type <- match.arg(type)
    pm <- ps <- pn <- pt <- px0 <- px1 <- 0
    mu.coef.table <- NULL
    sigma.coef.table <- NULL
    nu.coef.table <- NULL
    tau.coef.table <- NULL
    xi0.coef.table <- NULL
    xi1.coef.table <- NULL
    if (type == "vcov") {
        covmat <- try(suppressWarnings(vcov.gamlssinf0to1(object, type = "all", robust = robust, hessian.fun = hessian.fun)), 
            silent = TRUE)
        if (any(class(covmat) %in% "try-error" || any(is.na(covmat$se)))) {
            warning(paste("summary: vcov has failed, option qr is used instead\n"))
            type <- "qr"
        }
    }
    ifWarning <- rep(FALSE, length(object$parameters))
    if (type == "vcov") {
        coef <- covmat$coef
        se <- covmat$se
        tvalue <- coef/se
        pvalue <- 2 * pt(-abs(tvalue), object$df.res)
        coef.table <- cbind(coef, se, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        cat("*******************************************************************")
        cat("\nFamily: ", deparse(object$family), "\n")
        cat("\nCall: ", deparse(object$call), "\n", fill = TRUE)
        cat("Fitting method:", deparse(object$method), "\n\n")
        est.disp <- FALSE
        if ("mu" %in% object$parameters) {
            ifWarning[1] <- (!is.null(unlist(attr(terms(formula(object, "mu"), specials = .gamlss.sm.list), 
                "specials"))))
            if (object$mu.df != 0) {
                pm <- object$mu.qr$rank
                p1 <- 1:pm
                cat("-------------------------------------------------------------------\n")
                cat("Mu link function: ", object$mu.link)
                cat("\n")
                cat("Mu Coefficients:")
                if (is.character(co <- object$contrasts)) 
                  cat("  [contrasts: ", apply(cbind(names(co), co), 1, paste, collapse = "="), "]")
                cat("\n")
                printCoefmat(coef.table[p1, , drop = FALSE], digits = digits, signif.stars = TRUE)
                cat("\n")
            }
            else if (object$mu.fix == TRUE) {
                cat("-------------------------------------------------------------------\n")
                cat("Mu parameter is fixed \n")
                if (all(object$mu.fv == object$mu.fv[1])) 
                  cat("Mu = ", object$mu.fv[1], "\n")
                else cat("Mu is equal with the vector (", object$mu.fv[1], ",", object$mu.fv[2], ",", object$mu.fv[3], 
                  ",", object$mu.fv[4], ", ...) \n")
            }
        }
        if ("sigma" %in% object$parameters) {
            ifWarning[2] <- (!is.null(unlist(attr(terms(formula(object, "sigma"), specials = .gamlss.sm.list), 
                "specials"))))
            if (object$sigma.df != 0) {
                ps <- object$sigma.qr$rank
                p1 <- (pm + 1):(pm + ps)
                cat("-------------------------------------------------------------------\n")
                cat("Sigma link function: ", object$sigma.link)
                cat("\n")
                cat("Sigma Coefficients:")
                cat("\n")
                printCoefmat(coef.table[p1, , drop = FALSE], digits = digits, signif.stars = TRUE)
                cat("\n")
            }
            else if (object$sigma.fix == TRUE) {
                cat("-------------------------------------------------------------------\n")
                cat("Sigma parameter is fixed")
                cat("\n")
                if (all(object$sigma.fv == object$sigma.fv[1])) 
                  cat("Sigma = ", object$sigma.fv[1], "\n")
                else cat("Sigma is equal with the vector (", object$sigma.fv[1], ",", object$sigma.fv[2], 
                  ",", object$sigma.fv[3], ",", object$sigma.fv[4], ", ...) \n")
            }
        }
        if ("nu" %in% object$parameters) {
            ifWarning[3] <- (!is.null(unlist(attr(terms(formula(object, "nu"), specials = .gamlss.sm.list), 
                "specials"))))
            if (object$nu.df != 0) {
                pn <- object$nu.qr$rank
                p1 <- (pm + ps + 1):(pm + ps + pn)
                cat("-------------------------------------------------------------------\n")
                cat("Nu link function: ", object$nu.link, "\n")
                cat("Nu Coefficients:")
                cat("\n")
                printCoefmat(coef.table[p1, , drop = FALSE], digits = digits, signif.stars = TRUE)
                cat("\n")
            }
            else if (object$nu.fix == TRUE) {
                cat("-------------------------------------------------------------------\n")
                cat("Nu parameter is fixed")
                cat("\n")
                if (all(object$nu.fv == object$nu.fv[1])) 
                  cat("Nu = ", object$nu.fv[1], "\n")
                else cat("Nu is equal with the vector (", object$nu.fv[1], ",", object$nu.fv[2], ",", object$nu.fv[3], 
                  ",", object$nu.fv[4], ", ...) \n")
            }
        }
        if ("tau" %in% object$parameters) {
            ifWarning[4] <- (!is.null(unlist(attr(terms(formula(object, "tau"), specials = .gamlss.sm.list), 
                "specials"))))
            if (object$tau.df != 0) {
                pt <- object$tau.qr$rank
                p1 <- (pm + ps + pn + 1):(pm + ps + pn + pt)
                cat("-------------------------------------------------------------------\n")
                cat("Tau link function: ", object$tau.link, "\n")
                cat("Tau Coefficients:")
                cat("\n")
                printCoefmat(coef.table[p1, , drop = FALSE], digits = digits, signif.stars = TRUE)
                cat("\n")
            }
            else if (object$tau.fix == TRUE) {
                cat("-------------------------------------------------------------------\n")
                cat("Tau parameter is fixed")
                cat("\n")
                if (all(object$tau.fv == object$tau.fv[1])) 
                  cat("Tau = ", object$tau.fv[1], "\n")
                else cat("Tau is equal with the vector (", object$tau.fv[1], ",", object$tau.fv[2], ",", 
                  object$tau.fv[3], ",", object$tau.fv[4], ", ...) \n")
            }
        }
        if ("xi0" %in% object$parameters) {
            ifWarning[4] <- (!is.null(unlist(attr(terms(formula(object, "xi0"), specials = .gamlss.sm.list), 
                "specials"))))
            if (object$xi0.df != 0) {
                px0 <- object$xi0.qr$rank
                p1 <- (pm + ps + pn + pt + 1):(pm + ps + pn + pt + px0)
                cat("-------------------------------------------------------------------\n")
                cat("xi0 link function: ", object$xi0.link, "\n")
                cat("xi0 Coefficients:")
                cat("\n")
                printCoefmat(coef.table[p1, , drop = FALSE], digits = digits, signif.stars = TRUE)
                cat("\n")
            }
            else if (object$xi0.fix == TRUE) {
                cat("-------------------------------------------------------------------\n")
                cat("xi0 parameter is fixed")
                cat("\n")
                if (all(object$xi0.fv == object$xi0.fv[1])) 
                  cat("xi0 = ", object$xi0.fv[1], "\n")
                else cat("xi0 is equal with the vector (", object$xi0.fv[1], ",", object$xi0.fv[2], ",", 
                  object$xi0.fv[3], ",", object$xi0.fv[4], ", ...) \n")
            }
        }
        else px0 <- 0
        if ("xi1" %in% object$parameters) {
            ifWarning[4] <- (!is.null(unlist(attr(terms(formula(object, "xi1"), specials = .gamlss.sm.list), 
                "specials"))))
            if (object$xi1.df != 0) {
                px1 <- object$xi1.qr$rank
                p1 <- (pm + ps + pn + pt + px0 + 1):(pm + ps + pn + pt + px0 + px1)
                cat("-------------------------------------------------------------------\n")
                cat("xi1 link function: ", object$xi1.link, "\n")
                cat("xi1 Coefficients:")
                cat("\n")
                printCoefmat(coef.table[p1, , drop = FALSE], digits = digits, signif.stars = TRUE)
                cat("\n")
            }
            else if (object$xi1.fix == TRUE) {
                cat("-------------------------------------------------------------------\n")
                cat("xi1 parameter is fixed")
                cat("\n")
                if (all(object$xi1.fv == object$xi1.fv[1])) 
                  cat("xi1 = ", object$xi1.fv[1], "\n")
                else cat("xi1 is equal with the vector (", object$xi1.fv[1], ",", object$xi1.fv[2], ",", 
                  object$xi1.fv[3], ",", object$xi1.fv[4], ", ...) \n")
            }
        }
        if (any(ifWarning)) {
            cat("-------------------------------------------------------------------\n")
            cat("NOTE: Additive smoothing terms exist in the formulas: \n")
            cat(" i) Std. Error for smoothers are for the linear effect only. \n")
            cat("ii) Std. Error for the linear terms maybe are not accurate. \n")
        }
        cat("-------------------------------------------------------------------\n")
        cat("No. of observations in the fit: ", object$noObs, "\n")
        cat("Degrees of Freedom for the fit: ", object$df.fit)
        cat("\n")
        cat("      Residual Deg. of Freedom: ", object$df.residual, "\n")
        cat("                      at cycle: ", object$iter, "\n \n")
        cat("Global Deviance:    ", object$G.deviance, "\n            AIC:    ", object$aic, "\n            SBC:    ", 
            object$sbc, "\n")
        cat("*******************************************************************")
        cat("\n")
    }
    if (type == "qr") {
        estimatesgamlss <- function(object, Qr, p1, coef.p, est.disp, df.r, digits = max(3, getOption("digits") - 
            3), covmat.unscaled, ...) {
            dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
            covmat <- covmat.unscaled
            var.cf <- diag(covmat)
            s.err <- sqrt(var.cf)
            tvalue <- coef.p/s.err
            dn <- c("Estimate", "Std. Error")
            if (!est.disp) {
                pvalue <- 2 * pnorm(-abs(tvalue))
                coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
                dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", "Pr(>|z|)"))
            }
            else if (df.r > 0) {
                pvalue <- 2 * pt(-abs(tvalue), df.r)
                coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
                dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", "Pr(>|t|)"))
            }
            else {
                coef.table <- cbind(coef.p, Inf)
                dimnames(coef.table) <- list(names(coef.p), dn)
            }
            return(coef.table)
        }
        dispersion <- NULL
        cat("*******************************************************************")
        cat("\nFamily: ", deparse(object$family), "\n")
        cat("\nCall: ", deparse(object$call), "\n", fill = TRUE)
        cat("Fitting method:", deparse(object$method), "\n\n")
        est.disp <- FALSE
        df.r <- object$noObs - object$mu.df
        if ("mu" %in% object$parameters) {
            ifWarning[1] <- (!is.null(unlist(attr(terms(formula(object, "mu"), specials = .gamlss.sm.list), 
                "specials"))))
            if (object$mu.df != 0) {
                Qr <- object$mu.qr
                df.r <- object$noObs - object$mu.df
                if (is.null(dispersion)) 
                  dispersion <- if (any(object$family == c("PO", "BI", "EX", "P1"))) 
                    1
                  else if (df.r > 0) {
                    est.disp <- TRUE
                    if (any(object$weights == 0)) 
                      warning(paste("observations with zero weight", "not used for calculating dispersion"))
                  }
                  else Inf
                p <- object$mu.df
                p1 <- 1:(p - object$mu.nl.df)
                covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
                mu.coef.table <- estimatesgamlss(object = object, Qr = object$mu.qr, p1 = p1, coef.p = object$mu.coefficients[Qr$pivot[p1]], 
                  est.disp = est.disp, df.r = df.r, covmat.unscaled = covmat.unscaled)
                cat("-------------------------------------------------------------------\n")
                cat("Mu link function: ", object$mu.link)
                cat("\n")
                cat("Mu Coefficients:")
                if (is.character(co <- object$contrasts)) 
                  cat("  [contrasts: ", apply(cbind(names(co), co), 1, paste, collapse = "="), "]")
                cat("\n")
                printCoefmat(mu.coef.table, digits = digits, signif.stars = TRUE)
                cat("\n")
            }
            else if (object$mu.fix == TRUE) {
                cat("-------------------------------------------------------------------\n")
                cat("Mu parameter is fixed")
                cat("\n")
                if (all(object$mu.fv == object$mu.fv[1])) 
                  cat("Mu = ", object$mu.fv[1], "\n")
                else cat("Mu is equal with the vector (", object$mu.fv[1], ",", object$mu.fv[2], ",", object$mu.fv[3], 
                  ",", object$mu.fv[4], ", ...) \n")
            }
            coef.table <- mu.coef.table
        }
        else {
            if (df.r > 0) {
                est.disp <- TRUE
                if (any(object$weights == 0)) 
                  warning(paste("observations with zero weight", "not used for calculating dispersion"))
            }
        }
        if ("sigma" %in% object$parameters) {
            ifWarning[2] <- (!is.null(unlist(attr(terms(formula(object, "sigma"), specials = .gamlss.sm.list), 
                "specials"))))
            if (object$sigma.df != 0) {
                Qr <- object$sigma.qr
                df.r <- object$noObs - object$sigma.df
                p <- object$sigma.df
                p1 <- 1:(p - object$sigma.nl.df)
                covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
                sigma.coef.table <- estimatesgamlss(object = object, Qr = object$sigma.qr, p1 = p1, coef.p = object$sigma.coefficients[Qr$pivot[p1]], 
                  est.disp = est.disp, df.r = df.r, covmat.unscaled = covmat.unscaled)
                cat("-------------------------------------------------------------------\n")
                cat("Sigma link function: ", object$sigma.link)
                cat("\n")
                cat("Sigma Coefficients:")
                cat("\n")
                printCoefmat(sigma.coef.table, digits = digits, signif.stars = TRUE)
                cat("\n")
            }
            else if (object$sigma.fix == TRUE) {
                cat("-------------------------------------------------------------------\n")
                cat("Sigma parameter is fixed")
                cat("\n")
                if (all(object$sigma.fv == object$sigma.fv[1])) 
                  cat("Sigma = ", object$sigma.fv[1], "\n")
                else cat("Sigma is equal with the vector (", object$sigma.fv[1], ",", object$sigma.fv[2], 
                  ",", object$sigma.fv[3], ",", object$sigma.fv[4], ", ...) \n")
            }
            coef.table <- rbind(mu.coef.table, sigma.coef.table)
        }
        if ("nu" %in% object$parameters) {
            ifWarning[3] <- (!is.null(unlist(attr(terms(formula(object, "nu"), specials = .gamlss.sm.list), 
                "specials"))))
            if (object$nu.df != 0) {
                Qr <- object$nu.qr
                df.r <- object$noObs - object$nu.df
                p <- object$nu.df
                p1 <- 1:(p - object$nu.nl.df)
                covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
                nu.coef.table <- estimatesgamlss(object = object, Qr = object$nu.qr, p1 = p1, coef.p = object$nu.coefficients[Qr$pivot[p1]], 
                  est.disp = est.disp, df.r = df.r, covmat.unscaled = covmat.unscaled)
                cat("-------------------------------------------------------------------\n")
                cat("Nu link function: ", object$nu.link, "\n")
                cat("Nu Coefficients:")
                cat("\n")
                printCoefmat(nu.coef.table, digits = digits, signif.stars = TRUE)
                cat("\n")
            }
            else if (object$nu.fix == TRUE) {
                cat("-------------------------------------------------------------------\n")
                cat("Nu parameter is fixed")
                cat("\n")
                if (all(object$nu.fv == object$nu.fv[1])) 
                  cat("Nu = ", object$nu.fv[1], "\n")
                else cat("Nu is equal with the vector (", object$nu.fv[1], ",", object$nu.fv[2], ",", object$nu.fv[3], 
                  ",", object$nu.fv[4], ", ...) \n")
            }
            coef.table <- rbind(mu.coef.table, sigma.coef.table, nu.coef.table)
        }
        if ("tau" %in% object$parameters) {
            ifWarning[4] <- (!is.null(unlist(attr(terms(formula(object, "tau"), specials = .gamlss.sm.list), 
                "specials"))))
            if (object$tau.df != 0) {
                Qr <- object$tau.qr
                df.r <- object$noObs - object$tau.df
                p <- object$tau.df
                p1 <- 1:(p - object$tau.nl.df)
                covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
                tau.coef.table <- estimatesgamlss(object = object, Qr = object$tau.qr, p1 = p1, coef.p = object$tau.coefficients[Qr$pivot[p1]], 
                  est.disp = est.disp, df.r = df.r, covmat.unscaled = covmat.unscaled)
                cat("-------------------------------------------------------------------\n")
                cat("Tau link function: ", object$tau.link, "\n")
                cat("Tau Coefficients:")
                cat("\n")
                printCoefmat(tau.coef.table, digits = digits, signif.stars = TRUE)
                cat("\n")
            }
            else if (object$tau.fix == TRUE) {
                cat("-------------------------------------------------------------------\n")
                cat("Tau parameter is fixed")
                cat("\n")
                if (all(object$tau.fv == object$tau.fv[1])) 
                  cat("Tau = ", object$tau.fv[1], "\n")
                else cat("Tau is equal with the vector (", object$tau.fv[1], ",", object$tau.fv[2], ",", 
                  object$tau.fv[3], ",", object$tau.fv[4], ", ...) \n")
            }
            coef.table <- rbind(mu.coef.table, sigma.coef.table, nu.coef.table, tau.coef.table)
        }
        if ("xi0" %in% object$parameters) {
            ifWarning[4] <- (!is.null(unlist(attr(terms(formula(object, "xi0"), specials = .gamlss.sm.list), 
                "specials"))))
            if (object$xi0.df != 0) {
                Qr <- object$xi0.qr
                df.r <- object$noObs - object$xi0.df
                p <- object$xi0.df
                p1 <- 1:(p - object$xi0.nl.df)
                covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
                xi0.coef.table <- estimatesgamlss(object = object, Qr = object$xi0.qr, p1 = p1, coef.p = object$xi0.coefficients[Qr$pivot[p1]], 
                  est.disp = est.disp, df.r = df.r, covmat.unscaled = covmat.unscaled)
                cat("-------------------------------------------------------------------\n")
                cat("xi0 link function: ", object$xi0.link, "\n")
                cat("xi0 Coefficients:")
                cat("\n")
                printCoefmat(xi0.coef.table, digits = digits, signif.stars = TRUE)
                cat("\n")
            }
            else if (object$xi0.fix == TRUE) {
                cat("-------------------------------------------------------------------\n")
                cat("xi0 parameter is fixed")
                cat("\n")
                if (all(object$xi0.fv == object$xi0.fv[1])) 
                  cat("xi0 = ", object$xi0.fv[1], "\n")
                else cat("xi0 is equal with the vector (", object$xi0.fv[1], ",", object$xi0.fv[2], ",", 
                  object$xi0.fv[3], ",", object$xi0.fv[4], ", ...) \n")
            }
            coef.table <- rbind(mu.coef.table, sigma.coef.table, nu.coef.table, tau.coef.table, xi0.coef.table)
        }
        if (!("xi0" %in% object$parameters)) 
            xi0.coef.table <- NULL
        if ("xi1" %in% object$parameters) {
            ifWarning[4] <- (!is.null(unlist(attr(terms(formula(object, "xi1"), specials = .gamlss.sm.list), 
                "specials"))))
            if (object$xi1.df != 0) {
                Qr <- object$xi1.qr
                df.r <- object$noObs - object$xi1.df
                p <- object$xi1.df
                p1 <- 1:(p - object$xi1.nl.df)
                covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
                xi1.coef.table <- estimatesgamlss(object = object, Qr = object$xi1.qr, p1 = p1, coef.p = object$xi1.coefficients[Qr$pivot[p1]], 
                  est.disp = est.disp, df.r = df.r, covmat.unscaled = covmat.unscaled)
                cat("-------------------------------------------------------------------\n")
                cat("xi1 link function: ", object$xi1.link, "\n")
                cat("xi1 Coefficients:")
                cat("\n")
                printCoefmat(xi1.coef.table, digits = digits, signif.stars = TRUE)
                cat("\n")
            }
            else if (object$xi1.fix == TRUE) {
                cat("-------------------------------------------------------------------\n")
                cat("xi1 parameter is fixed")
                cat("\n")
                if (all(object$xi1.fv == object$xi1.fv[1])) 
                  cat("xi1 = ", object$xi1.fv[1], "\n")
                else cat("xi1 is equal with the vector (", object$xi1.fv[1], ",", object$xi1.fv[2], ",", 
                  object$xi1.fv[3], ",", object$xi1.fv[4], ", ...) \n")
            }
            coef.table <- rbind(mu.coef.table, sigma.coef.table, nu.coef.table, tau.coef.table, xi0.coef.table, 
                xi1.coef.table)
        }
        if (!("xi0" %in% object$parameters)) 
            xi0.coef.table <- NULL
        if (any(ifWarning)) {
            cat("-------------------------------------------------------------------\n")
            cat("NOTE: Additive smoothing terms exist in the formulas: \n")
            cat(" i) Std. Error for smoothers are for the linear effect only. \n")
            cat("ii) Std. Error for the linear terms may not be reliable. \n")
        }
        cat("-------------------------------------------------------------------\n")
        cat("No. of observations in the fit: ", object$noObs, "\n")
        cat("Degrees of Freedom for the fit: ", object$df.fit)
        cat("\n")
        cat("      Residual Deg. of Freedom: ", object$df.residual, "\n")
        cat("                      at cycle: ", object$iter, "\n \n")
        cat("Global Deviance:    ", object$G.deviance, "\n            AIC:    ", object$aic, "\n            SBC:    ", 
            object$sbc, "\n")
        cat("*******************************************************************")
        cat("\n")
    }
    if (save == TRUE) {
        out <- as.list(environment())
        return(out)
    }
    invisible(coef.table)
}
formula.gamlssinf0to1 <- function(x, parameter = c("mu", "sigma", "nu", "tau", "xi0", "xi1"), ...) {
    parameter <- match.arg(parameter)
    if (!parameter %in% x$par) 
        stop(paste(parameter, "is not a parameter in the object", "\n"))
    fo <- x[[paste(parameter, "formula", sep = ".")]]
    if (length(fo) == 2 && "." %in% strsplit(as.character(fo), split = "")[[2]]) 
        fo <- formula(x[[paste(parameter, "terms", sep = ".")]])
    if (length(fo) == 3 && "." %in% strsplit(as.character(fo), split = "")[[3]]) 
        fo <- formula(x[[paste(parameter, "terms", sep = ".")]])
    fo
}
predict.gamlssinf0to1 <- function(object, parameter = c("mu", "sigma", "nu", "tau", "xi0", "xi1"), newdata = NULL, 
    type = c("link", "response", "terms"), terms = NULL, se.fit = FALSE, data = NULL, ...) {
    par <- object$par
    parameter <- match.arg(parameter)
    type <- match.arg(type)
    if (!(parameter %in% par)) 
        stop("the asking parameter is not fitted")
    if (is.null(newdata)) {
        whetherBIorMN3 <- any(is.na((match(c("xi0", "xi1"), par))))
        if (parameter == "mu") 
            thepar <- predict(object$dist, "mu", type = type, terms = terms, se.fit = se.fit, data = data)
        if (parameter == "sigma") 
            thepar <- predict(object$dist, "sigma", type = type, terms = terms, se.fit = se.fit, data = data)
        if (parameter == "nu") 
            thepar <- predict(object$dist, "nu", type = type, terms = terms, se.fit = se.fit, data = data)
        if (parameter == "tau") 
            thepar <- predict(object$dist, "tau", type = type, terms = terms, se.fit = se.fit, data = data)
        if (parameter == "xi0") 
            thepar <- predict(object$multin, "mu", type = type, terms = terms, se.fit = se.fit, data = data)
        if (parameter == "xi1") {
            if (whetherBIorMN3) {
                thepar <- predict(object$multin, "mu", type = type, terms = terms, se.fit = se.fit, data = data)
            }
            else {
                thepar <- predict(object$multin, "sigma", type = type, terms = terms, se.fit = se.fit, data = data)
            }
        }
        return(thepar)
    }
    else {
        whetherBIorMN3 <- any(is.na((match(c("xi0", "xi1"), par))))
        if (parameter == "mu") 
            thepar <- predict(object$dist, "mu", newdata = newdata, type = type, terms = terms, se.fit = se.fit, 
                data = data)
        if (parameter == "sigma") 
            thepar <- predict(object$dist, "sigma", newdata = newdata, type = type, terms = terms, se.fit = se.fit, 
                data = data)
        if (parameter == "nu") 
            thepar <- predict(object$dist, "nu", newdata = newdata, type = type, terms = terms, se.fit = se.fit, 
                data = data)
        if (parameter == "tau") 
            thepar <- predict(object$dist, "tau", newdata = newdata, type = type, terms = terms, se.fit = se.fit, 
                data = data)
        if (parameter == "xi0") 
            thepar <- predict(object$multin, "mu", newdata = newdata, type = type, terms = terms, se.fit = se.fit, 
                data = data)
        if (parameter == "xi1") {
            if (whetherBIorMN3) {
                thepar <- predict(object$multin, "mu", newdata = newdata, type = type, terms = terms, se.fit = se.fit, 
                  data = data)
            }
            else {
                thepar <- predict(object$multin, "sigma", newdata = newdata, type = type, terms = terms, 
                  se.fit = se.fit, data = data)
            }
        }
        return(thepar)
    }
}
Inf0to1.d <- function(family = "BE", type.of.Inflation = c("Zero&One", "Zero", "One"), ...) {
    xi0 <- xi1 <- mu <- sigma <- nu <- tau <- 1
    typeInf <- match.arg(type.of.Inflation)
    fname <- family
    if (mode(family) != "character" && mode(family) != "name") 
        fname <- as.character(substitute(family))
    distype <- eval(gamlss.family(family))$type
    if (!grepl("tr", fname)) {
        if (!body(eval(gamlss.family(family))$y.valid) == "all(y > 0 & y < 1)") 
            stop("the function is not defined on 0 to 1")
    }
    nopar <- eval(gamlss.family(family))$nopar
    if (!distype == "Continuous") 
        stop("the family should be continuous")
    FindOrig <- ifelse(grepl("logit", fname), 1, 0)
    FindOrig <- ifelse(grepl("tr", fname), 2, FindOrig)
    if (FindOrig == 1) {
        nameFam <- paste("d", sub("logit", "", fname), sep = "")
        orig.dFam <- get(nameFam)
    }
    if (FindOrig == 2) {
        nameFam <- paste("d", gsub("tr", "", fname), sep = "")
        orig.dFam <- get(nameFam)
    }
    dfun <- paste("d", fname, sep = "")
    pdf <- eval(parse(text = dfun))
    fun <- function(x, log = FALSE, ...) {
        if (typeInf == "Zero") {
            if (any(xi0 <= 0)) 
                stop(paste("xi0 must greated than 0", "\n", ""))
        }
        if (typeInf == "One") {
            if (any(xi1 <= 0)) 
                stop(paste("xi1 must greated than 0", "\n", ""))
        }
        if (typeInf == "Zero&One") {
            if (any(xi0 <= 0)) 
                stop(paste("xi0 must greated than 0", "\n", ""))
            if (any(xi1 <= 0)) 
                stop(paste("xi1 must greated than 0", "\n", ""))
        }
        logfy <- ifelse((x > 0 & x < 1), switch(nopar, pdf(ifelse((x == 0 | x == 1), 0.5, x), mu, log = TRUE), 
            pdf(ifelse((x == 0 | x == 1), 0.5, x), mu, sigma, log = TRUE), pdf(ifelse((x == 0 | x == 1), 
                0.5, x), mu, sigma, nu, log = TRUE), pdf(ifelse((x == 0 | x == 1), 0.5, x), mu, sigma, nu, 
                tau, log = TRUE)), 0)
        if (typeInf == "Zero") {
            logfy <- ifelse((x == 0), log(xi0), log(1 - xi0) + logfy)
            logfy <- ifelse((x < 0) | (x >= 1), NA, logfy)
        }
        if (typeInf == "One") {
            logfy <- ifelse((x == 1), log(xi1), log(1 - xi1) + logfy)
            logfy <- ifelse((x <= 0) | (x > 1), NA, logfy)
        }
        if (typeInf == "Zero&One") {
            logfy <- ifelse((x == 0), log(xi0), logfy)
            logfy <- ifelse((x == 1), log(xi1), logfy)
            logfy <- logfy - log(1 + xi0 + xi1)
            logfy <- ifelse((x < 0) | (x > 1), NA, logfy)
        }
        if (log == FALSE) 
            fy <- exp(logfy)
        else fy <- logfy
        fy
    }
    if (typeInf == "Zero") {
        if (FindOrig == 0) {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(x = NULL, mu = formals(pdf)[[2]], 
                xi0 = 0.2, log = FALSE), list(x = NULL, mu = formals(pdf)[[2]], sigma = formals(pdf)[[3]], 
                xi0 = 0.2, log = FALSE), list(x = NULL, mu = formals(pdf)[[2]], sigma = formals(pdf)[[3]], 
                nu = formals(pdf)[[4]], xi0 = 0.2, log = FALSE), list(x = NULL, mu = formals(pdf)[[2]], sigma = formals(pdf)[[3]], 
                nu = formals(pdf)[[4]], tau = formals(pdf)[[5]], xi0 = 0.2, log = FALSE))
        }
        else {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(x = NULL, mu = formals(orig.dFam)[[2]], 
                xi0 = 0.2, log = FALSE), list(x = NULL, mu = formals(orig.dFam)[[2]], sigma = formals(orig.dFam)[[3]], 
                xi0 = 0.2, log = FALSE), list(x = NULL, mu = formals(orig.dFam)[[2]], sigma = formals(orig.dFam)[[3]], 
                nu = formals(orig.dFam)[[4]], xi0 = 0.2, log = FALSE), list(x = NULL, mu = formals(orig.dFam)[[2]], 
                sigma = formals(orig.dFam)[[3]], nu = formals(orig.dFam)[[4]], tau = formals(orig.dFam)[[5]], 
                xi0 = 0.2, log = FALSE))
        }
    }
    if (typeInf == "One") {
        if (FindOrig == 0) {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(x = NULL, mu = formals(pdf)[[2]], 
                xi1 = 0.2, log = FALSE), list(x = NULL, mu = formals(pdf)[[2]], sigma = formals(pdf)[[3]], 
                xi1 = 0.2, log = FALSE), list(x = NULL, mu = formals(pdf)[[2]], sigma = formals(pdf)[[3]], 
                nu = formals(pdf)[[4]], xi1 = 0.2, log = FALSE), list(x = NULL, mu = formals(pdf)[[2]], sigma = formals(pdf)[[3]], 
                nu = formals(pdf)[[4]], tau = formals(pdf)[[5]], xi1 = 0.2, log = FALSE))
        }
        else {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(x = NULL, mu = formals(orig.dFam)[[2]], 
                xi1 = 0.2, log = FALSE), list(x = NULL, mu = formals(orig.dFam)[[2]], sigma = formals(orig.dFam)[[3]], 
                xi1 = 0.2, log = FALSE), list(x = NULL, mu = formals(orig.dFam)[[2]], sigma = formals(orig.dFam)[[3]], 
                nu = formals(orig.dFam)[[4]], xi1 = 0.2, log = FALSE), list(x = NULL, mu = formals(orig.dFam)[[2]], 
                sigma = formals(orig.dFam)[[3]], nu = formals(orig.dFam)[[4]], tau = formals(orig.dFam)[[5]], 
                xi1 = 0.2, log = FALSE))
        }
    }
    if (typeInf == "Zero&One") {
        if (FindOrig == 0) {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(x = NULL, mu = formals(pdf)[[2]], 
                xi0 = 0.1, xi1 = 0.1, log = FALSE), list(x = NULL, mu = formals(pdf)[[2]], sigma = formals(pdf)[[3]], 
                xi0 = 0.1, xi1 = 0.1, log = FALSE), list(x = NULL, mu = formals(pdf)[[2]], sigma = formals(pdf)[[3]], 
                nu = formals(pdf)[[4]], xi0 = 0.1, xi1 = 0.1, log = FALSE), list(x = NULL, mu = formals(pdf)[[2]], 
                sigma = formals(pdf)[[3]], nu = formals(pdf)[[4]], tau = formals(pdf)[[5]], xi0 = 0.1, xi1 = 0.1, 
                log = FALSE))
        }
        else {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(x = NULL, mu = formals(orig.dFam)[[2]], 
                xi0 = 0.1, xi1 = 0.1, log = FALSE), list(x = NULL, mu = formals(orig.dFam)[[2]], sigma = formals(orig.dFam)[[3]], 
                xi0 = 0.1, xi1 = 0.1, log = FALSE), list(x = NULL, mu = formals(orig.dFam)[[2]], sigma = formals(orig.dFam)[[3]], 
                nu = formals(orig.dFam)[[4]], xi0 = 0.1, xi1 = 0.1, log = FALSE), list(x = NULL, mu = formals(orig.dFam)[[2]], 
                sigma = formals(orig.dFam)[[3]], nu = formals(orig.dFam)[[4]], tau = formals(orig.dFam)[[5]], 
                xi0 = 0.1, xi1 = 0.1, log = FALSE))
        }
    }
    fun
}
Inf0to1.p <- function(family = "BE", type.of.Inflation = c("Zero&One", "Zero", "One"), ...) {
    xi0 = xi1 = mu = sigma = nu = tau = 1
    lower.tail = log.p = FALSE
    typeInf <- match.arg(type.of.Inflation)
    fname <- family
    if (mode(family) != "character" && mode(family) != "name") 
        fname <- as.character(substitute(family))
    distype <- eval(gamlss.family(family))$type
    if (!grepl("tr", fname)) {
        if (!body(eval(gamlss.family(family))$y.valid) == "all(y > 0 & y < 1)") 
            stop("the function is not defined on 0 to 1")
    }
    nopar <- eval(gamlss.family(family))$nopar
    if (!distype == "Continuous") 
        stop("the family should be continuous")
    FindOrig <- ifelse(grepl("logit", fname), 1, 0)
    FindOrig <- ifelse(grepl("tr", fname), 2, FindOrig)
    if (FindOrig == 1) {
        nameFam <- paste("p", sub("logit", "", fname), sep = "")
        orig.pFam <- get(nameFam)
    }
    if (FindOrig == 2) {
        nameFam <- paste("p", gsub("tr", "", fname), sep = "")
        orig.pFam <- get(nameFam)
    }
    pfun <- paste("p", fname, sep = "")
    cdf <- eval(parse(text = pfun))
    fun <- function(q, log = FALSE, ...) {
        if (typeInf == "Zero") {
            if (any((xi0 <= 0) | (xi0 >= 1))) 
                stop(paste("xi0 must between  0 and 1", "\n", ""))
        }
        if (typeInf == "One") {
            if (any((xi1 <= 0) | (xi1 >= 1))) 
                stop(paste("xi1 must between  0 and 1", "\n", ""))
        }
        if (typeInf == "Zero&One") {
            if (any(xi0 <= 0)) 
                stop(paste("xi0 must greated than 0", "\n", ""))
            if (any(xi1 <= 0)) 
                stop(paste("xi1 must greated than 0", "\n", ""))
        }
        cdfy <- ifelse((q > 0 & q < 1), switch(nopar, cdf(ifelse((q == 0 | q == 1), 0.5, q), mu, lower.tail = TRUE, 
            log.p = FALSE), cdf(ifelse((q == 0 | q == 1), 0.5, q), mu, sigma, lower.tail = TRUE, log.p = FALSE), 
            cdf(ifelse((q == 0 | q == 1), 0.5, q), mu, sigma, nu, lower.tail = TRUE, log.p = FALSE), cdf(ifelse((q == 
                0 | q == 1), 0.5, q), mu, sigma, nu, tau, lower.tail = TRUE, log.p = FALSE), ), 0)
        if (typeInf == "Zero") {
            cdfy <- ifelse((q == 0), xi0, xi0 + (1 - xi0) * cdfy)
            cdfy <- ifelse((q < 0) | (q >= 1), NA, cdfy)
        }
        if (typeInf == "One") {
            cdfy <- ifelse((q == 1), 1, (1 - xi1) * cdfy)
            cdfy <- ifelse((q <= 0) | (q > 1), NA, cdfy)
        }
        if (typeInf == "Zero&One") {
            cdfy <- ifelse((q == 0), xi0, ifelse((q == 1), 1 + xi0 + xi1, xi0 + cdfy))
            cdfy <- (cdfy)/(1 + xi0 + xi1)
            cdfy <- ifelse((q < 0) | (q > 1), NA, cdfy)
        }
        if (lower.tail == TRUE) 
            cdfy <- cdfy
        else cdfy <- 1 - cdfy
        if (log.p == FALSE) 
            cdfy <- cdfy
        else cdfy <- log(cdfy)
        cdfy
    }
    if (typeInf == "Zero") {
        if (FindOrig == 0) {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(q = NULL, mu = formals(cdf)[[2]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(cdf)[[2]], sigma = formals(cdf)[[3]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(cdf)[[2]], sigma = formals(cdf)[[3]], 
                nu = formals(pdf)[[4]], xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(cdf)[[2]], 
                sigma = formals(cdf)[[3]], nu = formals(cdf)[[4]], tau = formals(cdf)[[5]], xi0 = 0.1, lower.tail = TRUE, 
                log.p = FALSE))
        }
        else {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(q = NULL, mu = formals(orig.pFam)[[2]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, 
                mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], tau = formals(orig.pFam)[[5]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE))
        }
    }
    if (typeInf == "One") {
        if (FindOrig == 0) {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(q = NULL, mu = formals(cdf)[[2]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(cdf)[[2]], sigma = formals(cdf)[[3]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(cdf)[[2]], sigma = formals(cdf)[[3]], 
                nu = formals(pdf)[[4]], xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(cdf)[[2]], 
                sigma = formals(cdf)[[3]], nu = formals(cdf)[[4]], tau = formals(cdf)[[5]], xi1 = 0.1, lower.tail = TRUE, 
                log.p = FALSE))
        }
        else {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(q = NULL, mu = formals(orig.pFam)[[2]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, 
                mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], tau = formals(orig.pFam)[[5]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE))
        }
    }
    if (typeInf == "Zero&One") {
        if (FindOrig == 0) {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(q = NULL, mu = formals(cdf)[[2]], 
                xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(cdf)[[2]], 
                sigma = formals(cdf)[[3]], xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, 
                mu = formals(cdf)[[2]], sigma = formals(cdf)[[3]], nu = formals(pdf)[[4]], xi0 = 0.1, xi1 = 0.1, 
                lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(cdf)[[2]], sigma = formals(cdf)[[3]], 
                nu = formals(cdf)[[4]], tau = formals(cdf)[[5]], xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, 
                log.p = FALSE))
        }
        else {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(q = NULL, mu = formals(orig.pFam)[[2]], 
                xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), 
                list(q = NULL, mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], 
                  xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(q = NULL, mu = formals(orig.pFam)[[2]], 
                  sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], tau = formals(orig.pFam)[[5]], 
                  xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE))
        }
    }
    fun
}
Inf0to1.q <- function(family = "BE", type.of.Inflation = c("Zero&One", "Zero", "One"), ...) {
    xi0 = xi1 = mu = sigma = nu = tau = 1
    lower.tail = log.p = FALSE
    typeInf <- match.arg(type.of.Inflation)
    fname <- family
    if (mode(family) != "character" && mode(family) != "name") 
        fname <- as.character(substitute(family))
    distype <- eval(gamlss.family(family))$type
    if (!grepl("tr", fname)) {
        if (!body(eval(gamlss.family(family))$y.valid) == "all(y > 0 & y < 1)") 
            stop("the function is not defined on 0 to 1")
    }
    nopar <- eval(gamlss.family(family))$nopar
    if (!distype == "Continuous") 
        stop("the family should be continuous")
    FindOrig <- ifelse(grepl("logit", fname), 1, 0)
    FindOrig <- ifelse(grepl("tr", fname), 2, FindOrig)
    if (FindOrig == 1) {
        nameFam <- paste("q", sub("logit", "", fname), sep = "")
        orig.pFam <- get(nameFam)
    }
    if (FindOrig == 2) {
        nameFam <- paste("q", gsub("tr", "", fname), sep = "")
        orig.pFam <- get(nameFam)
    }
    qfun <- paste("q", fname, sep = "")
    inv <- eval(parse(text = qfun))
    fun <- function(p, log = FALSE, ...) {
        if (typeInf == "Zero") {
            if (any((xi0 <= 0) | (xi0 >= 1))) 
                stop(paste("xi0 must between  0 and 1", "\n", ""))
        }
        if (typeInf == "One") {
            if (any((xi1 <= 0) | (xi1 >= 1))) 
                stop(paste("xi1 must between  0 and 1", "\n", ""))
        }
        if (typeInf == "Zero&One") {
            if (any(xi0 <= 0)) 
                stop(paste("xi0 must greated than 0", "\n", ""))
            if (any(xi1 <= 0)) 
                stop(paste("xi1 must greated than 0", "\n", ""))
        }
        if (typeInf == "Zero") {
            p_xi0 <- ifelse((p - xi0)/(1 - xi0) <= 0, 0.5, (p - xi0)/(1 - xi0))
            q <- ifelse((p > xi0), switch(nopar, inv(p_xi0, mu, lower.tail = TRUE, log.p = FALSE), inv(p_xi0, 
                mu, sigma, lower.tail = TRUE, log.p = FALSE), inv(p_xi0, mu, sigma, nu, lower.tail = TRUE, 
                log.p = FALSE), inv(p_xi0, mu, sigma, nu, tau, lower.tail = TRUE, log.p = FALSE)), 0)
        }
        if (typeInf == "One") {
            p_xi1 <- ifelse(p > (1 - xi1), 0.5, p/(1 - xi1))
            q <- ifelse((p <= 1 - xi1), switch(nopar, inv(p_xi1, mu, lower.tail = TRUE, log.p = FALSE), inv(p_xi1, 
                mu, sigma, lower.tail = TRUE, log.p = FALSE), inv(p_xi1, mu, sigma, nu, lower.tail = TRUE, 
                log.p = FALSE), inv(p_xi1, mu, sigma, nu, tau, lower.tail = TRUE, log.p = FALSE)), 1)
        }
        if (typeInf == "Zero&One") {
            q <- rep(0, length(p))
            ratio <- xi0/(1 + xi0 + xi1)
            q[p > ratio] <- 1
            middle <- (p > ratio) & (p < ((1 + xi0)/(1 + xi0 + xi1)))
            pm <- (p - ratio)/(1/(1 + xi0 + xi1))
            if (!is.null(mu)) 
                .mu <- if (length(mu) == 1) 
                  mu
                else mu[middle]
            if (!is.null(sigma)) 
                .sigma <- if (length(sigma) == 1) 
                  sigma
                else sigma[middle]
            if (!is.null(nu)) 
                .nu <- if (length(nu) == 1) 
                  nu
                else nu[middle]
            if (!is.null(tau)) 
                .tau <- if (length(tau) == 1) 
                  tau
                else tau[middle]
            q[middle] <- switch(nopar, inv(pm[middle], .mu, lower.tail = TRUE, log.p = FALSE), inv(pm[middle], 
                .mu, .sigma, lower.tail = TRUE, log.p = FALSE), inv(pm[middle], .mu, .sigma, .nu, lower.tail = TRUE, 
                log.p = FALSE), inv(pm[middle], .mu, .sigma, .nu, .tau, lower.tail = TRUE, log.p = FALSE))
        }
        if (lower.tail == TRUE) 
            q <- q
        else q <- 1 - q
        if (log.p == FALSE) 
            q <- q
        else q <- log(q)
        q
    }
    if (typeInf == "Zero") {
        if (FindOrig == 0) {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(p = NULL, mu = formals(inv)[[2]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(inv)[[2]], sigma = formals(inv)[[3]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(inv)[[2]], sigma = formals(inv)[[3]], 
                nu = formals(inv)[[4]], xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(inv)[[2]], 
                sigma = formals(inv)[[3]], nu = formals(inv)[[4]], tau = formals(inv)[[5]], xi0 = 0.1, lower.tail = TRUE, 
                log.p = FALSE))
        }
        else {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(p = NULL, mu = formals(orig.pFam)[[2]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, 
                mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], tau = formals(orig.pFam)[[5]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE))
        }
    }
    if (typeInf == "One") {
        if (FindOrig == 0) {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(p = NULL, mu = formals(inv)[[2]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(inv)[[2]], sigma = formals(inv)[[3]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(inv)[[2]], sigma = formals(inv)[[3]], 
                nu = formals(inv)[[4]], xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(inv)[[2]], 
                sigma = formals(inv)[[3]], nu = formals(inv)[[4]], tau = formals(inv)[[5]], xi1 = 0.1, lower.tail = TRUE, 
                log.p = FALSE))
        }
        else {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(p = NULL, mu = formals(orig.pFam)[[2]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, 
                mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], tau = formals(orig.pFam)[[5]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE))
        }
    }
    if (typeInf == "Zero&One") {
        if (FindOrig == 0) {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(p = NULL, mu = formals(inv)[[2]], 
                xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(inv)[[2]], 
                sigma = formals(inv)[[3]], xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, 
                mu = formals(inv)[[2]], sigma = formals(inv)[[3]], nu = formals(inv)[[4]], xi0 = 0.1, xi1 = 0.1, 
                lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(inv)[[2]], sigma = formals(inv)[[3]], 
                nu = formals(inv)[[4]], tau = formals(inv)[[5]], xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, 
                log.p = FALSE))
        }
        else {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(p = NULL, mu = formals(orig.pFam)[[2]], 
                xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), 
                list(p = NULL, mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], 
                  xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(p = NULL, mu = formals(orig.pFam)[[2]], 
                  sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], tau = formals(orig.pFam)[[5]], 
                  xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE))
        }
    }
    fun
}
Inf0to1.r <- function(family = "BE", type.of.Inflation = c("Zero&One", "Zero", "One"), ...) {
    xi0 = xi1 = mu = sigma = nu = tau = 1
    typeInf <- match.arg(type.of.Inflation)
    fname <- family
    if (mode(family) != "character" && mode(family) != "name") 
        fname <- as.character(substitute(family))
    distype <- eval(gamlss.family(family))$type
    nopar <- eval(gamlss.family(family))$nopar
    if (!distype == "Continuous") 
        stop("the family should be continuous")
    FindOrig <- ifelse(grepl("logit", fname), 1, 0)
    FindOrig <- ifelse(grepl("tr", fname), 2, FindOrig)
    if (FindOrig == 1) {
        nameFam <- paste("r", sub("logit", "", fname), sep = "")
        orig.pFam <- get(nameFam)
    }
    if (FindOrig == 2) {
        nameFam <- paste("r", gsub("tr", "", fname), sep = "")
        orig.pFam <- get(nameFam)
    }
    rfun <- Inf0to1.q(family = family, type.of.Inflation = typeInf)
    fun <- function(x, log = FALSE, ...) {
        n <- ceiling(n)
        p <- runif(n)
        if (typeInf == "Zero") {
            r <- switch(nopar, rfun(p, mu, xi0, lower.tail = TRUE, log.p = FALSE), rfun(p, mu, sigma, xi0, 
                lower.tail = TRUE, log.p = FALSE), rfun(p, mu, sigma, nu, xi0, lower.tail = TRUE, log.p = FALSE), 
                rfun(p, mu, sigma, nu, tau, xi0, lower.tail = TRUE, log.p = FALSE))
        }
        if (typeInf == "One") {
            r <- switch(nopar, rfun(p, mu, xi1, lower.tail = TRUE, log.p = FALSE), rfun(p, mu, sigma, xi1, 
                lower.tail = TRUE, log.p = FALSE), rfun(p, mu, sigma, nu, xi1, lower.tail = TRUE, log.p = FALSE), 
                rfun(p, mu, sigma, nu, tau, xi1, lower.tail = TRUE, log.p = FALSE))
        }
        if (typeInf == "Zero&One") {
            r <- switch(nopar, rfun(p, mu, xi0, xi1, lower.tail = TRUE, log.p = FALSE), rfun(p, mu, sigma, 
                xi0, xi1, lower.tail = TRUE, log.p = FALSE), rfun(p, mu, sigma, nu, xi0, xi1, lower.tail = TRUE, 
                log.p = FALSE), rfun(p, mu, sigma, nu, tau, xi0, xi1, lower.tail = TRUE, log.p = FALSE))
        }
        r
    }
    if (typeInf == "Zero") {
        if (FindOrig == 0) {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(n = NULL, mu = formals(rfun)[[2]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(rfun)[[2]], sigma = formals(rfun)[[3]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(rfun)[[2]], sigma = formals(rfun)[[3]], 
                nu = formals(rfun)[[4]], xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(rfun)[[2]], 
                sigma = formals(rfun)[[3]], nu = formals(rfun)[[4]], tau = formals(rfun)[[5]], xi0 = 0.1, 
                lower.tail = TRUE, log.p = FALSE))
        }
        else {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(n = NULL, mu = formals(orig.pFam)[[2]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, 
                mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], tau = formals(orig.pFam)[[5]], 
                xi0 = 0.1, lower.tail = TRUE, log.p = FALSE))
        }
    }
    if (typeInf == "One") {
        if (FindOrig == 0) {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(n = NULL, mu = formals(rfun)[[2]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(rfun)[[2]], sigma = formals(rfun)[[3]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(rfun)[[2]], sigma = formals(rfun)[[3]], 
                nu = formals(rfun)[[4]], xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(rfun)[[2]], 
                sigma = formals(rfun)[[3]], nu = formals(rfun)[[4]], tau = formals(rfun)[[5]], xi1 = 0.1, 
                lower.tail = TRUE, log.p = FALSE))
        }
        else {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(n = NULL, mu = formals(orig.pFam)[[2]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, 
                mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], tau = formals(orig.pFam)[[5]], 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE))
        }
    }
    if (typeInf == "Zero&One") {
        if (FindOrig == 0) {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(n = NULL, mu = formals(rfun)[[2]], 
                xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(rfun)[[2]], 
                sigma = formals(rfun)[[3]], xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, 
                mu = formals(rfun)[[2]], sigma = formals(rfun)[[3]], nu = formals(rfun)[[4]], xi0 = 0.1, 
                xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(rfun)[[2]], sigma = formals(rfun)[[3]], 
                nu = formals(rfun)[[4]], tau = formals(rfun)[[5]], xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, 
                log.p = FALSE))
        }
        else {
            formals(fun, envir = environment(fun)) <- switch(nopar, list(n = NULL, mu = formals(orig.pFam)[[2]], 
                xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), 
                list(n = NULL, mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], 
                  xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE), list(n = NULL, mu = formals(orig.pFam)[[2]], 
                  sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], tau = formals(orig.pFam)[[5]], 
                  xi0 = 0.1, xi1 = 0.1, lower.tail = TRUE, log.p = FALSE))
        }
    }
    fun
}
gen.Inf0to1 <- function(family = "BE", type.of.Inflation = c("Zero&One", "Zero", "One"), ...) {
    xi0 <- xi1 <- mu <- sigma <- nu <- tau <- 1
    type <- switch(match.arg(type.of.Inflation), `Zero&One` = "0to1", Zero = "0", One = "1")
    fam <- as.gamlss.family(family)
    fname <- fam$family[[1]]
    dfun <- paste(paste("d", fname, "Inf", sep = ""), type, sep = "")
    pfun <- paste(paste("p", fname, "Inf", sep = ""), type, sep = "")
    qfun <- paste(paste("q", fname, "Inf", sep = ""), type, sep = "")
    rfun <- paste(paste("r", fname, "Inf", sep = ""), type, sep = "")
    plotfun <- paste(paste("plot", fname, "Inf", sep = ""), type, sep = "")
    alldislist <- c(dfun, pfun, qfun, rfun, plotfun)
    eval(dummy <- Inf0to1.d(family = family, type.of.Inflation = match.arg(type.of.Inflation), ...))
    eval(call("<-", as.name(dfun), dummy), envir = parent.frame(n = 1))
    eval(dummy <- Inf0to1.p(family = family, type.of.Inflation = match.arg(type.of.Inflation), ...))
    eval(call("<-", as.name(pfun), dummy), envir = parent.frame(n = 1))
    eval(dummy <- Inf0to1.q(family = family, type.of.Inflation = match.arg(type.of.Inflation), ...))
    eval(call("<-", as.name(qfun), dummy), envir = parent.frame(n = 1))
    eval(dummy <- Inf0to1.r(family = family, type.of.Inflation = match.arg(type.of.Inflation), ...))
    eval(call("<-", as.name(rfun), dummy), envir = parent.frame(n = 1))
    eval(dummy <- plotInf0to1(family = family, type.of.Inflation = match.arg(type.of.Inflation), ...))
    eval(call("<-", as.name(plotfun), dummy), envir = parent.frame(n = 1))
    cat("A ", type, "inflated", fname, "distribution has been generated \n", "and saved under the names: ", 
        "\n", paste(alldislist[1:4], sep = ","), "\n", paste(alldislist[5], sep = ","), "\n")
}
plotInf0to1 <- function(family = "BE", type.of.Inflation = c("Zero&One", "Zero", "One"), ...) {
    from <- to <- n <- xi0 <- xi1 <- mu <- sigma <- nu <- tau <- 1
    lower.tail <- log.p <- FALSE
    typeInf <- match.arg(type.of.Inflation)
    fname <- family
    if (mode(family) != "character" && mode(family) != "name") 
        fname <- as.character(substitute(family))
    distype <- eval(gamlss.family(family))$type
    nopar <- eval(gamlss.family(family))$nopar
    if (!distype == "Continuous") 
        stop("the family should be continuous")
    FindOrig <- ifelse(grepl("logit", fname), 1, 0)
    FindOrig <- ifelse(grepl("tr", fname), 2, FindOrig)
    if (FindOrig == 1) {
        nameFam <- paste("p", sub("logit", "", fname), sep = "")
        orig.pFam <- get(nameFam)
    }
    if (FindOrig == 2) {
        nameFam <- paste("p", gsub("tr", "", fname), sep = "")
        orig.pFam <- get(nameFam)
    }
    dfun <- Inf0to1.d(family = family, type.of.Inflation = typeInf)
    fun <- function(x, mu = 0.5, sigma = 0.5, nu = 1, tau = 1, xi0 = 0.1, xi1 = 0.1, from = 0.001, to = 0.999, 
        n = 101, log = FALSE, ...) {
        if (typeInf == "Zero") {
            fy <- switch(nopar, dfun(x = x, mu = mu, xi0 = xi0, log = log), dfun(x = x, mu = mu, sigma = sigma, 
                xi0 = xi0, log = log), dfun(x = x, mu = mu, sigma = sigma, nu = nu, xi0 = xi0, log = log), 
                dfun(x = x, mu = mu, sigma = sigma, nu = nu, tau = tau, xi0 = xi0, log = log))
        }
        if (typeInf == "One") {
            fy <- switch(nopar, dfun(x = x, mu = mu, xi1 = xi1, log = log), dfun(x = x, mu = mu, sigma = sigma, 
                xi1 = xi1, log = log), dfun(x = x, mu = mu, sigma = sigma, nu = nu, xi1 = xi1, log = log), 
                dfun(x = x, mu = mu, sigma = sigma, nu = nu, tau = tau, xi1 = xi1, log = log))
        }
        if (typeInf == "Zero&One") {
            fy <- switch(nopar, dfun(x = x, mu = mu, xi0 = xi0, xi1 = xi1, log = log), dfun(x = x, mu = mu, 
                sigma = sigma, xi0 = xi0, xi1 = xi1, log = log), dfun(x = x, mu = mu, sigma = sigma, nu = nu, 
                xi0 = xi0, xi1 = xi1, log = log), dfun(x = x, mu = mu, sigma = sigma, nu = nu, tau = tau, 
                xi0 = xi0, xi1 = xi1, log = log))
        }
        fy
    }
    plotfun <- function(...) {
        fy <- fun(x = seq(from, to, length = n), mu = mu, sigma = sigma, nu = nu, tau = tau, xi0 = xi0, xi1 = xi1, 
            log = log)
        maxfy <- max(fy)
        if (typeInf == "Zero") {
            pr <- fun(0, mu = mu, sigma = sigma, nu = nu, tau = tau, xi0 = xi0, xi1 = xi1, log = log)
            po <- 0
        }
        if (typeInf == "One") {
            pr <- fun(1, mu = mu, sigma = sigma, nu = nu, tau = tau, xi0 = xi0, xi1 = xi1, log = log)
            po <- 1
        }
        if (typeInf == "Zero&One") {
            pr <- c(fun(0, mu = mu, sigma = sigma, nu = nu, tau = tau, xi0 = xi0, xi1 = xi1, log = log), 
                fun(1, mu = mu, sigma = sigma, nu = nu, tau = tau, xi0 = xi0, xi1 = xi1, log = log))
            po <- c(0, 1)
        }
        allmax <- max(na.omit(c(pr, maxfy)))
        plot(function(y) fun(y, mu = mu, sigma = sigma, nu = nu, tau = tau, xi0 = xi0, xi1 = xi1, log = log), 
            from = from, to = to, n = n, ylim = c(0, allmax), ylab = "density", ...)
        points(po, pr, type = "h")
        points(po, pr, type = "p", col = "blue")
    }
    if (typeInf == "Zero") {
        if (FindOrig == 0) {
            formals(plotfun, envir = environment(plotfun)) <- switch(nopar, list(mu = formals(dfun)[[2]], 
                xi0 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE), list(mu = formals(dfun)[[2]], 
                sigma = formals(dfun)[[3]], xi0 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, 
                log = FALSE), list(mu = formals(dfun)[[2]], sigma = formals(dfun)[[3]], nu = formals(dfun)[[4]], 
                xi0 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE), list(mu = formals(dfun)[[2]], 
                sigma = formals(dfun)[[3]], nu = formals(dfun)[[4]], tau = formals(dfun)[[5]], xi0 = 0.1, 
                n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE))
        }
        else {
            formals(plotfun, envir = environment(plotfun)) <- switch(nopar, list(mu = formals(orig.pFam)[[2]], 
                xi0 = 0.1, n = 101, from = 0.001, to = 0.999, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, 
                log = FALSE), list(mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], xi0 = 0.1, 
                n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE), list(mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], xi0 = 0.1, n = 101, from = 0.001, 
                to = 0.999, lower.tail = TRUE, log = FALSE), list(mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], 
                nu = formals(orig.pFam)[[4]], tau = formals(orig.pFam)[[5]], xi0 = 0.1, n = 101, from = 0.001, 
                to = 0.999, lower.tail = TRUE, log = FALSE))
        }
    }
    if (typeInf == "One") {
        if (FindOrig == 0) {
            formals(plotfun, envir = environment(plotfun)) <- switch(nopar, list(mu = formals(dfun)[[2]], 
                xi1 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE), list(mu = formals(dfun)[[2]], 
                sigma = formals(dfun)[[3]], xi1 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, 
                log = FALSE), list(mu = formals(dfun)[[2]], sigma = formals(dfun)[[3]], nu = formals(dfun)[[4]], 
                xi1 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE), list(mu = formals(dfun)[[2]], 
                sigma = formals(dfun)[[3]], nu = formals(dfun)[[4]], tau = formals(dfun)[[5]], xi1 = 0.1, 
                n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE))
        }
        else {
            formals(plotfun, envir = environment(plotfun)) <- switch(nopar, list(mu = formals(orig.pFam)[[2]], 
                xi1 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE), list(mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], xi1 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, 
                log = FALSE), list(mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], 
                xi1 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE), list(mu = formals(orig.pFam)[[2]], 
                sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], tau = formals(orig.pFam)[[5]], 
                xi1 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE))
        }
    }
    if (typeInf == "Zero&One") {
        if (FindOrig == 0) {
            formals(plotfun, envir = environment(plotfun)) <- switch(nopar, list(mu = formals(dfun)[[2]], 
                xi0 = 0.1, xi1 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE), 
                list(mu = formals(dfun)[[2]], sigma = formals(dfun)[[3]], xi0 = 0.1, xi1 = 0.1, n = 101, 
                  from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE), list(mu = formals(dfun)[[2]], 
                  sigma = formals(dfun)[[3]], nu = formals(dfun)[[4]], xi0 = 0.1, xi1 = 0.1, n = 101, from = 0.001, 
                  to = 0.999, lower.tail = TRUE, log = FALSE), list(mu = formals(dfun)[[2]], sigma = formals(dfun)[[3]], 
                  nu = formals(dfun)[[4]], tau = formals(dfun)[[5]], xi0 = 0.1, xi1 = 0.1, n = 101, from = 0.001, 
                  to = 0.999, lower.tail = TRUE, log = FALSE))
        }
        else {
            formals(plotfun, envir = environment(plotfun)) <- switch(nopar, list(mu = formals(orig.pFam)[[2]], 
                xi0 = 0.1, xi1 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE), 
                list(mu = formals(orig.pFam)[[2]], sigma = formals(orig.pFam)[[3]], xi0 = 0.1, xi1 = 0.1, 
                  n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE), list(mu = formals(orig.pFam)[[2]], 
                  sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], xi0 = 0.1, xi1 = 0.1, n = 101, 
                  from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE), list(mu = formals(orig.pFam)[[2]], 
                  sigma = formals(orig.pFam)[[3]], nu = formals(orig.pFam)[[4]], tau = formals(orig.pFam)[[5]], 
                  xi0 = 0.1, xi1 = 0.1, n = 101, from = 0.001, to = 0.999, lower.tail = TRUE, log = FALSE))
        }
    }
    plotfun
}
meanbeInf <- function(obj) {
    if (obj$family[1] != "BEINF") 
        stop("the object do not have a BEINF distribution")
    meanofY <- (fitted(obj, "xi1") + fitted(obj, "mu"))/(1 + fitted(obj, "xi0") + fitted(obj, "xi1"))
    meanofY
}
centiles.Inf0to1 <- function(obj, xvar = NULL, cent = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6), legend = TRUE, 
    ylab = "y", xlab = "x", main = NULL, main.gsub = "@", xleg = min(xvar), yleg = max(obj$y), xlim = range(xvar), 
    ylim = range(obj$y), save = FALSE, plot = TRUE, points = TRUE, pch = 15, cex = 0.5, col = gray(0.7), 
    col.centiles = 1:length(cent) + 2, lty.centiles = 1, lwd.centiles = 1, ...) {
    if (!is(obj, "gamlssinf0to1")) 
        stop(paste("This is not an gamlss object", "\n", ""))
    if (is.null(xvar)) 
        stop(paste("The xvar argument is not specified", "\n", ""))
    fname <- obj$original.family[1]
    invcdf <- Inf0to1.q(fname, obj$typeInf)
    Title <- paste("Centile curves using", fname, sep = " ")
    main <- if (is.null(main)) 
        paste("Centile curves using", fname, sep = " ")
    else gsub(main.gsub, Title, main)
    oxvar <- xvar[order(xvar)]
    oyvar <- obj$y[order(xvar)]
    if (is.matrix(obj$y)) {
        oyvar <- obj$y[, 1][order(xvar)]
        ylim <- range(obj$y[, 1])
        yleg = max(obj$y[, 1])
    }
    if (plot) {
        lty.centiles <- rep(lty.centiles, length(cent))
        lwd.centiles <- rep(lwd.centiles, length(cent))
        col.centiles <- rep(col.centiles, length(cent))
        if (points == TRUE) {
            plot(oxvar, oyvar, type = "p", col = col, pch = pch, cex = cex, xlab = xlab, ylab = ylab, xlim = xlim, 
                ylim, ...)
        }
        else {
            plot(oxvar, oyvar, type = "n", col = col, pch = pch, xlab = xlab, ylab = ylab, xlim = xlim, ylim, 
                ...)
        }
        title(main)
    }
    col <- 3
    lpar <- eval(gamlss.family(obj$original.family))$nopar
    ii <- 0
    per <- rep(0, length(cent))
    for (var in cent) {
        if (obj$typeInf == "Zero") {
            if (lpar == 2) {
                newcall <- call("invcdf", var/100, mu = fitted(obj, "mu")[order(xvar)], sigma = fitted(obj, 
                  "sigma")[order(xvar)], xi0 = fitted(obj, "xi0")[order(xvar)])
            }
            else if (lpar == 3) {
                newcall <- call("invcdf", var/100, mu = fitted(obj, "mu")[order(xvar)], sigma = fitted(obj, 
                  "sigma")[order(xvar)], nu = fitted(obj, "nu")[order(xvar)], xi0 = fitted(obj, "xi0")[order(xvar)])
            }
            else if (lpar == 4) {
                newcall <- call("invcdf", var/100, mu = fitted(obj, "mu")[order(xvar)], sigma = fitted(obj, 
                  "sigma")[order(xvar)], nu = fitted(obj, "nu")[order(xvar)], tau = fitted(obj, "tau")[order(xvar)], 
                  xi0 = fitted(obj, "xi0")[order(xvar)])
            }
            else {
                stop("not the right parameters size")
            }
        }
        else if (obj$typeInf == "One") {
            if (lpar == 2) {
                newcall <- call("invcdf", var/100, mu = fitted(obj, "mu")[order(xvar)], sigma = fitted(obj, 
                  "sigma")[order(xvar)], xi1 = fitted(obj, "xi1")[order(xvar)])
            }
            else if (lpar == 3) {
                newcall <- call("invcdf", var/100, mu = fitted(obj, "mu")[order(xvar)], sigma = fitted(obj, 
                  "sigma")[order(xvar)], nu = fitted(obj, "nu")[order(xvar)], xi1 = fitted(obj, "xi1")[order(xvar)])
            }
            else if (lpar == 4) {
                newcall <- call("invcdf", var/100, mu = fitted(obj, "mu")[order(xvar)], sigma = fitted(obj, 
                  "sigma")[order(xvar)], nu = fitted(obj, "nu")[order(xvar)], tau = fitted(obj, "tau")[order(xvar)], 
                  xi1 = fitted(obj, "xi1")[order(xvar)])
            }
            else {
                stop("not the right parameters size")
            }
        }
        else {
            if (lpar == 2) {
                newcall <- call("invcdf", var/100, mu = fitted(obj, "mu")[order(xvar)], sigma = fitted(obj, 
                  "sigma")[order(xvar)], xi0 = fitted(obj, "xi0")[order(xvar)], xi1 = fitted(obj, "xi1")[order(xvar)])
            }
            else if (lpar == 3) {
                newcall <- call("invcdf", var/100, mu = fitted(obj, "mu")[order(xvar)], sigma = fitted(obj, 
                  "sigma")[order(xvar)], nu = fitted(obj, "nu")[order(xvar)], xi0 = fitted(obj, "xi0")[order(xvar)], 
                  xi1 = fitted(obj, "xi1")[order(xvar)])
            }
            else if (lpar == 4) {
                newcall <- call("invcdf", var/100, mu = fitted(obj, "mu")[order(xvar)], sigma = fitted(obj, 
                  "sigma")[order(xvar)], nu = fitted(obj, "nu")[order(xvar)], tau = fitted(obj, "tau")[order(xvar)], 
                  xi0 = fitted(obj, "xi0")[order(xvar)], xi1 = fitted(obj, "xi1")[order(xvar)])
            }
            else {
                stop("not the right parameters size")
            }
        }
        ii <- ii + 1
        ll <- eval(newcall)
        if (plot) {
            lines(oxvar, ll, col = col.centiles[ii], lty = lty.centiles[ii], lwd = lwd.centiles[ii], ...)
        }
        per[ii] <- (1 - sum(oyvar > ll)/length(oyvar)) * 100
        if (!save) 
            cat("% of cases below ", var, "centile is ", per[ii], "\n")
    }
    if (plot) {
        if (legend == TRUE) 
            legend(list(x = xleg, y = yleg), legend = cent, col = col.centiles, lty = lty.centiles, lwd = lwd.centiles, 
                ncol = 1, ...)
    }
    if (save) {
        return(cbind(cent, per))
    }
}
term.plotInf0to1 <- function(object, parameter = c("mu", "sigma", "nu", "tau", "xi0", "xi1"), ...) {
    par <- object$par
    parameter <- match.arg(parameter)
    if (!(parameter %in% par)) 
        stop("the asking parameter is not fitted")
    whetherBIorMN3 <- any(is.na((match(c("xi0", "xi1"), par))))
    invisible(switch(parameter, mu = term.plot(object$dist, "mu", ...), sigma = term.plot(object$dist, "sigma", 
        ...), nu = term.plot(object$dist, "nu", ...), tau = term.plot(object$dist, "tau", ...), xi0 = term.plot(object$multin, 
        "mu", ...), xi1 = ifelse(whetherBIorMN3, term.plot(object$multin, "mu", ...), term.plot(object$multin, 
        "sigma", ...))))
}
