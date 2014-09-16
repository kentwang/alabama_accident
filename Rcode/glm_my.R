#-------------------------------------------------------------------------------
#- Modified modeling function i.e. mf$drop.unused.levels <- FALSE
#
#- TODO: - modify glm.nb?
#-------------------------------------------------------------------------------

glm.my <- function(formula, family = gaussian, data, weights,
                subset, na.action, start = NULL,
                etastart, mustart, offset,
                control = list(...),
                model = TRUE, method = "glm.fit",
                x = FALSE, y = TRUE,
                contrasts = NULL, ...)
{
  call <- match.call()
  ## family
  if(is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  ## extract x, y, etc from the model formula and frame
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- FALSE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if(identical(method, "model.frame")) return(mf)
  
  if (!is.character(method) && !is.function(method))
    stop("invalid 'method' argument")
  ## for back-compatibility in return result
  if (identical(method, "glm.fit"))
    control <- do.call("glm.control", control)
  
  mt <- attr(mf, "terms") # allow model.frame to have updated it
  
  Y <- model.response(mf, "any") # e.g. factors are allowed
  ## avoid problems with 1D arrays, but keep names
  if(length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) names(Y) <- nm
  }
  ## null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(,NROW(Y), 0L)
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  weights <- as.vector(model.weights(mf))
  if(!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  ## check weights and offset
  if( !is.null(weights) && any(weights < 0) )
    stop("negative weights not allowed")
  
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
    if(length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
  }
  ## these allow starting values to be expressed in terms of other vars.
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  
  ## We want to set the name on this call and the one below for the
  ## sake of messages from the fitter function
  fit <- eval(call(if(is.function(method)) "method" else method,
                   x = X, y = Y, weights = weights, start = start,
                   etastart = etastart, mustart = mustart,
                   offset = offset, family = family, control = control,
                   intercept = attr(mt, "intercept") > 0L))
  
  ## This calculated the null deviance from the intercept-only model
  ## if there is one, otherwise from the offset-only model.
  ## We need to recalculate by a proper fit if there is intercept and
  ## offset.
  ##
  ## The glm.fit calculation could be wrong if the link depends on the
  ## observations, so we allow the null deviance to be forced to be
  ## re-calculated by setting an offset (provided there is an intercept).
  ## Prior to 2.4.0 this was only done for non-zero offsets.
  if(length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <-
      eval(call(if(is.function(method)) "method" else method,
                x = X[, "(Intercept)", drop=FALSE], y = Y,
                weights = weights, offset = offset, family = family,
                control = control, intercept = TRUE))
    ## That fit might not have converged ....
    if(!fit2$converged)
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    fit$null.deviance <- fit2$deviance
  }
  if(model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if(x) fit$x <- X
  if(!y) fit$y <- NULL
  fit <- c(fit, list(call = call, formula = formula,
                     terms = mt, data = data,
                     offset = offset, control = control, method = method,
                     contrasts = attr(X, "contrasts"),
                     xlevels = .getXlevels(mt, mf)))
  class(fit) <- c(fit$class, c("glm", "lm"))
  fit
}

glm.nb.my <- function (formula, data, weights, subset, na.action, start = NULL, 
                       etastart, mustart, control = glm.control(...), method = "glm.fit", 
                       model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ..., 
                       init.theta, link = log) 
{
  loglik <- function(n, th, mu, y, w) sum(w * (lgamma(th + 
                                                        y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y * 
                                                 log(mu + (y == 0)) - (th + y) * log(th + mu)))
  link <- substitute(link)
  fam0 <- if (missing(init.theta)) 
    do.call("poisson", list(link = link))
  else do.call("negative.binomial", list(theta = init.theta, 
                                         link = link))
  dots <- list(...)
  mf <- Call <- match.call()
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- FALSE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  Terms <- attr(mf, "terms")
  if (method == "model.frame") 
    return(mf)
  Y <- model.response(mf, "numeric")
  X <- if (!is.empty.model(Terms)) 
    model.matrix(Terms, mf, contrasts)
  else matrix(, NROW(Y), 0)
  w <- model.weights(mf)
  if (!length(w)) 
    w <- rep(1, nrow(mf))
  else if (any(w < 0)) 
    stop("negative weights not allowed")
  offset <- model.offset(mf)
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  n <- length(Y)
  if (!missing(method)) {
    if (!exists(method, mode = "function")) 
      stop(gettextf("unimplemented method: %s", sQuote(method)), 
           domain = NA)
    glm.fitter <- get(method)
  }
  else {
    method <- "glm.fit"
    glm.fitter <- stats::glm.fit
  }
  if (control$trace > 1) 
    message("Initial fit:")
  fit <- glm.fitter(x = X, y = Y, w = w, start = start, etastart = etastart, 
                    mustart = mustart, offset = offset, family = fam0, control = list(maxit = control$maxit, 
                                                                                      epsilon = control$epsilon, trace = control$trace > 
                                                                                        1), intercept = attr(Terms, "intercept") > 0)
  class(fit) <- c("glm", "lm")
  mu <- fit$fitted.values
  th <- as.vector(theta.ml(Y, mu, sum(w), w, limit = control$maxit, 
                           trace = control$trace > 2))
  if (control$trace > 1) 
    message(gettextf("Initial value for 'theta': %f", signif(th)), 
            domain = NA)
  fam <- do.call("negative.binomial", list(theta = th, link = link))
  iter <- 0
  d1 <- sqrt(2 * max(1, fit$df.residual))
  d2 <- del <- 1
  g <- fam$linkfun
  Lm <- loglik(n, th, mu, Y, w)
  Lm0 <- Lm + 2 * d1
  while ((iter <- iter + 1) <= control$maxit && (abs(Lm0 - 
                                                       Lm)/d1 + abs(del)/d2) > control$epsilon) {
    eta <- g(mu)
    fit <- glm.fitter(x = X, y = Y, w = w, etastart = eta, 
                      offset = offset, family = fam, control = list(maxit = control$maxit, 
                                                                    epsilon = control$epsilon, trace = control$trace > 
                                                                      1), intercept = attr(Terms, "intercept") > 
                        0)
    t0 <- th
    th <- theta.ml(Y, mu, sum(w), w, limit = control$maxit, 
                   trace = control$trace > 2)
    fam <- do.call("negative.binomial", list(theta = th, 
                                             link = link))
    mu <- fit$fitted.values
    del <- t0 - th
    Lm0 <- Lm
    Lm <- loglik(n, th, mu, Y, w)
    if (control$trace) {
      Ls <- loglik(n, th, Y, Y, w)
      Dev <- 2 * (Ls - Lm)
      message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f", 
                      iter, signif(th), signif(Dev)), domain = NA)
    }
  }
  if (!is.null(attr(th, "warn"))) 
    fit$th.warn <- attr(th, "warn")
  if (iter > control$maxit) {
    warning("alternation limit reached")
    fit$th.warn <- gettext("alternation limit reached")
  }
  if (length(offset) && attr(Terms, "intercept")) {
    null.deviance <- if (length(Terms)) 
      glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w, 
                 offset = offset, family = fam, control = list(maxit = control$maxit, 
                                                               epsilon = control$epsilon, trace = control$trace > 
                                                                 1), intercept = TRUE)$deviance
    else fit$deviance
    fit$null.deviance <- null.deviance
  }
  class(fit) <- c("negbin", "glm", "lm")
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  Call$init.theta <- signif(as.vector(th), 10)
  Call$link <- link
  fit$call <- Call
  if (model) 
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x) 
    fit$x <- X
  if (!y) 
    fit$y <- NULL
  fit$theta <- as.vector(th)
  fit$SE.theta <- attr(th, "SE")
  fit$twologlik <- as.vector(2 * Lm)
  fit$aic <- -fit$twologlik + 2 * fit$rank + 2
  fit$contrasts <- attr(X, "contrasts")
  fit$xlevels <- .getXlevels(Terms, mf)
  fit$method <- method
  fit$control <- control
  fit$offset <- offset
  fit
}
