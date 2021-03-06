#' Extract draws of the linear predictor and draw from the predictive
#' distribution of the projected submodel
#'
#' \code{proj_linpred} extracts draws of the linear predictor and
#' \code{proj_predict} draws from the predictive distribution of the projected
#' submodel or submodels. If the projection has not been performed, the
#' functions also perform the projection.
#'
#' @name proj-pred
#'
#' @param object Either an object returned by \link[=varsel]{varsel}, \link[=cv_varsel]{cv_varsel}
#' or \link[=init_refmodel]{init_refmodel}, or alternatively any object that can be converted to a reference model.
#' @param xnew The predictor values used in the prediction. If \code{vind} is
#' specified, then \code{xnew} should either be a dataframe containing column names
#' that correspond to \code{vind} or a matrix with the number and order of columns
#' corresponding to \code{vind}. If \code{vind} is unspecified, then \code{xnew} must
#' either be a dataframe containing all the column names as in the original data or a matrix
#' with the same columns at the same positions as in the original data.
#' @param ynew New (test) target variables. If given, then the log predictive density
#' for the new observations is computed.
#' @param offsetnew Offsets for the new observations. By default a vector of
#' zeros.
#' @param weightsnew Weights for the new observations. For binomial model,
#' corresponds to the number trials per observation. For \code{proj_linpred},
#' this argument matters only if \code{ynew} is specified. By default a vector
#' of ones.
#' @param transform Should the linear predictor be transformed using the
#' inverse-link function? Default is \code{FALSE}. For \code{proj_linpred} only.
#' @param integrated If \code{TRUE}, the output is averaged over the
#' parameters. Default is \code{FALSE}. For \code{proj_linpred} only.
#' @param nv Number of variables in the submodel (the variable combination is
#' taken from the variable selection information). If a vector with several values,
#' then results for all specified model sizes are returned. Ignored if \code{vind} is specified.
#' By default use the automatically suggested model size.
#' @param draws Number of draws to return from the predictive distribution of
#' the projection. The default is 1000.
#' For \code{proj_predict} only.
#' @param seed_samp An optional seed to use for drawing from the projection.
#' For \code{proj_predict} only.
#' @param ... Additional argument passed to \link{project} if \code{object}
#' is an object returned by \link{varsel} or \link{cv_varsel}.
#'
#' @return If the prediction is done for one submodel only (\code{nv} has length one
#' or \code{vind} is specified) and ynew is unspecified, a matrix or vector of
#' predictions (depending on the value of \code{integrated}). If \code{ynew} is specified,
#' returns a list with elements pred (predictions) and lpd (log predictive densities).
#' If the predictions are done for several submodel sizes, returns a list with one element
#' for each submodel.
#' 
#' @examples
#' \donttest{
#' if (requireNamespace('rstanarm', quietly=TRUE)) {
#'   ### Usage with stanreg objects
#'   n <- 30
#'   d <- 5
#'   x <- matrix(rnorm(n*d), nrow=n)
#'   y <- x[,1] + 0.5*rnorm(n)
#'   data <- data.frame(x,y)
#'   
#'   fit <- rstanarm::stan_glm(y~., gaussian(), data=data, chains=2, iter=500)
#'   vs <- varsel(fit)
#'   
#'   # compute predictions with 4 variables at the training points
#'   pred <- proj_linpred(vs, xnew=x, nv = 4)
#'   pred <- proj_predict(vs, xnew=x, nv = 4)
#' }
#' }
#'
NULL

# The 'helper' for proj_linpred and proj_predict, ie. does all the
# functionality that is common to them. It essentially checks all the arguments
# and sets them to their respective defaults and then loops over the
# projections. For each projection, it evaluates the fun-function, which
# calculates the linear predictor if called from proj_linpred and samples from
# the predictive distribution if called from proj_predict.
proj_helper <- function(object, xnew, offsetnew, weightsnew, nv, seed_samp,
                        fun, ...) {

  if (!inherits(xnew, c('data.frame', 'matrix')))
    stop('xnew must be a data.frame or a matrix. See ?proj-pred.')

  if (is.null(offsetnew)) offsetnew <- rep(0, nrow(xnew))
  if (is.null(weightsnew)) weightsnew <- rep(1, nrow(xnew))

  if ('projection' %in% class(object) ||
      (length(object)>0 && 'projection' %in% class(object[[1]]))) {
    proj <- object
  } else {
    # reference model or varsel object obtained, so run the projection
    proj <- project(object = object, nv = nv, ...)
  }

  if (!.is_proj_list(proj)) {
    proj <- list(proj)
  } else {
    # proj is not a projection object
    if(any(sapply(proj, function(x) !('family_kl' %in% names(x)))))
      stop('list contains objects not created by varsel, cv_varsel or project')
  }

  projected_sizes <- sapply(proj, function(x) NROW(x$beta))
  nv <- list(...)$nv %ORifNULL% projected_sizes

  if (!all(nv %in% projected_sizes))
    stop(paste0('Linear prediction requested for nv = ',
                paste(nv, collapse = ', '),
                ', but projection performed only for nv = ',
                paste(projected_sizes, collapse = ', '), '.'))

  projs <- Filter(function(x) NROW(x$beta) %in% nv, proj)
  names(projs) <- nv

  xnew_df <- is.data.frame(xnew)
  if (xnew_df) {
    terms <- unique(unlist(lapply(projs, function(x) names(x$vind))))
    xnew <- .df_to_model_mat(xnew, terms)
  }

  vind <- list(...)$vind
  if (!is.null(vind) && NCOL(xnew) != length(vind))
    stop(paste('The number of columns in xnew does not match the',
               'number of variable indices (vind).'))

  # set random seed but ensure the old RNG state is restored on exit
  if (exists('.Random.seed')) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed_samp)

  preds <- lapply(projs, function(proj) {
    if (xnew_df) {
      xtemp <- xnew[, min(1, length(proj$vind)):length(proj$vind), drop = F]
    } else if (!is.null(vind)) {
      # columns of xnew are assumed to match to the given variable indices
      xtemp <- xnew
    } else {
      # fetch the right columns from the feature matrix
      if (length(proj$vind) > 0 && max(proj$vind) > ncol(xnew))
        stop(paste('xnew has', ncol(xnew), 'columns, but vind expects',
                   max(proj$vind), 'columns.'))
      xtemp <- xnew[, proj$vind, drop = F]
    }
    mu <- proj$family_kl$mu_fun(xtemp, proj$alpha, proj$beta, offsetnew)

    fun(proj, mu, offsetnew, weightsnew)
  })

  .unlist_proj(preds)
}

#' @rdname proj-pred
#' @export
proj_linpred <- function(object, xnew, ynew = NULL, offsetnew = NULL,
                         weightsnew = NULL, nv = NULL, transform = FALSE,
                         integrated = FALSE, ...) {

  # function to perform to each projected submodel
  fun <- function(proj, mu, offset, weights) {
    pred <- t(mu)
    if (!transform) pred <- proj$family_kl$linkfun(pred)
    if (integrated) {
      # average over the parameters
      pred <- as.vector( proj$weights %*% pred )
    } else if (!is.null(dim(pred)) && dim(pred)[1]==1) {
      # return a vector if pred contains only one row
      pred <- as.vector(pred)
    }

    if (!is.null(ynew)) {
      # compute also the log-density
      target <- .get_standard_y(ynew, weights, proj$family_kl)
      ynew <- target$y
      weights <- target$weights
      lpd <- proj$family_kl$ll_fun(mu, proj$dis, ynew, weights)
      if (integrated && !is.null(dim(lpd))) {
        lpd <- as.vector(apply(lpd, 1, log_weighted_mean_exp, proj$weights))
      } else if (!is.null(dim(lpd))) {
        lpd <- t(lpd)
      }
      list(pred = pred, lpd = lpd)
    } else {
      pred
    }
  }

  # proj_helper lapplies fun to each projection in object
  proj_helper(object = object, xnew = xnew, offsetnew = offsetnew,
              weightsnew = weightsnew, nv = nv, seed_samp = NULL, fun = fun,
              ...)
}

#' @rdname proj-pred
#' @export
proj_predict <- function(object, xnew, offsetnew = NULL, weightsnew = NULL,
                         nv = NULL, draws = NULL, seed_samp = NULL, ...) {

  # function to perform to each projected submodel
  fun <- function(proj, mu, offset, weights) {
    if(is.null(draws)) draws <- 1000
    draw_inds <- sample(x = seq_along(proj$weights), size = draws,
                        replace = TRUE, prob = proj$weights)

    t(sapply(draw_inds, function(i) {
      proj$family_kl$ppd_fun(mu[,i], proj$dis[i], weights)
    }))
  }

  # proj_helper lapplies fun to each projection in object
  proj_helper(object = object, xnew = xnew, offsetnew = offsetnew,
              weightsnew = weightsnew, nv = nv, seed_samp = seed_samp,
              fun = fun, ...)
}

#' Plot or fetch summary statistics related to variable selection
#'
#' \code{varsel_stats} can be used to obtain summary statistics related to
#' variable selection. The same statistics can be plotted with
#' \code{varsel_plot}.
#'
#' @name varsel-statistics
#'
#' @param object The object returned by \link[=varsel]{varsel} or
#' \link[=cv_varsel]{cv_varsel}.
#' @param nv_max Maximum submodel size for which the statistics are calculated.
#' For \code{varsel_plot} it must be at least 1.
#' @param stats One or several strings determining which statistics to calculate. Available
#' statistics are: 
#' \itemize{
#'  \item{elpd:} {(Expected) sum of log predictive densities}
#'  \item{mlpd:} {Mean log predictive density, that is, elpd divided by the number of datapoints.}
#'  \item{mse:} {Mean squared error (gaussian family only)}
#'  \item{rmse:} {Root mean squared error (gaussian family only)}
#'  \item{acc/pctcorr:} {Classification accuracy (binomial family only)}
#'  \item{auc:} {Area under the ROC curve (binomial family only)}
#' }
#' Default is elpd.
#' @param type One or more items from 'mean', 'se', 'lower' and 'upper' indicating which of these to
#' compute (mean, standard error, and lower and upper credible bounds). The credible bounds are determined so
#' that \code{1-alpha} percent of the mass falls between them.
#' @param deltas If \code{TRUE}, the submodel statistics are estimated relative to the baseline model 
#' (see argument \code{baseline}) instead of estimating the actual values of the statistics.
#' Defaults to \code{FALSE}.
#' @param alpha A number indicating the desired coverage of the credible
#' intervals. For example \code{alpha=0.32} corresponds to 68\% probability mass
#' within the intervals, that is, one standard error intervals.
#' @param baseline Either 'ref' or 'best' indicating whether the baseline is the reference model or
#' the best submodel found. Default is 'ref' when the reference model exists, and 'best' otherwise.
#' @param ... Currently ignored.
#' 
#' 
#' @examples
#' \donttest{
#' if (requireNamespace('rstanarm', quietly=TRUE)) {
#'   ### Usage with stanreg objects
#'   n <- 30
#'   d <- 5
#'   x <- matrix(rnorm(n*d), nrow=n)
#'   y <- x[,1] + 0.5*rnorm(n)
#'   data <- data.frame(x,y)
#'   
#'   fit <- rstanarm::stan_glm(y~., gaussian(), data=data, chains=2, iter=500)
#'   vs <- cv_varsel(fit)
#'   varsel_plot(vs)
#'   
#'   # print out some stats
#'   varsel_stats(vs, stats=c('mse'), type = c('mean','se'))
#' }
#' }
#' 
NULL

#' @rdname varsel-statistics
#' @export
varsel_plot <- function(object, nv_max = NULL, stats = 'elpd', deltas = F, alpha = 0.32, baseline=NULL, ...) {

  .validate_vsel_object_stats(object, stats)
  baseline <- .validate_baseline(object$refmode, baseline, deltas)

  # compute all the statistics and fetch only those that were asked
  nfeat_baseline <- .get_nfeat_baseline(object, baseline, stats[1])
  tab <- rbind(.tabulate_stats(object, stats, alpha = alpha, nfeat_baseline=nfeat_baseline),
               .tabulate_stats(object, stats, alpha = alpha))
  stats_table <- subset(tab, tab$delta==deltas)
  stats_ref <- subset(stats_table, stats_table$size==Inf)
  stats_sub <- subset(stats_table, stats_table$size!=Inf)
  stats_bs <- subset(stats_table, stats_table$size == nfeat_baseline)
	
	
	if(NROW(stats_sub) == 0) {
	    stop(paste0(ifelse(length(stats)==1, 'Statistics ', 'Statistic '),
	                paste0(unique(stats), collapse=', '), ' not available.'))
	}
	
	if(is.null(nv_max))
	  nv_max <- max(stats_sub$size)
        else {
	  # don't exceed the maximum submodel size
	  nv_max <- min(nv_max, max(stats_sub$size))
	  if (nv_max < 1)
	    stop('nv_max must be at least 1')
	}
	ylab <- if(deltas) 'Difference to the baseline' else 'Value'
	
	# make sure that breaks on the x-axis are integers
	n_opts <- c(4,5,6)
	n_possible <- Filter(function(x) nv_max %% x == 0, n_opts)
	n_alt <- n_opts[which.min(n_opts - (nv_max %% n_opts))]
	nb <- ifelse(length(n_possible) > 0, min(n_possible), n_alt)
	by <- ceiling(nv_max/min(nv_max, nb))
	breaks <- seq(0, by*min(nv_max, nb), by)
	minor_breaks <- if(by%%2 == 0)
	    seq(by/2, by*min(nv_max, nb), by)
	else
	  NULL

	# plot submodel results
	pp <- ggplot(data = subset(stats_sub, stats_sub$size <= nv_max), mapping = aes_string(x = 'size')) +
		geom_linerange(aes_string(ymin = 'lq', ymax = 'uq', alpha=0.1)) +
    geom_line(aes_string(y = 'value')) +
    geom_point(aes_string(y = 'value'))
	
	if (!all(is.na(stats_ref$se)))
		# add reference model results if they exist
		pp <- pp + geom_hline(aes_string(yintercept = 'value'), data = stats_ref,
													color = 'darkred', linetype=2)
	if (baseline != 'ref')
	  # add the baseline result (if different from the reference model)
	  pp <- pp + geom_hline(aes_string(yintercept = 'value'), data = stats_bs,
	                        color = 'black', linetype=3)
	pp <- pp + 
		scale_x_continuous(breaks = breaks, minor_breaks = minor_breaks,
	                       limits = c(min(breaks), max(breaks))) +
	  labs(x = 'Number of variables in the submodel', y = ylab) +
	  theme(legend.position = 'none') +
	  facet_grid(statistic ~ ., scales = 'free_y') 
	pp
}

#' @rdname varsel-statistics
#' @export
varsel_stats <- function(object, nv_max = NULL, stats = 'elpd', type = c('mean','se'), 
                         deltas = F, alpha=0.32, baseline=NULL, ...) {

  .validate_vsel_object_stats(object, stats)
  baseline <- .validate_baseline(object$refmode, baseline, deltas)

  # fetch statistics
  if (deltas) {
    nfeat_baseline <- .get_nfeat_baseline(object, baseline, stats[1])
    tab <- .tabulate_stats(object, stats, alpha=alpha, nfeat_baseline=nfeat_baseline)
  } else {
    tab <- .tabulate_stats(object, stats, alpha=alpha)
  }
  stats_table <- subset(tab, tab$size != Inf)  
  
  
  # these are the corresponding names for mean, se, upper and lower in the stats_table, and their suffices
  # in the table to be returned 
  qty <- unname(sapply(type, function(t) switch(t, mean='value', upper='uq', lower='lq', se='se')))
  suffix <- unname(sapply(type, function(t) switch(t, mean='', upper='.upper', lower='.lower', se='.se')))

  # loop through all the required statistics 
  arr <- data.frame(size = unique(stats_table$size), vind = c(NA, object$vind))
  for (i in seq_along(stats)) {
    temp <- subset(stats_table, stats_table$statistic == stats[i], qty)
    newnames <- sapply(suffix, function(s) paste0(stats[i],s))
    colnames(temp) <- newnames
    arr <- cbind(arr, temp)
  }

  if(is.null(nv_max)) 
    nv_max <- max(stats_table$size)

  if('pctch' %in% names(object))
    arr$pctch <- c(NA, diag(object$pctch[,-1]))

  subset(arr, arr$size <= nv_max)
}

#' Print methods for vsel/cvsel objects
#'
#' The \code{print} methods for vsel/cvsel objects created by \code{\link{varsel}}
#' or \code{\link{cv_varsel}}) rely on \code{\link{varsel_stats}} to display
#' a summary of the results of the projection predictive variable selection.
#'
#' @name print-vsel
#'
#' @param x An object of class vsel/cvsel.
#' @param digits Number of decimal places to be reported (2 by default).
#' @param ... Further arguments passed to \code{\link{varsel_stats}}.
#'
#' @return Returns invisibly the data frame produced by \code{\link{varsel_stats}}.
#'
#' @export
#' @method print vsel
print.vsel <- function(x, digits=2, ...) {
  stats <- round(varsel_stats(x, ...), digits)
  print(stats[, -match("vind", colnames(stats))])
  invisible(stats)
}

#' @rdname print-vsel
#' @export
#' @method print cvsel
print.cvsel <- function(x, digits=2, ...) {
  stats <- round(varsel_stats(x, ...), digits)
  print(stats[, -match("vind", colnames(stats))])
  invisible(stats)
}


#' Suggest model size 
#'
#' This function can be used for suggesting an appropriate model size
#' based on a certain default rule. Notice that the decision rules are heuristic
#' and should be interpreted as guidelines. It is recommended that the user
#' studies the results via \code{varsel_plot} and/or \code{varsel_stats}
#' and makes the final decision based on what is most appropriate for the given
#' problem.
#'
#' @param object The object returned by \link[=varsel]{varsel} or
#' \link[=cv_varsel]{cv_varsel}.
#' @param stat Statistic used for the decision. Default is 'elpd'. See \code{varsel_stats} for
#' other possible choices. 
#' @param alpha A number indicating the desired coverage of the credible
#' intervals based on which the decision is made. E.g. \code{alpha=0.32} corresponds to
#' 68\% probability mass within the intervals (one standard error intervals). 
#' See details for more information.
#' @param pct Number indicating the relative proportion between baseline model and null model
#' utilities one is willing to sacrifice. See details for more information.
#' @param type Either 'upper' (default) or 'lower' determining whether the decisions are
#' based on the upper or lower credible bounds. See details for more information.
#' @param baseline Either 'ref' or 'best' indicating whether the baseline is the reference model or
#' the best submodel found. Default is 'ref' when the reference model exists, and 'best' otherwise.
#' @param warnings Whether to give warnings if automatic suggestion fails, mainly for internal use.
#' Default is TRUE, and usually no reason to set to FALSE.
#' @param ... Currently ignored.
#' 
#' @details The suggested model size is the smallest model for which
#' either the lower or upper (depending on argument \code{type}) credible bound 
#' of the submodel utility \eqn{u_k} with significance level \code{alpha} falls above
#'   \deqn{u_base - pct*(u_base - u_0)} 
#' Here \eqn{u_base} denotes the utility for the baseline model and \eqn{u_0} the null model utility.
#' The baseline is either the reference model or the best submodel found (see argument \code{baseline}).
#' The lower and upper bounds are defined to contain the submodel utility with 
#' probability 1-alpha (each tail has mass alpha/2).
#' 
#' By default \code{ratio=0}, \code{alpha=0.32} and \code{type='upper'} which means that we select the smallest
#' model for which the upper tail exceeds the baseline model level, that is, which is better than the baseline 
#' model with probability 0.16 (and consequently, worse with probability 0.84). In other words,
#' the estimated difference between the baseline model and submodel utilities is at most one standard error
#' away from zero, so the two utilities are considered to be close.
#' 
#' NOTE: Loss statistics like RMSE and MSE are converted to utilities by multiplying them by -1, so call
#' such as \code{suggest_size(object, stat='rmse', type='upper')} should be interpreted as finding
#' the smallest model whose upper credible bound of the \emph{negative} RMSE exceeds the cutoff level
#' (or equivalently has the lower credible bound of RMSE below the cutoff level). This is done to make
#' the interpretation of the argument \code{type} the same regardless of argument \code{stat}.
#' 
#' @examples
#' \donttest{
#' if (requireNamespace('rstanarm', quietly=TRUE)) {
#'   ### Usage with stanreg objects
#'   n <- 30
#'   d <- 5
#'   x <- matrix(rnorm(n*d), nrow=n)
#'   y <- x[,1] + 0.5*rnorm(n)
#'   data <- data.frame(x,y)
#'   
#'   fit <- rstanarm::stan_glm(y~., gaussian(), data=data, chains=2, iter=500)
#'   vs <- cv_varsel(fit)
#'   suggest_size(vs)
#' }
#' }
#' 

#' @export
suggest_size <- function(object, stat = 'elpd', alpha = 0.32, pct = 0.0, type='upper', 
                         baseline=NULL, warnings=TRUE, ...) {

  .validate_vsel_object_stats(object, stat)
  if (length(stat) > 1)
    stop('Only one statistic can be specified to suggest_size')
  
  if (.is_util(stat)) {
    sgn <- 1
  } else {
    sgn <- -1
    if (type == 'upper')
      type <- 'lower'
    else
      type <- 'upper'
  }
  bound <- paste0(stat,'.',type)
  stats <- varsel_stats(object, stats=stat, alpha=alpha, type=c('mean','upper','lower'),
                        baseline=baseline, deltas = T)
  util_null <- sgn*unlist(unname(subset(stats, stats$size==0, stat)))
  util_cutoff <- pct*util_null
  res <- subset(stats, sgn*stats[,bound] >= util_cutoff, 'size')

  if(nrow(res) == 0) {
    # no submodel satisfying the criterion found
    if (object$nv_max == object$nv_all)
      ssize <- object$nv_max
    else {
      ssize <- NA
      if (warnings)
        warning(paste('Could not suggest model size. Investigate varsel_plot to identify',
                      'if the search was terminated too early. If this is the case,',
                      'run variable selection with larger value for nv_max.'))
    }
    
  } else {
    ssize <- min(res)
  }
   
  ssize
}




#' @method as.matrix projection
#' @export
as.matrix.projection <- function(x, ...) {
  if (x$p_type) {
    warning(paste0('Note, that projection was performed using',
                   'clustering and the clusters might have different weights.'))
  }
  res <- t(x$beta)
  if (ncol(res) > 0) colnames(res) <- names(x$vind)
  if (x$intercept) res <- cbind('(Intercept)' = x$alpha, res)
  if (x$family_kl$family == 'gaussian') res <- cbind(res, sigma = x$dis)
  res
}





#' Create cross-validation indices
#'
#' Divide indices from 1 to \code{n} into subsets for \code{k}-fold cross validation.
#' These functions are potentially useful when creating the \code{cvfits} and \code{cvfun}
#' arguments for \link[=init_refmodel]{init_refmodel}. The returned value is different for
#' these two methods, see below for details.
#' 
#' @name cv-indices
#'
#' @param n Number of data points.
#' @param k Number of folds. Must be at least 2 and not exceed \code{n}.
#' @param out Format of the output, either 'foldwise' (default) or 'indices'. See below for details.
#' @param seed Random seed so that the same division could be obtained again if needed.
#'
#' @return \code{cvfolds} returns a vector of length \code{n} such that each element is an integer
#' between 1 and \code{k} denoting which fold the corresponding data point belongs to.
#' The returned value of \code{cvind} depends on the \code{out}-argument. If \code{out}='foldwise',
#' the returned value is a list with \code{k} elements, 
#' each having fields \code{tr} and \code{ts} which give the training and test indices, respectively,
#' for the corresponding fold. If \code{out}='indices', the returned value is a list with fields \code{tr}
#' and \code{ts}
#' each of which is a list with \code{k} elements giving the training and test indices for each fold.
#' @examples
#' \donttest{
#' ### compute sample means within each fold
#' n <- 100
#' y <- rnorm(n)
#' cv <- cvind(n, k=5)
#' cvmeans <- lapply(cv, function(fold) mean(y[fold$tr]))
#' }
#' 
NULL

#' @rdname cv-indices
#' @export
cvfolds <- function(n, k, seed=NULL) {

  .validate_num_folds(k, n)
  
  # set random seed but ensure the old RNG state is restored on exit
  if (exists('.Random.seed')) {
    rng_state_old <- .Random.seed
    on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
  }
  set.seed(seed)
  
  # create and shuffle the indices
  folds <- rep_len(1:k, length.out = n)
  folds <- sample(folds, n, replace=FALSE)
  
  return(folds)
  
}

#' @rdname cv-indices
#' @export
cvind <- function(n, k, out=c('foldwise', 'indices'), seed=NULL) {

	.validate_num_folds(k, n)
	out <- match.arg(out)

	# set random seed but ensure the old RNG state is restored on exit
	if (exists('.Random.seed')) {
	  rng_state_old <- .Random.seed
	  on.exit(assign(".Random.seed", rng_state_old, envir = .GlobalEnv))
	}
	set.seed(seed)
	
	# shuffle the indices
	ind <- sample(1:n, n, replace=FALSE)
	
	if (out == 'foldwise') {
		cv <- lapply(1:k, function(i) {
			ts <- sort(ind[seq(i,n,k)])  # test set
			tr <- setdiff(1:n, ts) # training set
			list(tr=tr,ts=ts)
		})
	}	else if (out == 'indices') {
		cv <- list()
		cv$tr <- list()
		cv$ts <- list()
		for (i in 1:k) {
			ts <- sort(ind[seq(i,n,k)]) # test set
			tr <- setdiff(1:n, ts) # training set
			cv$tr[[i]] <- tr
			cv$ts[[i]] <- ts
		}
	}
	
	return(cv)
}
