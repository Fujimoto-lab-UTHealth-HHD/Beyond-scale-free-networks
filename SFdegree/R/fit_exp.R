#' @title Fits the Discretized Exponential Model
#'
#' @description Standard numerical optimization routine to obtain the maximum likelihood estimate for the discrete exponential model.
#'
#' @param k a numeric vector of each observed degree
#' @param kmin a numeric scalar to indicate where the right-tail cutoff is
#' @param inits a named vector indicating the initial values for the numerical optimization
#' @param opt.controls control parameter for the optimization routine
#'
#' @return a list object containing the following elements:
#' `par` is the parameter estimate in their original scale (e.g., labmda > 0). `log.par` is the estimated parameter from the `optim(...)` scaled to be ranged from (-Inf, Inf).
#' `LL` is the overall log-likelihood value evaluated at the estimated MLE, `optim` is the output from `optim(...)`.
#' `data` is a numeric vector for the data that was analyzed (i.e., where k >= kmin), `kmin` is the user specified `kmin`, and `data.orig` is the original vector of input data
#'
#' @export
#'
#' @note For the discrete exponential model, L-BFGS-B is used as the numerical optimization method
#'
#' @examples
#' # a data vector
#' set.seed(1)
#' my_k = rpois(n = 100, lambda = 2)
#' my_kmin = 2
#'
#' # fit
#' my_opt = fit_exp(
#'   k = my_k, kmin = my_kmin,
#'   inits = c(log.lambda = 1)
#' )

fit_exp = function(k, kmin, inits = c(log.lambda = 1), opt.controls = list(fnscale = -sum(k >= kmin), factr = 1e-5, maxit = 1e5)){
  # subset
  k.orig = k
  k = k[k >= kmin]

  # run optim
  opt = optim(
    par = inits,
    fn = function(par){
      LLi = exp_LLi(k = k, kmin = kmin, params = list(lambda = exp(par['log.lambda'])),
                    log = TRUE)
      LLi[!is.finite(LLi)] = min(LLi) + log(1e-10)
      LLi[!is.finite(LLi)] = log(1e-30)
      sum(LLi)
    },
    method = "L-BFGS-B",
    control = opt.controls,
    hessian = T
  )

  # output
  out = list(
    par = c(lambda = unname(exp(opt$par))),
    log.par = opt$par, LL = opt$value, optim = opt,
    data = k,
    kmin = kmin, data.orig = k.orig
  )
  return(out)
}
