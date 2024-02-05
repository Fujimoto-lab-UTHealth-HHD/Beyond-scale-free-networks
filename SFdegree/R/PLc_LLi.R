#' @title Likelihood function of the discrete power law model with exponential cutoff
#'
#' @description Computes the likelihood of the discrete exponential model with a given `kmin`
#'
#' @param k a numeric vector of each observed degree
#' @param kmin right-tail cutoff that should be estimated via the `igraph` package with the `fit_power_law(...)` function
#' @param params a list object holding the `alpha` and `lambda` parameter for the discrete power law with exponential cutoff model
#' @param log a logical scalar. Defaults to TRUE which tells the function to return the log-likelihood instead of the likelihood.
#'
#' @return a numeric vector with length equal to the number of data point that is greater than or equal to `kmin`.
#' Each element is the likelihood value of the data evaluated at the given (`alpha`, `lambda`) and right-tail cutoff `kmin`.
#'
#' @export
#'
#' @importFrom VGAM lerch
#'
#' @note Note that the `kmin` is inclusive, as in, for a vector `k = c(2,1,3,5,4)`, if `kmin` is 3, then the data analyzed would be `c(3, 5, 4)`.
#'
#' @examples
#' # a data vector c(1,2,3,4,5) with kmin = 1, thus all data are analyzed
#' PLc_LLi(k = c(1,2,3,4,5), kmin = 1, params = list(alpha = 1.5, lambda = 0.5))
#'
#' # a data vector c(2,1,3,5,4) with kmin = 3, thus only c(3, 5, 4) are analyzed
#' # note that shifting the kmin changes the location of the distribution,
#' # hence k=3 here will not yield the same likelihood value compare to the previous example
#' PLc_LLi(k = c(2,1,3,5,4), kmin = 3, params = list(alpha = 1.5, lambda = 0.5))

PLc_LLi = function(k, kmin, params, log = T){
  # input vector, restricted to only >= kmin
  k = k[k >= kmin]

  # parameters ~ alpha > 1; lambda > 0
  alpha = params$alpha
  lambda = params$lambda

  # log-likelihood
  if(log){
    log_Li = lambda*(kmin-k) - (alpha)*log(k) - log(VGAM::lerch(x=exp(-lambda), s=alpha, v=kmin))
    return(log_Li)
  }

  # likelihood
  if(!log){
    Li = exp(lambda*(kmin-k)) / exp((alpha)*log(k)) / VGAM::lerch(x=exp(-lambda), s=alpha, v=kmin)
    return(Li)
  }
}
