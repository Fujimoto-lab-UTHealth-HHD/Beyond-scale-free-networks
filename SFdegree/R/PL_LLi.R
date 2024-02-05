#' @title Likelihood function of the power law model
#'
#' @description Computes the likelihood of the discrete power law model with a given `kmin`
#' @param k a numeric vector of each observed degree
#' @param kmin right-tail cutoff that should be estimated via the `igraph` package with the `fit_power_law(...)` function
#' @param params a list object holding the alpha parameter for the discrete power law distribution
#' @param log a logical scalar. Defaults to TRUE which tells the function to return the log-likelihood instead of the likelihood.
#'
#' @return a numeric vector with length equal to the number of data point that is greater than or equal to `kmin`.
#' Each element is the likelihood value of the data evaluated at the given `alpha` and right-tail cutoff `kmin`.
#'
#' @export
#'
#' @importFrom VGAM zeta
#' @importFrom igraph fit_power_law
#'
#' @note Note that the `kmin` is inclusive, as in, for a vector `k = c(2,1,3,5,4)`, if `kmin` is 3, then the data analyzed would be `c(3, 5, 4)`.
#'
#' @examples
#' # a data vector c(1,2,3,4,5) with kmin = 1, thus all data are analyzed
#' PL_LLi(k = c(1,2,3,4,5), kmin = 1, params = list(alpha = 2.5))
#'
#' # a data vector c(2,1,3,5,4) with kmin = 3, thus only c(3, 5, 4) are analyzed
#' # note that shifting the kmin changes the location of the distribution,
#' # hence k=3 here will not yield the same likelihood value compare to the previous example
#' PL_LLi(k = c(2,1,3,5,4), kmin = 3, params = list(alpha = 2.5))

PL_LLi = function(k, kmin, params, log = TRUE){
  # input vector, restricted to only >= kmin
  k = k[k >= kmin]

  # parameters ~ alpha > 1
  alpha = params$alpha

  # log-likelihood
  if(log){
    log_Li = -(log(VGAM::zeta(x = alpha, shift = kmin, deriv = 0)) + alpha*log(k))
    return(log_Li)
  }

  # likelihood
  if(!log){
    Li = 1/VGAM::zeta(x = alpha, shift = kmin, deriv = 0) * k^(-alpha)
    return(Li)
  }
}
