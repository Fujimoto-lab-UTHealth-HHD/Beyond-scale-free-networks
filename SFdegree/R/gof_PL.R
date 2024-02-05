#' @title Semiparametric GoF test for the power law model
#'
#' @description Performs the sempiparametric approach to goodness of fit test described in `https://arxiv.org/pdf/0706.1062.pdf` (pg 17; paragraph 2-3).
#'
#' @param k a numeric vector for the observed data
#' @param kmin a numeric scalar indicating the right-tail cutoff
#' @param alpha a numeric scalar for the estimated `alpha` parameter from the power law model (can be obtained via `igraph::fit_power_law`)
#' @param Nsim number of times the simulations should be ran to obtain theoretical KS value (the more simulation, the more accurate the result).
#' For accuracy up to 0.01 decimal point, use `Nsim = 2500`. See Clauset's paper:
#' `Clauset, A., Shalizi, C. R., & Newman, M. E. (2009). Power-law distributions in empirical data. SIAM review, 51(4), 661-703.`
#' @param ss an integer value for setting simulation seed.
#'
#' @return a matrix object with number of rows equal to `Nsim`, where each row stores the simulated statistic obtained by fitting the power law model to the synthetic data
#'
#' @export
#'
#' @importFrom poweRlaw rpldis
#' @importFrom igraph fit_power_law
#'
#' @examples
#' # generate some dataset
#' set.seed(1)
#' my_N = 500
#' my_k = poweRlaw::rpldis(n = my_N, xmin = 3, alpha = 5)
#'
#' # power law fit
#' my_PLmodel = igraph::fit_power_law(x = my_k)
#'
#' # goodness of fit
#' ss = gof_PL(k = my_k, kmin = my_PLmodel$xmin, my_PLmodel$alpha, Nsim = 1000, ss = 123)
#' mean(ss[, "KS.stat"] > my_PLmodel$KS.stat)

gof_PL = function(k, kmin, alpha, Nsim = 100, ss){
  ## set seed
  set.seed(ss)

  ## some constants
  k.PL = k[k >= kmin]
  k.nonPL = k[k < kmin]
  N = length(k)
  Ntail = length(k.PL)
  PL.prob = Ntail/N

  ## simulate and compute KS for synthetic data
  cat("####  Result Log  ##### \n")
  sim_result = matrix(NA, nrow = Nsim, ncol = 6)
  for(ns in 1:Nsim){
    # sample left tail
    rand.nonPL = NULL
    if(length(k.nonPL) > 0){
      rand.nonPL = sample(k.nonPL, size = N, replace = TRUE)
    }

    # sample right tail
    rand.PL = poweRlaw::rpldis(n = N, xmin = kmin, alpha = alpha)

    # with probability Ntail/N, use right tail random number, otherwise use left tail
    r = runif(n = N)
    rand.k = ifelse(r <= PL.prob, rand.PL, rand.nonPL)

    # compute KS stat
    pl_model = igraph::fit_power_law(x = rand.k)
    sim_result[ns, ] = unlist(pl_model)

    # msg
    if((ns %% 100) == 0) cat(" - Sim ", paste(rep(" ", times = nchar(Nsim) - nchar(ns)), ns, sep = ""), "/", Nsim, " completed \n", sep = "")
  }

  # return
  colnames(sim_result) = names(pl_model)
  return(sim_result)
}
