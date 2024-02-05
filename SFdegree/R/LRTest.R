#' @title Likelihood Ratio test to compare model performance
#'
#' @description Performs the likelihood ratio test as described in `https://www.nature.com/articles/s41467-019-08746-5?ref=https://githubhelp.com`
#'
#' @param LLi1 a numeric vector of the log-likelihood of a fitted model
#' @param LLi2 a numeric vector of the log-likelihood of an other fitted model
#'
#' @return a list object containing the test statistics "R" and p-value (among some other relevant statistics).
#' `Rstat` is the R statistic, use its sign to tell which model performed better. `Rstat.sd` is the SE of the R statistic
#' `pval` is the test p-value, `L1` is the overall log-likelihood of model 1, and `L2` is the overall log-likelihood value for model 2
#'
#' @export
#'
#' @note This function takes LLi1 minus LLi2 to obtain the R statistic.
#' Hence a positive R statistic may indicate model 1 is better than model 2 (if p-value agrees). Simply running a paired t-test on LLi1 and LLi2 usually yield the same result
#'
#' @examples
#' # generate data and some theoretical likelihood
#' set.seed(1)
#' my_data = rpois(n = 100, lambda = 5)
#' my_data.PL = PL_LLi(k = my_data, kmin = 4, params = list(alpha = 2.1))
#' my_data.weib = weib_LLi(k = my_data, kmin = 4, params = list(a = 1.5, b = 1.5))
#'
#' # perform LRT
#' LRTest(LLi1 = my_data.PL, LLi2 = my_data.weib) # R statistic is greater than 0, p-value is < 0.05. So we can conclude that the power law model in this case was a better fit than Weibull

LRTest = function(LLi1, LLi2){
  # some stats
  ntail = length(LLi1)
  Ri = LLi1 - LLi2
  R = sum(Ri)
  s = sd(Ri)*sqrt(ntail)

  # compute p
  p = pnorm(q = -abs(R), mean = 0, sd = s, lower.tail = T) + pnorm(q = abs(R), mean = 0, sd = s, lower.tail = F)

  # return
  out = list(
    Rstat=R, Rstat.sd = s, pval=p,
    L1=sum(LLi1), L2=sum(LLi2)
  )
  return(out)
}
