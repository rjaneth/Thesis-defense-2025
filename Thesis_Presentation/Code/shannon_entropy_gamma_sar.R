shannon_entropy_gamma_sar <- function(L, mu) {
  L - log(L) + lgamma(L) + (1 - L) * digamma(L) + log(mu)
}