tsallis_entropy_gamma_sar <- function(mu, L, lambda) {
  (1 - exp((1 - lambda)*log(mu) +
             (lambda - 1)*log(L) +
             lgamma(lambda*(L - 1) + 1) -
             lambda*lgamma(L) -
             (lambda*(L - 1) + 1)*log(lambda))) / (lambda - 1)
}