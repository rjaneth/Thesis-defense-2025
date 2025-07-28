tsallis_gi0 <- function(alpha, mu, L, lambda){ 
  (1 - exp((1 - lambda)*log(mu) +
             (lambda - 1)*log(L) +
             (1 - lambda)*log(-alpha - 1) +                     
             lgamma(lambda*(L - 1) + 1) -
             lambda*lgamma(L) +
             lambda*(lgamma(L - alpha) - lgamma(-alpha)) +      
             lgamma(lambda*(1 - alpha) - 1) -
             lgamma(lambda*(L - alpha)))) /(lambda - 1)
}

tsallis_gammasar_log <- function(mu, L, lambda){ 
  (1 - exp((1 - lambda)*log(mu) +
             (lambda - 1)*log(L) +
             lgamma(lambda*(L - 1) + 1) -
             lambda*lgamma(L) -
             (lambda*(L - 1) + 1)*log(lambda))) /
  (lambda - 1)

}

#tsallis_gammasar_log(1, 5, 0.9)
 #tsallis_gi0(-1.5, 1, 5, 0.9)
# 
# 
# tsallis_gammasar <- function(mu, L, lambda) {
#   # Validación de parámetros
#   if (mu <= 0 || L <= 0 || lambda <= 0) {
#     stop("Parámetros deben ser positivos: mu > 0, L > 0, lambda > 0")
#   }
#   
#   # Término exponencial-logarítmico
#   exponent <- (1 - lambda) * log(mu) +
#     (lambda - 1) * log(L) +
#     lgamma(lambda * (L - 1) + 1) -
#     lambda * lgamma(L) -
#     (lambda * (L - 1) + 1) * log(lambda)
#   
#   # Entropía de Tsallis
#   (1 - exp(exponent)) / (lambda - 1)
# }
# 
tsallis_gi01 <- function(alpha, mu, L, lambda) {
  # Validación de parámetros
  if (alpha >= -1) stop("alpha debe ser < -1")
  if (mu <= 0 || L <= 0 || lambda <= 0) {
    stop("Parámetros deben ser positivos: mu > 0, L > 0, lambda > 0")
  }

  # Calcular término base (Gamma SAR)
  T_base <- tsallis_gammasar(mu, L, lambda)

  # Calcular J0 (en espacio logarítmico)
  log_J0 <- (1 - lambda) * log(mu) +
    (lambda - 1) * log(L) +
    lgamma(lambda * (L - 1) + 1) -
    lambda * lgamma(L) -
    (lambda * (L - 1) + 1) * log(lambda)

  # Calcular Qα (en espacio logarítmico)
  b <- lambda * (1 - alpha) - 1
  log_Q_alpha <- (1 - lambda) * log(-alpha - 1) +
    lambda * (lgamma(L - alpha) - lgamma(-alpha)) +
    lgamma(b) -
    lgamma(lambda * (L - alpha)) +
    (lambda * (L - 1) + 1) * log(lambda)

  # Calcular Δα
  J0_val <- exp(log_J0)
  Q_alpha_val <- exp(log_Q_alpha)
  Delta_alpha <- (J0_val / (lambda - 1)) * (1 - Q_alpha_val)

  # Entropía total
  T_base + Delta_alpha
}

#tsallis_gi01(-1.5, 1, 5, 0.9)
