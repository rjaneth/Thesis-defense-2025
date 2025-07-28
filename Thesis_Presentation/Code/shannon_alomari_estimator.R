
shannon_alomari_estimator <- function(x,m = round(sqrt(length(x)) + 0.5)) {
  n <- length(x)
  #m <- round(sqrt(n) + 0.5)
  xs <- sort(x)
  sum_log <- 0
  
  for (i in 1:n) {
    if (i <= m) {
      omega_i <- 3/2
      diff <- xs[i + m] - xs[1]  # Sin recorte superior
    } else if (i >= n - m + 1) {
      omega_i <- 3/2
      diff <- xs[n] - xs[i - m]  # Sin recorte inferior
    } else {
      omega_i <- 2
      diff <- xs[i + m] - xs[i - m]
    }
    term <- log(n * diff / (omega_i * m))
    sum_log <- sum_log + term
  }
  sum_log / n
}


# #------------------------------------------------------------------
# #  Shannon–Al-Omari (ω-spacing)  —  estimador rápido y exacto
# #------------------------------------------------------------------
# shannon_alomari_estimator <- function(x, m = round(sqrt(length(x)) + 0.5))
# {
#   n  <- length(x)
#   if (n < 10)
#     warning("Muestra muy pequeña; el estimador puede ser inestable")
#   
#   xs <- sort(x, method = "quick")
#   i  <- seq_len(n)
#   
#   # índices i−m e i+m (recortados a [1,n])
#   li <- pmax(i - m, 1L)
#   ri <- pmin(i + m, n)
#   
#   diff  <- xs[ri] - xs[li]
#   
#   # ω_i
#   omega <- ifelse(i <= m | i >= n - m + 1, 3/2, 2)
#   
#   mean(log(n * diff / (omega * m)), na.rm = TRUE)
# }

# #------------------------------------------------------------------
# #  Bootstrap acelerado + corrección de sesgo
# #------------------------------------------------------------------
# bootstrap_shannon_alomari <- function(x, B = 200L,
#                                       m = floor(sqrt(length(x)) + 0.5),
#                                       parallel = FALSE)
# {
#   n  <- length(x)
#   
#   # estimador en los datos originales
#   t0 <- shannon_alomari_estimator(x, m)
#   
#   # matriz de índices bootstrap  (n × B)
#   idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = n)
#   
#   f_boot <- function(col) shannon_alomari_estimator(x[col], m)
#   
#   boot_vals <- if (parallel) {
#     if (!requireNamespace("future.apply", quietly = TRUE))
#       stop("Instala 'future.apply' o usa parallel = FALSE")
#     future.apply::future_apply(idx, 2, f_boot)
#   } else {
#     apply(idx, 2, f_boot)
#   }
#   
#   # Estimador bias-corrected   2·t0 − mean(bootstrap)
#   2 * t0 - mean(boot_vals, na.rm = TRUE)
# }
# library(future)
# plan(multisession)           # o multicore en Linux/macOS
# res <- bootstrap_shannon_alomari(x, B = 1000, parallel = TRUE)
# 
# library(invgamma)
# 
# gi0_sample <- function(mu, alpha, L, n) {
#   if (alpha >= -1) stop("alpha debe ser < -1")
#   X <- rinvgamma(n, shape = -alpha, rate = mu * (-alpha - 1))
#   Y <- rgamma(n, shape = L, rate = L)
#   X * Y
# }
# 
# set.seed(42)
# n <- 2000
# B <- 200
# x <- gi0_sample(mu = 3, alpha = -5, L = 4, n)
# 
# # estimación simple
# h_hat <- shannon_alomari_estimator(x)
# # estimación bias-corrected
# h_bc  <- bootstrap_shannon_alomari(x, B)
# 
# cat("Estimador original :", h_hat, "\n")
# cat("Bias-corrected     :", h_bc , "\n")



