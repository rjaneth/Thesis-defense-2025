al_omari_1_estimator <- function(data) {
  n <- length(data)
  m <- round(sqrt(n) + 0.5)  # m-spacing
  data_sorted <- sort(data)
  sum_term <- 0
  
  for (i in 1:n) {
    if (i <= m) {
      omega_i <- 3/2
      diff_term <- data_sorted[i + m] - data_sorted[1]
    } else if (i >= n - m + 1) {
      omega_i <- 3/2
      diff_term <- data_sorted[n] - data_sorted[i - m]
    } else {
      omega_i <- 2
      diff_term <- data_sorted[i + m] - data_sorted[i - m]
    }
    
    sum_term <- sum_term + log((n / (omega_i * m)) * diff_term)
  }
  
  return(sum_term / n)
}

# set.seed(123)
# test_data <- rnorm(10)
# result <- al_omari_1_estimator(test_data)
# print(result)
# 
# data <- c(1.2, 1.5,3.6, 2.0, 2.1, 2.5, 2.6, 3.0, 3.1, 3.5)
# result <- al_omari_1_estimator(data)
# cat("Al-Omari-1 Estimator:", result, "\n")
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

# estimación simple
#h_hat <- al_omari_1_estimator(x)
# estimación bias-corrected
#h_bc  <- bootstrap_shannon_alomari(x, B)

#cat("Estimador original :", h_hat, "\n")
#cat("Bias-corrected     :", h_bc , "\n")
