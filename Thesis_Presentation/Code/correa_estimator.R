correa_estimator <- function(data) {
  n <- length(data)
  m <- round(sqrt(n) + 0.5)
  X_sorted <- sort(data)
  C_mn <- 0
  
  for (i in 1:n) {
    i_minus_m <- max(1, i - m)
    i_plus_m <- min(n, i + m)
    
    X_i <- X_sorted[i]
    X_window <- X_sorted[i_minus_m:i_plus_m]
    X_bar_i <- mean(X_window)
    
    num <- sum((X_window - X_bar_i) * (i_minus_m:i_plus_m - i))
    den <- n * sum((X_window - X_bar_i)^2)
    
    if (den == 0) {
      b_i <- 0
    } else {
      b_i <- num / den
    }
    
    if (b_i > 0) {
      C_mn <- C_mn - log(b_i)
    }
  }
  
  C_mn <- C_mn / n
  
  return(C_mn)
}


# # # Example usage:
# data <- data <- c( 0.429283, 0.429283, 0.9444219, 0.9608317, 1.451404, 1.451404, 0.429283, 0.429283, 0.429283 )#c(1.2, 1.5,3.6, 2.0, 2.1, 2.5, 2.6, 3.0, 3.1, 3.5 )
# result <- correa_estimator(data)
# cat("Correa Estimator:", result, "\n")



