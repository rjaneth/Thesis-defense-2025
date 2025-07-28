



#source("../MainFunctions/correa_estimator.R")


bootstrap_correa_estimator <-  function(x, B){
  
  v.Bootstrap <- rep(0, B)
  
  for(b in 1:B) {
    same_values <- TRUE
    while (same_values) {
      bx <- sample(x, replace = TRUE)
      if (!all(bx == bx[1])) {
        same_values <- FALSE
        entropy_result <- correa_estimator(bx)
        if (is.finite(entropy_result)) {
          v.Bootstrap[b] <- entropy_result
          #cat("Bootstrap Sample", b, ": ", bx, "\n")
          #cat("Entropy for Bootstrap Sample", b, ": ", v.Bootstrap[b], "\n")
        }
      }
    }
  }
  
  
  t <- correa_estimator(x)
  estimated_mean <- mean(v.Bootstrap[!is.na(v.Bootstrap)])
  
  return(2*t - estimated_mean)
  
}

# data <- c(1.2, 1.5,3.6, 2.0, 2.1, 2.5, 2.6, 3.0, 3.1, 3.5)
# result <- bootstrap_correa_estimator(data, 10)
# cat("Al-Omari-1 Estimator:", result, "\n")



