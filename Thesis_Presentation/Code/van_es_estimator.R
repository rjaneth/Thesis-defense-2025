# Non-Parametric Shannon Entropy Estimator
# Description: This function calculates Shannon entropy using a non-parametric estimator,
# known as the van Es estimator. 
# It is based on spacings between data points 
#
# Input:
#   - data: sample or dataset 
#
# Output:
#   - Estimated Shannon entropy 

van_es_estimator <- function(data) {
  n <- length(data)
  m <- round(sqrt(n)+0.5)# round(sqrt(n))# m-spacing # round(sqrt(n)+0.5)
  data_sorted <- sort(data)
  #cat("sort data ", data_sorted ,   "\n")
  sum_term1 <- 0
  sum_term2 <- 0
  
  for (i in 1:(n - m)) {
    term1 <-  log(((n + 1) / m) * (data_sorted[i + m] - data_sorted[i]))
    sum_term1 <-  sum_term1 + term1
    
  }
 
  for (k in m:n) {
    sum_term2 <- sum_term2 + 1 / k
   
  }
  term3 <- log(m / (n + 1))
 
  return((sum_term1 / (n - m)) + sum_term2+ term3)
  
}

# # Example usage:
# data <- c(1.229641, 1.451404, 1.229641, 0.9444219, 1.229641, 0.429283, 0.6902405, 1.904591, 0.4554712)
# # data <- c( 0.429283, 0.429283, 0.9444219, 0.9608317, 1.451404, 1.451404, 0.429283, 0.429283, 0.429283 ) # -inf
#  result <- van_es_estimator(data)
#  cat("van es :", result, "\n")