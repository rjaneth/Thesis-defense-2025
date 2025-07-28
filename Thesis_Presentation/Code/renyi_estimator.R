# Al-Labadi 
renyi_estimator <- function(x, lambda, m = round(sqrt(length(x)) + 0.5))#m = floor(sqrt(length(x)) + 0.5)
{
  n <- length(x)
  x <- sort(x, method = "quick")
  i <- seq_len(n)
  
  li <- pmax(i - m, 1L)
  ri <- pmin(i + m, n)
  diff <- x[ri] - x[li]
  
  ci <- ifelse(i <= m, (m + i - 1)/m,
               ifelse(i >= n - m + 1, (n + m - i)/m, 2))
  
  r <- n * diff / (ci * m)
  r[r <= .Machine$double.eps] <- NA
  s <- mean(r^(1 - lambda), na.rm = TRUE)
  
  (1 / (1 - lambda)) * log(s)
}
