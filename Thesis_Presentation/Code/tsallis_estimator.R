#------------------------------------------------------------
# 1.  Estimador base (spacing) de Tsallis  –  vectorizado
#------------------------------------------------------------
tsallis_estimator <- function(x, alpha,
                                    m = floor(sqrt(length(x)) + 0.5))
{
  if (alpha <= 0) stop("alpha debe ser > 0 (Tsallis)")
  n <- length(x)
  if (n < 10) warning("Muestra muy pequeña; resultados inestables")
  
  xs <- sort(x, method = "quick")
  i  <- seq_len(n)
  
  li <- pmax(i - m, 1L)         # índices i-m (cortados a 1)
  ri <- pmin(i + m, n)          # índices i+m (cortados a n)
  
  diff <- xs[ri] - xs[li]
  
  ci <- ifelse(i <= m,             (m + i - 1)/m,
               ifelse(i >= n - m + 1,     (n + m - i)/m, 2))
  
  r <- n * diff / (ci * m)
  r[r <= .Machine$double.eps] <- NA
  
  s <- sum(r^(1 - alpha), na.rm = TRUE)  # ∑ r_i^{1-α}
  
  (1 / (alpha - 1)) * (1 - s / n)
}