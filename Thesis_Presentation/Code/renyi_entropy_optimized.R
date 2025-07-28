renyi_entropy_optimized <- function(x, lambda, m = NULL) {
  n <- length(x)
  if (n < 5) m <- 1 else if (is.null(m)) m <- round(sqrt(n) + 0.5)
  m <- max(1, min(m, floor((n - 1)/2)))  # Asegurar 1 ≤ m ≤ (n-1)/2
  
  xs <- sort(x)
  i <- 1:n
  li <- pmax(i - m, 1)     # Índice izquierdo (ajustado a 1)
  ri <- pmin(i + m, n)     # Índice derecho (ajustado a n)
  diff <- xs[ri] - xs[li]  # Espaciamiento ajustado
  
  # Coeficientes c_i (definición exacta del paper)
  ci <- ifelse(
    i <= m, (m + i - 1)/m,
    ifelse(i >= n - m + 1, (n + m - i)/m, 2)
  )
  
  # Términos del estimador (manejo de espaciamientos cero)
  r <- n * diff / (ci * m)
  valid <- diff > .Machine$double.eps & !is.na(diff)
  if (sum(valid) == 0) return(NA)  # Caso todos espaciamientos inválidos
  
  s <- mean(r[valid]^(1 - lambda), na.rm = TRUE)
  (1 / (1 - lambda)) * log(s)
}