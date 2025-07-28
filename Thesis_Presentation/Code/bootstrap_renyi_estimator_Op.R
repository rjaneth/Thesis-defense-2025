bootstrap_renyi_estimator_Op <- function(x, B = 100L, lambda,
                                      m = round(sqrt(length(x)) + 0.5),
                                      parallel = FALSE)
{
  n <- length(x)
  
  # estimador en los datos originales
  t0 <- renyi_entropy_optimized (x, lambda, m)
  
  # matriz de índices bootstrap (n × B)
  idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = n)
  
  # función auxiliar (usa la misma m para todas las réplicas)
  f_boot <- function(col) renyi_entropy_optimized (x[col], lambda, m)
  
  # vector de estimaciones bootstrap
  if (parallel) {
    boot_vals <- future.apply::future_apply(idx, 2, f_boot)
  } else {
    boot_vals <- apply(idx, 2, f_boot)
  }
  
  # sesgo estimado  = mean(bootstrap) - t0
  bias_hat <- mean(boot_vals, na.rm = TRUE) - t0
  
  # versión “bias-corrected”   2·t0 − mean(bootstrap)
  2 * t0 - mean(boot_vals, na.rm = TRUE)
}