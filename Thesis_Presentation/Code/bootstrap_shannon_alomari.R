#------------------------------------------------------------------
#  Bootstrap acelerado + corrección de sesgo
#------------------------------------------------------------------
bootstrap_shannon_alomari <- function(x, B = 100L,
                                      m = round(sqrt(length(x)) + 0.5),
                                      parallel = FALSE)
{
  n  <- length(x)
  
  # estimador en los datos originales
  t0 <- shannon_alomari_estimator(x, m)
  
  # matriz de índices bootstrap  (n × B)
  idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = n)
  
  f_boot <- function(col) shannon_alomari_estimator(x[col], m)
  
  boot_vals <- if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE))
      stop("Instala 'future.apply' o usa parallel = FALSE")
    future.apply::future_apply(idx, 2, f_boot)
  } else {
    apply(idx, 2, f_boot)
  }
  
  # Estimador bias-corrected   2·t0 − mean(bootstrap)
  2 * t0 - mean(boot_vals, na.rm = TRUE)
}