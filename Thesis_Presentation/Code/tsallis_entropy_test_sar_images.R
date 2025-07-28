# ──────────────────────────────────────────────────────────────
#  0.  Load libraries and helper functions
# ──────────────────────────────────────────────────────────────
library(future.apply)          # optional for parallelization
plan(sequential)               # switch to multisession for multicore CPU
source("./Code/read_ENVI_images.R")

#  Estimators and utilities
source("./Code/gamma_sar_sample.R")                       # if used elsewhere
source("./Code/tsallis_estimator_optimized.R")
source("./Code/bootstrap_tsallis_entropy_optimized.R")

# ──────────────────────────────────────────────────────────────
#  1.  Experiment parameters
# ──────────────────────────────────────────────────────────────
img_path   <- "./Code/SAR/envi_dublin_1100_HH/Intensity_HH.img"
hdr_path   <- "./Code/SAR/envi_dublin_1100_HH/Intensity_HH.hdr"

L          <- 16        # number of looks
lambda     <- 0.9       # Tsallis order
B          <- 100       # bootstrap replicas
window_size<- 7         # window size (p×p)
parallel   <- TRUE      # TRUE = use all available cores

# ──────────────────────────────────────────────────────────────
# 2.  Read SAR image (matrix rows × cols)
# ──────────────────────────────────────────────────────────────
x <- myread.ENVI(file = img_path, headerfile = hdr_path)

rows <- nrow(x);  cols <- ncol(x)
out_rows <- rows - window_size + 1
out_cols <- cols - window_size + 1
difference_values <- matrix(NA_real_, nrow = out_rows, ncol = out_cols)

# ──────────────────────────────────────────────────────────────
# 3.  Function to compute the statistic over a window
# ──────────────────────────────────────────────────────────────
stat_tsallis_window <- function(mat_win) {
  z        <- as.vector(mat_win)
  mu_hat   <- mean(z)
  
  est  <- bootstrap_tsallis_entropy_optimized(z, B = B,
                                              lambda = lambda, parallel = FALSE)
  
  theo <- (1 - exp((1 - lambda) * log(mu_hat) +
                     (lambda - 1) * log(L) +
                     lgamma(lambda * (L - 1) + 1) -
                     lambda * lgamma(L) -
                     (lambda * (L - 1) + 1) * log(lambda))) /
    (lambda - 1)
  
  est - theo
}

# ──────────────────────────────────────────────────────────────
# 4.  Sliding window loop (can be parallelized with future.apply)
# ──────────────────────────────────────────────────────────────
start_time <- Sys.time()

# index of all top-left corners of each window
idx <- expand.grid(i = seq_len(out_rows), j = seq_len(out_cols))

compute_one <- function(k) {
  i <- idx$i[k];  j <- idx$j[k]
  win <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
  stat_tsallis_window(win)
}

if (parallel) {
  res <- future_sapply(seq_len(nrow(idx)), compute_one)
} else {
  res <- sapply(seq_len(nrow(idx)), compute_one)
}

difference_values[] <- matrix(res, nrow = out_rows, byrow = FALSE)

end_time <- Sys.time()
cat("Total time:", round(end_time - start_time, 2), "sec.\n")

# ──────────────────────────────────────────────────────────────
# 5.  Save results and the image itself
# ──────────────────────────────────────────────────────────────
save(difference_values, x,
     file = "./Data/dublin_tsallis_w7_lambda09_L16_b100.Rdata")


# # ──────────────────────────────────────────────────────────────
# #  0.  Carga de librerías y funciones auxiliares
# # ──────────────────────────────────────────────────────────────
# library(future.apply)          # opcional para paralelizar
# plan(sequential)               # cambia a multisession para CPU multinúcleo
# source("./Code/read_ENVI_images.R")
# 
# #  Estimadores y utilidades
# source("./Code/gamma_sar_sample.R")                       # si lo usas en otras partes
# source("./Code/tsallis_estimator_optimized.R")
# source("./Code/bootstrap_tsallis_entropy_optimized.R")
# 
# # ──────────────────────────────────────────────────────────────
# #  1.  Parámetros del experimento
# # ──────────────────────────────────────────────────────────────
# img_path   <- "./Code/SAR/envi_dublin_1100_HH/Intensity_HH.img"
# hdr_path   <- "./Code/SAR/envi_dublin_1100_HH/Intensity_HH.hdr"
# 
# L          <- 16        # número de looks
# lambda     <- 0.9       # orden de Tsallis
# B          <- 100       # réplicas bootstrap
# window_size<- 7         # tamaño de la ventana (p×p)
# parallel   <- TRUE      # TRUE = usa todos los núcleos disponibles
# 
# # ──────────────────────────────────────────────────────────────
# # 2.  Lectura de la imagen SAR (matriz rows × cols)
# # ──────────────────────────────────────────────────────────────
# x <- myread.ENVI(file = img_path, headerfile = hdr_path)
# 
# rows <- nrow(x);  cols <- ncol(x)
# out_rows <- rows - window_size + 1
# out_cols <- cols - window_size + 1
# difference_values <- matrix(NA_real_, nrow = out_rows, ncol = out_cols)
# 
# # ──────────────────────────────────────────────────────────────
# # 3.  Función que calcula el estadístico en una ventana
# # ──────────────────────────────────────────────────────────────
# stat_tsallis_window <- function(mat_win) {
#   z        <- as.vector(mat_win)
#   mu_hat   <- mean(z)
#   
#   est  <- bootstrap_tsallis_entropy_optimized(z, B = B,
#                                               lambda = lambda, parallel = FALSE)
#   
#   theo <- (1 - exp((1 - lambda) * log(mu_hat) +
#                      (lambda - 1) * log(L) +
#                      lgamma(lambda * (L - 1) + 1) -
#                      lambda * lgamma(L) -
#                      (lambda * (L - 1) + 1) * log(lambda))) /
#     (lambda - 1)
#   
#   est - theo
# }
# 
# # ──────────────────────────────────────────────────────────────
# # 4.  Bucle deslizante  (se puede paralelizar con future.apply)
# # ──────────────────────────────────────────────────────────────
# start_time <- Sys.time()
# 
# # índice de todas las esquinas superiores de cada ventana
# idx <- expand.grid(i = seq_len(out_rows), j = seq_len(out_cols))
# 
# compute_one <- function(k) {
#   i <- idx$i[k];  j <- idx$j[k]
#   win <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
#   stat_tsallis_window(win)
# }
# 
# if (parallel) {
#   res <- future_sapply(seq_len(nrow(idx)), compute_one)
# } else {
#   res <- sapply(seq_len(nrow(idx)), compute_one)
# }
# 
# difference_values[] <- matrix(res, nrow = out_rows, byrow = FALSE)
# 
# end_time <- Sys.time()
# cat("Tiempo total:", round(end_time - start_time, 2), "seg.\n")
# 
# # ──────────────────────────────────────────────────────────────
# # 5.  Guarda resultados y la propia imagen
# # ──────────────────────────────────────────────────────────────
# save(difference_values, x,
#      file = "./Data/dublin_tsallis_w7_lambda09_L16_b100.Rdata")
