# Start the timer
start_time <- Sys.time()

# Estimators
source("./Code/al_omari_1_estimator.R")
source("./Code/bootstrap_al_omari_1_estimator.R")

source("./Code/renyi_entropy_estimator_v1.R")
source("./Code/bootstrap_renyi_entropy_estimator_v1.R")

# Tsallis estimators
source("./Code/tsallis_estimator_optimized.R")
source("./Code/bootstrap_tsallis_entropy_optimized.R")

source("./Code/read_ENVI_images.R") # to read the images in ENVI format

# Load and read data from SAR images
x <- myread.ENVI(
  file = './Code/SAR/envi_dublin_1100_HH/Intensity_HH.img', 
  headerfile = './Code/SAR/envi_dublin_1100_HH/Intensity_HH.hdr'
)

rows <- nrow(x)
cols <- ncol(x)

# Parameters
L <- 16           # Number of looks
B <- 100          # Bootstrap replications
lambda <- 0.9     # Tsallis entropy order
window_size <- 7  # Sliding window size

difference_values <- matrix(
  NA, 
  nrow = rows - window_size + 1, 
  ncol = cols - window_size + 1
)

# Iterate over sliding windows
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Select local window
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
    mu_hat <- mean(window_data)
    
    difference_values[i, j] <- bootstrap_tsallis_entropy_optimized(window_data, B = B, lambda = lambda) - 
      ((1 - exp((1 - lambda) * log(mu_hat) +
                  (lambda - 1) * log(L) +
                  lgamma(lambda * (L - 1) + 1) -
                  lambda * lgamma(L) -
                  (lambda * (L - 1) + 1) * log(lambda))) / (lambda - 1))
  }
}

# Save the results
save(difference_values, x, file = "./Data/results_dublin_tsallis_w7_09_L16_b100.Rdata")

# Stop the timer
end_time <- Sys.time()
execution_time <- end_time - start_time
execution_time
