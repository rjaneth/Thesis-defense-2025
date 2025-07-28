# Function for Entropy Calculation of the Gamma SAR Distribution
# Description: This function calculates the entropy of the Gamma SAR distribution
# based on provided parameters L and mu.
#
# Input:
#   - L: Shape parameter of the Gamma distribution.
#   - mu: Mean parameter of the distribution.
#
# Output:
#   - Entropy value of the Gamma SAR distribution.

entropy_gamma_sar <- function(L, mu) {
  
  return(L - log(L) + lgamma(L) + (1 - L) * digamma(L) + log(mu))
  
}

#entropy_gamma_sar(5,10)