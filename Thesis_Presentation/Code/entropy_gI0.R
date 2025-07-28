# Function for Entropy Calculation of the gI0 Distribution
# Description: This function computes the entropy of the GI0 distribution
# based on provided parameters mu, alpha, and L.
#
# Input:
#   - mu: Mean parameter of the distribution.
#   - alpha: Shape parameter.
#   - L: Shape parameter of the distribution.
#
# Output:
#   - Entropy value of the gI0 distribution.
#
#shannon
entropy_gI0 <- function(mu, alpha, L) {
  
  term1 <- L - log(L) + lgamma(L) + (1 - L) * digamma(L) + log(mu) 
  term2 <- -L - lgamma(L-alpha) + (L-alpha)*(digamma(L - alpha))- (1-alpha)*digamma(- alpha)+log(-1 - alpha)+lgamma(-alpha)
  
  entropy <- term1 + term2 
  return(entropy)
}

# entropy_gI0( 10,-2.5,  5)
# (L - log(L) + lgamma(L) + (1 - L) * digamma(L) + log(mu)-L - lgamma(L-alpha) + (L-alpha)*(digamma(L - alpha))- (1-alpha)*digamma(- alpha)+log(-1 - alpha)+lgamma(-alpha))