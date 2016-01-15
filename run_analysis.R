
# This script is to:
# (under given environmental conditions: ca, k, MAP (see below for explanations))
# 1. identify the correct value for Abar in Eq. 14 by minimizing the difference
#    between the input, q (the assumed value of Abar) with the output, the resulting Abar
# 2. solve Eq. 14
# 3. derive average soil water level (w), transpiration rate (E), g1, ci/ca ratio, marginal water cost for
#    photosynthesis (lambda), transpiration rainfall ratio (EP), and theprobability at 0 soil water content (p0)

# Functions
source("1. Average A (Vcmax).r")
source("2. Derived variables.r")

#environmental conditions
ca1 <- c(400)  # Atmospheric CO2 concentration (ppm)
k1 <- c(0.025, 0.05, 0.1) # Rainfall frequency (per day)
MAP <- seq(0.5, 10, by=0.5)*365 # MAP=MDP*365; MAP: mean annual precipitation; MDP: mean daily precipitation
env1 <- as.vector(expand.grid(ca1, k1, MAP))
ca2 <- c(800)
k2 <- c(0.05)
env2 <- as.vector(expand.grid(ca2, k2, MAP))
env <- rbind(env1, env2)

# Initialize
dvs <- matrix(nrow=nrow(env), ncol=8)

# Run every parameter combination
for(i in 1:nrow(env)){
  begin <- proc.time()
  
  # Identify the correct value for Abar under given environmental conditions
  dvs[i, 1] <- qopt(env[i, 1], env[i, 2], env[i, 3])
  
  # Solve Eq. 14 and derive variables of interest
  dvs[i, 2:8] <- dvsf(env[i, 1], env[i, 2], env[i, 3], dvs[i, 1])
  end <- proc.time()
  message(sprintf("%s/%s completed in %.2f min",i, nrow(env), (end[3]-begin[3])/60))
}

# Collect results
res <- cbind(env, dvs)
colnames(res) <- c("ca", "k", "MAP", "A", "w", "E", "g1", "cica", "p0", "lambda", "EP") 

write.csv(res, "Derived variables.csv", row.names = FALSE)
