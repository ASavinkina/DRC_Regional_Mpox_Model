# Code to run the mpox regional model multiple times
# 


setwd("/home/as3963/project/DRC_Regional_Mpox_Model")

# Call in model code:
source("Stochastic_Model_Code_2.26.25.R")
reactions <- unlist(reactions)


n=10

# Replace the sequential lapply with parallel processing
library(parallel)
# Create a "cluster" of CPU cores to run the model on
cl <- makeCluster(n)
# Export necessary objects and functions to the CPU-cluster
clusterExport(cl, c(ls(), "ssa", "ssa.btl"))
# Load required packages on each worker
clusterEvalQ(cl, {
  library(GillespieSSA)
})

# Run the simulations in parallel
out <- parLapply(cl = cl, X = 1:n, fun = function(x) {
  # Set a different random seed for each run based on the run number
  # This ensures that each run is both independent and reproducible
  set.seed(12345 + x)
  ssa(y_init, reactions, t(nu), params, tf = 100, method = ssa.btl(f = 10))$data
})

# Stop the cluster when done
stopCluster(cl)

# Add new column to each list in "out" numbering the run

Run <- c(1:n)

out2 <- mapply(cbind, out, "Run"=Run, SIMPLIFY=F)

# Create one long dataset of all output from all runs

output_all <- as.data.frame(do.call(rbind, out2))


write.csv(output_all, file="Full_output_test_32125.csv")