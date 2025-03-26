# Code to run the mpox regional model multiple times
# 

# Call in model code:
source("Stochastic_Model_Code_2.26.25.R")

n=1
#set.seed(NULL)


# Run model N times and save output to a list

out <- lapply(X=1:n,FUN=function(x) ssa(y_init, reactions2, nu2, params, tf = 365, method=ssa.btl(f=5))$data)


#out <- ssa(y_init, reactions2, nu2, params, tf = 365, method=ssa.btl(f=5))

# Add new column to each list in "out" numbering the run

Run <- c(1:n)

out2 <- mapply(cbind, out, "Run"=Run, SIMPLIFY=F)

# Create one long dataset of all output from all runs

output_all <- as.data.frame(do.call(rbind, out2))

write.csv(output_all, file="Output/Full_output_test_32125.csv")