# Calculate R0 for countries with Mpox Clade 1b only: this should be the R0 of Clade 1b.
# 
# 
# Using Burundi: R0 for sexual transmission of Clade 1b would be ~  4.144812[ 3.599115 , 4.773768 ],
# lower end: Reproduction number estimate using  Attack Rate  method.
#R :  1.002538[ 1.002454 , 1.002623 ]
#Upped end:
#Reproduction number estimate using  Maximum Likelihood  method.
#R :  9.917741[ 9.225208 , 10.64454 ]

library(tidyverse)
library(R0)

# Read in epidemic curve data for countries:
# 

data <- read.csv("Data/weekly AFR cases by country as of 16 February 2025.csv")

#Also read in demo data for population
#
#
pop <- read.csv("Data/Demographic_data_countries_2_27_25.csv")

data_Burundi <- data[which(data$country=="Burundi"),]

# Convert date
data_Burundi$week_end_date <- as.Date(data_Burundi$week_end_date)

# Order data by date
data_Burundi <- data_Burundi[order(data_Burundi$week_end_date), ]

# Define generational interval for mpox (assumed mean 11.3 days, SD 3 days based on 2022 outbreak estimates)
si_mean <- 11.3
si_sd <- 3

# population estimate
# 
pop_c <- pop[which(pop$Country=="Burundi"),"U5"] +pop[which(pop$Country=="Burundi"),"U15"] +pop[which(pop$Country=="Burundi"),"O15"] 

cases <- data_Burundi$new_confirmed_cases

cases <- cases[2:31]

# Ensure generation time is properly defined
mGT <- generation.time("gamma", c(si_mean, si_sd))

# Estimate R0 using multiple methods
estR0 <- estimate.R(
  cases,
  mGT,
  # begin = 2,
  # end= 31,
  methods = c("EG", "ML", "TD", "AR", "SB"),
  pop.size = (pop_c*0.05),
  nsim = 1000
)

# Print estimated R0
print(estR0)

# Uganda
# 
data_Uganda <- data[which(data$country=="Uganda"),]

# Convert date
data_Uganda$week_end_date <- as.Date(data_Uganda$week_end_date)

# Order data by date
data_Uganda <- data_Uganda[order(data_Uganda$week_end_date), ]


# population estimate
# 
pop_c <- pop[which(pop$Country=="Uganda"),"U5"] +pop[which(pop$Country=="Uganda"),"U15"] +pop[which(pop$Country=="Uganda"),"O15"] 

cases <- data_Uganda$new_confirmed_cases

cases <- cases[2:33]

# Ensure generation time is properly defined
mGT <- generation.time("gamma", c(si_mean, si_sd))

# Estimate R0 using multiple methods
estR0 <- estimate.R(
  cases,
  mGT,
  # begin = 2,
  # end= 31,
  methods = c("EG", "ML", "TD", "AR", "SB"),
  pop.size = pop_c,
  nsim = 1000
)

# Print estimated R0
print(estR0)



