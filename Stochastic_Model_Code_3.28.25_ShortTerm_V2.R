library(GillespieSSA)
library(ggplot2)
library(tidyverse)
library(gdata)
library(reshape2)

# This code re-codes model into GillespieSSA2 package, which speeds up processing significantly.
# Currently running for max 1 year, as births/deaths and aging is turned off for short-term analyses
# Ie Uganda/DRC vaccine model
# Also, currently only running DRC and Uganda


# Read in demographics etc data

demog = read.csv("Data/Demographic_data_countries_3_18_25.csv")

demog$location <- demog$Country
demog$location <- ifelse(demog$Country=="DRC",demog$Province, demog$Country)

# Read in transmission movement matrix, and clean row/column names

movement = read.csv("Data/transmission_movement_matrix_core.csv", row.names = 1)

rownames(movement)[rownames(movement) == "South Sudan"] = "SouthSudan"
colnames(movement)[colnames(movement) == "South.Sudan"] = "SouthSudan"

rownames(movement)[rownames(movement) == "Central African Republic"] = "CAR"
colnames(movement)[colnames(movement) == "Central.African.Republic"] = "CAR"

rownames(movement)[rownames(movement) == "South Sudan"] = "SouthSudan"
colnames(movement)[colnames(movement) == "South_Sudan"] = "SouthSudan"

names(movement) <- gsub(x = names(movement), pattern = "\\_", replacement = "")

rownames(movement) <- gsub(x = rownames(movement), pattern = "\\_", replacement = "")

demog$location <-gsub(x = demog$location, pattern = "\\-", replacement = "")


### Transform movement from matrix to a named list so it can be used within 
### parameters list
### 

movement4 <- as.data.frame(movement)

movement4$Country1 <- rownames(movement4)

movement5 <- gather(movement4, key="Country", value="Value", -Country1)

movement5$Country2 <- paste0(movement5$Country1,"_",movement5$Country)

movement5$Value <- as.numeric(movement5$Value)

movement6 <- setNames(as.list(movement5$Value), nm = movement5$Country2)

movement7 <- unlist(movement6)

# Add columns to demographics for pop size of high risk and low risk adult age groups
# This data comes from demographics file, elevated for areas of DRC (and for Uganda) where
# sexual transmission is more likely

demog$O15HR <- demog$O15 * demog$percHR
demog$O15LR <- demog$O15 * (1-demog$percHR)


# Take average of cases from 2021 and 2022 (pre-epidemic) to estimate number of sylvatic infections for DRC
# 

demog$avg_sylvatic <- (demog$Cases_2021 + demog$Cases_2022)/2

# Proportion of sylvatic infections in each age group, calibrated in DRC paper

prop_sylvatic_adult = 0.3 #proportion of sylvatic infections in adults
prop_sylvatic_u15 = 0.45 #0.6 #proportion of sylvatic infections in those 5-15

# New dataframe of only data needed for sylvatic transmission (exogenous shock) calculation

exog_shock <- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("Location", "exogshock_u5", "exogshock_u15","exogshock_o15HR", "exogshock_o15LR")
colnames(exog_shock) <- x

# Calculate rate of exogenous shocks in each age group: this takes the average number of sylvatic infections
# from pre-epidemic years, multiplies it by the proportion attributed to each age group, and then divides by population
# by age. Then divides it all by 365 to give a daily rate.

for (i in 1:length(demog$location)) {
  
  exog_shock[i,1] <- demog[i,19]
  
  # this assumes that half of all sylvatic transmission happens in adults, and 1/4 each in children's age groups
  exog_shock[i,2]= (((demog[i,22])*(1-prop_sylvatic_adult-prop_sylvatic_u15))/demog[i,4]) / 365
  exog_shock[i,3]= (((demog[i,22])*prop_sylvatic_u15)/demog[i,5]) / 365
  exog_shock[i,4]= ((((demog[i,22])*prop_sylvatic_adult)/demog[i,6])) / 365
  exog_shock[i,5]= ((((demog[i,22])*prop_sylvatic_adult)/demog[i,6])) / 365
  
  # sets exogenous shock rate to 0 for countries with no data
  
  exog_shock[is.na(exog_shock)] <- 0
  
}

# Transform into lists which can be appended to parameters list later

exog_shock_u5 <- setNames(as.list(exog_shock[,2]), nm=paste0(exog_shock$Location,'_','exogshock_U5'))
exog_shock_u15 <- setNames(as.list(exog_shock[,3]), nm=paste0(exog_shock$Location,'_','exogshock_U15'))
exog_shock_o15HR <- setNames(as.list(exog_shock[,4]), nm=paste0(exog_shock$Location,'_','exogshock_O15HR'))
exog_shock_o15LR <- setNames(as.list(exog_shock[,5]), nm=paste0(exog_shock$Location,'_','exogshock_O15LR'))

exog_shock_u5 <- unlist(exog_shock_u5)
exog_shock_u15 <- unlist(exog_shock_u15)
exog_shock_o15HR <- unlist(exog_shock_o15HR)
exog_shock_o15LR <- unlist(exog_shock_o15LR)

# Multiplier to attribute infections to clade 1a vs clade 1b- needs to be calibrated

exog_shock_multi_1a <- 0.95
exog_shock_multi_1b <- 0.05

# For when we use births/deaths rate- taken from country birth rate, meant to balance.

#births <- demog[,c(17,18)]

#births$birth_rate <- births$birth_per_1000/1000/365

# Additional inputs from past paper.
# Transmission parameters are calculated by considering number of contacts by and proportion of contacts with each other age,
# as well as secondary attack rate.

recoveryrate=1/21

contacts_u5 = 5.44 # number of daily physical contacts, under 5
contacts_u15 = 6.88 # number of daily physical contacts, under 15
contacts_o15 = 4.788 # number of daily physical contacts, over 15
contacts_sex = 10/30 # sexual contacts for high risk, daily (average 15 a month- made up)

contactprop_u5u5 = 0.235 # proportion of physical contacts for under 5 with under 5
contactprop_u5u15 = 0.235 # proportion of physical contacts for under 5 with under 15
contactprop_u5o15 = 1-contactprop_u5u5-contactprop_u5u15 # proportion of physical contacts for under 5 with over 15

contactprop_u15u5 = 0.235# proportion of physical contacts for under 15 with under 5
contactprop_u15u15 = 0.235 # proportion of physical contacts for under 15 with under 15
contactprop_u15o15 = 1-contactprop_u15u5-contactprop_u15u15 # proportion of physical contacts for under 15 with over 15

contactprop_o15u5 = 0.135 # proportion of physical contacts for over 15 with under 5
contactprop_o15u15 = 0.135 # proportion of physical contacts for over 15 with under 15
contactprop_o15o15 = 1-contactprop_o15u5-contactprop_o15u15 # proportion of physical contacts for over 15 with over 15

SA_u5 = 0.1168/3 # secondary attack rate when infected person is under 5, halved to account for data coming from household transmission
SA_u15 = 0.0803/4 # secondary attack rate when infected person is under 15, halved to account for data coming from household transmission
SA_o15 = 0.0133/2 # secondary attack rate when infected person is over 15, halved to account for data coming from household transmission
SA_sex = 0.005 # secondary attack rate from sex, made up for now

# Mpox mortality by age
mortality_u5 = 0.10 # mortality rate under 5
mortality_u15 = 0.08 # mortality rate under 15
mortality_o15 = 0.035 # mortality rate over 15


# Create parameters list for model runs

params <- c(
  #beta = 0.3,       # Transmission rate
  rho = 1/21,      # Recovery rate
  mu = 1/50,        # Birth and natural death rate (assuming lifespan of 50 years)
  aging1= 1/(5*365),
  aging2 = 1/(10*365), #c(0,0),#c(1/(5*365),1/(10*365)),#c(1/20, 1/30),  # Aging rates from young to adult and adult to senior
  #demog= demog,
  #movement = movement,
  #exog_shock = exog_shock,
  #percHR = percHR,
  beta_U5U5= SA_u5*contacts_u5*contactprop_u5u5, #beta for transmission between u5 and u5
  beta_U5U15= SA_u5*contacts_u5*contactprop_u5u15, #beta for transmission between u5 and u15
  beta_U5O15LR= SA_u5*contacts_u5*contactprop_u5o15,#beta for transmission between u5 and o15
  beta_U5O15HR= SA_u5*contacts_u5*contactprop_u5o15,#beta for transmission between u5 and o15
  
  beta_U15U5= SA_u15*contacts_u15*contactprop_u15u5, #beta for transmission between u15 and u5
  beta_U15U15= SA_u15*contacts_u15*contactprop_u15u15, #beta for transmission between u15 and u15
  beta_U15O15LR= SA_u15*contacts_u15*contactprop_u15o15, #beta for transmission between u15 and o15
  beta_U15O15HR= SA_u15*contacts_u15*contactprop_u15o15, #beta for transmission between u15 and o15
  
  beta_O15LRU5= SA_o15*contacts_o15*contactprop_o15u5, #beta for transmission between o15 and u5
  beta_O15LRU15= SA_o15*contacts_o15*contactprop_o15u15, #beta for transmission between o15 and u15
  beta_O15LRO15LR = SA_o15*contacts_o15*contactprop_o15o15*.95, #beta for transmission between o15 and u15
  beta_O15LRO15HR = SA_o15*contacts_o15*contactprop_o15o15*.05, #beta for transmission between o15 and u15
  
  beta_O15HRU5= SA_o15*contacts_o15*contactprop_o15u5, #beta for transmission between o15 and u5
  beta_O15HRU15= SA_o15*contacts_o15*contactprop_o15u15, #beta for transmission between o15 and u15
  beta_O15HRO15LR = SA_o15*contacts_o15*contactprop_o15o15*.95, #beta for transmission between o15 and u15
  beta_O15HRO15HR = SA_o15*contacts_o15*contactprop_o15o15*.05, #beta for transmission between o15 and u15
  
  beta_sex = 2.4*(recoveryrate),
  
  # (u5_pop/(population)*(beta_u5u5+beta_u5u15+beta_u5o15)/3 + u15_pop/(population)*(beta_u15u5+beta_u15u15+beta_u15o15)/3 +
  #  o15_pop/(population)*(beta_o15u5+beta_o15u15+beta_o15o15)/3)*21 # Estimate of R0
  #  
 # rho= recoveryrate, #recovery rate
  mu_U5 = mortality_u5/21, #death rate, under 5
  mu_U15 = mortality_u15/21, #death rate, under 15
  mu_O15LR = mortality_o15/21, #death rate, over 15
  mu_O15HR = mortality_o15/21, #death rate, over 15
  #theta_U5= exograte_u5, #rate of exogenous shocks
  #theta_U15= exograte_u15, #rate of exogenous shocks
  #theta_O15HR= exograte_o15*propHR, #rate of exogenous shocks
  #theta_O15LR = exograte_o15*(1-propHR),
  omega_U5= 0, #vaccination rate, under 5
  omega_U15= 0, #vaccination rate, under 15
  omega_O15LR= 0, #vaccination rate, over 15 low risk
  omega_O15HR= 0,  #vaccination rate, over 15 high risk 
  #theta=0.01,
  omega=0,
  exog_shock_multi_1a <- exog_shock_multi_1a
)


# Append all parameter lists into a new parameters list for code
# 
params2 <- append(params,movement7)
params3 <- append(params2,exog_shock_u5)
params4 <- append(params3,exog_shock_u15)
params5 <- append(params4,exog_shock_o15HR)
params6 <- append(params5,exog_shock_o15LR)


# Define locations and age groups

demog2 <- demog[which(demog$Country == "DRC"|demog$Country == "Uganda"),]
locations <- demog2$location

age_groups <- c("U5", "U15", "O15HR", "O15LR")

# Create dataset of initial values

y_init <- c()
for (c in 1:length(locations)) {
  for (b in 1:length(age_groups)) {
    
    age <- age_groups[b]
    location <- locations[c]
    
    y_init[paste0("S_", age, "_", location)] <- demog[which(demog$location==location), age]
      #params$demog[which(demog$location==location), age]
    y_init[paste0("Ia_", age, "_", location)] <- 0
    y_init[paste0("Ib_", age, "_", location)] <- 0
    y_init[paste0("V_", age, "_", location)] <- 0
    y_init[paste0("R_", age, "_", location)] <- 0
    y_init[paste0("D_", age, "_", location)] <- 0
  }
}

# Create reactions and effects datasets

reactions= list()
nu <- matrix(0, nrow = 0, ncol = length(y_init))
colnames(nu) <- names(y_init)


# Create all reactions and effects that happen within the model

for (c in 1:length(locations)) {
  for (b in 1:length(age_groups)) {
    
    age<- age_groups[b]
    location <- locations[c]
    
    
    # # Births
    # 
    # if (b==1) {
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("S_", age, "_", location, " -> I_", age, "_", location),
    #     paste0("births[births$location== '",location, "',", 3,"]* (S_O15LR","_", location, '+ S_O15HR','_', location, " +V_O15LR","_", location, '+ V_O15HR','_', location, " +R_O15LR","_", location, '+ R_O15HR','_', location,')'))))
    #   
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("S_", age, "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    # }
    # 
    # # Natural deaths- only in adults for now, to balance births
    # # 
    # 
    # if (b==3) {
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("S_", age, "_", location, " -> I_", age, "_", location),
    #     paste0("births[births$location== '",location, "',", 3,"]* (S_O15HR","_", location,')'))))
    #   
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    #   nu <- rbind(nu, new_row)
    #   
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("S_", age, "_", location, " -> I_", age, "_", location),
    #     paste0("births[births$location== '",location, "',", 3,"]* (V_O15HR","_", location,')'))))
    #   
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("V_", age, "_", location), names(y_init))] <- -1
    #   nu <- rbind(nu, new_row)
    #   
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("S_", age, "_", location, " -> I_", age, "_", location),
    #     paste0("births[births$location== '",location, "',", 3,"]* (R_O15HR","_", location,')'))))
    #   
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("R_", age, "_", location), names(y_init))] <- -1
    #   nu <- rbind(nu, new_row)
    #   
    # }
    # 
    # if (b==4) {
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("S_", age, "_", location, " -> I_", age, "_", location),
    #     paste0("births[births$location== '",location, "',", 3,"]* (S_O15LR","_", location,')'))))
    #   
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    #   nu <- rbind(nu, new_row)
    #   
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("S_", age, "_", location, " -> I_", age, "_", location),
    #     paste0("births[births$location== '",location, "',", 3,"]* (V_O15LR","_", location,')'))))
    #   
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("V_", age, "_", location), names(y_init))] <- -1
    #   nu <- rbind(nu, new_row)
    #   
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("S_", age, "_", location, " -> I_", age, "_", location),
    #     paste0("births[births$location== '",location, "',", 3,"]* (R_O15LR","_", location,')'))))
    #   
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("R_", age, "_", location), names(y_init))] <- -1
    #   nu <- rbind(nu, new_row)
    #   
    # }
    
    # Human-to-human transmission: clade 1a
    
    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> I_", age, "_", location),
      paste0("(beta_",age,age_groups[1], "* S_", age, "_", location, " * Ia_", age_groups[1], "_", location, 
             " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))",
             "+(beta_", age, age_groups[2], "* S_", age, "_", location, " * Ia_", age_groups[2], "_", location,  " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))",
             "+ (beta_", age, age_groups[3], "* S_", age, "_", location, " * Ia_", age_groups[3], "_", location,  " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))",
             "+ (beta_", age, age_groups[4], "* S_", age, "_", location, " * Ia_", age_groups[4], "_", location,  " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))"))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    # Human-to-human transmission: clade 1b
    
    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> I_", age, "_", location),
      paste0("(beta_",age,age_groups[1], "* S_", age, "_", location, " * Ib_", age_groups[1], "_", location, 
             " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))",
             "+(beta_", age, age_groups[2], "* S_", age, "_", location, " * Ib_", age_groups[2], "_", location,  " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))",
             "+ (beta_", age, age_groups[3], "* S_", age, "_", location, " * Ib_", age_groups[3], "_", location,  " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))",
             "+ (beta_", age, age_groups[4], "* S_", age, "_", location, " * Ib_", age_groups[4], "_", location,  " / (S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
             "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
             "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
             "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
             "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,"))"))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    # 
    # reactions <- append(reactions, list(c(
    #   #paste0("S_", age, "_", location, " -> I_", age, "_", location),
    #   paste0("0.95*",location, "_exogshock_", age,"* S_", age, "_", location,"/ (S_", age, "_", location, "+Ia_", age, "_", location,
    #          "+Ib_", age, "_", location,"+R_", age, "_", location,
    #          "+V_", age, "_", location,")"))))
    
    
    # Sylvatic transmission, clade 1a
    
    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> I_", age, "_", location),
      paste0("(0.95*",location, "_exogshock_", age,"* S_", age, "_", location,") - (0.95*",location, "_exogshock_", age,"* (V_", age, "_", location,"/ (S_", age, "_", location, "+Ia_", age, "_", location,
                                                                                  "+Ib_", age, "_", location,"+R_", age, "_", location,
                                                                                  "+V_", age, "_", location,")))"))))
  
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    # reactions <- append(reactions, list(c(
    #   #paste0("S_", age, "_", location, " -> I_", age, "_", location),
    #   paste0("0.05*",location, "_exogshock_", age,"* S_", age, "_", location,"/ (S_", age, "_", location, "+Ia_", age, "_", location,
    #          "+Ib_", age, "_", location,"+R_", age, "_", location,
    #          "+V_", age, "_", location,")"))))
    #          
    #          

    # Sylvatic transmission, clade 1b
    
    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> I_", age, "_", location),
      paste0("(0.05*",location, "_exogshock_", age,"* S_", age, "_", location,") - (0.05*",location, "_exogshock_", age,"* (V_", age, "_", location,"/ (S_", age, "_", location, "+Ia_", age, "_", location,
             "+Ib_", age, "_", location,"+R_", age, "_", location,
             "+V_", age, "_", location,")))"))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    # Sexual transmission among high-risk adults
    
    if (b==3) {
      
      reactions <- append(reactions, list(c(
        #paste0("S_", age, "_", location, " -> I_", age, "_", location),
        paste0("beta_sex* S_", age, "_", location, "*Ib_",age,"_",location,"/ (S_", age, "_", location, "+Ia_", age, "_", location,
               "+Ib_", age, "_", location,"+R_", age, "_", location,
               "+V_", age, "_", location,")"))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
    }
    
    # Transmission between different locations
    
    
    for (i in 1:length(locations)) {

      other_country <- locations[i]
      
      reactions <- append(reactions, list(c(
        #paste0("S_", age, "_", location, " -> I_", age, "_", location),
        paste0(location,"_", other_country, "* S_",age,"_", location, "*(Ia_", age_groups[1], "_", other_country,
               "+Ia_", age_groups[2], "_", other_country,"+ Ia_", age_groups[3], "_", other_country,"+ Ia_", age_groups[4], "_", other_country,')'))))
      
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("S_", age, "_", location, " -> I_", age, "_", location),
        paste0(location,"_", other_country, "* S_",age,"_", location, "*(Ib_", age_groups[1], "_", other_country,
               "+Ia_", age_groups[2], "_", other_country,"+ Ib_", age_groups[3], "_", other_country,"+ Ib_", age_groups[4], "_", other_country,')'))))
      
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
    }
    
    
    # Vaccination
    
    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> V", age, "_", location),
      paste0("omega * S_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("V_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    
    # Recovery, clade 1a
    
    reactions <- append(reactions, list(c(
      #paste0("I_", age, "_", location, " -> R", age, "_", location),
      paste0("rho * Ia_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("R_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
  
  # Recovery, clade 1b
    
    reactions <- append(reactions, list(c(
      #paste0("I_", age, "_", location, " -> R", age, "_", location),
      paste0("rho * Ib_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("R_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
  
  # Death, clade 1a
    
    reactions <- append(reactions, list(c(
      #paste0("I_", age, "_", location, " -> D", age, "_", location),
      paste0("mu * Ia_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("D_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
# Death, clade 1b
    
     reactions <- append(reactions, list(c(
      #paste0("I_", age, "_", location, " -> D", age, "_", location),
      paste0("mu * Ib_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("D_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    
# Aging: 
# 
    # if (b == 1) {
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("S_", age, "_", location, " -> S_", age_groups[2], "_", location),
    #     paste0("aging1 * S_", age, "_", location))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("S_", age_groups[2], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("I_", age, "_", location, " -> I_", age_groups[2], "_", location),
    #     paste0("aging1 * Ia_", age, "_", location))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("Ia_", age_groups[2], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("I_", age, "_", location, " -> I_", age_groups[2], "_", location),
    #     paste0("aging1 * Ib_", age, "_", location))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("Ib_", age_groups[2], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("V_", age, "_", location, " -> V_", age_groups[2], "_", location),
    #     paste0("aging1 * V_", age, "_", location))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("V_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("V_", age_groups[2], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("R_", age, "_", location, " -> R_", age_groups[2], "_", location),
    #     paste0("aging1 * R_", age, "_", location))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("R_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("R_", age_groups[2], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   
    # }
    # 
    # 
    # if (b == 2) {
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("S_", age, "_", location, " -> S_", age_groups[3], "_", location),
    #     paste0("aging[", b , "] * S_", age, "_", location,'* percHR'))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("S_", age_groups[3], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("S_", age, "_", location, " -> S_", age_groups[4], "_", location),
    #     paste0("aging[", b, "] * S_", age, "_", location,'* (1-percHR)'))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("S_", age_groups[4], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("I_", age, "_", location, " -> I_", age_groups[3], "_", location),
    #     paste0("aging2 * Ia_", age, "_", location,'* percHR'))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("Ia_", age_groups[3], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("I_", age, "_", location, " -> I_", age_groups[3], "_", location),
    #     paste0("aging2 * Ib_", age, "_", location,'* percHR'))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("Ib_", age_groups[3], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("I_", age, "_", location, " -> I_", age_groups[4], "_", location),
    #     paste0("aging[", b, "] * Ia_", age, "_", location,'* (1-percHR)'))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("Ia_", age_groups[4], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("I_", age, "_", location, " -> I_", age_groups[4], "_", location),
    #     paste0("aging2 * Ib_", age, "_", location,'* (1-percHR)'))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("Ib_", age_groups[4], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("V_", age, "_", location, " -> V_", age_groups[3], "_", location),
    #     paste0("aging2 * V_", age, "_", location,'* percHR'))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("V_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("V_", age_groups[3], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("V_", age, "_", location, " -> V_", age_groups[4], "_", location),
    #     paste0("aging2 * V_", age, "_", location,'* (1-percHR)'))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("V_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("V_", age_groups[4], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("R_", age, "_", location, " -> R_", age_groups[3], "_", location),
    #     paste0("aging2 * R_", age, "_", location,'* (percHR)'))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("R_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("R_", age_groups[3], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    #   
    #   reactions <- append(reactions, list(c(
    #     #paste0("R_", age, "_", location, " -> R_", age_groups[4], "_", location),
    #     paste0("aging2 * R_", age, "_", location,'* (1-percHR)'))))
    #   
    #   new_row <- rep(0, length(y_init))
    #   names(new_row) <- names(y_init)
    #   new_row[match(paste0("R_", age, "_", location), names(y_init))] <- -1
    #   new_row[match(paste0("R_", age_groups[4], "_", location), names(y_init))] <- 1
    #   nu <- rbind(nu, new_row)
    # }
  }
}

reactions2 <- as.character(reactions)
nu2 <- t(nu)

out <- 
  GillespieSSA2::ssa(
    initial_state = y_init,
    reactions = port_reactions(x0= y_init, a = reactions2, nu = nu2),
    params = params6,
    method = ssa_btl(),
    final_time = 365,
    census_interval = .001,
    verbose = FALSE
  )


out2 <- out$state

time <- out$time

out2 <- cbind(out2, time)

  

# Run Gillespie simulation
#Binomial tau-leap method implementation of the SSA as described by Chatterjee et al. (2005). 
#Should be passed as method argument for ssa().
# system.time(out <- ssa(y_init, reactions2, nu2, params, tf = 365, method=ssa.btl(f=10)))
# 
#  # Convert to dataframe
# output_df <- as.data.frame(out$data)
# 
# 
# output_df2 <- output_df %>% 
#   pivot_longer(
#     cols = S_U5_BasUele:D_O15LR_HautLomami,
#     names_to = c("Status", "Age", "Location"), 
#     names_pattern = "(.*)_(.*)_(.*)",
#     values_to = "Count"
#   )
# 
# output_df2_noage <- output_df2 %>%
#   group_by(t, Status, Location) %>%
#   summarise(Count=sum(Count))
# 
# 
# output_df2_noage <- output_df2_noage[which(output_df2_noage$t<365),]
# 
# 
# output_df2_noage_noS <- output_df2_noage[which(output_df2_noage$Status!="S" &output_df2_noage$Status!="R"),]
# 
# 
# # Plot results for selected countries
# ggplot(output_df2_noage_noS) + geom_line(aes(y =Count, x=t, color = Status)) + facet_wrap(~Location, scales = "free")
#   
# 
#   labs(title = 'Stochastic SIR Model with Age and Country Structure', x = 'Time', y = 'Infected Population') +
#   scale_color_manual('', values = c('red', 'blue', 'green'))