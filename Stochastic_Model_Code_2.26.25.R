library(GillespieSSA)
library(ggplot2)
library(adaptivetau)
library(tidyverse)

demog = read.csv("Data/Demographic_data_countries_2_27_25.csv")

demog$location <- demog$Country
demog$location <- ifelse(demog$Country=="DRC",demog$Province, demog$Country)

#movement = read.csv("Data/distance_matrix.csv", row.names = 1)
#
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

percHR <- 0.05

propHR=percHR

demog$O15HR <- demog$O15 * percHR
demog$O15LR <- demog$O15 * (1-percHR)


# Take average of cases from 2021 and 2022 (pre-epidemic) to estimate number of sylvatic infections
# 

demog$avg_sylvatic <- (demog$Cases_2021 + demog$Cases_2022)/2

prop_sylvatic_adult = 0.3 #proportion of sylvatic infections in adults
prop_sylvatic_u15 = 0.45 #0.6 #proportion of sylvatic infections in those 5-15

exog_shock <- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("Location", "exogshock_u5", "exogshock_u15","exogshock_o15HR", "exogshock_o15LR")
colnames(exog_shock) <- x
  
for (i in 1:length(demog$location)) {

  exog_shock[i,1] <- demog[i,17]
  
# this assumes that half of all sylvatic transmission happens in adults, and 1/4 each in children's age groups
  exog_shock[i,2]= (((demog[i,20])*(1-prop_sylvatic_adult-prop_sylvatic_u15))/demog[i,4]) / 365
  exog_shock[i,3]= (((demog[i,20])*prop_sylvatic_u15)/demog[i,5]) / 365
  exog_shock[i,4]= ((((demog[i,20])*prop_sylvatic_adult)/demog[i,6])) / 365
  exog_shock[i,5]= ((((demog[i,20])*prop_sylvatic_adult)/demog[i,6])) / 365

  exog_shock[is.na(exog_shock)] <- 0

}

exog_shock_multi_1a <- 0.95
exog_shock_multi_1b <- 0.05


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


mortality_u5 = 0.10 # mortality rate under 5
mortality_u15 = 0.08 # mortality rate under 15
mortality_o15 = 0.035 # mortality rate over 15


# exograte_u5 <- 
# exograte_u15 <- 0
# exograte_o15 <- 0
# exograte_o15HR <- 0
# exograte_o15LR <- 0

params <- list(
  #beta = 0.3,       # Transmission rate
  rho = 1/21,      # Recovery rate
  mu = 1/50,        # Birth and natural death rate (assuming lifespan of 50 years)
  aging = c(0,0),#c(1/5,1/10),#c(1/20, 1/30),  # Aging rates from young to adult and adult to senior
  demog= demog,
  movement = movement,
  exog_shock = exog_shock,
  percHR = percHR,
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
  beta_O15LRO15LR = SA_o15*contacts_o15*contactprop_o15o15, #beta for transmission between o15 and u15
  beta_O15LRO15HR = SA_o15*contacts_o15*contactprop_o15o15, #beta for transmission between o15 and u15
  
  beta_O15HRU5= SA_o15*contacts_o15*contactprop_o15u5, #beta for transmission between o15 and u5
  beta_O15HRU15= SA_o15*contacts_o15*contactprop_o15u15, #beta for transmission between o15 and u15
  beta_O15HRO15LR = SA_o15*contacts_o15*contactprop_o15o15*.95, #beta for transmission between o15 and u15
  beta_O15HRO15HR = SA_o15*contacts_o15*contactprop_o15o15*.05 + SA_sex * contacts_sex, #beta for transmission between o15 and u15

  beta_sex = 2.4*(recoveryrate),
  
  # (u5_pop/(population)*(beta_u5u5+beta_u5u15+beta_u5o15)/3 + u15_pop/(population)*(beta_u15u5+beta_u15u15+beta_u15o15)/3 +
  #  o15_pop/(population)*(beta_o15u5+beta_o15u15+beta_o15o15)/3)*21 # Estimate of R0
  #  
  rho= recoveryrate, #recovery rate
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
  omega=0
)

# Define parameters

locations <- demog$location
locations <- locations[1:4]

age_groups <- c("U5", "U15", "O15HR", "O15LR")

y_init <- c()
for (c in 1:length(locations)) {
  for (b in 1:length(age_groups)) {
    
    age <- age_groups[b]
    location <- locations[c]
    
    y_init[paste0("S_", age, "_", location)] <- params$demog[which(demog$location==location), age]
    y_init[paste0("Ia_", age, "_", location)] <- 0
    y_init[paste0("Ib_", age, "_", location)] <- 0
    y_init[paste0("V_", age, "_", location)] <- 0
    y_init[paste0("R_", age, "_", location)] <- 0
    y_init[paste0("D_", age, "_", location)] <- 0
  }
}


reactions= list()
nu <- matrix(0, nrow = 0, ncol = length(y_init))
colnames(nu) <- names(y_init)


for (c in 1:length(locations)) {
  for (b in 1:length(age_groups)) {
    
    age<- age_groups[b]
    location <- locations[c]
  

    N <- paste0("(S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
                "+Ia_", age_groups[1], "_", location, "+ Ia_", age_groups[2], "_", location, "+Ia_", age_groups[3], "_", location, "+Ia_", age_groups[4], "_", location,
                "+Ib_", age_groups[1], "_", location, "+ Ib_", age_groups[2], "_", location, "+Ib_", age_groups[3], "_", location, "+Ib_", age_groups[4], "_", location,
                "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+V_", age_groups[4], "_", location,
                "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+R_", age_groups[4], "_", location,")")

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

    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> I_", age, "_", location),
                                          paste0("exog_shock_multi_1a* exog_shock[exog_shock$Location== '",location, "',", b+1,"]* S_", age, "_", location, "-exog_shock[exog_shock$Location== '",location, "',", b+1,"] / sum(S_", age, "_", location, "+Ia_", age, "_", location,
                                                 "+Ib_", age, "_", location,"+R_", age, "_", location,
                                                 "+V_", age, "_", location,")"))))
                                          
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> I_", age, "_", location),
      paste0("exog_shock_multi_1b* exog_shock[exog_shock$Location== '",location, "',", b+1,"]* S_", age, "_", location, "-exog_shock[exog_shock$Location== '",location, "',", b+1,"] / sum(S_", age, "_", location, "+Ia_", age, "_", location,
             "+Ib_", age, "_", location,"+R_", age, "_", location,
             "+V_", age, "_", location,")"))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    if (b==3) {
    
    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> I_", age, "_", location),
      paste0("beta_sex* S_", age, "_", location, "*Ib_",age,"_",location,"/ sum(S_", age, "_", location, "+Ia_", age, "_", location,
             "+Ib_", age, "_", location,"+R_", age, "_", location,
             "+V_", age, "_", location,")"))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    }
    
    for (i in 1:length(locations)) {

      other_country <- locations[i]

      reactions <- append(reactions, list(c(
        #paste0("S_", age, "_", location, " -> I_", age, "_", location),
                                            paste0("movement['", location,"','", other_country, "'] * S_",age,"_", location, "*Ia_", age, "_", other_country))))
      
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("S_", age, "_", location, " -> I_", age, "_", location),
        paste0("movement['", location,"','", other_country, "'] * S_",age,"_", location, "*Ib_", age, "_", other_country))))
      
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
    }


    reactions <- append(reactions, list(c(
      #paste0("S_", age, "_", location, " -> V", age, "_", location),
                                          paste0("omega * S_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("V_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)

    reactions <- append(reactions, list(c(
      #paste0("I_", age, "_", location, " -> R", age, "_", location),
                                          paste0("rho * Ia_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("R_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    reactions <- append(reactions, list(c(
      #paste0("I_", age, "_", location, " -> R", age, "_", location),
      paste0("rho * Ib_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("R_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)

    reactions <- append(reactions, list(c(
      #paste0("I_", age, "_", location, " -> D", age, "_", location),
                                          paste0("mu * Ia_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("D_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)
    
    reactions <- append(reactions, list(c(
      #paste0("I_", age, "_", location, " -> D", age, "_", location),
      paste0("mu * Ib_", age, "_", location))))
    
    new_row <- rep(0, length(y_init))
    names(new_row) <- names(y_init)
    new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
    new_row[match(paste0("D_", age, "_", location), names(y_init))] <- 1
    nu <- rbind(nu, new_row)


    if (b == 1) {

      reactions <- append(reactions, list(c(
        #paste0("S_", age, "_", location, " -> S_", age_groups[2], "_", location),
                                            paste0("aging[", b, "] * S_", age, "_", location))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("S_", age_groups[2], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      
      reactions <- append(reactions, list(c(
        #paste0("I_", age, "_", location, " -> I_", age_groups[2], "_", location),
                                            paste0("aging[", b, "] * Ia_", age, "_", location))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ia_", age_groups[2], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("I_", age, "_", location, " -> I_", age_groups[2], "_", location),
        paste0("aging[", b, "] * Ib_", age, "_", location))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ib_", age_groups[2], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("V_", age, "_", location, " -> V_", age_groups[2], "_", location),
                                            paste0("aging[", b, "] * V_", age, "_", location))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("V_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("V_", age_groups[2], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("R_", age, "_", location, " -> R_", age_groups[2], "_", location),
                                            paste0("aging[", b, "] * R_", age, "_", location))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("R_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("R_", age_groups[2], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      
    }


    if (b == 2) {

      reactions <- append(reactions, list(c(
        #paste0("S_", age, "_", location, " -> S_", age_groups[3], "_", location),
                                            paste0("aging[", b , "] * S_", age, "_", location,'* percHR'))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("S_", age_groups[3], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("S_", age, "_", location, " -> S_", age_groups[4], "_", location),
                                            paste0("aging[", b, "] * S_", age, "_", location,'* (1-percHR)'))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("S_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("S_", age_groups[4], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("I_", age, "_", location, " -> I_", age_groups[3], "_", location),
                                            paste0("aging[", b, "] * Ia_", age, "_", location,'* percHR'))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ia_", age_groups[3], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("I_", age, "_", location, " -> I_", age_groups[3], "_", location),
        paste0("aging[", b, "] * Ib_", age, "_", location,'* percHR'))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ib_", age_groups[3], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("I_", age, "_", location, " -> I_", age_groups[4], "_", location),
                                            paste0("aging[", b, "] * Ia_", age, "_", location,'* (1-percHR)'))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("Ia_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ia_", age_groups[4], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("I_", age, "_", location, " -> I_", age_groups[4], "_", location),
        paste0("aging[", b, "] * Ib_", age, "_", location,'* (1-percHR)'))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("Ib_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("Ib_", age_groups[4], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      
      
      reactions <- append(reactions, list(c(
        #paste0("V_", age, "_", location, " -> V_", age_groups[3], "_", location),
                                            paste0("aging[", b, "] * V_", age, "_", location,'* percHR'))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("V_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("V_", age_groups[3], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("V_", age, "_", location, " -> V_", age_groups[4], "_", location),
                                            paste0("aging[", b, "] * V_", age, "_", location,'* (1-percHR)'))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("V_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("V_", age_groups[4], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("R_", age, "_", location, " -> R_", age_groups[3], "_", location),
                                            paste0("aging[", b, "] * R_", age, "_", location,'* (percHR)'))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("R_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("R_", age_groups[3], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
      
      reactions <- append(reactions, list(c(
        #paste0("R_", age, "_", location, " -> R_", age_groups[4], "_", location),
                                            paste0("aging[", b, "] * R_", age, "_", location,'* (1-percHR)'))))
      
      new_row <- rep(0, length(y_init))
      names(new_row) <- names(y_init)
      new_row[match(paste0("R_", age, "_", location), names(y_init))] <- -1
      new_row[match(paste0("R_", age_groups[4], "_", location), names(y_init))] <- 1
      nu <- rbind(nu, new_row)
    }
  }
}

reactions2 <- as.character(reactions)
nu2 <- t(nu)


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