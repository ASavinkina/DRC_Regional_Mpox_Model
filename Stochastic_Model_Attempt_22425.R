#################
#################
#################

# Code for modelling estimated vaccine needs for mpox in central Africa, stochastically
# Allowing for aging between age groups and movement between countries
# 
# Alexandra Savinkina
# 2.24.25
# 

# 
## Import libraries

library(adaptivetau)
library(ggplot2)
library(tidyverse)
library(scales)
library(gridExtra)
library(DT)

# read in country demographics data
# 

demographics <- read.csv("Data/Demographic_data_countries_2_24_25.csv")

# For locations, DRC locations and other countries will each be treated as a location
demographics$location <- ifelse(demographics$Country=="DRC",demographics$location, demographics$Country)

# read in distance matrix
# 
dist_mat <- read.csv("Data/distance_matrix.csv")

location_list <- demographics$location

location_list[37] <- "All"

# for the model, general inputs
# 

runs=10

set.seed(1)

# create datasets to capture data in model
# 
# dataset for all infected over time
location_all_infected_bytime <- list()
# dataset for total infections over time
location_totalinfections_bytime <- list()
#dataset for final infection statistics
total_infections_location <- data.frame(matrix(0, nrow=0, ncol=8))
colnames(total_infections_location) <- c("location", "n_min","Q1", "Q95l","median", "Q3","Q95u","max")
#dataset for final deaths statistics
total_deaths_location <- data.frame(matrix(0, nrow=0, ncol=8))
colnames(total_deaths_location) <- c("location", "n_min","Q1", "Q95l","median", "Q3","Q95u","max")
#dataset for final infections statistics, bye age (medial only)
total_infections_age_location <- data.frame(matrix(0, nrow=0, ncol=6))
colnames(total_infections_age_location) <- c("location", "total_cases_u5_median","total_cases_u15_median",  "total_cases_o15LR_median","total_cases_o15HR_median", "total_cases_o15_median")
#dataset for final deaths statistics, bye age (medial only)
total_deaths_age_location <- data.frame(matrix(0, nrow=0, ncol=5))
colnames(total_deaths_age_location) <- c("location", "D_u5_median","D_u15_median",  "D_o15LR_median",
                                         "D_o15HR_median")  


#dataset for total vaccinated and vaccine doses
number_vaccinated_location <- data.frame(matrix(0, nrow=0, ncol=3))
colnames(number_vaccinated_location) <- c("location", "n_vaccinated", "n_doses")

# Set initial conditions for model

initial_inf_u5 = 0
initial_inf_u15 = 0
initial_inf_o15LR = 0
initial_inf_o15HR = 0

# Set vaccination strategy
# 
vax_effect <- 0.75
vax_S_u5 <- vax_effect*0.8
vax_S_u15 <- vax_effect*0
vax_S_o15LR <- vax_effect*0
vax_S_o15HR <- vax_effect*0

#proportion of infections that are sylvatic

prop_sylvatic = 0.7


for (i in 1:nrow(demographics)) {
  
  location = demographics[i,1]
  
  # this assumes that half of all sylvatic transmission happens in adults, and 1/4 each in children's age groups
  exogshock_u5= (demographics[i,4]*prop_sylvatic)*.2
  exogshock_u15= (demographics[i,5]*prop_sylvatic)*.6
  exogshock_o15= (demographics[i,6]*prop_sylvatic)*.2
  
  population = as.numeric(demographics[i,2]) + as.numeric(demographics[i,3]) + as.numeric(demographics[i,4]) #population of region
  exograte_u5 = exogshock_u5/365  # rate of exogenous shocks, under 5
  exograte_u15 = exogshock_u15/365  # rate of exogenous shocks, under 15
  exograte_o15 = exogshock_o15/365  # rate of exogenous shocks, over 15
  recoveryrate = 1/21 # recovery rate
  
  u5_pop = as.numeric(demographics[i,2]) # proportion of pop under 5
  o15_pop = as.numeric(demographics[i,3]) # proportion of pop over  15
  u15_pop = as.numeric(demographics[i,4]) # proportion of pop over 5 and under 15
  propHR = as.numeric(demographics[i,5])
  
  contacts_u5 = 5.44 # number of daily physical contacts, under 5
  contacts_u15 = 6.88 # number of daily physical contacts, under 15
  contacts_o15 = 4.788 # number of daily physical contacts, over 15
  contacts_sex = 10/30 # sexual contacts for high risk, daily (average 10 a month- made up)
  
  contactprop_u5u5 = 0.235 # proportion of physical contacts for under 5 with under 5
  contactprop_u5u15 = 0.235 # proportion of physical contacts for under 5 with under 15
  contactprop_u5o15 = 1-contactprop_u5u5-contactprop_u5u15 # proportion of physical contacts for under 5 with over 15
  
  contactprop_u15u5 = 0.235# proportion of physical contacts for under 15 with under 5
  contactprop_u15u15 = 0.235 # proportion of physical contacts for under 15 with under 15
  contactprop_u15o15 = 1-contactprop_u15u5-contactprop_u15u15 # proportion of physical contacts for under 15 with over 15
  
  contactprop_o15u5 = 0.135 # proportion of physical contacts for over 15 with under 5
  contactprop_o15u15 = 0.135 # proportion of physical contacts for over 15 with under 15
  contactprop_o15o15 = 1-contactprop_o15u5-contactprop_o15u15 # proportion of physical contacts for over 15 with over 15
  
  SA_u5 = 0.034#0.04#0.11/3#0.1168/3 #0.1168/4 # secondary attack rate when infected person is under 5, halved to account for data coming from household transmission
  SA_u15 = 0.008#0.011#0.0133#0.0803/7 #0.0803/4 # secondary attack rate when infected person is under 15, halved to account for data coming from household transmission
  SA_o15 = 0.005#0.011#0.0133 #0.0133/3 # secondary attack rate when infected person is over 15, halved to account for data coming from household transmission
  SA_sex = 0.15 #0.15 # secondary attack rate from sex, made up for now
  
  
  mortality_u5 = 0.10 # mortality rate under 5
  mortality_u15 = 0.08 # mortality rate under 15
  mortality_o15 = 0.035 # mortality rate over 15
  
  
  init.values = c(
    S_u5 = as.integer(u5_pop-initial_inf_u5-u5_pop*vax_S_u5),  
    I_u5 = as.integer(initial_inf_u5), 
    V_u5=as.integer(u5_pop*vax_S_u5), 
    D_u5=0, 
    R_u5=0,
    S_u15 = as.integer(u15_pop-initial_inf_u15-u15_pop*vax_S_u15),
    I_u15 = as.integer(initial_inf_u15),
    V_u15=as.integer(u15_pop*vax_S_u15),
    D_u15=0,
    R_u15=0,
    S_o15LR = as.integer(o15_pop-initial_inf_o15LR-o15_pop*vax_S_o15LR), 
    I_o15LR = as.integer(initial_inf_o15LR), 
    V_o15LR= as.integer(o15_pop*vax_S_o15LR), 
    D_o15LR=0, 
    R_o15LR=0,
    S_o15HR = as.integer(o15_pop*propHR-initial_inf_o15HR-(o15_pop*propHR*vax_S_o15HR+o15_pop*propHR*vax_S_o15LR)), 
    I_o15HR = as.integer(initial_inf_o15HR), 
    V_o15HR=as.integer(o15_pop*propHR*vax_S_o15HR+o15_pop*propHR*vax_S_o15LR), 
    D_o15HR=0 , 
    R_o15HR=0
  )
  
  
  # Specify all transitions
  
  transitions = list(
    # under 5
    c(S_u5 = -1, I_u5 = +1), # movement from susceptible to infected, endogenous infection
    c(S_u5 = -1, I_u5 = +1), # movement from susceptible to infected, exogenous infection
    c(S_u5= -1, I_u5 = +1), # movement from susceptible to infected, from infection from other province
    c(S_u5= -1, V_u5= +1), #movement from susceptible to vaccinated
    c(I_u5 = -1, D_u5 = +1), #movement from infected to dead
    c(I_u5 = -1, R_u5 = +1), #movement from infected to recovered
    
    # u5 aging
    # 
    
    c(S_u5= -1, S_u15 = +1), #aging among susceptibles
    c(I_u5= -1, I_u15 = +1), #aging among susceptibles
    c(V_u5= -1, V_u15 = +1), #aging among susceptibles
    c(R_u5= -1, R_u15 = +1), #aging among susceptibles
    
    # 5-15
    c(S_u15 = -1, I_u15 = +1), # movement from susceptible to infected, endogenous infection
    c(S_u15 = -1, I_u15 = +1), # movement from susceptible to infected, exogenous infection 
    c(S_u15= -1, I_u15 = +1), # movement from susceptible to infected, from infection from other province
    c(S_u15= -1, V_u15= +1), #movement from susceptible to vaccinated
    c(I_u15 = -1, D_u15 = +1), #movement from infected to dead
    c(I_u15 = -1, R_u15 = +1), #movement from infected to recovered
    
    # u15 aging, high risk and low risk
    # 
    
    c(S_u15= -1, S_o15HR = +1), #aging among susceptibles
    c(I_u15= -1, I_o15HR = +1), #aging among infected
    c(V_u15= -1, V_o15HR = +1), #aging among vaccinated
    c(R_u15= -1, R_o15HR = +1), #aging among recovered
    
    c(S_u15= -1, S_o15LR = +1), #aging among susceptibles
    c(I_u15= -1, I_o15LR = +1), #aging among infected
    c(V_u15= -1, V_o15LR = +1), #aging among vaccinated
    c(R_u15= -1, R_o15LR = +1), #aging among recovered
    
    
    # 15+, LR
    c(S_o15LR = -1, I_o15LR = +1), # movement from susceptible to infected, endogenous infection
    c(S_o15LR = -1, I_o15LR = +1), # movement from susceptible to infected, exogenous infection
    c(S_o15LR= -1, I_o15LR = +1), # movement from susceptible to infected, from infection from other province
    c(S_o15LR= -1, V_o15LR= +1), #movement from susceptible to vaccinated
    c(I_o15LR = -1, D_o15LR = +1), #movement from infected to dead
    c(I_o15LR = -1, R_o15LR = +1), #movement from infected to recovered
    
    # 15+, HR
    c(S_o15HR = -1, I_o15HR = +1), # movement from susceptible to infected, endogenous infection
    c(S_o15HR = -1, I_o15HR = +1), # movement from susceptible to infected, exogenous infection 
    c(S_o15HR= -1, I_o15HR = +1), # movement from susceptible to infected, from infection from other province
    c(S_o15HR= -1, V_o15HR= +1), #movement from susceptible to vaccinated
    c(I_o15HR = -1, D_o15HR = +1), #movement from infected to dead
    c(I_o15HR = -1, R_o15HR = +1) #movement from infected to recovered
    
  )
  
  # Rates for all transitions
  # (In same order as "transitions")
  RateF <- function(x, pars, times) {
    
    #under 5
    return(c(
      pars$beta_u5u5*x["S_u5"]*x["I_u5"]/(x["S_u5"] + x["I_u5"] + x["R_u5"] + x["V_u5"] + x["S_u15"] + x["I_u15"] + x["R_u15"] + x["V_u15"]+
                                            x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"]+
                                            x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"])+
        pars$beta_u5u15*x["S_u5"]*x["I_u15"]/(x["S_u5"]  + x["I_u5"]  + x["R_u5"] + x["V_u5"]+ x["S_u15"]  + x["I_u15"] + x["R_u15"] + x["V_u15"]+
                                                x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"]+
                                                x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"])+
        pars$beta_u5o15*x["S_u5"]*x["I_o15LR"]/(x["S_u5"]  + x["I_u5"]  + x["R_u5"] + x["V_u5"]+ x["S_u15"]  + x["I_u15"] + x["R_u15"] +  x["V_u15"] +
                                                  x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"] +
                                                  x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"]) +
        pars$beta_u5o15*x["S_u5"]*x["I_o15HR"]/(x["S_u5"]  + x["I_u5"]  + x["R_u5"] + x["V_u5"]+ x["S_u15"]  + x["I_u15"] + x["R_u15"] +  x["V_u15"] +
                                                  x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"] +
                                                  x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"]), # movement from susceptible to pre-symptomatic, due to endogenous infection
      #for exogenous shocks- we assume a standard amount annually, broken up into daily rates
      #this is assumed to happen by age, an homogenously across vaccination status, though only those
      #unvaccinated will be infected
      pars$theta_u5-pars$theta_u5*x["V_u5"]/(x["S_u5"]  + x["I_u5"] + x["R_u5"] + x["V_u5"]), # movement from susceptible to infected, exogenous infection,
      
  
      

      
      pars$omega_u5*x["S_u5"], #movement from susceptible to vaccinated
      pars$mu_u5*x["I_u5"], #movement from infected to dead
      pars$rho*x["I_u5"], #movement from infected to recovered
      
      
      # 5- 15
      pars$beta_u15u5*x["S_u15"]*x["I_u5"]/(x["S_u5"] + x["I_u5"] + x["R_u5"] + x["V_u5"] + x["S_u15"] + x["I_u15"] + x["R_u15"] + x["V_u15"]+
                                              x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"]+
                                              x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"])+
        pars$beta_u15u15*x["S_u15"]*x["I_u15"]/(x["S_u5"]  + x["I_u5"]  + x["R_u5"] + x["V_u5"]+ x["S_u15"]  + x["I_u15"] + x["R_u15"] + x["V_u15"]+
                                                  x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"]+
                                                  x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"])+
        pars$beta_u15o15*x["S_u15"]*x["I_o15LR"]/(x["S_u5"]  + x["I_u5"]  + x["R_u5"] + x["V_u5"]+ x["S_u15"]  + x["I_u15"] + x["R_u15"] +  x["V_u15"] +
                                                    x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"] +
                                                    x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"]) +
        pars$beta_u15o15*x["S_u15"]*x["I_o15HR"]/(x["S_u5"]  + x["I_u5"]  + x["R_u5"] + x["V_u5"]+ x["S_u15"]  + x["I_u15"] + x["R_u15"] +  x["V_u15"] +
                                                    x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"] +
                                                    x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"]), # movement from susceptible to pre-symptomatic, due to endogenous infection
      pars$theta_u15-pars$theta_u15*x["V_u15"]/(x["S_u15"]  + x["I_u15"] + x["R_u15"] + x["V_u15"]), # movement from susceptible to infected, exogenous infection
      pars$omega_u15*x["S_u15"], #movement from susceptible to vaccinated
      pars$mu_u15*x["I_u15"], #movement from infected to dead
      pars$rho*x["I_u15"], #movement from infected to recovered
      
      # 15+ LR
      pars$beta_o15u5*x["S_o15LR"]*x["I_u5"]/(x["S_u5"] + x["I_u5"] + x["R_u5"] + x["V_u5"] + x["S_u15"] + x["I_u15"] + x["R_u15"] + x["V_u15"]+
                                                x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"]+
                                                x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"])+
        pars$beta_o15u15*x["S_o15LR"]*x["I_u15"]/(x["S_u5"]  + x["I_u5"]  + x["R_u5"] + x["V_u5"]+ x["S_u15"]  + x["I_u15"] + x["R_u15"] + x["V_u15"]+
                                                    x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"]+
                                                    x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"])+
        pars$beta_o15o15*x["S_o15LR"]*x["I_o15LR"]/(x["S_u5"]  + x["I_u5"]  + x["R_u5"] + x["V_u5"]+ x["S_u15"]  + x["I_u15"] + x["R_u15"] +  x["V_u15"] +
                                                      x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"] +
                                                      x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"]) +
        pars$beta_o15o15*x["S_o15LR"]*x["I_o15HR"]/(x["S_u5"]  + x["I_u5"]  + x["R_u5"] + x["V_u5"]+ x["S_u15"]  + x["I_u15"] + x["R_u15"] +  x["V_u15"] +
                                                      x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"] +
                                                      x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"]), # movement from susceptible to pre-symptomatic, due to endogenous infection
      pars$theta_o15LR-pars$theta_o15LR*x["V_o15LR"]/(x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"]), # movement from susceptible to infected, exogenous infection
      pars$omega_o15*x["S_o15LR"], #movement from susceptible to vaccinated
      pars$mu_o15*x["I_o15LR"], #movement from infected to dead
      pars$rho*x["I_o15LR"], #movement from infected to recovered
      
      # 15+ HR
      pars$beta_o15u5*x["S_o15HR"]*x["I_u5"]/(x["S_u5"] + x["I_u5"] + x["R_u5"] + x["V_u5"] + x["S_u15"] + x["I_u15"] + x["R_u15"] + x["V_u15"]+
                                                x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"]+
                                                x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"])+
        pars$beta_o15u15*x["S_o15HR"]*x["I_u15"]/(x["S_u5"]  + x["I_u5"]  + x["R_u5"] + x["V_u5"]+ x["S_u15"]  + x["I_u15"] + x["R_u15"] + x["V_u15"]+
                                                    x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"]+
                                                    x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"])+
        pars$beta_o15o15*x["S_o15HR"]*x["I_o15LR"]/(x["S_u5"]  + x["I_u5"]  + x["R_u5"] + x["V_u5"]+ x["S_u15"]  + x["I_u15"] + x["R_u15"] +  x["V_u15"] +
                                                      x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"] +
                                                      x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"]) +
        pars$beta_o15o15*x["S_o15HR"]*x["I_o15HR"]/(x["S_u5"]  + x["I_u5"]  + x["R_u5"] + x["V_u5"]+ x["S_u15"]  + x["I_u15"] + x["R_u15"] +  x["V_u15"] +
                                                      x["S_o15LR"]  + x["I_o15LR"] + x["R_o15LR"] + x["V_o15LR"] +
                                                      x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"])+
        pars$beta_o15o15HR*x["S_o15HR"]*x["I_o15HR"]/ (x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"]), # movement from susceptible to pre-symptomatic, due to endogenous infection
      pars$theta_o15HR - pars$theta_o15HR*x["V_o15HR"]/(x["S_o15HR"]  + x["I_o15HR"] + x["R_o15HR"] + x["V_o15HR"]), # movement from susceptible to infected, exogenous infection
      pars$omega_o15HR*x["S_o15HR"] +  pars$omega_o15*x["S_o15HR"], #movement from susceptible to vaccinated
      pars$mu_o15*x["I_o15HR"], #movement from infected to dead
      pars$rho*x["I_o15HR"] #movement from infected to recovered
    ))
  }
  
  
  # Setting parameters
  pars = list(
    beta_u5u5= SA_u5*contacts_u5*contactprop_u5u5, #beta for transmission between u5 and u5
    beta_u5u15= SA_u5*contacts_u5*contactprop_u5u15, #beta for transmission between u5 and u15
    beta_u5o15= SA_u5*contacts_u5*contactprop_u5o15,#beta for transmission between u5 and o15
    
    beta_u15u5= SA_u15*contacts_u15*contactprop_u15u5, #beta for transmission between u15 and u5
    beta_u15u15= SA_u15*contacts_u15*contactprop_u15u15, #beta for transmission between u15 and u15
    beta_u15o15= SA_u15*contacts_u15*contactprop_u15o15, #beta for transmission between u15 and o15
    
    beta_o15u5= SA_o15*contacts_o15*contactprop_o15u5, #beta for transmission between o15 and u5
    beta_o15u15= SA_o15*contacts_o15*contactprop_o15u15, #beta for transmission between o15 and u15
    beta_o15o15= SA_o15*contacts_o15*contactprop_o15o15, #beta for transmission between o15 and o15
    
    beta_o15o15HR = SA_sex * contacts_sex,
    
    # (u5_pop/(population)*(beta_u5u5+beta_u5u15+beta_u5o15)/3 + u15_pop/(population)*(beta_u15u5+beta_u15u15+beta_u15o15)/3 +
    #  o15_pop/(population)*(beta_o15u5+beta_o15u15+beta_o15o15)/3)*21 # Estimate of R0
    #  
    rho= recoveryrate, #recovery rate
    mu_u5 = mortality_u5/21, #death rate, under 5
    mu_u15 = mortality_u15/21, #death rate, under 15
    mu_o15 = mortality_o15/21, #death rate, over 15
    theta_u5= exograte_u5, #rate of exogenous shocks
    theta_u15= exograte_u15, #rate of exogenous shocks
    theta_o15HR= exograte_o15*propHR, #rate of exogenous shocks
    theta_o15LR = exograte_o15*(1-propHR),
    omega_u5= 0, #vaccination rate, under 5
    omega_u15= 0, #vaccination rate, under 15
    omega_o15= 0, #vaccination rate, over 15 low risk
    omega_o15HR= 0  #vaccination rate, over 15 high risk 
    
  )
  
  
  #  Create dataset
  
  results_all_susceptible <- data.frame(matrix(0, nrow=0, ncol=6))
  colnames(results_all_susceptible) <- c("time", "S_u5","S_u15", "S_o15LR","S_o15HR", "run")
  
  results_all_vaccinated <- data.frame(matrix(0, nrow=0, ncol=6))
  colnames(results_all_vaccinated) <- c("time", "V_u5","V_u15", "V_o15LR", "V_o15HR", "run")
  
  results_all_infected <- data.frame(matrix(0, nrow=0, ncol=6))
  colnames(results_all_infected) <- c("time", "I_u5","I_u15", "I_o15LR","I_o15HR","run")
  
  results_all_recovered <- data.frame(matrix(0, nrow=0, ncol=6))
  colnames(results_all_recovered) <- c("time", "R_u5","R_u15", "R_o15LR", "R_o15HR","run")
  
  results_all_dead<- data.frame(matrix(0, nrow=0, ncol=6))
  colnames(results_all_dead) <- c("time", "D_u5", "D_u15", "D_o15LR","D_o15HR", "run")
  
  
  
  