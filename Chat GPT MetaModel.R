library(GillespieSSA)
library(ggplot2)
library(adaptivetau)

demog = read.csv("Data/Demographic_data_countries_2_24_25.csv")

demog$location <- demog$Country
demog$location <- ifelse(demog$Country=="DRC",demog$Province, demog$Country)

percHR <- 0.05

demog$O15HR <- demog$O15 * percHR
demog$O15LR <- demog$O15 * (1-percHR)

movement = read.csv("Data/distance_matrix.csv", row.names = 1)

rownames(movement)[rownames(movement) == "South Sudan"] = "SouthSudan"
colnames(movement)[colnames(movement) == "South.Sudan"] = "SouthSudan"

rownames(movement)[rownames(movement) == "Central African Republic"] = "CAR"
colnames(movement)[colnames(movement) == "Central.African.Republic"] = "CAR"

rownames(movement)[rownames(movement) == "South Sudan"] = "SouthSudan"
colnames(movement)[colnames(movement) == "South_Sudan"] = "SouthSudan"

rownames(movement)[rownames(movement) == "Équateur"] = "Equateur"
colnames(movement)[colnames(movement) == "Équateur"] = "Equateur"

names(movement) <- gsub(x = names(movement), pattern = "\\.", replacement = "_")

rownames(movement) <- gsub(x = rownames(movement), pattern = "\\-", replacement = "_")

demog$location <-gsub(x = demog$location, pattern = "\\-", replacement = "_")

params <- list(
  beta = 0.3,       # Transmission rate
  rho = 0.1,      # Recovery rate
  mu = 1/50,        # Birth and natural death rate (assuming lifespan of 50 years)
  aging = c(1/20, 1/30),  # Aging rates from young to adult and adult to senior
  demog= demog,
  movement = movement,
  percHR = percHR
)

# Define parameters

locations <- demog$location
age_groups <- c("U5", "U15", "O15HR", "O15LR")

y_init <- c()
for (c in 1:length(locations)) {
  for (b in 1:length(age_groups)) {

    age <- age_groups[b]
    location <- locations[c]
    
    y_init[paste0("S_", age, "_", location)] <- params$demog[which(demog$location==location), age]
    y_init[paste0("I_", age, "_", location)] <- 0
    y_init[paste0("V_", age, "_", location)] <- 0
    y_init[paste0("R_", age, "_", location)] <- 0
  }
}

transition = list()
    
for (c in 1:length(locations)) {
  for (b in 1:length(age_groups)) {
      
    
      age <- age_groups[b]
      location <- locations[c]
      
      c((as.character(paste0("S_", age, "_", location))=-1),paste0("I_", age, "_", location)=1))
      c(as.character(paste0("S_", age, "_", location))=1, Y=2)
      
      c(y_init[1]=1)
      
      transition= append(c((paste0("S_", age, "_", location)=1), (paste0("I_", age, "_", location)=1)))
      
    }}

      
      transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> I_", age, "_", location))))

      for (i in 1:length(locations)) {
        
        other_country <- locations[i]
        
      transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> I_", age, "_", location))))
      
      transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> V", age, "_", location))))
      
      transition <- append(transition, list(c(paste0("I_", age, "_", location, " -> R", age, "_", location))))
      
      transition <- append(transition, list(c(paste0("I_", age, "_", location, " -> D", age, "_", location))))
      
      }
      
      if (b == 1) {
        
        transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> S_", age_groups[2], "_", location))))
        
        transition <- append(transition, list(c(paste0("I_", age, "_", location, " -> I_", age_groups[2], "_", location))))
        
         
        transition <- append(transition, list(c(paste0("R_", age, "_", location, " -> R_", age_groups[2], "_", location))))
        
        
      }
      
      
      if (b == 2) {
        
        transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> S_", age_groups[3], "_", location))))
        
      
        transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> S_", age_groups[4], "_", location))))
        
            
        transition <- append(transition, list(c(paste0("I_", age, "_", location, " -> I_", age_groups[3], "_", location))))
        
      
        
        transition <- append(transition, list(c(paste0("I_", age, "_", location, " -> I_", age_groups[4], "_", location))))
        
        
        transition <- append(transition, list(c(paste0("R_", age, "_", location, " -> R_", age_groups[3], "_", location))))
        
        transition <- append(transition, list(c(paste0("R_", age, "_", location, " -> R_", age_groups[4], "_", location))))
        
      
    }
  }
  }  
) 





# Define reaction set dynamically
transition <- list()
rate <- list()

for (c in 1:length(locations)) {
  for (b in 1:length(age_groups)) {
    
    age <- age_groups[b]
    location <- locations[c]
    
    
    N <- paste0("S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
                "+I_", age_groups[1], "_", location, "+ I_", age_groups[2], "_", location, "+I_", age_groups[3], "_", location, "+S_", age_groups[4], "_", location,
                "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+S_", age_groups[4], "_", location,
                "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+S_", age_groups[4], "_", location)
    
    transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> I_", age, "_", location))))
                                          
    rate <-  append(rate, list(c(paste0("(beta",age,age_groups[1], "* S_", age, "_", location, " * I_", age_groups[1], "_", location,  " / N)",
                    "+ (beta", age, age_groups[2], "* S_", age, "_", location, " * I_", age_groups[2], "_", location,  " / N)",
                     "+ (beta", age, age_groups[3], "* S_", age, "_", location, " * I_", age_groups[3], "_", location,  " / N)",
                      "+ (beta", age, age_groups[4], "* S_", age, "_", location, " * I_", age_groups[4], "_", location,  " / N)"))))
    
    transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> I_", age, "_", location))))
    
    rate <- append(rate, list(c(paste0("theta_", age,"_",location, "* S_", age, "_", location, "-theta_", age,"_",location," * V_", age, "_", location, 
                                          " / sum(S_", age, "_", location, "+I", age, "_", location,"+R", age, "_", location))))
    
    for (i in 1:length(locations)) {
      
    other_country <- locations[i]

    transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> I_", age, "_", location))))
    
    rate <- append(rate, list(c(paste0("movement['", location,"','", other_country, "'] * S_",age,"_", location, "*I_", age, "_", other_country))))
    }  

    
    transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> V", age, "_", location))))
    
    rate <- append(rate, list(c(paste0("omega * S_", age, "_", location))))
    
    transition <- append(transition, list(c(paste0("I_", age, "_", location, " -> R", age, "_", location))))
    
    rate <- append(rate, list(c(paste0("rho * I_", age, "_", location))))
    
    transition <- append(transition, list(c(paste0("I_", age, "_", location, " -> D", age, "_", location))))
    
    rate <- append(rate, list(c(paste0("mu * I_", age, "_", location))))
    
    
    if (b == 1) {
      
      transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> S_", age_groups[2], "_", location))))
      rate <- append(rate, list(c(paste0("aging[", b, "] * S_", age, "_", location))))


      transition <- append(transition, list(c(paste0("I_", age, "_", location, " -> I_", age_groups[2], "_", location))))
      
      rate <- append(rate, list(c(paste0("aging[", b, "] * I_", age, "_", location))))
      
      transition <- append(transition, list(c(paste0("R_", age, "_", location, " -> R_", age_groups[2], "_", location))))
      
      rate <- append(rate, list(c(paste0("aging[", b, "] * R_", age, "_", location))))
    }
    
    
    if (b == 2) {
      
      transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> S_", age_groups[3], "_", location))))
      
      rate <- append(rate, list(c(paste0("aging[", b , "] * S_", age, "_", location,'* percHR'))))
      
      
      transition <- append(transition, list(c(paste0("S_", age, "_", location, " -> S_", age_groups[4], "_", location))))
      
      rate <- append(rate, list(c(paste0("aging[", b, "] * S_", age, "_", location,'* (1-percHR)'))))
      
      
      
      transition <- append(transition, list(c(paste0("I_", age, "_", location, " -> I_", age_groups[3], "_", location))))
      

      rate <- append(rate, list(c(paste0("aging[", b, "] * I_", age, "_", location,'* percHR'))))
      
      transition <- append(transition, list(c(paste0("I_", age, "_", location, " -> I_", age_groups[4], "_", location))))
      
      rate <- append(rate, list(c(paste0("aging[", b, "] * I_", age, "_", location,'* (1-percHR)'))))

      transition <- append(transition, list(c(paste0("R_", age, "_", location, " -> R_", age_groups[3], "_", location))))
      
      rate <- append(rate, list(c(paste0("aging[", b, "] * R_", age, "_", location,'* percHR'))))
      
      transition <- append(transition, list(c(paste0("R_", age, "_", location, " -> R_", age_groups[4], "_", location))))
      
      rate <- append(rate, list(c(paste0("aging[", b, "] * R_", age, "_", location,'* (1-percHR)'))))
    }
  }
}

# for (c in 1:length(locations)) {
#   for (b in 1:length(age_groups)) {
#     
#     age <- age_groups[b]
#     location <- locations[c]
#     
#     
#     N <- paste0("S_", age_groups[1], "_", location, "+ S_", age_groups[2], "_", location, "+S_", age_groups[3], "_", location,"+S_", age_groups[4], "_", location,
#                 "+I_", age_groups[1], "_", location, "+ I_", age_groups[2], "_", location, "+I_", age_groups[3], "_", location, "+S_", age_groups[4], "_", location,
#                 "+V_", age_groups[1], "_", location, "+ V_", age_groups[2], "_", location, "+V_", age_groups[3], "_", location, "+S_", age_groups[4], "_", location,
#                 "+R_", age_groups[1], "_", location, "+ R_", age_groups[2], "_", location, "+R_", age_groups[3], "_", location, "+S_", age_groups[4], "_", location)
#     
#     reactions <- append(reactions, list(c(paste0("S_", age, "_", location, " -> I_", age, "_", location), 
#                                           paste0("(beta",age,age_groups[1], "* S_", age, "_", location, " * I_", age_groups[1], "_", location,  " / N)",
#                                                  "+ (beta", age, age_groups[2], "* S_", age, "_", location, " * I_", age_groups[2], "_", location,  " / N)",
#                                                  "+ (beta", age, age_groups[3], "* S_", age, "_", location, " * I_", age_groups[3], "_", location,  " / N)",
#                                                  "+ (beta", age, age_groups[4], "* S_", age, "_", location, " * I_", age_groups[4], "_", location,  " / N)"))))
#     
#     reactions <- append(reactions, list(c(paste0("S_", age, "_", location, " -> I_", age, "_", location), 
#                                           paste0("theta_", age,"_",location, "* S_", age, "_", location, "-theta_", age,"_",location," * V_", age, "_", location, 
#                                                  " / sum(S_", age, "_", location, "+I", age, "_", location,"+R", age, "_", location))))
#     
#     for (i in 1:length(locations)) {
#       
#       other_country <- locations[i]
#       
#       reactions <- append(reactions, list(c(paste0("S_", age, "_", location, " -> I_", age, "_", location), 
#                                             paste0("movement['", location,"','", other_country, "'] * S_",age,"_", location, "*I_", age, "_", other_country))))
#     }  
#     
#     
#     reactions <- append(reactions, list(c(paste0("S_", age, "_", location, " -> V", age, "_", location), 
#                                           paste0("omega * S_", age, "_", location))))
#     
#     reactions <- append(reactions, list(c(paste0("I_", age, "_", location, " -> R", age, "_", location), 
#                                           paste0("rho * I_", age, "_", location))))
#     
#     reactions <- append(reactions, list(c(paste0("I_", age, "_", location, " -> D", age, "_", location), 
#                                           paste0("mu * I_", age, "_", location))))
#     
#     
#     if (b == 1) {
#       
#       reactions <- append(reactions, list(c(paste0("S_", age, "_", location, " -> S_", age_groups[2], "_", location), 
#                                             paste0("aging[", b, "] * S_", age, "_", location))))
#       reactions <- append(reactions, list(c(paste0("I_", age, "_", location, " -> I_", age_groups[2], "_", location), 
#                                             paste0("aging[", b, "] * I_", age, "_", location))))
#       reactions <- append(reactions, list(c(paste0("R_", age, "_", location, " -> R_", age_groups[2], "_", location), 
#                                             paste0("aging[", b, "] * R_", age, "_", location))))
#     }
#     
#     
#     if (b == 2) {
#       
#       reactions <- append(reactions, list(c(paste0("S_", age, "_", location, " -> S_", age_groups[3], "_", location), 
#                                             paste0("aging[", b , "] * S_", age, "_", location,'* percHR'))))
#       reactions <- append(reactions, list(c(paste0("S_", age, "_", location, " -> S_", age_groups[4], "_", location), 
#                                             paste0("aging[", b, "] * S_", age, "_", location,'* (1-percHR)'))))
#       reactions <- append(reactions, list(c(paste0("I_", age, "_", location, " -> I_", age_groups[3], "_", location), 
#                                             paste0("aging[", b, "] * I_", age, "_", location,'* percHR'))))
#       reactions <- append(reactions, list(c(paste0("I_", age, "_", location, " -> I_", age_groups[4], "_", location), 
#                                             paste0("aging[", b, "] * I_", age, "_", location,'* (1-percHR)'))))
#       reactions <- append(reactions, list(c(paste0("R_", age, "_", location, " -> R_", age_groups[3], "_", location), 
#                                             paste0("aging[", b, "] * R_", age, "_", location,'* percHR'))))
#       reactions <- append(reactions, list(c(paste0("R_", age, "_", location, " -> R_", age_groups[4], "_", location), 
#                                             paste0("aging[", b, "] * R_", age, "_", location,'* (1-percHR)'))))
#     }
#   }
# }


transition <- as.character(transition)

rate[1]
# Run Gillespie simulation
time <- seq(0, 10, by = 1)
out <- ssa(y_init, transition, rate, params, tf = max(time))

results = as.data.frame(ssa.adaptivetau(y_init, 
                                        transition, 
                                        rate, 
                                        params, 
                                        tf=10))


# Convert to dataframe
output_df <- as.data.frame(out)

# Plot results for selected countries
ggplot(output_df, aes(x = time)) +
  geom_line(aes(y = I1C1, color = 'I1 Country 1')) +
  geom_line(aes(y = I1C2, color = 'I1 Country 2')) +
  geom_line(aes(y = I1C3, color = 'I1 Country 3')) +
  labs(title = 'Stochastic SIR Model with Age and Country Structure', x = 'Time', y = 'Infected Population') +
  scale_color_manual('', values = c('red', 'blue', 'green'))
