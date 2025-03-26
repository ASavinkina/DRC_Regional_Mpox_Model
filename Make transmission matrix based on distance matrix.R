

# 3.4.25
# Code to create transmission matrix due to movement between countries/provinces. Using data from Fraser et al:
# 
# # Using Kramer et al: crossing core borders: 
# 1/ 1+ e^(B0 + B1*(dij/(pi*pj)^B2)
# If crossing a core border: B3 * (1/ 1+ e^(B0 + B1*(dij/(pi*pj)^B2))
# Where: B0= 5.792, B1 = 105.7, B2=0.186, B3= 0.15
# Can be found in Kramer et al Supplement
# 

# Using Kramer et al: crossing country borders: 
# 1/ 1+ e^(B0 + B1*(dij/(pi*pj)^B2)
# If crossing a core border: B3 * (1/ 1+ e^(B0 + B1*(dij/(pi*pj)^B2))
# Where: B0= 5.166, B1 = 157.1, B2=0.189, B3= 0.507
# Can be found in Kramer et al Supplement
# 

# Read in movement data
movement = read.csv("Data/distance_matrix.csv", row.names = 1)


# Read in demographic data
demog = read.csv("Data/Demographic_data_countries_2_27_25.csv")

# Add in new column for "location" which will match locations in movement matrix
demog$location <- demog$Country
demog$location <- ifelse(demog$Country=="DRC",demog$Province, demog$Country)

# Add element with total population
# 
demog$pop <- demog$U5 + demog$U15 + demog$O15

# For movement across borders, add a marker for same country (DRC)
# 

demog$diff_country <- ifelse(demog$Province=="", 1,0)

#
#Clean up movement and demographic country names
#
rownames(movement)[rownames(movement) == "South Sudan"] = "SouthSudan"
colnames(movement)[colnames(movement) == "South.Sudan"] = "SouthSudan"

rownames(movement)[rownames(movement) == "Central African Republic"] = "CAR"
colnames(movement)[colnames(movement) == "Central.African.Republic"] = "CAR"

rownames(movement)[rownames(movement) == "South Sudan"] = "SouthSudan"
colnames(movement)[colnames(movement) == "South_Sudan"] = "SouthSudan"


names(movement) <- gsub(x = names(movement), pattern = "\\.", replacement = "_")

rownames(movement) <- gsub(x = rownames(movement), pattern = "\\-", replacement = "_")

demog$location <-gsub(x = demog$location, pattern = "\\-", replacement = "_")

# beta parameters from kramer et al
# 
# Core borders
B0_core= 5.792
B1_core = 105.7
B2_core=0.186
B3_core= 0.15

# Country borders
B0_count= 5.166
B1_count = 157.1
B2_count= 0.189
B3_count = 0.507

locations <- demog$location

transmission_movement_matrix <- movement

transmission_movement_matrix[transmission_movement_matrix>0] <- 0

for (i in 1:length(locations)) {
  for (j in 1:length(locations)) {
    
    B3 <- ifelse(demog[i,"diff_country"]|demog[j,"diff_country"]==1,B3_core, 1 )
    
    transmission_movement_matrix[i,j] <- B3 * (1/ (1+ exp(B0_core + B1_core*((movement[i,j]/1000)/((demog[i,"pop"]*demog[j,"pop"])^B2_core)))))
                                                    
       
  }
 
}


diag(transmission_movement_matrix) <- 0

write.csv(transmission_movement_matrix, file="Data/transmission_movement_matrix_core.csv")

transmission_movement_matrix[transmission_movement_matrix>0] <- 0

for (i in 1:length(locations)) {
  for (j in 1:length(locations)) {
    
    B3 <- ifelse(demog[i,"diff_country"]|demog[j,"diff_country"]==1,B3_count, 1 )
    
    transmission_movement_matrix[i,j] <- B3 * (1/ (1+ exp(B0_count + B1_count*((movement[i,j]/1000)/((demog[i,"pop"]*demog[j,"pop"])^B2_count)))))
    
    
  }
  
}

diag(transmission_movement_matrix) <- 0

write.csv(transmission_movement_matrix, file="Data/transmission_movement_matrix_country.csv")

# Core borders
B0_core= 5.792
B1_core = 105.7
B2_core=0.186
B3_core= 0



transmission_movement_matrix[transmission_movement_matrix>0] <- 0

for (i in 1:length(locations)) {
  for (j in 1:length(locations)) {
    
    B3 <- ifelse(demog[i,"diff_country"]|demog[j,"diff_country"]==1,B3_count, 1 )
    
    transmission_movement_matrix[i,j] <- B3 * (1/ (1+ exp(B0_count + B1_count*((movement[i,j]/1000)/((demog[i,"pop"]*demog[j,"pop"])^B2_count)))))
    
    
  }
  
}

diag(transmission_movement_matrix) <- 0

write.csv(transmission_movement_matrix, file="Data/transmission_movement_matrix_core_noB3.csv")



