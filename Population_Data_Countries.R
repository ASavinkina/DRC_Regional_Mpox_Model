# Get by-age populations for each country of interest: Angola, Burundi, Kenya, Rwanda, Uganda, Zambia, Zimbabwe, CAR, Sudan,
# South Sudan
# 
library(tidyverse)

#Read in demographic data to fill
demo_data <- read.csv("Data/Demographic_data_2.20.25.csv")

#List of countries 
#
Countries <- unique(demo_data$Country)

# Exclude DRC (already have pop) and Zambia (couldn't download pop)
# 
Countries <- Countries[!Countries %in% c("DRC")]

demo_data$birth_per_1000 <- 0

for (i in 1:length(Countries)) {
  
  i=Countries[i]
  
data <- read.csv(paste0("Data/Populations/World_Bank_Group_",i,".csv"), skip=4)

data_pop <- data %>% filter(grepl(c('Population'), Indicator.Name))

data_mort <- data %>% filter(grepl(c('Mortality'), Indicator.Name))

data_birth <- data %>% filter(grepl(c('Birth rate, crude'), Indicator.Name))

data_pop2 <- data_pop %>% filter(Indicator.Name %in% c("Population, male", "Population, female", "Population ages 15-64, female",
                                        "Population ages 15-64, male","Population ages 0-14, female",
                                        "Population ages 0-14, male", "Population ages 00-04, male (% of male population)",
                                        "Population ages 00-04, female (% of female population)", "Population ages 65 and above, male",
                                        "Population ages 65 and above, female", "")) 

data_pop2 <- data_pop2 %>% select(Country.Name, Indicator.Name,X2023)

data_pop2_wide <- spread(data_pop2, Indicator.Name, X2023)

data_pop2_wide$O15 <- data_pop2_wide$`Population ages 15-64, female` + data_pop2_wide$`Population ages 15-64, male` +
  data_pop2_wide$`Population ages 65 and above, female` + data_pop2_wide$`Population ages 65 and above, male`

data_pop2_wide$U5 <- (data_pop2_wide$`Population ages 00-04, female (% of female population)`/100 * data_pop2_wide$`Population, female`) +
 (data_pop2_wide$`Population ages 00-04, male (% of male population)`/100 * data_pop2_wide$`Population, male`) 

data_pop2_wide$U15 <- data_pop2_wide$`Population ages 0-14, female` + data_pop2_wide$`Population ages 0-14, male` -
  data_pop2_wide$U5

#Check
sum(data_pop2_wide$`Population, female`, data_pop2_wide$`Population, male`)
sum(data_pop2_wide$U5, data_pop2_wide$O15, data_pop2_wide$U15)


demo_data$U5 <- ifelse(demo_data$Country==i, data_pop2_wide$U5,demo_data$U5 )
demo_data$U15 <- ifelse(demo_data$Country==i, data_pop2_wide$U15,demo_data$U15 )
demo_data$O15 <- ifelse(demo_data$Country==i, data_pop2_wide$O15,demo_data$O15 )

demo_data$birth_per_1000 <- ifelse(demo_data$Country==i, data_birth$X2022, demo_data$birth_per_1000)

}

demo_data$birth_per_1000 <- ifelse(demo_data$Country=="DRC", 40.08, demo_data$birth_per_1000)

write.csv(demo_data, file="Data/Demographic_data_countries_3_18_25.csv")
