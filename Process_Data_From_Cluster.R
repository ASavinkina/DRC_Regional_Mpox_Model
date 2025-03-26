# Process output from cluster
# 

library(tidyverse)


output_all <- read.csv("Output/Full_output_test_32425_st.csv")

# For summary statistics (and for calibration), keep only last time period
# This way, can calculate total cases (including deaths, recovered, and currently infected)

output_all_end <- output_all %>%
  group_by(Run) %>%
  top_n(1, t)

output_all_end <- unique(output_all_end)

# Keep only infected, recovered, dead (for now, dropping vaccinated and susceptible)

df2 <- output_all_end %>%
  select(starts_with("I")|starts_with("D")|starts_with("R"))

# Pivot dataset longer to make summary stats easier

df3 <- df2 %>%
  pivot_longer(
    cols= 1:(length(df2)-1),
    names_to = "Group",
    values_to = "Count"
  )

# Separate out longer group names into separate columns of group, age and location

df4 <- df3 %>% 
  separate(Group,c("Group", "Age", "Location"), '_') 

# FOR ALL INFECTED:
# First sum across all groups by run

df5_allinfected <- df4 %>%
  group_by(Age, Location,Run) %>%
  summarise(Count=sum(Count))

# Then summarize across runs

df6_allinfected <- df5_allinfected %>%
  group_by(Age, Location) %>%
  summarise_at(vars(Count), 
               list(min=min, median=median, max=max,
                    mean=mean, sd=sd)
  )

# Drop age to have overall case counts

df6_allinfected_noage <- df5_allinfected %>%
  group_by(Location) %>%
  summarise_at(vars(Count), 
               list(min=min, median=median, max=max,
                    mean=mean, sd=sd)
  )

# FOR ALL DEATHS ONLY
# First drop all other than deaths

df_alldeaths <- df4[which(df4$Group=="D"),]


# Collapse across runs to summarize
#
df6_alldeaths <- df_alldeaths %>%
  group_by(Age, Location) %>%
  summarise_at(vars(Count), 
               list(min=min, median=median, max=max,
                    mean=mean, sd=sd)
  )

# Collapse across ages, as well

df6_alldeaths_noage <- df_alldeaths %>%
  group_by(Location) %>%
  summarise_at(vars(Count), 
               list(min=min, median=median, max=max,
                    mean=mean, sd=sd)
  )


# FOR PLOTS:
# 
# Pivot from wide to long, with the entire dataset
# 


fulldata_long <- output_all %>%
  pivot_longer(
    cols= 2:(length(output_all)-1),
    names_to = "Group",
    values_to = "Count"
  )



# Separate out longer group names into separate columns of group, age and location

fulldata_long1 <- fulldata_long %>% 
  separate(Group,c("Group", "Age", "Location"), '_') 

# Limit to only those currently infected

fulldata_infected <- fulldata_long1[which(fulldata_long1$Group=="Ia"|fulldata_long1$Group=="Ib"),]

# Collapse across groups to get all cases at each time period
# 

fulldata_infected1 <- fulldata_infected %>%
  group_by(t, Run, Group, Location) %>%
  summarise(Count=sum(Count))

# Graph of current infections

ggplot(fulldata_infected1) + geom_line(aes(x=t, y=Count, group=interaction(Run,Group), color=Group)) + 
  facet_wrap(~Location, scales='free') + theme_classic()


