#############################################
# Pumpkins Dataset Script                   #
# A script to analyse the Pumpkins dataset  #
# Author : Hannah Byrne                     #
# Date : 12/11/25                           #
#############################################

#Read in pumpkins.csv
pumpkins <- read.csv("C:/Users/hanna/OneDrive/Documents/Computational Bio/Coding Challenge/pumpkins_03.csv")
library(tidyverse)

#Identify the heaviest pumpkin (variety, location, year)
#Combine 'city', 'state_prov' and 'country' columns to get single location column
pumpkins_2 <- pumpkins %>%
  mutate(location=paste(city, state_prov, country, sep = ", ")) 
#Create a tibble of just variety, location, year(id) and weight_lbs columns
pumpkins_heavy <- pumpkins_2 %>%
  select(variety, location, id, weight_lbs) %>%
  slice_max(weight_lbs)


#Convert lbs -> kg and create a new column ('weight_kg')
pounds_to_kilos <- (kilos = pumpkins_2$weight_lbs * 0.453592) #Weight conversion function
pumpkins_3 <- pumpkins_2 %>%
    mutate(weight_kg=pounds_to_kilos, .after=weight_lbs) #Add new column 'weight_kg"

#Add a ' weight_class' column (light, medium, heavy)
pumpkins_3 %>%
    summarise(mean_weight = mean(weight_lbs, na.rm = TRUE),
              max_weight = max(weight_lbs, na.rm = TRUE),
              min_weight = min(weight_lbs, na.rm = TRUE)) #Determine parameters for weight class based on column summary 
pumpkins_4 <- subset(pumpkins_3, est_weight != 0.0000 & est_weight < 9000)#Removes all estimated weights = 0.0000 and any above 9000 lbs
pumpkins_wc <- pumpkins_4 %>%
    mutate(weight_class = case_when(
      weight_lbs <= 750 ~ 'Light (<750 lbs)',
      weight_lbs <=1500 ~ 'Medium (750-1500 lbs)',
      TRUE ~ 'Heavy (>1500 lbs)'
    )) #Creates a new column and categorises 'weight_lbs' into different weight classes

#Plot a graph of estimated weight vs. actual weight (coloured by class)
est_vs_actual_plot <- ggplot(pumpkins_wc, aes(x = est_weight, y = weight_lbs, colour=weight_class))+
  geom_point(size = 0.5)+
  scale_colour_manual(values = c('Light (<750 lbs)' = 'hotpink2',
                               'Medium (750-1500 lbs)' = 'deeppink2',
                               'Heavy (>1500 lbs)' = 'deeppink4'))+
  labs(title = "Estimated vs. Actual Weight (lbs) Plot", x = "Estimated Weight (lbs)", y = "Actual Weight (lbs)", colour = "Weight Class") 
#Save plot
ggsave("estimated_vs_actual_plot.png", est_vs_actual_plot, path = "C:/Users/hanna/OneDrive/Documents/Computational Bio/Coding Challenge/R coursework/" )

#Filter for 3 countries and save 'pumpkins_filtered.csv'
pumpkins_countries <- pumpkins %>%
  count(country) %>%
  arrange(desc(n))#Determine countries in dataset
pumpkins_filtered <- pumpkins %>%
  filter(country == "United States" | country == "Canada" | country == "Germany") #Filter data to 3 countries
write.csv(pumpkins_filtered, file='C:/Users/hanna/OneDrive/Documents/Computational Bio/Coding Challenge/R coursework/pumpkins_filtered.csv', na = '')

#Summarise mean weights by country and variety
#Mean weight in lbs by the 3 filtered countries
pumpkins_filtered_country <- pumpkins_filtered %>%
  group_by(country) %>%
  summarise(mean_weight_lbs = mean(weight_lbs, na.rm = TRUE))
pumpkins_max_weight <- pumpkins_filtered_country %>%
  slice_max(mean_weight_lbs)


#Mean weight in lbs by country and variety
pumpkins_filtered2 <- pumpkins_filtered %>%
  group_by(country, variety) %>%
  summarise(mean_weight_lbs = mean(weight_lbs, na.rm = TRUE))
pumpkins_min_weight <- pumpkins_filtered2 %>%
  slice_min(mean_weight_lbs)

#Boxplot of weight for the three countries
filtered_weight_boxplot <- ggplot(aes(x = country,y = weight_lbs, fill = country), data = pumpkins_filtered)+
  geom_boxplot()+
  labs(x = "Country", y = "Weight (lbs)", title = "Boxplot of Weight Distribution per Country")+
  scale_fill_manual(values = c('Canada' = 'hotpink2',
                               'Germany' = 'olivedrab3',
                               'United States' = 'steelblue3'))
#Save boxplot
ggsave("country_weight_boxplot.png", filtered_weight_boxplot, path = "C:/Users/hanna/OneDrive/Documents/Computational Bio/Coding Challenge/R coursework/")

#Facet plot for each variety 
facet_plot <- ggplot(aes(country, weight_lbs, fill = country), data = pumpkins_filtered)+
  geom_boxplot()+
  facet_wrap(~variety)+
  labs(x = "Country", y = "Weight (lbs)", title = "Facet Plot of Weight Distribution per Country for each Variety")+
  theme(
    axis.text.x = element_text(angle = 90)
  )+
  scale_fill_manual(values = c('Canada' = 'hotpink2',
                               'Germany' = 'olivedrab3',
                               'United States' = 'steelblue3'))
#Save facet plot
ggsave("facet_variety&country_plot.png", facet_plot, path = "C:/Users/hanna/OneDrive/Documents/Computational Bio/Coding Challenge/R coursework/")
