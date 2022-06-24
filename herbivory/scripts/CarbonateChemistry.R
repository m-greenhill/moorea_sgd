# See how bad carb chem data is, need to repeat hell weeek?

# bahaha

# loading libraries
library(tidyverse)
library(here)
library(vegan)


# read in files


carb <- read_csv(here('herbivory', 'data', 'carbonate_chemistry.csv')) %>%
  filter(Sample!="SandSGD3") %>% 
  filter(Group == "ReefAmbient" | Group =="SandAmbient")

carb_plot <-  carb %>%
 # filter(Day_Night=="night") %>%
  filter(Sample!="SandSGD3") %>%
  ggplot(aes(x = Tide, y=Salinity_Gump )) +
  geom_boxplot() +
  facet_wrap(~Group)

carb_pH <- read_csv(here('herbivory', 'data', 'carbonate_chemistry.csv')) %>%
  # filter(Day_Night=="night") %>%
  #filter(Sample!="SandSGD3") %>%
  unite(col = "Tide_Time", Tide, Day_Night,sep = "_", remove = F) %>%
  filter(Tide_Time != "low_day") %>%
  ggplot(aes(x = Tide, y=pH_Gump )) +
  geom_boxplot() +
  facet_wrap(~Group)


aov <- lm(Salinity_Gump~Tide, data = carb)

summary(aov)
