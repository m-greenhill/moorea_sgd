
# Description: Analysis for groundwater reef fish surveys & herbivory project
# Author: Maya Zeff
# created: 2023-2-27
# modified: 2023-2-27

######## Load libraries ########

library(tidyverse)
library(here)
library(vegan)
library(ggplot2)

######## File paths ########

# set working directory 
setwd("~/Documents/R_sessions/silbiger_lab/moorea_sgd/herbivory/data")


# read in functional group data
functional.group <- read.csv("species_functional.csv") 

# read in herbivory cam data

all.cam <- read_csv("cam_surveys_prelim.csv") %>%
  # add name descriptions
  mutate(SGD = str_replace(SGD, "high", "High SGD"),
         SGD = str_replace(SGD, "low", "Low SGD"),
         tide = str_replace(tide, "high", "High Tide"),
         tide = str_replace(tide, "low", "Low Tide")) %>%
  mutate(SGD = as.factor(SGD),
         tide = as.factor(tide),
         rep = as.factor (rep),
         site = as.factor(site),
         habitat = as.factor(habitat)) %>%
  left_join(functional.group, by = 'code')





######## NMDS plots ########
###### Presence Absence NMDS #######

presence <- all.cam %>%
            # choose only presence observations
            filter(activity == "nat") %>%
            # group down to individual videos
            group_by(habitat, SGD, tide, site, rep) %>%
            # choose only unique species
            distinct(code) %>%
            # turn presence/absence into 1s and 0s
            mutate(present = 1) %>%
            # spread presence values into nmds format
            spread(key = code, value = present, fill = 0) %>%
            ungroup() 
           

# turn data frame into NMDS matrix
presence.matrix <- presence %>%
                dplyr::select(6:last_col()) %>%
                as.matrix()

# conduct NMDS
presence.results <- metaMDS(presence.matrix, distance = "bray")


# extract NMDS scores (x and y coordinates)
presence.scores <- as.data.frame(scores(presence.results)$sites)

# add columns to data frame 
presence.scores$habitat = presence$habitat
presence.scores$SGD = presence$SGD
presence.scores$tide = presence$tide
presence.scores$site = presence$site



# plot NMDS
presence.plot <- presence.scores %>%
                ggplot(aes(x = NMDS1, y = NMDS2)) + 
                  geom_point(size = 4, 
                             aes(shape = factor(SGD), 
                                 colour = factor(site), 
                                 fill = factor(tide), 
                                 stroke = 1.5)) + 
                  scale_shape_manual(values = c(21, 24)) +
                  #add6 colors 
                  scale_color_manual(values = c('gold2',
                                                'gold3',
                                                'darkorange2',
                                                'darkslategray2',
                                                'darkslategray3',
                                                'darkslategray4',
                                                'darkslategray')) +
                  scale_fill_manual(values = c("darkolivegreen3", 
                                               "darkgreen"))+
                  guides(fill = guide_legend(override.aes=list(shape=21))) + 
                  facet_wrap(~habitat)+
                  theme(axis.text.y = element_text(colour = "black", 
                                                   size = 12), 
                        axis.text.x = element_text(colour = "black", 
                                                   size = 12), 
                        legend.text = element_text(size = 12, 
                                                   colour ="black"), 
                        legend.position = "right", 
                        axis.title.y = element_text(size = 14), 
                        axis.title.x = element_text(size = 14, 
                                                    colour = "black"), 
                        legend.title = element_text(size = 14, 
                                                    colour = "black", 
                                                    face = "bold"), 
                        panel.background = element_blank(), 
                        panel.border = element_rect(colour = "black", 
                                                    fill = NA, 
                                                    size = 1.2),
                        legend.key = element_blank()) +
                  labs(x = "NMDS1", 
                       colour = "Site", 
                       y = "NMDS2", 
                       shape = "SGD",
                       fill = "Tide",
                       title = "Fish Species Presence by Habitat") 




####### Species Abundance NMDS ###### 
abundance.species <- all.cam %>%
  # choose only presence observations
  filter(activity == "nat") %>%
  # group down to individual videos
  group_by(habitat, SGD, tide, site, rep, code) %>%
  summarise(abundance = n()) %>%
  na.omit() %>%
  spread(key = code, value = abundance, fill = 0) %>%
  ungroup() %>%
  mutate(SGD = str_replace(SGD, "high", "High SGD"),
         SGD = str_replace(SGD, "low", "Low SGD"),
         tide = str_replace(tide, "high", "High Tide"),
         tide = str_replace(tide, "low", "Low Tide"))

# turn data frame into NMDS matrix
abundance.species.matrix <- abundance.species %>%
  dplyr::select(6:last_col()) %>%
  as.matrix()

# conduct NMDS
abundance.sp.results <- metaMDS(abundance.species.matrix, distance = "bray")


# extract NMDS scores (x and y coordinates)
abundance.sp.scores <- as.data.frame(scores(abundance.sp.results)$sites)

# add columns to data frame 
abundance.sp.scores$habitat = abundance.species$habitat
abundance.sp.scores$SGD = abundance.species$SGD
abundance.sp.scores$tide = abundance.species$tide
abundance.sp.scores$site = abundance.species$site


# plot NMDS
abundance.sp.plot <- abundance.sp.scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, 
             aes(shape = factor(SGD), 
                 colour = factor(site), 
                 fill = factor(tide), 
                 stroke = 1.5)) + 
  scale_shape_manual(values = c(21, 24)) + 
  scale_color_manual(values = c('gold2',
                                'gold3',
                                'darkorange2',
                                'darkslategray2',
                                'darkslategray3',
                                'darkslategray4')) +
  scale_fill_manual(values = c("darkolivegreen3", 
                               "darkgreen"))+
  guides(fill = guide_legend(override.aes=list(shape=21))) + 
  facet_wrap(~habitat)+
  theme(axis.text.y = element_text(colour = "black", 
                                   size = 12), 
        axis.text.x = element_text(colour = "black", 
                                   size = 12), 
        legend.text = element_text(size = 12, 
                                   colour ="black"), 
        legend.position = "right", 
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14, 
                                    colour = "black"), 
        legend.title = element_text(size = 14, 
                                    colour = "black", 
                                    face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    size = 1.2),
        legend.key = element_blank()) +
  labs(x = "NMDS1", 
       colour = "Site", 
       y = "NMDS2", 
       shape = "SGD",
       fill = "Tide",
       title = "Fish Species Abundance by Habitat") 


####### Functional Abundance NMDS ###### 
                  
abundanace.functional <- all.cam %>%
  # choose only presence observations
  filter(activity == "nat") %>%
  # group down to individual videos
  group_by(habitat, SGD, tide, site, rep, functional) %>%
  summarise(abundance = n()) %>%
  na.omit() %>%
  spread(key = functional, value = abundance, fill = 0) %>%
  ungroup() %>%
  mutate(SGD = str_replace(SGD, "high", "High SGD"),
         SGD = str_replace(SGD, "low", "Low SGD"),
         tide = str_replace(tide, "high", "High Tide"),
         tide = str_replace(tide, "low", "Low Tide"))

# turn data frame into NMDS matrix
abundance.matrix <- abundance.functional %>%
  dplyr::select(6:last_col()) %>%
  as.matrix()

# conduct NMDS
abundance.results <- metaMDS(abundance.matrix, distance = "bray")


# extract NMDS scores (x and y coordinates)
abundance.scores <- as.data.frame(scores(abundance.results)$sites)

# add columns to data frame 
abundance.scores$habitat = abundance$habitat
abundance.scores$SGD = abundance$SGD
abundance.scores$tide = abundance$tide
abundance.scores$site = abundance$site


# plot NMDS
abundance.plot <- abundance.scores %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, 
             aes(shape = factor(SGD), 
                 colour = factor(site), 
                 fill = factor(tide), 
                 stroke = 1.5)) + 
  scale_shape_manual(values = c(21, 24)) + 
  scale_color_manual(values = c('gold2',
                                'gold3',
                                'darkorange2',
                                'darkslategray2',
                                'darkslategray3',
                                'darkslategray4')) +
  scale_fill_manual(values = c("darkolivegreen3", 
                               "darkgreen"))+
  guides(fill = guide_legend(override.aes=list(shape=21))) + 
  facet_wrap(~habitat)+
  theme(axis.text.y = element_text(colour = "black", 
                                   size = 12), 
        axis.text.x = element_text(colour = "black", 
                                   size = 12), 
        legend.text = element_text(size = 12, 
                                   colour ="black"), 
        legend.position = "right", 
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14, 
                                    colour = "black"), 
        legend.title = element_text(size = 14, 
                                    colour = "black", 
                                    face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    size = 1.2),
        legend.key = element_blank()) +
  labs(x = "NMDS1", 
       colour = "Site", 
       y = "NMDS2", 
       shape = "SGD",
       fill = "Tide",
       title = "Fish Functional Abundance by Habitat") 

######## Herbivory plots ########

###### boxplot #######

rates <- all.cam %>%
  filter(activity == "feed") %>%
  group_by(habitat, SGD, tide, rep) %>%
  summarise(total_bites = sum(number_bites)) %>%
  group_by(habitat, SGD, tide) %>%
  summarise(average_total_bites = mean(total_bites)) %>%
  ggplot(aes(x = SGD, 
             y = average_total_bites)) +
  geom_bar(aes(fill = tide), 
           position = "dodge", 
           stat = "identity") +
  facet_wrap(~habitat)


rates <- all.cam %>%
  filter(activity == "feed") %>%
  group_by(habitat, SGD, tide, rep, functional) %>%
  summarise(total_bites = sum(number_bites)) %>%
  group_by(habitat, SGD, tide) %>%
  mutate(average_total_bites = mean(total_bites)) %>%
  ggplot(aes(x = tide, 
             y = average_total_bites)) +
  geom_bar(aes(fill = functional), 
           position = "dodge", 
           stat = "identity") +
  facet_wrap(SGD~habitat)
    
percent.rates <- all.cam %>%
  filter(activity == "feed",
         rep == "1") %>%
  group_by(habitat, SGD, tide, site, functional) %>%
  summarise(substrate_bites = sum(number_bites)) %>%
  na.omit()%>%    
  group_by(habitat, SGD, tide) %>%
  mutate(total_bites = sum(substrate_bites)) %>%
  ungroup()%>%
  mutate(percent_bite = (substrate_bites / total_bites) * 100) %>%
  na.omit()%>%
  ggplot(aes(x = tide, 
             y = percent_bite)) +
  geom_bar(aes(fill = functional), 
           stat = "identity") +
  facet_wrap(SGD~habitat)

 
  
  
          
  
  
  
  
          
  
          

