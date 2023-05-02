## Script for analysis of herbivory assays 

# loading libraries
library(tidyverse)
library(here)
library(vegan)
library(car)
library(lme4) # mixed effects models
library(lmerTest)
library(performance)

colors <- c("#9bb076", '#7aa3c9')

# read in file
nutrient_ug <- read_csv(here('herbivory', 'data', 'padina_nutrient.csv')) %>%
            group_by(SGD) %>%
  mutate(ug_N = ug_N / weight_mg,
         del15N = del15N / weight_mg) %>%
  ggplot(aes(x = SGD, y=ug_N)) +
  geom_boxplot(aes(fill = SGD), show.legend = FALSE, alpha = .6) +
  theme(legend.position ="none") +
  scale_fill_manual(values=colors) +
  #facet_wrap(~SGD, scales = "free" )+
  theme_bw() +
  labs(x = "SGD Algal Source",
       y = "ug_N",
       title = "Nitrogen Content (ug)") 

quartz()
ggsave(here('herbivory', 'outputs', 'nutrient_ug.jpg'))

nutrient_delN <- read_csv(here('herbivory', 'data', 'padina_nutrient.csv')) %>%
  group_by(SGD) %>%
  mutate(ug_N = ug_N / weight_mg,
         del15N = del15N / weight_mg) %>%
  ggplot(aes(x = SGD, y=del15N)) +
  geom_boxplot(aes(fill = SGD), show.legend = FALSE, alpha = .6) +
  theme(legend.position ="none") +
  scale_fill_manual(values=colors) +
  #facet_wrap(~SGD, scales = "free" )+
  theme_bw() +
  labs(x = "SGD Algal Source",
       y = "del15N",
       title = "Nitrogen Content (ug)") 



