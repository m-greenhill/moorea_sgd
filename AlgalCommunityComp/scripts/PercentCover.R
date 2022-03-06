library(dplyr)
library(magrittr)
library(lubridate)
library(broom)
library(readr)
library(ggplot2)
library(tidyr)
library(forcats)
library(readxl)
library(viridis)
library(vegan)
library(vegan3d)
library(rgl)
library(devtools)
library(pairwiseAdonis)
library(ggrepel)

setwd("~/Documents/R Sessions/Silbiger Lab/moorea_sgd/AlgalCommunityComp")

# clean up  
percent_cover <- read.csv("PercentCover.csv") 

cover <- percent_cover %>% 
          filter(TopBottom=="Top", site =="Varari") %>%
          mutate(TileID = factor(TileID, levels = c("V13", "V1", "V2", "V3", "V4", "V6", "V5", "V16", "V15", "V7", "V18", "V11", "V12", "V10", "V19", "V20", "V14"))) %>%
          select(site, TileID, CowTagID, TopBottom, placement, point, code) %>%
          group_by(site, TileID, code) %>%
          summarise(point = n()) %>%
          mutate(proportion = point/12)
  
plot <- ggplot(cover, aes(TileID, proportion)) +
        geom_point() +
        facet_wrap(vars(code)) 
       


