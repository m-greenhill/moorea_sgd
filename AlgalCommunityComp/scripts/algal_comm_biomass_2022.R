##
##    Project: Algal plate percent cover  analysis
##
##
##    Date Created: 2022-03-03
##    
##    Last Updated: 2022-03-03
.

# Set up --------------------------------------------------------

# load libraries
library(tidyverse)
library(here)
library(devtools)
library(geosphere)
library(vegan)

#load data 
biomass_2021 <- read_csv(here("AlgalCommunityComp", "data", "biomass_2021.csv"))
locations <- read_csv(here("AlgalCommunityComp", "data", "tile_locations.csv"))

# Seep distance -------------------------------------------------
  
  seep_dist <- locations %>%
  select(location=Location, CowTagID = Top_Plate_ID, lat, lon) %>%
  mutate(seep_lat = case_when(location == "Varari" ~ -17.54030,
                              location == "Cabral" ~ -17.51638),
         seep_lon = case_when(location == "Varari" ~ -149.8993,
                              location == "Cabral" ~ -149.9125)) %>%
  # find Haversine distance
  mutate(seep_dist_m = distHaversine(cbind(lon, lat), cbind(seep_lon, seep_lat))) %>%
  select(CowTagID, seep_dist_m)


# Biomass ####
biomass <- biomass_2021 %>%
  mutate(weight = total_weight-boat_weight) %>%
  select(site, CowTagID, PlateID, TopBottom = Top_Bottom, code = group, weight)%>%
  filter(code!="TOTAL_BIOMASS")%>%
  group_by(site, CowTagID, PlateID, TopBottom, code) %>%
  summarise(sqrt_weight = sqrt(sum(weight))) %>%
  pivot_wider(names_from = code, 
              values_from = sqrt_weight) %>% # pivot_wider() makes your implicit missing values explicit 
  pivot_longer(GFIL:SAND_MIX, names_to = "code", 
               values_to = "sqrt_weight") %>% # Turn to your desired format (long)
  mutate(sqrt_weight = replace_na(sqrt_weight, 0)) %>%
  pivot_wider(names_from = code, 
              values_from = sqrt_weight) %>% 
  mutate(total = sum(GFIL + RFIL + CCA + RCOR + CYO + RMUC + SAND + SAND_MIX),
        GFIL = sum(GFIL +SAND+SAND_MIX)) %>%
  pivot_longer(GFIL:total, names_to = "code", 
               values_to = "sqrt_weight") %>%
  left_join(seep_dist, by = 'CowTagID') %>%
  mutate(TopBottom = factor(TopBottom, levels = c('Top', 'Bottom'))) %>%
  filter(code!=c("SAND", "SAND_MIX"))

# 1. Cabral ####

# top bottom scatter
cabral_biomass_scatter <-  biomass %>%
  filter(site == "Cabral") %>%
  ggplot(aes(x = seep_dist_m, y = sqrt_weight)) +
  geom_point(show.legend = FALSE, 
             size = 1) +
  geom_smooth(method= "lm" , alpha = 0.3)+
  theme_bw() +
  facet_grid(TopBottom~code) +
  labs(x = "Distance from seep (m)",
       y = "√ Biomass (g)",
       title = "Algal biomass across SGD gradient - Cabral") +
  xlim(14, 154) +
  theme(axis.title.y=element_text(angle=0, vjust=0.5))

ggsave(here("AlgalCommunityComp", "outputs", "CabralBiomass.png"))


# 2. Varari #####

vararibiomass_scatter <-  biomass %>%
  filter(site == "Varari") %>%
  ggplot(aes(x = seep_dist_m, y = sqrt_weight)) +
  geom_point(show.legend = FALSE, 
             size = 1) +
  geom_smooth(method= "lm" , alpha = 0.3, color = "purple")+
  theme_bw() +
  facet_grid(TopBottom~code) +
  labs(x = "Distance from seep (m)",
       y = "√ Biomass (g)",
       title = "Algal biomass across SGD gradient - Varari") +
  xlim(14, 154) +
  theme(axis.title.y=element_text(angle=0, vjust=0.5))

ggsave(here("AlgalCommunityComp", "outputs", "VarariBiomass.png"))


# 3. Total

# gradient
total_biomass_seep <- biomass %>%
  filter(code == "total") %>%
  mutate(weight = sqrt_weight^2) %>%
  ggplot(aes(x = seep_dist_m, y = weight, color = site, linetype = TopBottom)) +
  geom_point(aes( shape = TopBottom),
             size = 1,
             show.legend=FALSE) +
  geom_smooth(method = "lm", alpha = 0.3) +
  theme_bw() +
  facet_grid(~site) +
  theme(axis.title.y=element_text(angle=0, vjust=0.5)) +
  labs(x = "Distance from seep (m)",
       y = "Biomass (g) ",
       title = "Total algal biomass")
  
ggsave(here("AlgalCommunityComp", "outputs", "TotalBiomassSeep.png"))

#boxplot

total_biomass <- biomass %>%
  filter(code == "total") %>%
  mutate(weight = sqrt_weight^2) %>%
  group_by(site, TopBottom) %>%
  summarise(mean_weight = mean(weight),
            sd = sd(weight)) %>%
  ggplot(aes(x = site, y = mean_weight, fill = TopBottom))+
  geom_bar(stat= "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin=mean_weight-sd, ymax=mean_weight+sd), width=.2,position=position_dodge(.9)) +
  scale_fill_brewer(palette="Blues") +
  theme_bw() +
  labs (x = "",
        y = "Mean biomass (g)",
        title = "Mean tile algal biomass across sites")+
  theme(axis.title.y=element_text(angle=0, vjust=0.5))

ggsave(here("AlgalCommunityComp", "outputs", "TotalBiomass.png"))

total_biomass_site <- biomass %>%
  filter(code == "total") %>%
  group_by(site) %>%
  summarise(mean_weight = mean(weight),
            sd = sd(weight)) %>%
  ggplot(aes(x = site, y = mean_weight, fill = site))+
  geom_bar(stat= "identity") +
  geom_errorbar(aes(ymin=mean_weight-sd, ymax=mean_weight+sd), width=.2,position=position_dodge(.9)) +
  scale_fill_brewer(palette="Blues") +
  theme_bw() 
  
ggsave(here("AlgalCommunityComp", "outputs", "TotalBiomassSite.png"))


# biomass nmds ####
pc.tb.wide <- percent_cover_dist %>%
  mutate(sqrt = sqrt(percent_cover)) %>%
  select(site, TileID, TopBottom, code, sqrt) %>%
  pivot_wider(names_from = code,
              values_from = sqrt) %>%
  ungroup()

pc.tb.matrix <-  pc.tb.wide %>%
  select(-site, -TileID, -TopBottom)

pc.tb.mds <- metaMDS(pc.tb.matrix, k = 2, trymax = 200)


pc.tb.data.scores <- as.data.frame(scores(pc.tb.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
pc.tb.data.scores$TileID <- pc.tb.wide$TileID
# create a column of site names, from the rownames of data.scores
pc.tb.data.scores$site <- pc.tb.wide$site  #  add the grp variable created earlier
pc.tb.data.scores$TopBottom <- pc.tb.wide$TopBottom



pc_tb_mds_plot <- ggplot(pc.tb.data.scores) +
  geom_point(stat ="identity", size = 2,
             aes(NMDS1, NMDS2, colour = TopBottom, shape=site)) +
  theme_bw()
