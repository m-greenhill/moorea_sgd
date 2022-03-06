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
library(viridis)
library(mgcv)
library(directlabels)
library(ggrepel)

#load data 
percent_cover <- read_csv(here("AlgalCommunityComp", "data", "percent_cover_2021.csv"))

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
                           


# Percent cover -------------------------------------------------

percent_cover_dist <- percent_cover %>%
  select(site, CowTagID, TileID, TopBottom, point, code) %>%
  group_by(site, CowTagID, TileID, TopBottom, code) %>%
  summarise(points = n()) %>%
  pivot_wider(names_from = code, 
              values_from = points) %>% # pivot_wider() makes your implicit missing values explicit 
  pivot_longer(GFIL:GHOL, names_to = "code", 
               values_to = "points") %>% # Turn to your desired format (long)
  mutate(points = replace_na(points, 0)) %>%  # Replace missing values (NA) with 0s
  group_by(site, CowTagID, TileID, TopBottom, code) %>%
  summarise(percent_cover = (points/12) * 100) %>%
  left_join(seep_dist, by = 'CowTagID') %>%
  mutate(sqrt_percent_cover = sqrt(percent_cover)) %>%
  ungroup() %>%
  filter(code!= "SAND") %>%
  mutate(code = recode(code,
                       'X' = "Blank")) %>%
  mutate(TopBottom = factor(TopBottom, levels = c('Top', 'Bottom')))


# 1. Cabral ####

# density   
cabral_density <- percent_cover_dist %>%
  filter(site == "Cabral") %>%
  ggplot(aes(x = percent_cover)) +
  geom_density() +
  facet_grid(TopBottom~code)

# Top bottom scatter
cabral_percent_scatter <-  percent_cover_dist %>%
  filter(site == "Cabral") %>%
  ggplot(aes(x = seep_dist_m, y = sqrt_percent_cover)) +
         geom_point(show.legend = FALSE, 
                    size = 1) +
        geom_smooth(method= "lm" , alpha = 0.3)+
        theme_bw() +
        facet_grid(TopBottom~code) +
        labs(x = "Distance from seep (m)",
             y = "√ % Cover",
            title = "Algal % cover across SGD gradient - Cabral") +
  xlim(14, 154) +
  theme(axis.title.y=element_text(angle=0, vjust=0.5))
  
  
ggsave(here("AlgalCommunityComp", "outputs", "CabralPercentScatter.png"))
 

# 2. Varari ####      

# Top bottom scatter
varari_percent_scatter <-  percent_cover_dist %>%
  filter(site == "Varari") %>%
  ggplot(aes(x = seep_dist_m, y = sqrt_percent_cover)) +
  geom_point(show.legend = FALSE, 
             size = 1) +
  geom_smooth(method= "lm" , alpha = 0.3, color = "purple")+
  theme_bw() +
  facet_grid(TopBottom~code) +
  labs(x = "Distance from seep (m)",
       y = "√ % Cover",
       title = "Algal % cover across SGD gradient - Varari") +
  theme(axis.title.y=element_text(angle=0, vjust=0.5)) +
  xlim(14, 154)

ggsave(here("AlgalCommunityComp", "outputs", "VarariPercentScatter.png"))

# Cyanobacteria  ------------------------------------------------

# 1. Cabral #### 
cyo_cabral <- percent_cover %>%
  mutate(TopBottom = factor(TopBottom, levels = c('Top', 'Bottom'))) %>%
  filter(site == "Cabral") %>%
  select(CowTagID, TopBottom, cyo) %>%
  group_by(CowTagID, TopBottom) %>%
  summarise(cyo = mean(cyo)) %>%
  left_join(seep_dist, by = 'CowTagID') %>%
  ggplot(aes(x = seep_dist_m, y = cyo)) +
    geom_point(aes(color=seep_dist_m),
                show.legend = FALSE,
                size = 2) +
    scale_color_gradient(low = "black",
                         high = "cyan") +
    theme_bw() +
    facet_grid(TopBottom~.) +
    labs(x = "Distance from seep (m)",
         y = "",
         title = "Cyanobacterial mat presence - Cabral") +
  theme(axis.title.y=element_text(angle=0, vjust=0.5),
        axis.text.y=element_blank(),)
    
ggsave(here("AlgalCommunityComp", "outputs", "CabralCyanobacteria.png")) 

# 2. Varari #### 
cyo_varari <- percent_cover %>%
  mutate(TopBottom = factor(TopBottom, levels = c('Top', 'Bottom'))) %>%
  filter(site == "Varari") %>%
  select(CowTagID, TopBottom, cyo) %>%
  group_by(CowTagID, TopBottom) %>%
  summarise(cyo = mean(cyo)) %>%
  left_join(seep_dist, by = 'CowTagID') %>%
  ggplot(aes(x = seep_dist_m, y = cyo)) +
  geom_point(aes(color=seep_dist_m),
             show.legend = FALSE,
             size = 2) +
  scale_color_gradient(low = "black",
                       high = "cyan") +
  theme_bw() +
  facet_grid(TopBottom~.) +
  labs(x = "Distance from seep (m)",
       y = "",
       title = "Cyanobacterial mat presence - Varari") +
  theme(axis.title.y=element_text(angle=0, vjust=0.5),
        axis.text.y=element_blank(),)

ggsave(here("AlgalCommunityComp", "outputs", "VarariCyanobacteria.png")) 


# MDS ####
# 1. Gradient nmds ---------------------------------------------------------


c.wide <- percent_cover_dist %>%
  mutate(sqrt = sqrt(percent_cover)) %>%
  select(site, TileID, TopBottom, code, sqrt, seep_dist_m) %>%
  pivot_wider(names_from = code,
              values_from = sqrt) %>%
  ungroup() %>%
  filter(site== "Varari")


c.matrix <- c.wide %>%
  select(-site, -TileID, -TopBottom, -seep_dist_m)
  
c.mds <- metaMDS(c.matrix, k = 2, trymax = 100)


c.scores <- as.data.frame(scores(c.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
c.scores$TileID <- c.wide$TileID
# create a column of site names, from the rownames of data.scores
c.scores$site <- c.wide$site  #  add the grp variable created earlier
c.scores$TopBottom <- c.wide$TopBottom

names(c.scores)[c(1, 2)] <- c("x", "y")
c.scores$z <- NA

c.sf <- ordisurf(c.mds ~ c.wide$seep_dist_m, plot = FALSE, scaling = 3)

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}

c.contour.vals <- extract.xyz(obj = c.sf)

c.plot.curve <- ggplot(data = c.contour.vals, aes(x, y, z = z)) + 
  stat_contour(aes(colour = ..level..)) + 
  coord_cartesian(xlim = c(-2, 2), ylim = c(-1, 1.5)) + 
  theme_bw() 

c.plot <- c.plot.curve +
  geom_text(data = c.scores, 
            aes(x = x, y = y, label = TileID), 
            colour = "black",
            size = 3,
            position=position_jitter(width=.000000000000000000001,height=1)) + 
  xlim(-1.2,1.5) +
  ylim(-1,1) +
  coord_equal() + 
  theme_bw() + 
  labs(x = "NMDS1", 
       y = "NMDS2") +
  theme(panel.border = element_rect(fill = NA), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) 



ggsave(here("AlgalCommunityComp", "outputs", "v.sf.png"))




# 2. % Cover NMDS  --------------------------------------------------------


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


ggsave(here("AlgalCommunityComp", "outputs", "pc_tb_mds.png"))

quartz()
