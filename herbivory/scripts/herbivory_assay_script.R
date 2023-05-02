## Script for analysis of herbivory assays 

# loading libraries
library(tidyverse)
library(here)
library(vegan)
library(car)
library(lme4) # mixed effects models
library(lmerTest)
library(performance)
library(multcomp)

colors <- c("#9bb076", '#7aa3c9')

# read in file
assay <- read_csv(here('herbivory', 'data', 'herbivory_assay.csv')) %>%
  # calculate amount consumed
  mutate(diff = Pre_algae-Post_algae ) %>%
  # select working columns
  dplyr:: select(Date, Sample, Group,HighLow, Source, Location, Place, diff) %>%
  # remove rockwall sites
  filter(Group!="RockwallSGD")%>%
  # pivot wider to control for flow
  pivot_wider(names_from = 'Place', values_from = 'diff') %>%
  # subtract caged controls 
  mutate(consumed = open - cage) %>%
  # correcting for water weight with minimum weight value
  mutate(consumed = consumed - min(consumed))%>%
  # change characters to factors
  mutate(Date = as.factor(Date),
        Group = as.factor(Group),
         Source = as.factor(Source),
         Sample = as.factor (Sample),
         HighLow = as.factor(HighLow),
         Location = as.factor(Location)) %>%
  # rename factors
  mutate(Group = recode_factor(Group, 
                        ReefAmbient = 'Low SGD Reef',
                        ReefSGD = 'High SGD Reef',
                        SandAmbient = 'Low SGD Sand Patch',
                        SandSGD = 'High SGD Sand Patch')) %>%
  mutate(Source = recode_factor(Source,
                         SGD = 'High SGD',
                         ambient = 'Low SGD')) %>%
 
  # rename columns for clarity
  rename(Habitat = Location,
         Algal_Source = Source,
         High_Low_SGD_Habitat = HighLow) 


# visualize data

assay_plot_location <- assay %>%
  filter(Habitat == "Reef") %>% 
  ggplot(aes(x = Algal_Source, y=consumed)) +
  geom_boxplot(aes(fill = Algal_Source), show.legend = FALSE,alpha = .6) +
  theme(legend.position ="none") +
  scale_fill_manual(values=colors) +
  #facet_wrap(~Group, scales = "free" )+
  facet_wrap(~Group) +
  theme_bw() +
  labs(x = "Algae Collection Location",
       y = "Herbivory rate (g / 6 hours)",
       title = "") 

ggsave(here('herbivory', 'outputs', 'herbivory_reef.jpg'))

  quartz()

# 3 way ANOVA (date as random effect)
mod1<-lmer(consumed ~ Algal_Source*High_Low_SGD_Habitat*Habitat +(1|Date), data = assay)
summary(mod1)
anova(mod1)

check_model(mod1)  

aov.mod1<-aov(consumed ~ Algal_Source*High_Low_SGD_Habitat*Habitat, data = assay)
TukeyHSD(aov.mod1)
summary(aov.mod1)
anova(aov.mod1)

modelparm(mod1)


#summary
assay %>% group_by(High_Low_SGD_Habitat, Habitat, Algal_Source) %>% summarise(mean = mean(consumed))

## plot residual vs fitted

mod1_frame <- dplyr::select(assay,Date,Algal_Source,Habitat, High_Low_SGD_Habitat)  %>%
              dplyr::mutate(fits=fitted(mod1),
                    resids=resid(mod1),
                    sresids=rstudent(mod1))

    #by date
    mod1_plot_habitat <- ggplot(data=mod1_frame, mapping=aes(x=fits,y=resids)) +
                geom_point(aes(color=factor(Date))) +
                geom_hline(yintercept=0,linetype="dashed")

    #by habitat
    mod1_plot <- ggplot(data=mod1_frame, mapping=aes(x=fits,y=resids)) +
                geom_point(aes(color=factor(Habitat), shape = Algal_Source)) +
                geom_hline(yintercept=0,linetype="dashed")

## check with sqrt transformed data
    sqrt_assay <-  assay %>%
                  mutate(consumed_sqrt = sqrt(consumed)) 
      
    mod1_sqrt<-lmer(consumed_sqrt ~ Algal_Source*Habitat*High_Low_SGD_Habitat +(1|Date), data = sqrt_assay)
    
    
    check_model(mod1_sqrt)  
                  
    summary(mod1_sqrt)
    anova(mod1_sqrt)             
    
   mod <-  aov(consumed ~ Algal_Source*High_Low_SGD_Habitat, data = assay)
    summary(mod)
    