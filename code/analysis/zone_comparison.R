## ZONE COMPARISONS ##

# The first step of my analysis should be to compare the differences
# in density and frequency of presence between the different transect zones
# for each species to determine if they are significantly different.
# Because our data do not meet the assumption of normality for an ANOVA,
# I will apply a Kruskall-Wallis test instead.


# setup -------------------------------------------------------------------

library(tidyverse)

mbm_data <- 
  read_csv('data/clean/mbm_master.csv') %>% 
  rename(species = `Species_code`)


# species comparisons -----------------------------------------------------

## glaucous gull -----------------------------------------------------------

GL.kw <- mbm_data %>% filter(species == "GL")

# run an Kruskall-Wallis test on the data
print(kruskal.test(Density ~ Zone, data = GL.kw)) # there are highly significant differences between the zones
# post-hoc testing
pairwise.wilcox.test(GL.kw$Density, GL.kw$Zone,
                     p.adjust.method = "bonferroni" # apply Bonferroni's correction
)

# visualize the data
GL_zone_comp <- 
  GL.kw %>% 
  mutate(species = 'Glaucous-winged Gull') %>% 
  ggplot(aes(x = factor(Zone), y = Density, color = factor(Zone))) +
  geom_jitter(alpha = 0.5) +
  geom_boxplot(color = "black", fill = NA) +
  annotate(geom = "text", x = 1, y = 80, label = "a", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 2, y = 80, label = "a", color = "grey22", fontface = "italic") +
  geom_rect(aes(xmin = 0.5, xmax = 2.45, ymin = 0, ymax = 90), fill = NA, color = "grey", alpha = 0.6) +
  annotate(geom = "text", x = 3, y = 80, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 4, y = 80, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 5, y = 80, label = "c", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 6, y = 80, label = "b", color = "grey22", fontface = "italic") +
  geom_rect(aes(xmin = 4.55, xmax = 5.45, ymin = 0, ymax = 90), fill = NA, color = "grey", alpha = 0.6) +
  coord_flip() +
  xlab('Transect Zone') +
  ylab("Density (indv./km2)") +
  facet_grid(species~.) +
  theme_classic() + 
  theme(legend.position = "bottom") 


## common murre ------------------------------------------------------------

CM.kw <- 
  mbm_data %>% filter(species == "CoMu") %>% 
  mutate(Zone = factor(Zone))

# run an anova on the data
kruskal.test(Density ~ Zone, data = CM.kw) # there are highly significant differences between the zones
# post-hoc testing
pairwise.wilcox.test(CM.kw$Density, CM.kw$Zone,
                     p.adjust.method = "bonferroni" # apply Bonferroni's correction
)

# visualize the data
CM_zone_comp <- 
  CM.kw %>%
  mutate(species = 'Common Murre') %>% 
  ggplot(aes(x = Zone, y = Density, color = Zone)) +
  geom_jitter(alpha = 0.5) +
  geom_boxplot(color = "black", fill = NA) +
  annotate(geom = "text", x = 1, y = 170, label = "a", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 2, y = 170, label = "a", color = "grey22", fontface = "italic") +
  geom_rect(aes(xmin = 0.5, xmax = 2.45, ymin = 0, ymax = 180), fill = NA, color = "grey", alpha = 0.6) +
  annotate(geom = "text", x = 3, y = 170, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 4, y = 170, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 5, y = 170, label = "b", color = "grey22", fontface = "italic") +
  # geom_rect(aes(xmin = 2.55, xmax = 5.45, ymin = 0, ymax = 180), fill = NA, color = "grey", alpha = 0.6) +
  annotate(geom = "text", x = 6, y = 170, label = "c", color = "grey22", fontface = "italic")+
  geom_rect(aes(xmin = 5.55, xmax = 6.5, ymin = 0, ymax = 180), fill = NA, color = "grey", alpha = 0.6) +
  coord_flip() +
  xlab('Transect Zone') +
  ylab("Density (indv./km2)") +
  facet_grid(species~.) +
  theme_classic() + 
  theme(legend.position = "bottom")

## Harbor seal -----------------------------------------------------------

HSeal <- 
  mbm_data %>% 
  filter(species == "HSeal") %>% 
  mutate(Zone = factor(Zone))

# run an anova on the data 
kruskal.test(Density ~ Zone, data = HSeal) # there are highly significant differences between the zones
# post-hoc testing
pairwise.wilcox.test(HSeal$Density, HSeal$Zone,
                     p.adjust.method = "bonferroni" # apply Bonferroni's correction
)

# visualize the data
HS_zone_comp <- 
  HSeal %>% 
  mutate(species = 'Harbor Seal') %>% 
  ggplot(aes(x = Zone, y = Density, color = Zone)) +
  geom_jitter(alpha = 0.5) +
  geom_boxplot(color = "black", fill = NA) +
  annotate(geom = "text", x = 1, y = 15, label = "a", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 2, y = 15, label = "a", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 3, y = 15, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 4, y = 15, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 5, y = 15, label = "b", color = "grey22", fontface = "italic") +
  geom_rect(aes(xmin = 2.55, xmax = 5.45, ymin = 0, ymax = 20), fill = NA, color = "grey", alpha = 0.6) +
  annotate(geom = "text", x = 6, y = 15, label = "a", color = "grey22", fontface = "italic") +
  coord_flip() +
  xlab('Transect Zone') +
  ylab("Density (indv./km2)") +
  facet_grid(species~.) +
  theme_classic() + 
  theme(legend.position = "bottom") 

# harbor porpoise ---------------------------------------------------------

HPorp <- mbm_data %>% filter(species == "HPorp") %>% 
  mutate(Zone = factor(Zone))

# run an anova on the data
kruskal.test(Density ~ Zone, data = HPorp) # there are highly significant differences between the zones

# post-hoc testing
pairwise.wilcox.test(HPorp$Density, HPorp$Zone,
                     p.adjust.method = "bonferroni" # apply Bonferroni's correction
)
# visualize the data
HP_zone_comp <-
  HPorp %>% 
  mutate(species = 'Harbor Porpoise') %>% 
  ggplot(aes(x = Zone, y = Density, color = Zone)) +
  geom_jitter(alpha = 0.5) +
  geom_boxplot(color = "black", fill = NA) +
  annotate(geom = "text", x = 1, y = 06, label = "a", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 2, y = 06, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 3, y = 06, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 4, y = 06, label = "b", color = "grey22", fontface = "italic") +
  geom_rect(aes(xmin = 1.55, xmax = 4.45, ymin = 0, ymax = 08), fill = NA, color = "grey", alpha = 0.6) +
  annotate(geom = "text", x = 5, y = 06, label = "a", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 6, y = 06, label = "a", color = "grey22", fontface = "italic") +
  coord_flip() +
  xlab('Transect Zone') +
  ylab("Density (indv./km2)") +
  facet_grid(species~.) +
  theme_classic() + 
  theme(legend.position = "bottom")


# Save figure 2 -----------------------------------------------------------

zone_comp_plots <-
  list(
    GL_zone_comp,
    CM_zone_comp,
    HS_zone_comp,
    HP_zone_comp) %>% 
  set_names(
    c(
      'GL_zone_comp',
      'CM_zone_comp',
      'HS_zone_comp',
      'HP_zone_comp'))

for (x in names(zone_comp_plots)) {
  fname <-
    paste0('products/figure2/raw/', (x), '.jpg')
  ggsave(fname,
         zone_comp_plots[[x]],
         width = 4.48,
         height = 3.36,
         units = 'in')
  print(paste('Done with', x))}


# year ANOVA --------------------------------------------------------------


## glaucous-winged gull ----------------------------------------------------
GL.kw <- 
  GL.kw %>%
  mutate(year = lubridate::year(Date)) # %>% 
  filter(year >= 2017)

# this chunk will compare differences in species abundance between years
print(kruskal.test(Density ~ year, data = GL.kw)) # there are highly significant differences between the zones
# post-hoc testing
pairwise.wilcox.test(GL.kw$Density, GL.kw$year,
                     p.adjust.method = "bonferroni" # apply Bonferroni's correction
)

## common murre ------------------------------------------------------------

CM.kw <- 
  CM.kw %>% 
  mutate(year = lubridate::year(Date)) # %>% 
  filter(year >= 2017)

# this chunk will compare differences in species abundance between years
print(kruskal.test(Density ~ year, data = CM.kw)) # there are highly significant differences between the zones
# post-hoc testing
# no significant difference between years


## harbor seal -------------------------------------------------------------
HSeal <- 
  HSeal %>%
  mutate(year = lubridate::year(Date)) #%>% 
  filter(year >= 2017)

# this chunk will compare differences in species abundance between years
print(kruskal.test(Density ~ year, data = HSeal)) # there are highly significant differences between the zones
# post-hoc testing
pairwise.wilcox.test(HSeal$Density, HSeal$year,
                     p.adjust.method = "bonferroni" # apply Bonferroni's correction
)


## harbor porpoise ---------------------------------------------------------
HPorp <- 
  HPorp %>%
  mutate(year = lubridate::year(Date)) #%>% 
  filter(year >= 2017)

# this chunk will compare differences in species abundance between years
print(kruskal.test(Density ~ year, data = HPorp)) # there are highly significant differences between the zones
# post-hoc testing
pairwise.wilcox.test(HPorp$Density, HPorp$year,
                     p.adjust.method = "bonferroni" # apply Bonferroni's correction
)


# week ANOVA --------------------------------------------------------------

## glaucous-winged gull ---------------------------------------------------------

GL.kw %>% 
  mutate(week = 
           lubridate::week(Date)) %>% 
  # this chunk will compare differences in species abundance between week
  kruskal.test(Density ~ week, data = .) %>% 
  print() # p = 0.159

## common murre ---------------------------------------------------------

CM.kw %>% 
  mutate(week = 
           lubridate::week(Date)) %>% 
  # this chunk will compare differences in species abundance between week
  kruskal.test(Density ~ week, data = .) %>% 
  print() # p = 0.000759

CM.kw <-
  CM.kw %>% 
  mutate(week = 
           lubridate::week(Date))

pairwise.wilcox.test(CM.kw$Density, CM.kw$week,
                     p.adjust.method = "bonferroni" # apply Bonferroni's correction
)

## harbor seal ---------------------------------------------------------

HSeal %>% 
  mutate(week = 
           lubridate::week(Date)) %>% 
  # this chunk will compare differences in species abundance between week
  kruskal.test(Density ~ week, data = .) %>% 
  print() # p = 0.0005277

## harbor porpoise ---------------------------------------------------------

HPorp %>% 
  mutate(week = 
           lubridate::week(Date)) %>% 
  # this chunk will compare differences in species abundance between weeks
  kruskal.test(Density ~ week, data = .) %>% 
  print() # there are no significant differences between the zones
0

