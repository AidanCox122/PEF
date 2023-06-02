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
Gl.kw %>% 
  ggplot(aes(x = factor(Zone), y = Density, color = factor(Zone))) +
  geom_jitter(alpha = 0.5) +
  geom_boxplot(color = "black", fill = NA) +
  annotate(geom = "text", x = 1, y = 100, label = "a", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 2, y = 100, label = "a", color = "grey22", fontface = "italic") +
  geom_rect(aes(xmin = 0.5, xmax = 2.45, ymin = 0, ymax = 120), fill = NA, color = "grey", alpha = 0.6) +
  annotate(geom = "text", x = 3, y = 100, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 4, y = 100, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 5, y = 100, label = "c", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 6, y = 160, label = "b", color = "grey22", fontface = "italic") +
  geom_rect(aes(xmin = 4.55, xmax = 5.45, ymin = 0, ymax = 170), fill = NA, color = "grey", alpha = 0.6) +
  ylab("Density (indv./km2)") +
  facet_grid(species~., switch = "y") +
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
CM.kw %>% 
  ggplot(aes(x = Zone, y = Density, color = Zone)) +
  geom_jitter(alpha = 0.5) +
  geom_boxplot(color = "black", fill = NA) +
  annotate(geom = "text", x = 1, y = 250, label = "a", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 2, y = 250, label = "a", color = "grey22", fontface = "italic") +
  geom_rect(aes(xmin = 0.5, xmax = 2.45, ymin = 0, ymax = 270), fill = NA, color = "grey", alpha = 0.6) +
  annotate(geom = "text", x = 3, y = 750, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 4, y = 750, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 5, y = 750, label = "b", color = "grey22", fontface = "italic") +
  geom_rect(aes(xmin = 2.55, xmax = 5.45, ymin = 0, ymax = 780), fill = NA, color = "grey", alpha = 0.6) +
  annotate(geom = "text", x = 6, y = 1000, label = "c", color = "grey22", fontface = "italic")+
  geom_rect(aes(xmin = 5.55, xmax = 6.5, ymin = 0, ymax = 1020), fill = NA, color = "grey", alpha = 0.6) +
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
ggplot(data = HSeal, aes(x = Zone, y = Density, color = Zone)) +
  geom_boxplot(color = "black", fill = NA) +
  geom_jitter(alpha = 0.5) +
  annotate(geom = "text", x = 1, y = 15, label = "a", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 2, y = 15, label = "a", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 3, y = 30, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 4, y = 20, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 5, y = 20, label = "b", color = "grey22", fontface = "italic") +
  geom_rect(aes(xmin = 2.55, xmax = 5.45, ymin = 0, ymax = 35), fill = NA, color = "grey", alpha = 0.6) +
  annotate(geom = "text", x = 6, y = 15, label = "a", color = "grey22", fontface = "italic") +
  ylab("Density (indv./km2)") +
  facet_grid(species~., switch = "y") +
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
ggplot(data = HPorp, aes(x = Zone, y = Density, color = Zone)) +
  geom_jitter(alpha = 0.5) +
  geom_boxplot(color = "black", fill = NA) +
  annotate(geom = "text", x = 1, y = 10, label = "a", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 2, y = 10, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 3, y = 20, label = "b", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 4, y = 20, label = "b", color = "grey22", fontface = "italic") +
  geom_rect(aes(xmin = 1.55, xmax = 4.45, ymin = 0, ymax = 25), fill = NA, color = "grey", alpha = 0.6) +
  annotate(geom = "text", x = 5, y = 10, label = "a", color = "grey22", fontface = "italic") +
  annotate(geom = "text", x = 6, y = 10, label = "a", color = "grey22", fontface = "italic") +
  ylab("Density (indv./km2)") +
  facet_grid(species~.) +
  theme_classic() + 
  theme(legend.position = "bottom")



