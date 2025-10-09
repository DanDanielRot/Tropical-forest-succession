library(lmerTest)
library(emmeans)
library(ggpubr)
library(fitdistrplus)
library(ggh4x)
library(car)
library(readxl)
library(tidyverse)
library(vegan)
library(cowplot)
library(goeveg)
library(corrplot)
library(ggpubr)


environmental <- read_xlsx("~/Downloads/EM_Tropical_forest_succession(1).xlsx", sheet="Site_environmental_data")

vegetation <- read_xlsx("~/Downloads/EM_Tropical_forest_succession(1).xlsx", sheet="Site_vegetation_data")


vegetation <- vegetation %>%
  mutate(
    Forest_type = factor(Forest_type, levels = c("Dry", "Wet")),
    Age = as.factor(Age),
    Site = as.factor(Site)
  )

data_full <- vegetation %>%
  left_join(environmental, by = c("Forest_type", "Site"))

glimpse(data_full)


##shannon diversity plotted over the years, grouped by forest type, 
ggplot(data_full, aes(x = Age, y = Shannon_diversity, group = Forest_type, color = Forest_type)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  labs(title = "Mean Shannon Diversity Trends",
       x = "Year", y = "Shannon Diversity Index") +
  theme_pubr()


##species richness plotted over the years, grouped by forest type, 
ggplot(data_full, aes(x = Age, y = Species_richness, group = Forest_type, color = Forest_type)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  labs(title = "Mean species richness",
       x = "Year", y = "species richness") +
  theme_pubr()


##leaf area index plotted over the years, grouped by forest type, (boxplot)
ggplot(data_full, aes(x = Age, y = Leaf_area_index, color = Forest_type)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.7) +
  labs(title = "Leaf Area Index by Age and Forest Type",
       y = "Leaf Area Index", x = "Year") +
  theme_pubr()



##total touch plotted over the years, grouped by forest type, (boxplot)
ggplot(data_full, aes(x = Age, y = Total_touch, color = Forest_type)) +
  geom_boxplot() +
  labs(title = "Vegetation Density (Total Touch) across Succession",
       y = "Total Touch Count") +
  theme_pubr()


##Transpose dataset
growth_forms <- data_full %>%
  pivot_longer(cols = c(Herb, Grass, Fern, Vine, Liana, Shrub, Tree),
               names_to = "Growth_form", values_to = "Touches")

# Relative abundance per plot
growth_forms <- growth_forms %>%
  group_by(Forest_type, Age, Site) %>%
  mutate(Proportion = Touches / sum(Touches))



##changes in growth form composition over time 
ggplot(growth_forms, aes(x = Age, y = Proportion, fill = Growth_form)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Forest_type) +
  labs(title = "Changes in Growth Form Composition over Time",
       x = "Successional Age", y = "Proportion of Total Touches") +
  scale_fill_brewer(palette = "Set2") 


##
env_numeric <- forest_environmental %>%
  select(Soil_pH, Organic_matter, Nitrogen, Phosphorus, Potassium,
         starts_with("Surrounding_forest_cover"))

##correlation plot
corrplot(cor(env_numeric, use = "complete.obs"), method = "color", tl.col = "black")


##Effect of Soil Organic Matter on Species Richness 
ggplot(data_full, aes(x = Organic_matter, y = Species_richness, color = Forest_type)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Effect of Soil Organic Matter on Species Richness",
       x = "Soil Organic Matter (%)", y = "Species Richness")


############################### linear mixed models

###shannon diversity 

lmm_shannon <- lmer(Shannon_diversity ~ Age * Forest_type + (1|Site), data = veg)

# Compute EMMs
em_shannon <- emmeans(lmm_shannon, ~ Age | Forest_type)

# Pairwise comparisons with Tukey adjustment
pair_shannon <- contrast(em_shannon, method = "pairwise", adjust = "tukey")

# View results
summary(pair_shannon)

anova(lmm_shannon)



###total touch lmm

lmm_touch <- lmer(Total_touch ~ Age * Forest_type + (1|Site), data = veg)

# Compute EMMs
em_touch <- emmeans(lmm_touch, ~ Age | Forest_type)

# Pairwise comparisons with Tukey adjustment
pair_touch <- contrast(em_touch, method = "pairwise", adjust = "tukey")

# View results
summary(pair_touch)

anova(lmm_touch)




###---lmer leaf area index

lmm_lai <- lmer(Leaf_area_index ~ Age * Forest_type + (1|Site), data = veg)

# Compute EMMs
em_lai <- emmeans(lmm_lai, ~ Age | Forest_type)

# Pairwise comparisons with Tukey adjustment
pair_lai <- contrast(em_lai, method = "pairwise", adjust = "tukey")

# View results
summary(pair_lai)

anova(lmm_lai)













