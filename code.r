#Load library
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

#Read files
environmental <- read_xlsx("~/Downloads/EM_Tropical_forest_succession(1).xlsx", sheet="Site_environmental_data")
vegetation <- read_xlsx("~/Downloads/EM_Tropical_forest_succession(1).xlsx", sheet="Site_vegetation_data")

#make columns factors
vegetation <- vegetation %>%
  mutate(
    Forest_type = factor(Forest_type, levels = c("Dry", "Wet")),
    Age = as.factor(Age),
    Site = as.factor(Site)
  )
###Create different datasets---

##merge datasets
data_full <- vegetation %>%
  left_join(environmental, by = c("Forest_type", "Site"))

glimpse(data_full)

##Pivot dataset with only growthform touches
growth_forms <- data_full %>%
  pivot_longer(cols = c(Herb, Grass, Fern, Vine, Liana, Shrub, Tree),
               names_to = "Growth_form", values_to = "Touches")

#Relative abundance per plot
growth_forms <- growth_forms %>%
  group_by(Forest_type, Age, Site) %>%
  mutate(Proportion = Touches / sum(Touches))

##Environmental data
env_numeric <- environmental %>%
  select(Soil_pH, Organic_matter, Nitrogen, Phosphorus, Potassium,
         starts_with("Surrounding_forest_cover"))

###Line graphs ----

##Effect of Soil Organic Matter on Species Richness 
ggplot(data_full, aes(x = Organic_matter, y = Species_richness, color = Forest_type)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Effect of Soil Organic Matter on Species Richness",
       x = "Soil Organic Matter (%)", y = "Species Richness") +
  theme_pubr()

##shannon diversity plotted over the years, grouped by forest type, 
ggplot(data_full, aes(x = Age, y = Shannon_diversity, group = Forest_type, color = Forest_type)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  labs(title = "Mean Shannon Diversity Trends",
       x = "Year", y = "Shannon Diversity Index") +
  theme_pubr()

##vegetation richness plotted over the years, grouped by forest type, 
ggplot(data_full, aes(x = Age, y = Species_richness, group = Forest_type, color = Forest_type)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  labs(title = "Mean vegetation richness",
       x = "Year", y = "vegetation richness") +
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

##changes in growth form composition over time 
ggplot(growth_forms, aes(x = Age, y = Proportion, fill = Growth_form)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Forest_type) +
  labs(title = "Changes in Growth Form Composition over Time",
       x = "Successional Age", y = "Proportion of Total Touches") +
  scale_fill_brewer(palette = "Set2") +
  theme_pubr()

#Herbs
ggplot(data = vegetation, mapping = aes(x = Age, y = Herb, colour = Forest_type)) + 
  geom_boxplot()+ 
  labs(y = "Number of leaf touches", x = "Years since abandonment")+
  ggtitle("Herbs over time") +
  theme_pubr()

#Grass
ggplot(data = vegetation, mapping = aes(x = Age, y = Grass, colour = Forest_type)) + 
  geom_boxplot()+ 
  labs(y = "Number of leaf touches", x = "Years since abandonment")+
  ggtitle("Grass over time") +
  theme_pubr()

#Ferns
ggplot(data = vegetation, mapping = aes(x = Age, y = Fern, colour = Forest_type)) + 
  geom_boxplot() + 
  labs(y = "Number of leaf touches", x = "Years since abandonment")+
  ggtitle("Ferns over time") +
  theme_pubr()

#Vines
ggplot(data = vegetation, mapping = aes(x = Age, y = Vine, colour = Forest_type)) + 
  geom_boxplot()+ 
  labs(y = "Number of leaf touches", x = "Years since abandonment")+
  ggtitle("Vines over time") +
  theme_pubr()

#Lianas
ggplot(data = vegetation, mapping = aes(x = Age, y = Liana, colour = Forest_type)) + 
  geom_boxplot()+ 
  labs(y = "Number of leaf touches", x = "Years since abandonment")+
  ggtitle("Lianas over time") +
  theme_pubr()

#Shrubs
ggplot(data = vegetation, mapping = aes(x = Age, y = Shrub, colour = Forest_type)) + 
  geom_boxplot()+ 
  labs(y = "Number of leaf touches", x = "Years since abandonment")+
  ggtitle("Shrubs over time") +
  theme_pubr()

#Trees
ggplot(data = vegetation, mapping = aes(x = Age, y = Tree, colour = Forest_type)) + 
  geom_boxplot()+ 
  labs(y = "Number of leaf touches", x = "Years since abandonment")+
  ggtitle("Trees over time") +
  theme_pubr()

###Boxplots ----
#Soil pH
ggplot(environmental, aes(x = Forest_type, y = Soil_pH, fill = Forest_type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +   # boxplot
  geom_jitter(width = 0.15, alpha = 0.7) +          # add points
  labs(x = "Forest type", y = "Soil pH") +
  theme_pubr() +
  theme(legend.position = "none")

#Organic matter
ggplot(environmental, aes(x = Forest_type, y = Organic_matter, fill = Forest_type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +   # boxplot
  geom_jitter(width = 0.15, alpha = 0.7) +          # add points
  labs(x = "Forest type", y = "Organic matter in soil (%)") +
  theme_pubr() +
  theme(legend.position = "none")

#Nitrogen
ggplot(environmental, aes(x = Forest_type, y = Nitrogen, fill = Forest_type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +   # boxplot
  geom_jitter(width = 0.15, alpha = 0.7) +          # add points
  labs(x = "Forest type", y = "Nitrogen (mg/g)") +
  theme_pubr() +
  theme(legend.position = "none")

#Phosphorus
ggplot(environmental, aes(x = Forest_type, y = Phosphorus, fill = Forest_type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +   # boxplot
  geom_jitter(width = 0.15, alpha = 0.7) +          # add points
  labs(x = "Forest type", y = "Phosphorus (μg/g)") +
  theme_pubr() +
  theme(legend.position = "none")

#Potassium
ggplot(environmental, aes(x = Forest_type, y = Potassium, fill = Forest_type)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +   # boxplot
  geom_jitter(width = 0.15, alpha = 0.7) +          # add points
  labs(x = "Forest type", y = "Potassium (meq/100g)") +
  theme_minimal() +
  theme(legend.position = "none")


##correlation plot
corrplot(cor(env_numeric, use = "complete.obs"), method = "color", tl.col = "black")





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

# 4.1 Species richness: check distribution & overdispersion for Poisson
hist(data_full$Species_richness, main = "Species richness distribution")
# Fit Poisson GLMM (simple)
pois_glmm <- glmer(Species_richness ~ Age * Forest_type + (1|Site), data = data_full, family = poisson)

pois_glmm

# Fit Negative Binomial GLMM (preferred if overdispersion > ~1.3)
nb_glmm <- glmer.nb(Species_richness ~ Age * Forest_type + (1|Site), data = data_full)

# Summarize NB GLMM
summary(nb_glmm)
# Likelihood ratio tests for fixed effects (Type III)
anova(nb_glmm)











