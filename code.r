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

###Functions ----
###Function tpaired.randtest ----
tpaired.krandtest <- function(mat1, mat2, nperm=999, list.all=FALSE)
{
  epsilon <- .Machine$double.eps
  n = nrow(mat1)
  p = ncol(mat1)
  if(nrow(mat2) != n) stop("Unequal number of rows")
  if(ncol(mat2) != p) stop("Unequal number of species")
  
  # Check species names in files mat1 and mat2
  sp.diff <- which(colnames(mat1) != colnames(mat2))
  if(length(sp.diff)>0) 
    cat("Warning – The following species names differ between mat1 & mat2:\n",sp.diff,"\n\n")
  
  # Select the species that have variances > epsilon. Discard the other species.
  comm.gr = rbind(mat1, mat2)
  tmp = apply(comm.gr, 2, var)
  sel.sp = which(tmp > epsilon)
  pp <- length(sel.sp)
  cat((p-pp),"species were eliminated because they did not vary in the combined data set\n")
  cat(pp,"species retained:\n")
  mat1.sel = mat1[,sel.sp]
  mat2.sel = mat2[,sel.sp]
  
  # 
  res = matrix(0,pp,4)
  sp.names <- rownames(res) <- colnames(mat1.sel)
  
  #
  Star <- NA
  Keep <- NA
  Reject <- NA
  for(j in 1:pp) {
    vec1 = mat1[,j]
    vec2 = mat2[,j]
    # if(any(abs(vec1-vec2) > 0)) Keep = c(Keep,j) else Reject = c(Reject,j)
    cat(j," ")      # Print ID of species being processed
    tmp <- t.test(vec1, vec2, paired=TRUE)
    if(tmp$estimate == 0) {
      Reject = c(Reject,j) 
      res[j,3:4] = c(-999,-999)
      Star <- c(Star, " ")
    } else {
      Keep = c(Keep,j) 
      if(sign(tmp$estimate)<0) tail="less" else if(sign(tmp$estimate)>=0) tail="greater"
      t.res <- tpaired.randtest(vec1, vec2, nperm=nperm, alternative=tail, silent=TRUE)
      res[j,1:4] <- c(t.res$estim, t.res$t.ref, t.res$p.param, t.res$p.perm)
      Star <- c(Star, ifelse(t.res$p.perm>0.05, " ","*"))
    }
  }
  cat("\n\n")
  Keep <- Keep[-1]
  Reject <- Reject[-1]
  tmp <- length(Reject)
  if(tmp>0) cat(tmp,"species not tested because t.stat = 0. See 'No_test' output list\n\n")
  res = data.frame(res,Star[-1],sign(res[,1]))
  colnames(res) = c("mean(T1-T2)","t.stat","p.param","p.perm","p<=0.05","Sign(T1-T2)")
  if(list.all) 
    out <- list(t.tests=res, Tested=sp.names[Keep], No_test=sp.names[Reject])
  else out <- list(t.tests=res[Keep,],Tested=sp.names[Keep],No_test=sp.names[Reject])
  out
}

###data preperation ----
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

##Matrix with growth form types
vegetation$Site <- paste(vegetation$Site, vegetation$Age)
matrix <- vegetation %>% column_to_rownames("Site")
matrix <- select(matrix, "Herb", "Grass", "Fern", "Vine", "Liana", "Shrub", "Tree")
#assign rows to their respective years
year1 <- c(1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,61,64,67,70,73,76,79,82,85,88)
year2 <- c(2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89)
year3 <- c(3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90)

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

###Other plots ----

##correlation plot
corrplot(cor(env_numeric, use = "complete.obs"), method = "color", tl.col = "black")

###RQ1---
#How does vegetation structure, diversity, and composition change in the early succession of dry and wet forests?

##differences in Shannon diversity over the years between the forest types
lmm_shannon <- lmer(Shannon_diversity ~ Age * Forest_type + (1|Site), data = vegetation)
# Compute EMMs
em_shannon <- emmeans(lmm_shannon, ~ Age | Forest_type)
# Pairwise comparisons with Tukey adjustment
pair_shannon <- contrast(em_shannon, method = "pairwise", adjust = "tukey")
# View results
summary(pair_shannon)
anova(lmm_shannon)

##differences in Leaf area index (coverage) over the years between the forest types
lmm_lai <- lmer(Leaf_area_index ~ Age * Forest_type + (1|Site), data = vegetation)
# Compute EMMs
em_lai <- emmeans(lmm_lai, ~ Age | Forest_type)
# Pairwise comparisons with Tukey adjustment
pair_lai <- contrast(em_lai, method = "pairwise", adjust = "tukey")
# View results
summary(pair_lai)
anova(lmm_lai)

##differences in growth type coverage between the years
#Temporal Beta-diversity Indices (TBI) with Bray-Curtis distances for the different years
TBI1 <-TBI(matrix[year1, ], matrix[year3, ], method = "%difference")
TBI2 <-TBI(matrix[year1, ], matrix[year2, ], method = "%difference")
TBI3 <-TBI(matrix[year2, ], matrix[year3, ], method = "%difference")

#Create dataset with TBI data
Sites <- as.data.frame(c(01:30))
Sites$Forest_type <- c(rep("Dry", 16), rep("Wet", 14))
Sites$Number <- c(01:16, 01:14)
Sites$Number <- sprintf("%02d", Sites$Number)
Sites$Site <- paste(Sites$Forest_type, Sites$Number, sep = "")
Sites$TBI1_3 <- TBI1$TBI
Sites$TBI1_2 <- TBI2$TBI
Sites$TBI2_3 <- TBI3$TBI

##T-test for differences between TBI's between the years
#year 1-2
tapply(X = Sites$TBI1_2, INDEX = Sites$Forest_type, FUN = shapiro.test)
leveneTest(TBI1_2~Forest_type, data = Sites, center = median)
t.test(TBI1_2 ~ Forest_type, data = Sites, var.equal = TRUE)

#year 1-3
tapply(X = Sites$TBI1_3, INDEX = Sites$Forest_type, FUN = shapiro.test)
leveneTest(TBI1_3~Forest_type, data = Sites, center = median)
t.test(TBI1_3 ~ Forest_type, data = Sites, var.equal = TRUE)

#year 2-3
tapply(X = Sites$TBI2_3, INDEX = Sites$Forest_type, FUN = shapiro.test)
leveneTest(TBI2_3~Forest_type, data = Sites, center = median)
t.test(TBI2_3 ~ Forest_type, data = Sites, var.equal = TRUE)

##Create dataset with TBI data and environmental data
distance <- merge(Sites, environmental, by = c("Site"))
distance$Site <- as.factor(distance$Site)
distance$Forest_type.x <- as.factor(distance$Forest_type.x)

##GLM's for TBI's of the different years against environmental data
#year 1-2
glm1 <- glm(log(TBI1_2)~ Forest_type.x + scale(Surrounding_forest_cover_200) + scale(Organic_matter) +scale(Nitrogen) + scale(Phosphorus) + scale(Potassium), data= distance)
shapiro.test(glm1$residuals)
vif(glm1)
summary(glm1)            
glm1 <- step(glm1, direction = "backward")
summary(glm1)

#year 1-3
glm2 <- glm(log(TBI1_3)~ Forest_type.x + scale(Surrounding_forest_cover_800) + scale(Organic_matter) +scale(Nitrogen) + scale(Phosphorus) + scale(Potassium), data= distance)
shapiro.test(glm2$residuals)
vif(glm2)
summary(glm2)            
glm2 <- step(glm2, direction = "backward")
summary(glm2)

#year 2-3
glm3 <- glm(log(TBI2_3)~ Forest_type.x + scale(Surrounding_forest_cover_800) + scale(Organic_matter) +scale(Nitrogen) + scale(Phosphorus) + scale(Potassium), data= distance)
shapiro.test(glm3$residuals)
vif(glm3)
summary(glm3)            
glm3 <- step(glm3, direction = "backward")
summary(glm3)
