# -------------------------------------------------------------------------
#PERMANOVA: to test global effects
#How community composition (species assemblages divided into life-history traits) varies with SO₂ and year?
#Are the species identities and abundances different between time periods?
# -------------------------------------------------------------------------

library(vegan)
library(dplyr)
library(tidyr)
library(readxl)
library(permute)

### Step 1: Species matrix ###

# Read species matrix: rows = samples, columns = species, values = total abundance
sp_matrix <- read_excel("Dataframe_macrotemporal.xlsx", sheet = "sp_matrix")

### Step 2: Read environmental metadata ###
env_data <- read_excel("Dataframe_macrotemporal.xlsx", sheet = "env_data")

# Ensure grouping variables are factors
env_data$Site <- as.factor(env_data$Site)
env_data$Policy_period <- as.factor(env_data$Policy_period)
env_data$Woody.species <- as.factor(env_data$Woody.species)

### Step 3: Multivariate Dispersion (PERMDISP) ###
# Calculate Bray-Curtis dissimilarity
bray_dist <- vegdist(sqrt(sp_df), method = "bray")

# Check for homogeneity of multivariate dispersions
dispersion <- betadisper(bray_dist, env_data$Policy_period)
anova(dispersion)
pairwise.perm.test <- TukeyHSD(dispersion)
print(pairwise.perm.test)

### Step 4: PERMANOVA ###
# Model testing marginal effects of predictors, stratified by nine sites (four tree species)
# Note: Temperature was not used to copy GAMM model for consistency.

perm_control <- how(plots = Plots(strata = env_data$Site), nperm = 9999)
set.seed(123)
adonis_full <- adonis2(
  bray_dist ~ Pollution + Policy_period + Precipitation + Wind, 
  data = env_data, 
  permutations = perm_control,
  by = "margin"
)
print(adonis_full)
