#PERMANOVA: to test global effects
#How community composition (species assemblages divided into life-history traits) varies with SO₂ and year?
#Are the species identities and abundances different between time periods?

library(vegan)
library(dplyr)
library(tidyr)
library(readxl)
library(permute)

format1 <- read_excel("Dataframe_final.xlsx")

### Step 1: Create Species Matrix ###
# Create a unique sample ID for grouping variable
format1 <- format1 %>%
  mutate(SampleID = paste(Year, Site, sep = "_"))

# Create species matrix: rows = samples, columns = species, values = total abundance
sp_matrix <- format1 %>%
  group_by(SampleID, Species) %>%
  summarise(Abundance = sum(Number), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0)

# Convert to data frame and set rownames
sp_df <- as.data.frame(sp_matrix)
rownames(sp_df) <- sp_df$SampleID
sp_df$SampleID <- NULL

### Step 2: Summarize environmental metadata ###
env_data <- format1 %>%
  group_by(SampleID) %>%
  summarise(
    Site = first(Site),
    Policy_period = first(Policy_period),
    Woody.species = first(Woody.species),
    Pollution = mean(Pollution, na.rm = TRUE),
    T = mean(T, na.rm = TRUE),
    Precipitation = mean(Precipitation, na.rm = TRUE),
    Wind = mean(Wind, na.rm = TRUE),
    .groups = "drop"
  )

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
