# -------------------------------------------------------------------------
#CWMs and FD Rao's Q: to test Time series (Time.period) with Site (n=9)
#Will environmental filtering reduce or increase increase overall functional 
# diversity with driving severe directional shifts in isolated traits?
# -------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readxl)
library(writexl)
library(scales)
library(mgcv)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(ggplot2)

### Read dataset ###
cwm_df <- read_excel("Dataframe_microtemporal.xlsx",sheet = "Sheet1")

### Step 1: Convert grouping variables to factors for mgcv ###
cwm_df$Woody.species <- as.factor(cwm_df$Woody.species)
cwm_df$Year_factor   <- as.factor(cwm_df$Year)
cwm_df$Month         <- as.numeric(cwm_df$Month)
cwm_df$Site          <- as.numeric(cwm_df$Site)

### Step 2: Fit GAMM Models for each trait ###
# Note: Temperate was dropped, because it was confounded with Month variable.
# Observed collinearity, model fit was not well behaving.

# Trophic strategy
mod_diet <- gam(
  CWM_Diet ~ Woody.species + 
    s(Pollution, k = 15) +       
    s(Precipitation, k = 15) +   
    s(Wind, k = 15) +
    s(Month, k = 5, bs = "tp") +
    s(Site, bs = "re") + 
    s(Year_factor, bs = "re"),
  data = cwm_df,
  family = quasibinomial(link = "logit"),
  method = "REML",
  knots  = list(Month = c(3, 11))
)
summary(mod_diet)

# Dispersal ability
mod_wing <- gam(
  CWM_Wing ~ Woody.species + 
    s(Pollution, k = 15) +       
    s(Precipitation, k = 15) +   
    s(Wind, k = 15) + 
    s(Month, k = 5, bs = "tp") +
    s(Site, bs = "re") + 
    s(Year_factor, bs = "re"),
  data = cwm_df,
  family = quasibinomial(link = "logit"),
  method = "REML",
  knots  = list(Month = c(3, 11))
)
summary(mod_wing)

# Body size
mod_size <- gam(
  CWM_Size ~ Woody.species + 
    s(Pollution, k = 15) +       
    s(Precipitation, k = 15) +   
    s(Wind, k = 15) +
    s(Month, k = 5, bs = "tp") +
    s(Site, bs = "re") + 
    s(Year_factor, bs = "re"),
  data = cwm_df,
  family = gaussian(link = "log"),
  method = "REML",
  knots  = list(Month = c(3, 11))
)
summary(mod_size)

# Check for residuals
par(mfrow = c(2, 2))
gam.check(mod_wing)
concurvity(mod_wing, full = TRUE)
gratia::draw(mod_wing)

### Check loadings and collinearity of CWM traits ###
### Step 1: PCA of robust traits ###

cwm_mat_robust <- cwm_df %>%
  dplyr::select(
    `Dispersal ability` = CWM_Wing,
    `Trophic strategy`  = CWM_Diet,
    `Body size`         = CWM_Size
  ) %>%
  na.omit()

res_pca <- PCA(cwm_mat_robust, scale.unit = TRUE, graph = FALSE)

tiff("PCA_robust.tiff", width = 5, height = 5, units = "in", res = 300, compression = "lzw") 
fviz_pca_var(res_pca, repel = TRUE, col.var = "black")
dev.off()

### Step 2: Collinearity of CWM traits ###
# Spearman correlation matrix and visualization

cwm_traits <- cwm_df[, c("CWM_Diet", "CWM_Wing", "CWM_Size")]
cor_mat <- cor(cwm_traits, method = "spearman", use = "pairwise.complete.obs")
colnames(cor_mat) <- rownames(cor_mat) <- c(
  "Trophic strategy",
  "Dispersal ability",
  "Body size"
)
corrplot(
  cor_mat,
  method = "number",
  type = "upper",
  tl.col = "black",
  tl.cex = 1.1,
  tl.srt = 45,
  addCoef.col = "black",
  number.cex = 1.1
)

####################################################################################################################

### Functional Diversity Rao's Entropy ###

library(FD)
library(vegan)
library(mgcv)
library(tibble)

# Rebuild sp_df to match the SampleID format
sp_df <- format1 %>%
  mutate(SampleID = paste(Time.period, Site, sep = "_")) %>%
  group_by(SampleID, Species) %>%
  summarise(Number = sum(Number, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = Number, values_fill = list(Number = 0)) %>%
  column_to_rownames("SampleID")

# Extract unique species and their traits
trait_df <- format1 %>%
  distinct(Species, Dietary, Wing.m, Size) %>%
  mutate(
    Dietary = as.factor(trimws(Dietary)),
    Wing.m = ordered(trimws(Wing.m), levels = c("B", "M/B", "M")),
    Size = as.numeric(as.character(Size))
  ) %>%
  column_to_rownames("Species")

trait_df <- trait_df[colnames(sp_df), , drop = FALSE]

# Calculate Gower's distance handling mixed variables
trait_dist <- gowdis(trait_df)

# Calculate Functional Diversity metrics using dbFD
fd_res <- dbFD(
  x = trait_dist, 
  a = sp_df, 
  calc.FRic = FALSE, 
  calc.FDiv = FALSE, 
  corr = "cailliez"
)

# Extract the Rao's Q values into new dataframe
rao_df <- data.frame(
  SampleID = rownames(sp_df),
  RaoQ = fd_res$RaoQ
)

# Merge Rao's Q directly into cwm_df
env_data_rao <- cwm_df %>%
  left_join(rao_df, by = "SampleID")
print(paste("Number of NAs in RaoQ:", sum(is.na(env_data_rao$RaoQ))))
env_data_rao <- env_data_rao %>% arrange(Time.period, as.numeric(Site))

env_data_rao$Year_factor   <- as.factor(env_data_rao$Year)
# Fit the GAMM with a Tweedie distribution and log link
mod_rao <- gam(
  RaoQ ~
    Woody.species +                 
    s(Pollution, k = 15) +          
    s(Precipitation) +              
    s(Wind) +                       
    s(Month, bs = "tp", k = 5) +    
    s(Site, bs = "re") +            
    s(Year_factor, bs = "re"),      
  data   = env_data_rao,
  family = tw(link = "log"),
  method = "REML",
  knots  = list(Month = c(3, 11))
)

# Check the results and diagnostic plots
# Note: positive temporal autocorrelation, 
# not an insufficient basis dimension.

summary(mod_rao)
par(mfrow = c(2, 2))
gam.check(mod_rao)
concurvity(mod_rao, full = TRUE)
gratia::draw(mod_rao)
