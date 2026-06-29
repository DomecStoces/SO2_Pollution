library(readxl)
library(dplyr)
library(tidyr)
library(mgcv)
library(ggplot2)
library(patchwork)

### Read dataset ###
env_data_pd <- read_excel("Dataframe_macrotemporal.xlsx", sheet = "Phylogenetic_diversity")

### Step 1: Summaries and statistical tests ###

# Wilcoxon tests
wilcox.test(SESpd_richness ~ Period_2002, data = env_data_pd)
wilcox.test(SESpd_swap ~ Period_2002, data = env_data_pd)

### Step 1.1: Merge with environemtnal data and model fit ###
# Prepare data for GAMMs

env_data_pd <- env_data_pd %>%
  mutate(
    Year_factor  = factor(Year),
    Site         = factor(Site),
    Woody.species = factor(Woody.species)
  )

### Step 2: Run the hypothesis GAMMs ###

# Model A: Testing lineage survival (richness null)

mod_pd_richness <- gam(
  SESpd_richness ~ Woody.species +
    s(Pollution, k = 15) +       
    s(Precipitation, k = 15) +   
    s(Wind, k = 15) +      
    s(Site, bs = "re") + 
    s(Year_factor, bs = "re"),
  data = env_data_pd, 
  method = "REML"
)

# Model B: Testing specific species dominance (independent swap null)

mod_pd_swap <- gam(
  SESpd_swap ~ Woody.species +
    s(Pollution, k = 15) +       
    s(Precipitation, k = 15) +   
    s(Wind, k = 15) +          
    s(Site, bs = "re") + 
    s(Year_factor, bs = "re"),
  data = env_data_pd,
  method = "REML"
)

# View results

summary(mod_pd_richness)
summary(mod_pd_swap)

####################################################################################################################

### Graphical representation ###

## Step 1: Build the prediction grid ###
# Smooth sequence of pollution, holding the tree species constant (birch)

pollution_seq <- seq(min(env_data_pd$Pollution, na.rm = TRUE),
                     max(env_data_pd$Pollution, na.rm = TRUE), length.out = 100)

new_data <- data.frame(
  Pollution = pollution_seq,
  Woody.species = "Birch",
  Policy.period = env_data_pd$Policy.period[1], # <--- ADDED THIS LINE
  Precipitation = mean(env_data_pd$Precipitation, na.rm = TRUE), 
  Wind = mean(env_data_pd$Wind, na.rm = TRUE),                    
  Year_factor = env_data_pd$Year_factor[1],
  Site = env_data_pd$Site[1]                  
)

### Step 2: Predict for both models ###
# Richness Predictions

pred_rich <- predict(
  mod_pd_richness, 
  newdata = new_data, 
  exclude = c("s(Year_factor)", "s(Site)", "Policy.period"), 
  se.fit = TRUE
)
# Swap Predictions

pred_swap <- predict(
  mod_pd_swap, 
  newdata = new_data, 
  exclude = c("s(Year_factor)", "s(Site)", "Policy.period"), 
  se.fit = TRUE
)

### Step 3: Format data for ggplot ###

plot_df <- bind_rows(
  new_data %>% mutate(
    fit = pred_rich$fit,
    se = pred_rich$se.fit,
    lower = fit - (1.96 * se),
    upper = fit + (1.96 * se),
    Model = "Richness null model" 
  ),
  new_data %>% mutate(
    fit = pred_swap$fit,
    se = pred_swap$se.fit,
    lower = fit - (1.96 * se),
    upper = fit + (1.96 * se),
    Model = "Independent-swap null model" 
  )
)

# Format the raw data points to match facets
raw_rich <- env_data_pd %>% 
  select(Pollution, SESpd = SESpd_richness) %>% 
  mutate(Model = "Richness null model") 

raw_swap <- env_data_pd %>% 
  select(Pollution, SESpd = SESpd_swap) %>% 
  mutate(Model = "Independent-swap null model") 

raw_points <- bind_rows(raw_rich, raw_swap)

# Set factor levels to force the facet order
plot_df <- plot_df %>%
  mutate(Model = factor(Model, levels = c("Richness null model", "Independent-swap null model")))

raw_points <- raw_points %>%
  mutate(Model = factor(Model, levels = c("Richness null model", "Independent-swap null model")))

# Paired test to compare the two null model outputs
test_nulls <- wilcox.test(
  env_data_pd$SESpd_richness, 
  env_data_pd$SESpd_swap, 
  paired = TRUE
)
print(test_nulls)

### Step 4: Plotting ###

### Step 4.1: Side-by-side (linear dependency of SESpd and SO2 concentrations) ###

p_phylo_combined <- ggplot() +
  
  # Raw observations
  geom_jitter(
    data = raw_points,
    aes(x = Pollution, y = SESpd),
    width = 0.5,
    height = 0,
    alpha = 0.25,
    size = 1.2,
    colour = "black"
  ) +
  
  # 95% CI
  geom_ribbon(
    data = plot_df,
    aes(x = Pollution, ymin = lower, ymax = upper),
    fill = "grey70",
    alpha = 0.4
  ) +
  
  # GAM smooth
  geom_line(
    data = plot_df,
    aes(x = Pollution, y = fit),
    linewidth = 1.1,
    colour = "black"
  ) +
  
  # Null expectation
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    colour = "black",
    linewidth = 0.7
  ) +  
  
  geom_hline(
    yintercept = c(-1.96, 1.96),
    linetype = "dashed",
    colour = "grey30",
    linewidth = 0.6
  ) +
  
  facet_wrap(~Model, nrow = 1) +
  
  labs(
    x = expression("Mean annual SO"[2]~"concentration ["*mu*"g · m"^-3*"]"),
    y = "Faith's phylogenetic diversity (SESpd)"
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 14
    ),
    
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    panel.border = element_rect()
  )
print(p_phylo_combined)

### Step 4.2: Side-by-side (boxplot of SESpd and Air quality policy period) ###

# Format the data into "long" format for faceting
box_df <- bind_rows(
  env_data_pd %>% 
    select(Period_2002, SESpd = SESpd_richness) %>% 
    mutate(Model = "Richness null model"),
  
  env_data_pd %>% 
    select(Period_2002, SESpd = SESpd_swap) %>% 
    mutate(Model = "Independent-swap null model")
) %>%
  mutate(Model = factor(Model, levels = c("Richness null model", "Independent-swap null model")))

# Run the Wilcoxon tests to report in your results
print("Richness null: before vs. after")
wilcox.test(SESpd ~ Period_2002, data = filter(box_df, Model == "Richness null model"))

print("Swap null: before vs. after")
wilcox.test(SESpd ~ Period_2002, data = filter(box_df, Model == "Independent-swap null model"))

# Create the dual-panel boxplot
p_box_combined <- ggplot(
  box_df,
  aes(
    x = factor(Period_2002, levels = c("Before_2002", "After_2002")),
    y = SESpd
  )
) +
  geom_boxplot(fill = "grey85", colour = "black", width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.2, colour = "grey40") +
  
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    colour = "black",
    linewidth = 0.7
  ) +  
  
  geom_hline(
    yintercept = c(-1.96, 1.96),
    linetype = "dashed",
    colour = "grey30",
    linewidth = 0.6
  ) +
  
  facet_wrap(~Model, nrow = 1) +
  
  scale_x_discrete(
    labels = c(
      "Before_2002" = "Before 2002",
      "After_2002"  = "After 2002"
    )
  ) +
  labs(
    x = "Air quality policy period",
    y = "Faith's phylogenetic diversity (SESpd)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 14),
    
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(),
    panel.border = element_rect()
  )
print(p_box_combined)

### Step 4.3: Combination panel of linear and box plots ###

final_figure <- p_box_combined / p_phylo_combined  + 
  plot_annotation(tag_levels = 'A')
print(final_figure)

pdf("Phylogenetic_diversity.pdf", width = 8, height = 10)
print(final_figure)
dev.off()
