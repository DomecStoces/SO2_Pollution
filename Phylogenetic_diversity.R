library(readxl)
library(writexl)
library(dplyr)
library(picante)
library(ape)
library(rotl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mgcv)
library(patchwork)

### Step 1: Load and format data ###

dat <- read_excel("Dataframe_final.xlsx", sheet = "Sheet1") %>%
  mutate(
    Date = as.Date(Date),
    Species = as.character(Species),
    Site = as.character(Site),
    Woody.species = as.character(Woody.species),
    Year = as.numeric(Year),
    Number = as.numeric(Number),
    Time.period = as.numeric(Time.period)
  )

## Step 2: Filter data for two periods ###

dat <- dat %>%
  mutate(
    Period_2002 = case_when(
      Date >= as.Date("1989-04-15") & Date < as.Date("2001-10-27") ~ "Before_2002",
      Date >= as.Date("2002-04-04") & Date < as.Date("2015-10-24") ~ "After_2002",
      TRUE ~ NA_character_
    ),
    SampleID_PD = paste(Year, Site, sep = "_")
  ) %>%
  filter(!is.na(Period_2002))

### Step 3: Community matrix ###

comm <- xtabs(Number ~ SampleID_PD + Species, data = dat)
comm <- as.matrix(comm)
comm <- comm[, colSums(comm) > 0, drop = FALSE]
sample_meta <- dat %>%
  distinct(SampleID_PD, Period_2002)
bad_ids <- sample_meta %>%
  count(SampleID_PD) %>% filter(n > 1)
if (nrow(bad_ids) > 0) {
  warning("Some SampleID have multiple Period_2002 assignments. Check dates within SampleID.")
}
bad_ids <- dat %>%
  distinct(SampleID_PD, Period_2002) %>%
  count(SampleID_PD) %>%
  filter(n > 1)

bad_ids

sample_meta <- sample_meta %>% distinct(SampleID_PD, .keep_all = TRUE)

### Step 4: Build a phylogeny for species from OpenTree of rotl ###

clean_to_binomial <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  # remove OpenTree suffix, e.g. "_ott12345"
  x <- sub("_ott[0-9]+$", "", x)
  
  # replace underscores with spaces (OpenTree format)
  x <- gsub("_", " ", x)
  
  # collapse whitespace
  x <- gsub("\\s+", " ", x)
  
  # keep only first two tokens (Genus species)
  x <- sub("^([A-Z][a-z]+)\\s+([a-z-]+).*", "\\1 \\2", x)
  
  x
}

# Extract your species list (assuming they are the column names from your community matrix)
raw_species <- colnames(comm)

# Apply your cleaning function to get valid binomials
cleaned_species <- clean_to_binomial(raw_species)

# Query the Open Tree of Life to create the 'matches' dataframe

matches <- tnrs_match_names(names = cleaned_species)

### Step 4.1: Filter ###

matches_ok <- matches %>% dplyr::filter(!is.na(ott_id))
ott_ids <- matches_ok$ott_id

### Step 4.2: Induce subtree from OpenTree ###

ott_info <- rotl::taxonomy_taxon_info(ott_ids = ott_ids)
ott_df <- dplyr::bind_rows(lapply(ott_info, function(rec) {
  tibble::tibble(
    ott_id = rec$ott_id,
    name = rec$name,
    unique_name = rec$unique_name,
    rank = rec$rank,
    source = rec$source,
    suppressed = isTRUE(rec$is_suppressed) | isTRUE(rec$is_suppressed_from_synth)
  )
}))

message("OTT matches: ", nrow(ott_df),
        " | suppressed: ", sum(ott_df$suppressed),
        " | unique taxa: ", dplyr::n_distinct(ott_df$unique_name))

# Flag ambiguous/non-species matches
suspect <- ott_df %>% dplyr::filter(rank != "species" | grepl("\\(", unique_name))
if (nrow(suspect) > 0) {
  message("Potentially ambiguous matches (inspect):")
  print(suspect, n = Inf)
}

# Step 4.3: Induce subtree from OpenTree

tree <- rotl::tol_induced_subtree(ott_ids = ott_ids)
tree$tip.label <- clean_to_binomial(tree$tip.label)
colnames(comm) <- clean_to_binomial(colnames(comm))
if (is.null(tree$edge.length)) {
  tree <- ape::compute.brlen(tree, method = "Grafen")
}

# Step 4.4: Remove any duplicated tip labels

if (any(duplicated(tree$tip.label))) {
  dup <- unique(tree$tip.label[duplicated(tree$tip.label)])
  warning("Duplicate tip labels after cleaning; dropping duplicates: ",
          paste(dup, collapse = ", "))
  tree <- ape::drop.tip(tree, tree$tip.label[duplicated(tree$tip.label)])
}

### Step 4.5: Resolve polytomies ###

tree <- ape::multi2di(tree)

### Step 4.6: Align species between community matrix and phylogeny ###
keep_spp <- intersect(colnames(comm), tree$tip.label)
if (length(keep_spp) < 5) warning("Very few species overlap between community matrix and tree. Check naming.")

comm2 <- comm[, keep_spp, drop = FALSE]
tree2 <- ape::drop.tip(tree, setdiff(tree$tip.label, keep_spp))

### Step 5: Compute PD and SESpd (Faith’s PD; richness-controlled null) ###

sample_meta <- dat %>%
  distinct(SampleID_PD, Period_2002) %>%
  filter(!is.na(Period_2002))

sample_meta <- sample_meta %>%
  filter(SampleID_PD %in% rownames(comm2))

sample_meta <- sample_meta %>%
  mutate(SampleID_PD = factor(SampleID_PD, levels = rownames(comm2))) %>%
  arrange(SampleID_PD)

stopifnot(identical(as.character(sample_meta$SampleID_PD), rownames(comm2)))
comm_pa <- (comm2 > 0) * 1
keep_rows <- rowSums(comm_pa) > 0
comm_pa <- comm_pa[keep_rows, , drop = FALSE]
sample_meta <- sample_meta[keep_rows, , drop = FALSE]
if (is.null(tree2$edge.length)) stop("tree2 has no branch lengths. Add compute.brlen() in Step 4.")
pd_obs <- picante::pd(comm_pa, tree2, include.root = FALSE)

# Null Model 1: Richness
set.seed(1)
ses_rich <- picante::ses.pd(
  samp = comm_pa,
  tree = tree2,
  null.model = "richness",
  runs = 999,
  iterations = 1000
)

# Null Model 2: Independent Swap
set.seed(1)
ses_swap <- picante::ses.pd(
  samp = comm_pa,
  tree = tree2,
  null.model = "independentswap",
  runs = 999,
  iterations = 1000
)

# Create pd_res
pd_res <- tibble::tibble(
  SampleID_PD = rownames(comm_pa),
  Period_2002 = sample_meta$Period_2002,
  Richness = rowSums(comm_pa),
  SESpd_richness = ses_rich$pd.obs.z,      
  SESpd_swap = ses_swap$pd.obs.z       
)

### Step 6: Summaries and statistical tests ###

# Summary stats (using SESpd_richness as the primary metric for the summary)
pd_res %>%
  group_by(Period_2002) %>%
  summarise(
    n_samples = n(),
    mean_richness = mean(Richness),
    mean_SESpd_richness = mean(SESpd_richness),
    .groups = "drop"
  )

# Wilcoxon tests
wilcox.test(SESpd_richness ~ Period_2002, data = pd_res)
wilcox.test(SESpd_swap ~ Period_2002, data = pd_res)

### Step 6.1: Merge with environemtnal data and model fit ###
# Prepare data for GAM
# Collapse the original data to one row per SampleID_PD to get environmental variables

env_data_pd <- dat %>%
  group_by(SampleID_PD) %>%
  summarise(
    Pollution = mean(Pollution, na.rm = TRUE),
    Precipitation = mean(Precipitation, na.rm = TRUE), 
    Wind = mean(Wind, na.rm = TRUE),                   
    Year = first(Year),
    Policy.period = first(Policy_period),
    Woody.species = first(Woody.species),
    Site = first(Site),             
    .groups = "drop"
  ) %>%
  inner_join(pd_res, by = "SampleID_PD") %>%
  mutate(
    Year_factor = as.factor(Year),
    Site = as.factor(Site), 
    Woody.species = as.factor(Woody.species),
    Policy.period = factor(Policy.period, levels = c("Before2002", "After2002"))
  )

### Step 7: Run the hypothesis GAMMs ###

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
