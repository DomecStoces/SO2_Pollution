#vegan::NMDS based on community matrix of vegan::PERMANOVA -> firstly do PERMANOVA then NMDS

library(vegan)
library(ggplot2)
library(dplyr)
library(labdsv)
library(ggforce)
library(concaveman)
library(ggrepel)
library(ggnewscale)
library(shadowtext)
library(indicspecies)
library(ggtext)

### Step 1: Run NMDS ###
set.seed(123)
nmds_result <- metaMDS(sqrt(sp_df), distance = "bray", k = 2, trymax = 50, autotransform = FALSE)

# Check stress
nmds_result$stress

### Step 2: Extract NMDS Scores and Prepare Data ###
# Define zoom parameters so they can be used for trimming
x_zoom <- c(-3.0, 1)
y_zoom <- c(-2.0, 2)

# Site scores function
site_scores <- as.data.frame(vegan::scores(nmds_result, display = "sites"))
site_scores$SampleID <- rownames(site_scores)
nmds_plot_data <- left_join(site_scores, env_data, by = "SampleID")
colnames(nmds_plot_data)[1:2] <- c("NMDS1", "NMDS2")

# Trim the NMDS plot data for outliers and zoom limits
nmds_plot_trimmed <- nmds_plot_data %>% filter(NMDS1 < 10)
nmds_trim <- nmds_plot_trimmed %>% filter(between(NMDS1, x_zoom[1], x_zoom[2]) & between(NMDS2, y_zoom[1], y_zoom[2]))

# Convex hulls and centroids
hulls <- nmds_trim %>% group_by(Policy_period) %>% slice(chull(NMDS1, NMDS2))
centroids <- nmds_trim %>% group_by(Policy_period) %>% summarize(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups="drop")
centroids <- centroids %>%
  mutate(
    Abbrev = dplyr::recode(Policy_period, "Before2002" = "B", "After2002" = "A")
  )

# Connect each sample to centroid
spider_lines <- left_join(nmds_trim, centroids, by = "Policy_period", suffix = c("", ".centroid"))

### Step 3: Fit Environmental Variables ###
ef <- envfit(nmds_result, env_data[, c("Pollution")], permutations = 999)
ef_scores <- as.data.frame(vegan::scores(ef, display = "vectors"))
ef_scores$Variable <- rownames(ef_scores)
ef_scores$pval <- ef$vectors$pvals
print(round(ef_scores$pval, 3))

ef_sig <- ef_scores %>% filter(pval < 0.05)

arrow_multiplier <- ordiArrowMul(ef)
ef_sig$NMDS1 <- ef_sig$NMDS1 * arrow_multiplier
ef_sig$NMDS2 <- ef_sig$NMDS2 * arrow_multiplier

### Step 4: Indicator Species Analysis (indicspecies::multipatt()) ###
set.seed(123)
inv <- multipatt(sp_df, env_data$Policy_period, 
                 func = "IndVal.g", 
                 control = how(nperm=9999))

# Extract the results table from multipatt
inv_res <- inv$sign
inv_res$Species <- rownames(inv_res)

# multipatt returns the square root of IndVal in the 'stat' column. Square it to get true IndVal.
inv_res$IndVal <- inv_res$stat^2 

# Identify the names of your groups
group_names <- levels(env_data$Policy_period)

# Filter for strong (IndVal >= 0.3) and significant (p <= 0.01) indicators
sig_indicators <- inv_res %>%
  filter(p.value < 0.05 & IndVal > 0.3)
group_cols <- setdiff(colnames(inv_res), c("index", "stat", "p.value", "Species", "IndVal"))
# Determine which group each species indicates (extracting the column name where value is 1)
sig_indicators$Group <- apply(sig_indicators[, group_cols, drop = FALSE], 1, function(row) {
  raw_name <- names(row)[row == 1][1]
  gsub("^s\\.", "", raw_name)
})
print(sig_indicators[, c("Species", "Group", "IndVal", "p.value")])
### Step 5: NMDS Species Scores ###
species_scores <- as.data.frame(vegan::scores(nmds_result, display = "species"))
species_scores$Species <- rownames(species_scores)

# Join the NMDS species coordinates with the significant indicator dataframe
species_scores_sig <- species_scores %>% 
  inner_join(sig_indicators %>% select(Species, Group), by = "Species")

# Prepare abbreviated labels
species_scores_sig$Label <- abbreviate(species_scores_sig$Species, minlength = 4)
species_scores_sig$Group <- factor(species_scores_sig$Group, levels = c("After2002", "Before2002"))

### Step 6: Plotting ###
final <- ggplot() +
  # Sites (Using full data)
  geom_point(data = nmds_plot_data, aes(x = NMDS1, y = NMDS2), shape = 1, size = 1.5, color = "black") +
  geom_segment(data = spider_lines,
               aes(x = NMDS1, y = NMDS2, xend = NMDS1.centroid, yend = NMDS2.centroid, color = Policy_period),
               linewidth = 0.4) +
  
  # Convex hulls
  geom_polygon(data = hulls,
               aes(x = NMDS1, y = NMDS2, color = Policy_period, group = Policy_period), fill=NA,
               alpha = 0.2, linewidth = 0.8) +   
  scale_color_brewer(palette = "Pastel1") + ggnewscale::new_scale_fill() + 
  
  # Species points
  geom_point(data = species_scores_sig,
             aes(x = NMDS1, y = NMDS2, fill = Group),
             shape = 21, size = 3, color = "black", stroke = 0.5) +
  scale_fill_manual(values = c("Before2002" = "white", "After2002" = "black")) +
  
  # Species labels
  ggrepel::geom_text_repel(data = species_scores_sig,
                           aes(x = NMDS1, y = NMDS2, label = Label),
                           color = "black", fontface = "bold", size = 4,
                           box.padding = 0.5, point.padding = 0.3, 
                           segment.size = 0.1, max.overlaps = Inf) +
  
  # Centroid text
  geom_label(data = centroids,
             aes(x = NMDS1, y = NMDS2, label = Abbrev),
             fontface = "bold", size = 5, fill="white", linewidth = 0.1) +
  # Arrow of SO2 pollution
  geom_segment(data = ef_sig,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 1.2) +
  
  # The label (placed slightly past the tip of the arrow)
  geom_richtext(data = ef_sig,
                aes(x = NMDS1 * 1.15, y = NMDS2 * 1.15, 
                    label = "SO<sub>2</sub> concentration [\u03bcg\u00b7m<sup>-3</sup>]"),
                color = "black", size = 5, 
                fill = NA, label.color = NA) +
  
  # Expanding xlim entirely encompasses the left outlier.
  coord_fixed(ratio = 1, xlim = c(-3.0, 1.0), ylim = c(-1.5, 1.2), clip = "on") +
  
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.line = element_blank(), 
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 10)) 
  )
print(final)

tiff("NMDS.tiff", 
     width = 15, height = 10,     
     units = "in",                  
     res = 600,                     
     compression = "lzw")           
print(final)
dev.off()

