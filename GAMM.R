# -------------------------------------------------------------------------
# Long-term Activity-Density GAM Analysis
# Research Question: How does carabid activity and species richness change 
# under declining SO₂ levels while accounting for 
# environmental variables and succession?
# -------------------------------------------------------------------------

library(mgcv)
library(gratia)
library(ggplot2)
library(dplyr)
library(patchwork)
library(broom)
library(readxl)
library(tidyr)
library(ggtext)
library(corrplot)
library(car)

### Read dataset ###
abundance_richness_df <- read_excel("Dataframe_microtemporal.xlsx")

### Step 1: Set correct factor levels for modelling ###
abundance_richness_df$Woody.species <- as.factor(abundance_richness_df$Woody.species)
abundance_richness_df$Year_factor <- as.factor(abundance_richness_df$Year)
abundance_richness_df$Month <- as.numeric(abundance_richness_df$Month)
abundance_richness_df$Date <- as.Date(abundance_richness_df$Date, format = "%Y-%m-%d")

### Step 2: Multicollinearity checks ###
# 2.1) Spearman's rank correlation for non-linear predictors
cor_matrix <- cor(
  abundance_richness_df[, c("Pollution", "T", "Precipitation", "Wind")],
  use = "complete.obs",
  method = "spearman"
)

colnames(cor_matrix) <- rownames(cor_matrix) <- c(
  "SO2 concentration", "Temperature", "Precipitation", "Wind speed"
)

# Export correlation plot
tiff("Spearman_Correlation.tiff", width = 9, height = 6, units = "in", res = 300, compression = "lzw")           
corrplot(cor_matrix, method = "number", type = "upper", number.cex = 1.5, number.font = 2,
         number.digits = 2, addCoef.col = "black", tl.cex = 1.3, tl.col = "black", tl.srt = 45,
         cl.cex = 1.3, cl.ratio = 0.2, cl.align.text = "l", cl.offset = 0.5,
         col = colorRampPalette(c("red", "white", "blue"))(200), mar = c(0, 0, 1, 2))
dev.off()

# 2.2) Variance Inflation Factor (VIF)
vif_model <- lm(Pollution ~ T + Precipitation + Wind, data = abundance_richness_df)
print(vif(vif_model)) 

####################################################################################################################

### Activity-density ###
# Tweedie GAM with a log link function.
# Includes a parametric term for dominant tree species (Woody.species), independent environmental smooths 
# (SO2, Precipitation, Wind), a cyclic smooth for seasonality (Month), 
# and random intercepts for Site and Year.
# Note: Temperature was collinear with Month, hence dropped for edf.

### Step 1: Model fit of number of individuals mgcv::gam() ###
mod_abundance <- gam(
  Abundance ~ Woody.species + 
    s(Pollution, k = 15) +       
    s(Precipitation, k = 15) +   
    s(Wind, k = 15) +            
    s(Month, bs = "cc", k = 8) + 
    s(Site, bs = "re") +         
    s(Year_factor, bs = "re"),
  data = abundance_richness_df,
  family = tw(link = "log"),     
  method = "REML",
  knots = list(Month = c(3, 11)) 
)
summary(mod_abundance)

### Step 2: Model diagnostics ###
par(mfrow = c(2, 2))
gam.check(mod_abundance)
concurvity(mod_abundance, full = TRUE)
gratia::draw(mod_abundance)

### Step 3: Vizualization for activity-density trend ###
# Fit new GAM for trend
gam_abundance <- gam(
  Abundance ~ s(as.numeric(Date), bs = "cs"), 
  data   = abundance_richness_df,
  family = tw(link = "log"),  
  method = "REML"
)

# Extract predictions and CI
pred_abund <- predict(gam_abundance, newdata = abundance_richness_df, se.fit = TRUE, type = "response")

# Extract fitted values from the main abundance model (mod_abundance)
used_rows_abund <- as.numeric(rownames(model.frame(mod_abundance)))
abundance_richness_df$Fitted_Abund <- NA_real_
abundance_richness_df$Fitted_Abund[used_rows_abund] <- fitted(mod_abundance, type = "response")

# Re-update filtered_data object so it includes the new column
filtered_data <- abundance_richness_df

abundance_richness_df$Abund_pred  <- pred_abund$fit
abundance_richness_df$Abund_se    <- pred_abund$se.fit
abundance_richness_df$Abund_lower <- abundance_richness_df$Abund_pred - (1.96 * abundance_richness_df$Abund_se)
abundance_richness_df$Abund_upper <- abundance_richness_df$Abund_pred + (1.96 * abundance_richness_df$Abund_se)

# Calculate a dynamic scale factor for activity-density
# We find the maximum upper bound of abundance to scale it accurately against the 160 SO2 limit
max_abund <- ceiling(max(abundance_richness_df$Abund_upper, na.rm = TRUE))
scale_factor_abund <- 160 / max_abund

# Scale activity-density metrics to the 0-160 SO2 axis
abundance_richness_df <- abundance_richness_df %>%
  mutate(
    Abund_scaled       = Abund_pred * scale_factor_abund,
    Abund_lower_scaled = Abund_lower * scale_factor_abund,
    Abund_upper_scaled = Abund_upper * scale_factor_abund
  )

# Build the abundance plot
plot_abundance <- ggplot() +
  # 1. SO2 values
  geom_point(
    data  = filtered_data,
    aes(x = Date, y = Pollution, size = Fitted_Abund), 
    shape = 21, fill = "gray70", colour = "black", alpha = 0.7
  ) +
  
  # 2. SO2 Pollution Trend
  # FIX: Moved color inside aes() and mapped it to "pollution"
  geom_smooth(
    data = abundance_richness_df, 
    aes(x = Date, y = Pollution, linetype = "pollution", fill = "pollution", color = "pollution"),
    method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.3, se = TRUE
  ) +
  
  # 3. Predicted ABUNDANCE CI Ribbon
  geom_ribbon(
    data = abundance_richness_df, 
    aes(x = Date, ymin = Abund_lower_scaled, ymax = Abund_upper_scaled, fill = "predicted"),
    alpha = 0.3
  ) +
  
  # 4. Predicted ABUNDANCE Trend
  # FIX: Moved color inside aes() and mapped it to "predicted"
  geom_line(
    data = abundance_richness_df, 
    aes(x = Date, y = Abund_scaled, linetype = "predicted", color = "predicted"),
    linewidth = 1.2
  ) +
  
  # 5. Policy Line
  geom_vline(xintercept = policy_date, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  scale_size_continuous(
    name = expression("Predicted activity-density for SO"[2] * " values"),
    range = c(0.5, 8),
    limits = c(0, 30),
    breaks = c(5, 10, 15, 20)
  ) + 
  
  scale_fill_manual(
    name = NULL,
    values = c(predicted = "grey60", pollution = "grey80"),
    labels = c(predicted = "Predicted activity-density", pollution = expression("SO"[2] * " concentration trend"))
  ) +
  
  scale_color_manual(
    name = NULL,
    values = c(predicted = "black", pollution = "black"),
    labels = c(predicted = "Predicted activity-density", pollution = expression("SO"[2] * " concentration trend"))
  ) +
  
  scale_linetype_manual(
    name = NULL,
    values = c(predicted = "dashed", pollution = "solid"),
    labels = c(predicted = "Predicted activity-density", pollution = expression("SO"[2] * " concentration trend"))
  ) +
  
  # FIX: Simplified guides. Because everything is mapped in aes() now, 
  # ggplot merges them perfectly without needing messy overrides.
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 1)
  ) +
  
  # Y-Axis Formatting
  scale_y_continuous(
    limits = c(0, 160), 
    breaks = seq(0, 160, by = 20),
    expand = expansion(mult = c(0, 0.05)), 
    sec.axis = sec_axis(
      transform = ~ . / scale_factor_abund, 
      name      = "Predicted activity-density",
      breaks    = seq(0, 30, by = 5)
    )
  ) +
  
  # X-Axis and Theme
  scale_x_date(
    date_breaks = "3 years", date_labels = "%Y",
    expand = expansion(add = c(250, 0)),
    limits = range(abundance_richness_df$Date, na.rm = TRUE)
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
    axis.title.y.right = element_text(angle = 90, margin = margin(l = 15)),
    axis.title.y.left = ggtext::element_markdown()
  ) +
  labs(x = "Year", y = "SO<sub>2</sub> concentration [&mu;g&middot;m<sup>-3</sup>]")
plot_abundance
####################################################################################################################

### Species richness ###
# Quasipoisson distribution and log link to account for underdispersion in the count data
# Includes a parametric term for dominant tree species (Woody.species), independent environmental smooths 
# (SO2, Precipitation, Wind), a cyclic cubic regression spline for Month (constrained between March and November), 
# and random intercepts for Site and Year to account for spatial and inter-annual dependencies.

### Step 1: Species richness mgcv::gam() ###
# Data are underdispersed (0.637)
mod_richness <- gam(
  Richness ~ Woody.species + 
    s(Pollution, k = 15) +       
    s(Precipitation, k = 15) +   
    s(Wind, k = 15) +            
    s(Month, bs = "cc", k = 8) + 
    s(Site, bs = "re") +         
    s(Year_factor, bs = "re"),
  data = abundance_richness_df,
  family = quasipoisson(link = "log"), 
  method = "REML",
  knots = list(Month = c(3, 11))
)

### Step 2: Model diagnostics ###
summary(mod_richness)
par(mfrow = c(2, 2))
gam.check(mod_richness)
concurvity(mod_richness, full = TRUE)
gratia::draw(mod_richness)

### Step 3: Plot SO2 vs. predicted richness over time ###
# Fit a secondary GAM to isolate the overall temporal trend of species richness
gam_richness <- gam(
  Richness ~ s(as.numeric(Date), bs = "cs"), 
  data   = abundance_richness_df,
  family = quasipoisson(link = "log"),
  method = "REML"
)

# Extract predictions and confidence intervals
pred_gam <- predict(gam_richness, newdata = abundance_richness_df, se.fit = TRUE, type = "response")

abundance_richness_df$Richness_pred  <- pred_gam$fit
abundance_richness_df$Richness_se    <- pred_gam$se.fit
abundance_richness_df$Richness_lower <- abundance_richness_df$Richness_pred - (1.96 * abundance_richness_df$Richness_se)
abundance_richness_df$Richness_upper <- abundance_richness_df$Richness_pred + (1.96 * abundance_richness_df$Richness_se)

# Rescale richness predictions to match the primary Pollution (SO2) y-axis
range_Pollution <- range(abundance_richness_df$Pollution, na.rm = TRUE)
range_richness  <- range(abundance_richness_df$Richness_pred, na.rm = TRUE)

scale_factor <- 160/6

abundance_richness_df <- abundance_richness_df %>%
  mutate(
    Richness_scaled       = Richness_pred * scale_factor,
    Richness_lower_scaled = Richness_lower * scale_factor,
    Richness_upper_scaled = Richness_upper * scale_factor
  )

# Extract fitted values from the main model (mod_richness) to size the bubbles
used_rows <- as.numeric(rownames(model.frame(mod_richness)))
abundance_richness_df$Fitted <- NA_real_
abundance_richness_df$Fitted[used_rows] <- fitted(mod_richness, type = "response")

# Filter points for clean visualization (adjusted for lower richness counts)
# You may want to lower this to Fitted >= 1 if too many points disappear
filtered_data <- abundance_richness_df

policy_date <- as.Date("2002-06-01")

# Determine y-axis boundaries safely
ci_range <- range(abundance_richness_df$Richness_lower_scaled, abundance_richness_df$Richness_upper_scaled, na.rm = TRUE)
padding  <- 0.05 * diff(ci_range)
y_min    <- floor(min(0, ci_range[1] - padding))
y_max    <- ceiling(max(150, ci_range[2] + padding))

### Step 4: Build the dual-axis plot ###
plot_richness <- ggplot() +
  # 1. SO2 values sized by main model predictions
  geom_point(
    data  = filtered_data,
    aes(x = Date, y = Pollution, size = Fitted),
    shape = 21, fill = "gray70", colour = "black", alpha = 0.7
  ) +
  
  # 2. SO2 Pollution Trend (Solid Line)
  geom_smooth(
    data = abundance_richness_df, aes(x = Date, y = Pollution, linetype = "pollution", fill = "pollution"),
    color = "black", method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.3, se = TRUE
  ) +
  
  # 3. Predicted Richness CI Ribbon
  geom_ribbon(
    data = abundance_richness_df, aes(x = Date, ymin = Richness_lower_scaled, ymax = Richness_upper_scaled, fill = "predicted"),
    alpha = 0.3
  ) +
  
  # 4. Predicted Richness Trend (Dotted Line)
  geom_line(
    data = abundance_richness_df, aes(x = Date, y = Richness_scaled, linetype = "predicted"),
    color = "black", linewidth = 1.2
  ) +
  
  # 5. Policy Update Vertical Line
  geom_vline(xintercept = policy_date, linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate("text", x = policy_date, y = 150, label = "Air policy update (01/06/2002)", 
           hjust = 0.8, angle = 90, vjust = -1.3, size = 4.5, color = "black") +
  
  # Y-Axis Formatting (Primary and Secondary)
  scale_y_continuous(
    limits = c(0, 160), 
    breaks = seq(0, 160, by = 20),
    expand = expansion(mult = c(0, 0.05)), 
    sec.axis = sec_axis(
      transform = ~ . / scale_factor, 
      name      = "Predicted species richness",
      breaks    = seq(0, 6, by = 2)
    )
  ) +
  
  # X-Axis Formatting
  scale_x_date(
    date_breaks = "3 years", date_labels = "%Y",
    expand = expansion(add = c(250, 0)),
    limits = range(abundance_richness_df$Date, na.rm = TRUE)
  ) +
  
  # Legend Mapping
  scale_fill_manual(
    name = NULL,
    values = c(predicted = "grey60", pollution = "grey80"),
    breaks = c("predicted", "pollution"),
    labels = c("Predicted species richness", expression("SO"[2] * " pollution")),
    guide = "none"
  ) + 
  scale_size_continuous(
    name = expression("Predicted richness for SO"[2] * " values"),
    range = c(0.5, 5), 
    limits = c(0.5, 6),      
    breaks = c(1, 3, 5) 
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(predicted = "dotted", pollution = "solid"),
    breaks = c("predicted", "pollution"),
    labels = c("Predicted species richness", expression("SO"[2] * " concentration trend"))
  ) +
  guides(
    linetype = guide_legend(order = 1, override.aes = list(color = "black", linewidth = 1.1))
  ) +
  
  # Final Theme
  theme_classic(base_size = 14) +
  theme(
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(5, "pt"),
    axis.title.y.right = element_text(angle = 90, margin = margin(l = 15)),
    axis.title.y.left = ggtext::element_markdown(),
    legend.position = "right"
  ) +
  labs(
    x = "Year",
    y = "SO<sub>2</sub> concentration [&mu;g&middot;m<sup>-3</sup>]"
  )
plot_richness

####################################################################################################################

### Combine both plots using patchwork ###
combined_plot <- (plot_abundance / plot_richness) + 
  plot_layout(guides = "collect") + 
  theme(legend.position = "right") + plot_annotation(tag_levels = 'A')
combined_plot

pdf("SO2_abun_rich.pdf", width = 12, height = 13)
print(combined_plot)
dev.off()