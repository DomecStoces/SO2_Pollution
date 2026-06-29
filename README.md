# Community Reorganization of Carabid Beetles Following Long-Term Sulfur Dioxide Pollution Decline and Forest Succession

This repository contains the data and R scripts used to perform the statistical analyses for the manuscript: **"Community reorganization of carabid beetles following long-term sulfur dioxide pollution decline and forest succession."** The analytical workflow evaluates environmental drivers across taxonomic, functional, and phylogenetic metrics using Generalized Additive Mixed Models (GAMMs), multivariate ordinations, and null-model frameworks. 

## Repository Structure & Data

The raw observations of epigeic arthropod communities—sampled across nine sites located within a strict no-logging zone (collected using a 4% formaldehyde solution preservative)—are contained in the master dataset. This master file was processed to create two temporally distinct dataframes to capture different ecological processes:

* **`Dataframe_FINAL.xlsx`**: The master raw observation dataset.
* **`Dataframe_microtemporal.xlsx`**: Intra-annual resolution data. Trap catches are aggregated by site and sampling period to capture within-season dynamics (phenological noise preserved).
* **`Dataframe_macrotemporal.xlsx`**: Inter-annual resolution data. Trap catches are pooled into annual assemblages per site to isolate baseline macro-ecological state shifts.

---

## Methods & Analytical Pipeline

All statistical analyses were performed in **R 4.5.1**. The analysis is divided into four main parts corresponding to the methods described in the manuscript.

### Part 1: Taxon-Based Taxonomic Diversity (Intra-annual Scale)
**Script:** `GAMM.R`
**Data:** `Dataframe_microtemporal.xlsx`

This script models how carabid activity-density and species richness change under declining SO₂ levels, accounting for environmental variables and forest succession.
* **Activity-Density:** Modeled using a Tweedie distribution (log link) to handle zero-inflation and severe overdispersion. 
* **Species Richness:** Modeled using a quasipoisson distribution (log link) to account for underdispersion.
* **Structure:** GAMMs share a core predictor structure, including parametric fixed effects for dominant tree species (4 levels), independent non-linear penalized regression splines for SO₂, precipitation, and wind speed, and a cyclic cubic regression spline for the sampling month to account for seasonality.

### Part 2: Functional Trait Diversity (Intra-annual Scale)
**Script:** `CWM_FD.R`
**Data:** `Dataframe_microtemporal.xlsx`

This script examines shifts in the functional composition of beetle communities.
* **Community-Weighted Means (CWM):** Calculated for trophic strategy, dispersal ability (wing morphology), and body size. Collinearity among traits is evaluated via PCA.
* **Functional Diversity:** Quantified as Rao's quadratic entropy (Q) utilizing Gower trait dissimilarities with a Cailliez correction.
* **GAMMs:** Fits quasibinomial models for rescaled trophic/dispersal CWMs, a Gaussian model for body size, and a Tweedie model for Rao's Q. 

### Part 3: Community Composition Shifts (Inter-annual Scale)
**Scripts:** `PERMANOVA.R` and `NMDS.R`
**Data:** `Dataframe_macrotemporal.xlsx`

These scripts test for and visualize total compositional shifts between historically high-pollution (pre-2002) and recovery (post-2002) air policy periods.
* **PERMANOVA (`PERMANOVA.R`):** Tests marginal effects (Type III SS) on square-root transformed Bray-Curtis dissimilarities. Homogeneity of multivariate dispersions is confirmed prior to analysis.
* **NMDS & Indicator Species (`NMDS.R`):** Visualizes the two-dimensional community structure, fits an SO₂ concentration vector to the ordination, and utilizes the `IndVal.g` method to identify robust indicator taxa characteristic of each policy period.

### Part 4: Phylogenetic Diversity (Inter-annual Scale)
**Scripts:** `Phylogeny_from_OpenTree.R` and `Phylogenetic_diversity.R`
**Data:** `Dataframe_FINAL.xlsx` (generates `Phylogenetic_diversity.xlsx`)

This step isolates the effects of environmental filtering from regional species dominance by mapping the community to the Open Tree of Life.
* **Phylogeny Construction (`Phylogeny_from_OpenTree.R`):** Induces a subtree for the community, estimates branch lengths (Grafen’s method), and randomly resolves polytomies.
* **Diversity Modeling (`Phylogenetic_diversity.R`):** Computes the Standardized Effect Size of Faith’s PD (SESpd) against 999 null communities using two approaches: a richness-constrained model and an independent-swap algorithm. The outputs are compared via Wilcoxon rank-sum tests and modeled against continuous SO₂ using GAMMs.

---

## Dependencies

To replicate the analysis, ensure the following R packages are installed:

**Data Manipulation & General:**
`dplyr`, `tidyr`, `readxl`, `writexl`, `tibble`, `broom`, `scales`

**Modeling & Statistics:**
`mgcv`, `car`, `vegan`, `permute`, `FD`, `picante`, `indicspecies`

**Phylogeny:**
`ape`, `rotl`

**Multivariate & Trait Analysis:**
`FactoMineR`, `factoextra`, `labdsv`

**Visualization:**
`ggplot2`, `gratia`, `patchwork`, `ggtext`, `corrplot`, `ggforce`, `concaveman`, `ggrepel`, `ggnewscale`, `shadowtext`
