library(readxl)
library(dplyr)
library(picante)
library(ape)
library(rotl)

### Step 1: Load and format data ###

dat <- read_excel("Dataframe_FINAL.xlsx", sheet = "Sheet1") %>%
  mutate(
    Date = as.Date(Date),
    Species = as.character(Species),
    Site = as.character(Site),
    Woody.species = as.character(Woody.species),
    Year = as.numeric(Year),
    Number = as.numeric(Number),
    Time.period = as.numeric(Time.period)
  )

### Step 2: Filter data for two periods ###

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