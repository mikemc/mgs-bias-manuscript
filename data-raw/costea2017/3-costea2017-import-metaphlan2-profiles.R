library(tidyverse)
dotenv::load_dot_env(here::here("data-raw", ".env"))

# Import metaphlan2 profiles --------------------------------

# Read in the merged metaphlan2 profiles
profiles_fn <- file.path(Sys.getenv("DATA_PATH"), "costea2017", "metaphlan2",
    "costea2017_metaphlan2_custom_profiles.tsv")
tb <- read_tsv(profiles_fn)
tb %>% glimpse
# Improve column names
tb <- tb %>% 
    rename(Clade = ID) %>%
    rename_at(vars(-Clade), ~str_extract(., "ERR[0-9]+"))
tb %>% names
# Save for loading with the `data` function
costea2017_metaphlan2_profiles <- tb
usethis::use_data(costea2017_metaphlan2_profiles)

# Import Metaphlan2 phylogeny -------------------------------

# Download tree from https://github.com/waldronlab/curatedMetagenomicData
dl_path <- file.path(Sys.getenv("DATA_PATH"), "costea2017")
tree_url <- "https://github.com/waldronlab/curatedMetagenomicData/raw/master/inst/extdata/metaphlan2_selected.tree.reroot.nwk.bz2"
tree_fn <- file.path(dl_path, basename(tree_url))
download.file(tree_url, tree_fn)
tree <- ape::read.tree(tree_fn)
tree
# Tip labels are clade strings at the species level, ending in a GCF
# identifier, except for a few with no taxonomy info.
head(tree$tip.label)
tree$tip.label %>%
    {.[!str_detect(., "s__")]}
# To match the clade names in the above, we want to trim the GCF_.+ at the end
tree$tip.label <- tree$tip.label %>%
    str_extract(".+(?=\\|+GCF)")
tree <- ape::drop.tip(tree, NA_character_)
# Most bacterial Metaphlan species id's are in the tree, with two exceptions
tb$Clade %>%
    str_subset("s__") %>%
    str_subset("t__", negate = TRUE) %>%
    str_subset("k__Viruses", negate = TRUE) %>%
    str_subset("k__Eukaryota", negate = TRUE) %>%
    str_subset("unclassified", negate = TRUE) %>%
    setdiff(tree$tip.label)
#> [1] "k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae|g__Lactobacillus|s__Lactobacillus_casei_paracasei"        
#> [2] "k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae|g__Streptococcus|s__Streptococcus_mitis_oralis_pneumoniae"
# Let's manually adjust these to match
tree$tip.label <- tree$tip.label %>%
    str_replace("Streptococcus_mitis", 
        "Streptococcus_mitis_oralis_pneumoniae") %>%
    str_replace("Lactobacillus_casei", "Lactobacillus_casei_paracasei")
# Save for loading with the `data` function
costea2017_metaphlan2_tree <- tree
usethis::use_data(costea2017_metaphlan2_tree)
