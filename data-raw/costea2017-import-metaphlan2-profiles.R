library(tidyverse)

# Read in the merged metaphlan2 profiles
dotenv::load_dot_env(here::here("data-raw", ".env"))
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
