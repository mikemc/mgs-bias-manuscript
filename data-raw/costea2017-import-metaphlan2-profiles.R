# library(tidyverse)

## First we combine the individual sample profiles into a single profile

# Fill in the paths to the metaphlan2 profiles and for the utility script
# merge_metaphlan_tables.py included with metaphlan2 (in the "utils" folder)
out_path = "/folder/for/metaphlan2/output"
merge_metaphlan_path = "/path/to/merge_metaphlan_tables.py"

command <- paste(
    merge_metaphlan_path, 
    paste0(out_path, "/*_profiled_metagenome.tsv"),
    ">",
    paste0(out_path, "/costea2017_metaphlan2_profiles.tsv")
    )
system(command)

## Then, read into R and save for loading with the `data` function

tb <- readr::read_tsv(paste0(out_path, "/costea2017_metaphlan2_profiles.tsv"))
tb %>% glimpse
tb <- tb %>% 
    rename(Clade = ID) %>%
    rename_at(vars(-Clade), ~str_extract(., "ERR[0-9]+"))
tb %>% glimpse

costea2017_metaphlan2_profiles <- tb
usethis::use_data(costea2017_metaphlan2_profiles)
