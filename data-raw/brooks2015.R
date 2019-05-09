library(tidyverse)
library(here)

# Download needed files from the supplement of Brooks2015, available at
# https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-015-0351-6
#
# "Additional file 2 Experimental design. Table of the prescribed mixing
# proportions, plate, and barcode for the experiments mixing equal proportions of
# cells, DNA, and PCR product."
# https://static-content.springer.com/esm/art%3A10.1186%2Fs12866-015-0351-6/MediaObjects/12866_2015_351_MOESM2_ESM.csv
#
# "Additional file 10 Table of above-threshold counts. Above-threshold counts for
# each sample in the experiments mixing equal amounts of cells, DNA, and PCR
# product."
# https://static-content.springer.com/esm/art%3A10.1186%2Fs12866-015-0351-6/MediaObjects/12866_2015_351_MOESM10_ESM.csv
#
# "Additional file 11 Table of below-threshold counts. Below-threshold counts for
# each sample in the experiments mixing equal amounts of cells, DNA, and PCR
# product."
# https://static-content.springer.com/esm/art%3A10.1186%2Fs12866-015-0351-6/MediaObjects/12866_2015_351_MOESM11_ESM.csv

temp_path <- tempdir(check = TRUE)
file_numbers <- c(2, 10, 11)
urls <- paste0("https://static-content.springer.com/esm/",
    "art%3A10.1186%2Fs12866-015-0351-6/MediaObjects/",
    "12866_2015_351_MOESM", file_numbers, "_ESM.csv")
file_names <- file.path(temp_path,
    paste0("AdditionalFile", file_numbers, ".csv"))
walk2(urls, file_names, download.file)

##################

# Read the downloaded csv files. The counts files yield an empty column of all
# NAs due to trailing commas on each line, which we remove.
names(file_names) <- c("design", "counts_above", "counts_below")
design <- readr::read_csv(file_names["design"])
above <- readr::read_csv(file_names["counts_above"]) %>%
    select_if(~!all(is.na(.)))
below <- readr::read_csv(file_names["counts_below"]) %>%
    select_if(~!all(is.na(.)))

design
above[1:5, 1:5]
below[1:5, 1:5]

# The names corresponding to the mock taxa differ between the count and design
# tables, and also from the Genus_species format we ultimately want to use. 
# A table connecting the different names:
mock_taxa <- tibble(
    Taxon = c("Gardnerella_vaginalis", "Atopobium_vaginae",
        "Lactobacillus_crispatus", "Lactobacillus_iners",
        "Prevotella_bivia", "Sneathia_amnii", "Streptococcus_agalactiae"),
    Design_name = c("Gvaginalis", "Avaginae", "Lcrispatus", "Liners",
        "Pbivia", "Samnii", "GroupBStrep"),
    Count_name = c("Gardnerella vaginalis", "Atopobium vaginae",
        "Lactobacillus crispatus_cluster", "Lactobacillus iners",
        "Prevotella bivia", "Sneathia amnii", "Streptococcus agalactiae")
)

# Join the above and below threshold counts into a single dataframe (tibble) in
# "tidy" or tall format
above0 <- above %>%
    gather("Taxon", "Count", -Sample) 
below0 <- below %>%
    gather("Taxon", "Count", -Sample) %>%
    mutate(Taxon = str_extract(Taxon, ".+(?=BT)"))
tb <- list(above = above0, below = below0) %>%
    bind_rows(.id = "Table")
tb <- tb %>%
    mutate(Mock = (Taxon %in% mock_taxa$Count_name))
tb
# Lump all the non-mock counts into an "Other" category.
tbm <- tb %>%
    rename(Count_name = Taxon) %>%
    left_join(mock_taxa, by = "Count_name") %>%
    mutate(Taxon = ifelse(is.na(Taxon), "Other", Taxon)) %>%
    select(Sample, Taxon, Count, Table) %>%
    group_by(Sample, Taxon, Table) %>%
    summarize(Count = sum(Count)) %>%
    ungroup
tbm

# Check that the counts match what was reported in the Brooks2015 SI Rmd
tbm %>%
    group_by(Table) %>%
    summarize(Sum = sum(Count))
tbm %>%
    group_by(Table, Taxon == "Other") %>%
    summarize(Sum = sum(Count))

# The sample names in the count tables and design tables don't match and are
# cumbersome. We will rename the samples to the format `sP-B` where P is the
# plate and B is the barcode number. For the count sample names, we can parse
# the existing sample names to get the plate and barcode numbers, following the
# matching pattern used in the Brooks SI. Sample names have the format
# `TRUTH{P}_{B}_?-?` where the bracketed P and B are the plate and barcode
# numbers (the final ?-? seems to be another representation of the barcode).
tbm <- tbm %>%
    mutate(Plate_Barcode = str_extract(Sample, "([1-6])_([0-9]+)")) %>%
    separate(Plate_Barcode, c("Plate", "Barcode")) %>%
    mutate(Sample = paste0('s', Plate, '-', Barcode)) %>%
    arrange(Plate, Barcode)
tbm

# Next, get a more convenient form of the sample metadata from the design table
sam <- design %>%
    select(-SampleName) %>%
    gather("Design_name", "Proportion", mock_taxa$Design_name) %>%
    filter(Proportion > 1e-4) %>%
    left_join(mock_taxa, by = "Design_name") %>%
    group_by(Experiment, Plate, Barcode) %>%
    summarize(Num_species = length(Taxon), Species_list = list(Taxon)) %>%
    mutate(Species_list = map_chr(Species_list, paste, collapse = ";")) %>%
    ungroup
sam <- sam %>%
    mutate(Sample = paste0('s', Plate, '-', Barcode))
sam <- sam %>%
    mutate(Mixture_type = case_when(
            Experiment == "Extraction" ~ "Cells",
            Experiment == "PCR" ~ "DNA",
            Experiment == "Seq" ~ "PCR_product"
        )
    )
sam

# Clean up some more and save the counts and sample tibbles
brooks2015_counts <- tbm %>%
    select(Sample, Taxon, Table, Count)
brooks2015_sample_data <- sam %>%
    select(Sample, Plate, Barcode, Mixture_type, Num_species, Species_list)
brooks2015_counts 
brooks2015_sample_data 
usethis::use_data(brooks2015_counts, brooks2015_sample_data)


# The above tibbles are what is used in our analysis, with the actual
# composition extracted from the sample data "Species_list". For those who want
# the observed and actual abundances in "OTU table" format for their own
# analysis, you can run the below code, to create tables with samples as rows
# and taxa as columns. Replace "/tmp/" with your chosen destination folder.

# observed <- brooks2015_counts %>%
#     filter(Table == "above") %>% 
#     select(-Table) %>%
#     spread(Taxon, Count) %>%
#     # Move "Other" to last column
#     select(-Other, everything())
#
# actual <- crossing(Sample = observed$Sample, Taxon = mock_taxa$Taxon) %>%
#     left_join(brooks2015_sample_data, by = "Sample") %>%
#     mutate(Present = str_detect(Species_list, Taxon) %>% as.integer) %>%
#     select(Sample, Taxon, Present) %>%
#     spread(Taxon, Present)
#
# write_csv(observed, "/tmp/brooks2015-observed.csv")
# write_csv(actual, "/tmp/brooks2015-actual.csv")
# write_csv(brooks2015_sample_data, "/tmp/brooks2015-sample-data.csv")

# Download GTDB and rrnDB -----------------------------------------------------

# For our analysis of 16S copy numbers we will make use of the GTDB tree and
# metadata table, and the rrnDB table of 16S copy numbers linked to NCBI
# species names. Here we simply download these files, to be used later in in
# the Rmd document `analysis/brooks2015_species_info.Rmd`

dotenv::load_dot_env(here("data-raw", ".env"))
data_path <- Sys.getenv("DATA_PATH")

# GTDB
dir.create(file.path(data_path, "gtdb"))
urls <- c(
    "https://data.ace.uq.edu.au/public/gtdb/release86/bac120_r86.2.tree",
    "https://data.ace.uq.edu.au/public/gtdb/release86/bac_metadata_r86.tsv"
)
fns <- file.path(data_path, "gtdb", basename(urls))
walk2(urls, fns, download.file)

# rrnDB
dir.create(file.path(data_path, "rrndb"))
urls <- c("https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.5.tsv.zip")
fns <- file.path(data_path, "rrndb", basename(urls))
walk2(urls, fns, download.file)
