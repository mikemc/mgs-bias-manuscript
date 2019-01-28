# Setup -----------------------------------------------------------------------

library(tidyverse)

# Location for downloading files and storing intermediate files. Create a
# symlink if want to store large data (the sequence reads) in a different
# location)
dl_path <- file.path(here::here(), "data-raw", "costea2017-intermediate")
if (!dir.exists(dl_path)) {
    dir.create(dl_path)
}

# Sample information-----------------------------------------------------------

## Download the sample metadata from the supplemental file provided by the
## authors as Excel spreadsheets
# Supplementary Data 2 : Members and composition of mock community
# https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S5.xlsx
# Supplementary Data 3 : Sample description
# https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S6.xlsx
urls <- c("https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S5.xlsx",
    "https://media.nature.com/original/nature-assets/nbt/journal/v35/n11/extref/nbt.3960-S6.xlsx")
names(urls) <- c("mock_composition.xlsx", "sample_description.xlsx")
fns <- file.path(dl_path, names(urls))
walk2(urls, fns, download.file)
list.files(dl_path)

## Import the metadata for the Phase 3 experiment
# The metadata for the three experimental phases of the Costea2017 experiment
# is given on a single sheet in separate rectagular ranges. I'll extract just
# the Phase 3 metadata.
sam <- readxl::read_xlsx(fns[2], range = "H3:J32")
print(sam, n = Inf)
# The `Individual` identifier `Mock-only` can be unambiguously shortened to
# just `M` to give all `Individual`'s a uniform 1-character format.
sam <- sam %>%
    mutate(
        Individual = case_when(
            Individual == "Mock-only" ~ "M",
            TRUE ~ Individual
        )
    )
# We will use simpler sample names, keeping the original as SI_sample
sam <- sam %>%
    mutate(
        SI_sample = Sample,
        Sample = paste0(Protocol, Individual)
    )
print(sam, n = Inf)

## Sequence file metadata from the ENA
# https://www.ebi.ac.uk/ena/data/view/PRJEB14847
fn <- file.path(dl_path, "PRJEB14847.tsv")
download.file("https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB14847&result=read_run&download=txt", fn)
seqtb <- readr::read_tsv(fn)

## Link to ENA run accessions to connect samples to their sequence data
# The `SI_sample` in `sam` can be matched to the `library_name` in `seqtb`:
all(paste0(sam$SI_sample, "_DA") %in% seqtb$library_name)
sam <- sam %>%
    mutate(library_name = paste0(SI_sample, "_DA"))
sam <- left_join(sam, seqtb %>% select(run_accession, library_name),
    by = "library_name")
sam <- sam %>%
    select(-library_name) %>%
    rename(Run_accession = run_accession)
print(sam, n = Inf)

# This is the sample data we will use for our analysis, so save it to be loaded
# with `data`.
costea2017_sample_data <- sam
usethis::use_data(costea2017_sample_data)
# For analysis or inspection outside of this package, save to csv:
# write_csv(costea2017_sample_data, file.path(dl_path, "costea2017_sample_data.csv"))

##  Mock community composition
mock <- readxl::read_xlsx(fns[1], sheet=1) %>%
    rename(Taxon = `Bacterial species`) %>%
    mutate(Taxon = str_extract(Taxon, "\\S+ \\S+"),
        Taxon = str_replace(Taxon, " ", "_")) %>%
    select(-`NCBI tax ID`)
head(mock)
costea2017_mock_composition <- mock
usethis::use_data(costea2017_mock_composition)

# Download reads from the ENA--------------------------------------------------

## Setup
# Directory to download reads to
reads_path <- file.path(dl_path, "reads")
if (!dir.exists(reads_path)) {
    dir.create(reads_path)
}
# Paths to the aspera connect `ascp` program and the private key file to be
# used with the -i option
dotenv::load_dot_env(file.path(here::here(), "data-raw", ".env"))
Sys.getenv("ASPERA_ASCP")
Sys.getenv("ASPERA_KEY")

## Download with ascp (aspera connect command line tool)
# If don't have asprea, can download with wget (see below).
# The aspera urls are in the format "url/for/read1;url/for/read2" in `seqtb`
# and so we first split all out into a single list
aspera_urls <- seqtb$fastq_aspera %>% str_split(";", simplify=TRUE) %>% c
commands <- paste(
    Sys.getenv("ASPERA_ASCP"),
    "-QT -l 300m -P33001 -i", 
    Sys.getenv("ASPERA_KEY"),
    paste0("era-fasp@", aspera_urls),
    reads_path
    )
walk(commands, system)

# ## Alternately, download with wget:
# ftp_urls <- tb$fastq_ftp %>% str_split(";", simplify=TRUE) %>% c
# dir.create(file.path(data_path, "reads"), recursive = TRUE)
# commands <- paste("wget", "-P", file.path(data_path, "reads"), 
#     paste0("ftp://", ftp_urls)
# walk(commands, system)

## Check dowloaded files against md5sums
downloads <- aspera_urls %>% 
    str_extract("ERR[0-9]*_[1-2]\\.fastq\\.gz") %>%
    file.path(reads_path, .)
md5sums_expected <- seqtb$fastq_md5 %>% str_split(";", simplify=TRUE) %>% c
md5sums_actual <- tools::md5sum(downloads)
# Fraction of files that were successfully downloaded and gave an md5
mean(!is.na(md5sums_actual))
# Check that the downloaded files match the expected md5
all(md5sums_expected == md5sums_actual, na.rm = TRUE)
