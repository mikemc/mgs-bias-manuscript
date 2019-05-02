# Setup -----------------------------------------------------------------------

library(tidyverse)
library(here)

# Location for downloading files and storing intermediate files. Create a
# symlink if want to store large data (the sequence reads) in a different
# location)
dl_path <- here("data-raw", "costea2017-intermediate")
if (!dir.exists(dl_path)) {
    dir.create(dl_path)
}

# Download sample information -------------------------------------------------

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

## Sequence file metadata from the ENA
# https://www.ebi.ac.uk/ena/data/view/PRJEB14847
fn <- file.path(dl_path, "PRJEB14847.tsv")
download.file("https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB14847&result=read_run&download=txt", fn)

#  Mock community composition -------------------------------------------------
mock <- readxl::read_xlsx(file.path(dl_path, "mock_composition.xlsx"), 
    sheet = 1) %>%
    rename(Taxon = `Bacterial species`) %>%
    mutate(Taxon = str_extract(Taxon, "\\S+ \\S+"),
        Taxon = str_replace(Taxon, " ", "_")) %>%
    select(-`NCBI tax ID`)
head(mock)
costea2017_mock_composition <- mock
usethis::use_data(costea2017_mock_composition)

# Link sample metadata to sequence data ---------------------------------------

seqtb <- readr::read_tsv(file.path(dl_path, "PRJEB14847.tsv"))

# The metadata for the three experimental phases of the Costea2017 experiment
# is given on a single sheet in separate rectagular ranges.
fn <- file.path(dl_path, "sample_description.xlsx")
ranges <- c("1" = "A3:B192", "2" = "D3:F77", "3" = "H3:J32")
sam.all <- map_dfr(ranges, ~readxl::read_xlsx(fn, range = .),
    .id = "Phase")
sam.all %>%
    group_by(Phase) %>%
    count
# Keep copy of the sample names used in the SI spreadsheet so that we can
# change the sample names later
sam.all <- sam.all %>% 
    mutate(SI_sample = Sample)

# There are a few strings for each run that will help us in matching with the
# metadata
seqtb <- seqtb %>%
    mutate(
        String1 = str_extract(library_name, "(AWF)|(BYQ)"),
        String2 = str_extract(library_name, "(?<=(AWF|BYQ)_)[A-Z]{4,6}"),
        String3 = str_extract(sample_alias, 
            "([ABC][1I][_ -]{1,2}[:digit:]{3}([_ -][:digit:])?)|C[12]"),
    )
# seqtb$String1
# seqtb$String2
# seqtb$String3

# Note, there is considerable heterogeneity in the formatting of String3,
# including an apparent typo of AI that should be A1. The first letter [ABC]
# indicates (I think) the specimen. That is followed by a 1, an I (typo), or a
# 2, with the 2 only appearing for C2. 

# We are not interested in the C samples since these are not included in the
# sample metadata so let's drop these now. 
seqtb <- seqtb %>%
    filter(is.na(String3) | (str_sub(String3, 1, 1) != "C"))


# For the purposes of matching to the samples, let's fix up String3: fix the
# typo and make the separator formatting (of underscore, space, and dash)
# uniform 
seqtb <- seqtb %>%
    mutate(
        String3 = str_replace(String3, "AI", "A1"),
        String3 = str_replace_all(String3, "[_ -]+", "_"),
    )

# HERE. Challenge is to use the Strings to identify the phase

# Phase 3 samples are easy to identify and get the sample names of
seqtb <- seqtb %>%
    mutate(
        Phase = ifelse(String1 == "BYQ", "3", NA),
        SI_sample = ifelse(Phase == "3", paste("BYQ", String2, sep = "_"), NA),
    )
seqtb %>%
    filter(Phase == 3) %>%
    .$SI_sample

# But for Phase 1 and 2 we can't easily tell, and will need to match against
# the sample metadata. 

## To do this, first create the names we would see if these samples were from
## the particular phase. 
seqtb <- seqtb %>%
    mutate(
        Name_if_1 = paste(String3, String2, sep = "_"),
        Name_if_2 = String3,
    )
## Then, look for matches in the known sample names.
samples.phase1 <- sam.all %>% filter(Phase == 1) %>% .$SI_sample
samples.phase2 <- sam.all %>% filter(Phase == 2) %>% .$SI_sample
# Check what's missing in Phase 1
setdiff(samples.phase1, seqtb$Name_if_1)
# Both appear to be typos. We can find these two missing samples if we suppose:
# - Sample A1_081_A1OSW should be A1_081_AIOSW (typo of 1 <-> I)
# - Sample A1_054_HLOSW corresponds to String2 = HLOSW & String3 = A1_045; a
#   typo of (054 <- 045).

# Let's adjust Name_if_1 to match the sample metadata
seqtb <- seqtb %>%
    mutate(Name_if_1 = case_when(
        Name_if_1 == "A1_081_AIOSW" ~ "A1_081_A1OSW",
        Name_if_1 == "A1_045_HLOSW" ~ "A1_054_HLOSW",
        TRUE ~ Name_if_1
        )
    )
setdiff(samples.phase1, seqtb$Name_if_1)

# Check what's missing in Phase 2
setdiff(samples.phase2, seqtb$Name_if_2)
# all good!

seqtb <- seqtb %>%
    mutate(
        Phase = case_when(
            Name_if_1 %in% samples.phase1 ~ "1",
            Name_if_2 %in% samples.phase2 ~ "2",
            TRUE ~ Phase),
        SI_sample = case_when(
            Phase == 1 ~ Name_if_1,
            Phase == 2 ~ Name_if_2,
            Phase == 3 ~ SI_sample)
        ) %>%
    arrange(Phase)
# Check that the extracted names look good
seqtb %>%
    group_by(Phase) %>%
    top_n(5, SI_sample) %>%
    select(Phase, SI_sample)
# Check the counts by Phase and make sure no missing sample names
seqtb %>%
    group_by(Phase, is.na(SI_sample)) %>%
    count

# Note, we expect 189 samples for phase 1 and 74 samples for phase 2,
# indicating extra samples or duplicate samples.

dups <- seqtb %>%
    filter(duplicated(SI_sample) | 
            duplicated(SI_sample, fromLast = TRUE)
        ) %>%
    arrange(Phase, SI_sample, run_accession)
# The library_name and sample_alias are duplicated; while the run_alias and
# run_accessions differ
dups %>%
    select(Phase, SI_sample, library_name, sample_alias) %>%
    group_by(Phase) %>%
    top_n(4, SI_sample)
dups %>%
    select(Phase, SI_sample, run_accession, run_alias) %>%
    group_by(Phase) %>%
    top_n(4, SI_sample)
# (Differences appear in the HAYK8ADXX part of the run alias)
# For Phase 1, there can be differences in library layout and length
dups %>%
    select(Phase, SI_sample, library_layout, nominal_length) %>%
    group_by(Phase) %>%
    top_n(4, SI_sample)

# Unclear why there are multiple libraries for these samples.


# Let's get sample data for just Phases 2 and 3

# Get a minimal set of variables for joining with the sample metadata
seqtb0 <- seqtb %>%
    select(SI_sample, Phase, run_accession, 
        instrument_model, library_layout, nominal_length) %>%
    rename_at(vars(-SI_sample), str_to_sentence)
        # library_name, run_alias,
        # fastq_bytes, fastq_md5, fastq_ftp, fastq_aspera)

# Once we join with the sample data, we will have multiple rows for the
# SI_sample that are duplicated in the seqtb.
sam <- sam.all %>%
    filter(Phase %in% c(2, 3)) %>%
    left_join(seqtb0, by = c("Phase", "SI_sample")) %>%
    arrange(Phase, SI_sample, Run_accession) %>%
    group_by(SI_sample) %>%
    mutate(
        id = rank(Run_accession),
        n = n(),
        Sample = ifelse(n > 1,
            paste(SI_sample, id, sep = "_"), 
            SI_sample),
        ) %>%
    select(Sample, Phase, everything())

sam %>%
    # filter(duplicated(SI_sample) | duplicated(SI_sample, fromLast = TRUE)) %>%
    filter(n > 1) %>%
    select(Phase, Sample, SI_sample, Run_accession, id, n)

sam <- sam %>%
    mutate(
        Individual = case_when(
            Phase == 2 ~ str_sub(SI_sample, 1, 1),
            Phase == 3 ~ Individual),
        Sample = case_when(
            # Phase == 2 ~ paste0("Ph2", "L", Lab, Protocol, Individual),
            Phase == 2 ~ Sample,
            Phase == 3 ~ paste0(Protocol, Individual))
        )

write_csv(sam, here("data", "costea2017-phases23-sample-data.csv"))

# Download reads from the ENA for Phase 2 and Phase 3 -------------------------

## Setup

# Load .env file
dotenv::load_dot_env(here("data-raw", ".env"))
# Directory to download reads to
Sys.getenv("DATA_PATH")
# Paths to the aspera connect `ascp` program and the private key file to be
# used with the -i option
Sys.getenv("ASPERA_ASCP")
Sys.getenv("ASPERA_KEY")

# TODO: consider putting all the intermediate files into DATA_PATH as well

reads_path <- file.path(Sys.getenv("DATA_PATH"), "costea2017", "reads")
if (!dir.exists(reads_path)) {
    dir.create(reads_path)
}

## Download with ascp (aspera connect command line tool)
# If don't have aspera, can download with wget (see below).
# The aspera urls are in the format "url/for/read1;url/for/read2" in `seqtb`
# and so we first split all out into a single list
# aspera_urls <- seqtb$fastq_aspera %>% str_split(";", simplify=TRUE) %>% c
aspera_urls <- seqtb %>%
    filter(Phase %in% c(2, 3)) %>%
    .$fastq_aspera %>% 
    str_split(";", simplify=TRUE) %>%
    c
commands <- paste(
    Sys.getenv("ASPERA_ASCP"),
    "-QT -l 300m -P33001 -i", 
    Sys.getenv("ASPERA_KEY"),
    paste0("era-fasp@", aspera_urls),
    reads_path
    )
walk(commands, system)

# ## Alternately, download with wget:
# ftp_urls <- seqtb$fastq_ftp %>% str_split(";", simplify=TRUE) %>% c
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

# Submit all jobs to sbatch ---------------------------------------------------

# TODO: Create sbatch commands calling metaphlan or motus

# Check all mOTUS2 profiles were generated ------------------------------

accs <- seqtb %>%
    filter(Phase == 2) %>%
    .$run_accession

fns <- paste0(accs, "-counts.motus") %>%
    file.path("~/data/mgs_bias/costea2017/motus2", .)
names(fns) <- accs

all(file.exists(fns))

# Afterwards, merge all with
#> FLIST=$(echo ERR*-counts.motus | tr ' ' ,)
#> motus merge -i $FLIST -o costea2017-counts.motus

# Get commands for calling mOTUs2 map_snv

accs %>%
    paste("sbatch costea2017-motus2-map_snv.sh", .) %>%
    write_lines(here("data-raw", "sbatch-motus2-map_snv.sh"))
