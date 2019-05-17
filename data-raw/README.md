## Setup

Requirements: The `tidyverse`, `here`, and `dotenv` R packages; Aspera connect
client for faster sequence downloads, otherwise can download via FTP client
like wget; Bowtie2 and Metaphlan2 for profiling shotgun data.

#### `.env` file

Several paths need to be set in a file `data-raw/.env`. Copy the example
`.env_example` to `.env` and edit the paths. The Aspera Connect paths needed to
download the Costea2017 reads with `ascp`, but are not needed if downloading
via FTP. `DATA_PATH` must be set to to indicate where the intermediate files
will be downloaded to.

## Brooks2015

#### brooks2015.R

Download the sample information and read-count tables from the SI of
Brooks2015; format and import into the R package for loading with `data()`; and
download files from the GTDB and rrnDB for later use in estimating genome size
and 16S copy number.

## Costea2017

#### 0-costea2017.R

Download the reads, metadata, and mock composition; and format the metadata and
mock composition and save as Rdata objects in the `data/` folder to be accessed
with `data()` from within the package.

#### 1-costea2017-metaphlan2.sh

Run Metaphlan2 with default settings, saving the Bowtie2 output for the next
step.

To profile all accessions on a computer cluster running Slurm,
```sh
nproc=8
for i in {1971003..1971031}
do
    accession=ERR$i
    sbatch -c $nproc costea2017-metaphlan2.sh $accession $nproc
done
```

#### 2-costea2017-metaphlan2-custom.sh

Re-run Metaphlan2 with custom settings starting from Bowtie2 output. This is
fast because the reads have already been mapped.

Afterwards, merge the profiles into a single table using
`merge_metaphlan_tables.py` from the `utils` folder of the `metaphlan2`
package.

This script can used used to submit all accessions to our computing cluster:
```sh
nproc=2
for i in {1971003..1971031}
do
    accession=ERR$i
    sbatch -c $nproc costea2017-metaphlan2-custom.sh $accession $nproc
    # Or run locally with
    # bash costea2017-metaphlan2-custom.sh $accession $nproc
done
```

Afterwards, merge all accessions into a single table by running
```sh
/path/to/merge_metaphlan_tables.py *_profiled_metagenome_custom.tsv > costea2017_metaphlan2_custom_profiles.tsv
```

#### 3-costea2017-import-metaphlan2-profiles.R

Import the Metaphlan2 profiles with custom settings into the R package so that
it can be accessed through `data()`.
