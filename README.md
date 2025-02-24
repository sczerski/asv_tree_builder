# asv_tree_builder
Visualize and understand the evolutionary relationship among ASVs.

Note that due to the upload structure on github, after cloning/downloading, you will need to move the "V4" folder into the "all_required_tree_building_files_FL16SDB" folder so the path will look like this: "all_required_tree_building_files_FL16SDB/V4" . Please contact me if you have any questions or issues about this.

# Required Software
You must have the following software installed on your machine:
R >= version 4
Python >= version 3.12.2
mafft >= v7.505 (2022/Apr/10)
RAxML-NG- >= v 1.1 (29.11.2021)
**RECOMMENDED**
I highly recommend simply using the provided snake.yaml file to clone the environment which I use to run this pipeline.

# Inputs and Description of Pipeline
This script can be used following the successful generation of an ASV dataset using the DADA2 pipeline. The first input is the ASV dataset. The other input is the name of an organism, at the genus level. This string will be used for obtaining reads from the NCBI 16S database for said organism. All output files can be found in the output directory, output_files. The main output from this pipeline is a ORGANISM_tree.png file. A description of the pipeline is as follows: obtain reads from ncbi, filter/combine ASV and ncbi reads, make annotation table, align reads using Mafft einsi multiple sequence alignment algorithm, edit headers in alignment file, build tree with RAxML-ng algorithm, visualize tree using R, clean up output files.

# Note on NCBI data retrieval
The asv_build_tree pipeline gets sequences directly from NCBI database using esearch from the entrez edirect utility set. These don't seem to pull all of the sequences available on NCBI, likely they are not in sync. However, this script should also work with entries downloaded directly from NCBI, which I recommend, as all available sequences will be included. Just name the file exactly as it appears in the Snakefile and put it in the directory with the Snakefile. Snakemake will recognize this rule as being completed and continue with the next rule. IMPORTANT: When downloading seqs from NCBI web-version, do NOT select the box to show/include the GI number. These are no longer used by NCBI. 

Here is a link to search for your Genus of interest: https://www.ncbi.nlm.nih.gov/nuccore?term=33175%5BBioProject%5D+OR+33317%5BBioProject%5D

Simply search "Streptococcus" after the bioproject ID's in the searchbar, and select "Send to:", "Choose Destination": File, Format: FASTA, Sort by: Sequence Length (default). Then select "Create File"

The asv_build_tree_FL16SDB pipeline produces the same output, but is designed to deal with a specific version of the full-length 16S ref seq database. This is the database used to classify reads in MCSMRT (i.e. 16sMicrobial_lineage.fasta, or something like this). This file has been formatted for making a usearch database, which has an annoying header. There are two different scripts I needed to make to deal with this, as well as edit subsequent scripts to not break the pipeline.
