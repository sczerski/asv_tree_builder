#! /bin/bash

ORGANISM=$1
ref_db="/data/shared/homes/sam/apps/mcsmrt/all_required_mcsmrt_files/16sMicrobial_ncbi_lineage_reference_database.udb"

# FL
usearch -utax ${ORGANISM}_ncbi_fl.fasta -db ${ref_db} -utax_cutoff 0.8 -strand both -utaxout ${ORGANISM}_ncbi_fl.utax

# V4
usearch -utax ${ORGANISM}_ncbi_v4.fasta -db ${ref_db} -utax_cutoff 0.8 -strand both -utaxout ${ORGANISM}_ncbi_v4.utax

# ASV
usearch -utax ${ORGANISM}_asv_v4.fasta -db ${ref_db} -utax_cutoff 0.8 -strand both -utaxout ${ORGANISM}_asv_v4.utax
