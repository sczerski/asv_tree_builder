#! /bin/bash

ORGANISM=$1
fwd_primer="GTGCCAGCMGCCGCGGTAA"
rev_primer="GGACTACHVGGGTWTCTAAT"

#use this if you want to extract the V4 region from the inferred ASV-16S sequence. 

#usearch11 -search_pcr2 all_asv_fa/${ORGANISM}_asv.fasta -fwdprimer ${fwd_primer} -revprimer ${rev_primer} -strand both -threads 32 -fastaout ${ORGANISM}_asv_v4.fasta

#the command below will link my v4 ASV sequences that I have already processed. This is so my script doesn't break.
ln -s all_asv_fa/${ORGANISM}_asv.fasta ${ORGANISM}_asv_v4.fasta

usearch11 -search_pcr2 ${ORGANISM}_ncbi_fl.fa -fwdprimer ${fwd_primer} -revprimer ${rev_primer} -strand both -threads 32 -fastaout ${ORGANISM}_ncbi_v4.fa
