#! /bin/bash

#Author: Sam Czerski
#Description: gets user-defined spp reads from database. 'genera' is a filename of a file with a string of the genus you find to get ref seqs for.

ORGANISM=$1
FILE=genera.txt

#test to make sure organism was input
[[ -n $ORGANISM ]] && printf "\nYou are building a tree for $ORGANISM \n" || printf "\nArgument error. Please edit snakefile and assign a valid Genus name. \n"

#check to see if a file named, genera exists, otherwise, make it.
[[ -f $FILE ]] || printf ${ORGANISM} > genera.txt
      
#get sequences for genus from database
usearch --fastx_getseqs all_required_asv_tree_building_files_FL16SDB/ncbi_database/16sMicrobial_ncbi_lineage.fasta -labels genera.txt -fastaout ${ORGANISM}_ncbi_fl.fa -label_substr_match

#count number of unique sequences
UNIQ_SPP="grep '>' ${ORGANISM}_ncbi_fl.fa |cut -d':' -f8|cut -d';' -f1|sort|uniq -c|sort|wc -l"
TOTAL_SPP="grep '>' ${ORGANISM}_ncbi_fl.fa"

printf "\n${UNIQ_SPP} unique species found in ${TOTAL_SPP} sequences for $ORGANISM. \n" 
