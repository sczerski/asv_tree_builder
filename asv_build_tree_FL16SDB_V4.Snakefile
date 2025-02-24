#You must edit this line below to set the ORGANISM variable to organism of interest. This will be used for pulling FL16S reads from NCBI. This should be a GENUS.
ORGANISM = 'Streptococcus' 

rule target:
     input: expand("{orgn}_v4_tree.png", orgn=ORGANISM)
     params: organism=ORGANISM
     shell: "./all_required_asv_tree_building_files_FL16SDB/tidy_up.sh {params.organism}"

rule get_ncbi_seqs:
     output: expand("{orgn}_ncbi_fl.fa", orgn=ORGANISM)
     params: organism=ORGANISM
     threads: 1
     message: "Pulling NCBI full-length 16S reference sequences from database for organism of interest."
     shell: "./all_required_asv_tree_building_files_FL16SDB/get_16S_seqs_from_database.sh {params.organism}"

# I need to swap this rule with the following one... testing first. This requires changing the inputs/outputs
#rule edit_ncbi_seqs_headers:
#     input: expand("{orgn}_ncbi_fl.fa", orgn=ORGANISM)
#     output: expand("{orgn}_ncbi_fl.fasta", orgn=ORGANISM)
#     threads: 1
#     message: "Editing NCBI full-length 16S reference sequences headers."
#     params: organism=ORGANISM
#     shell: "./all_required_asv_tree_building_files_FL16SDB/reduce_ncbi_fasta_headers.sh {params.organism}"

rule extract_v4_reads:
     input: expand("{orgn}_ncbi_fl.fa", orgn=ORGANISM)
     output: expand("{orgn}_ncbi_v4.fa", orgn=ORGANISM)
     threads: 1
     message: "Extracting V4 region from NCBI FL16SDB (and ASVs, if specified)."
     params: organism=ORGANISM
     shell: "./all_required_asv_tree_building_files_FL16SDB/V4/extract_v4_reads.sh {params.organism}"

rule edit_ncbi_seqs_headers:
     input: expand("{orgn}_ncbi_v4.fa", orgn=ORGANISM)
     output: expand("{orgn}_ncbi_v4.fasta", orgn=ORGANISM)
     threads: 1
     message: "Editing NCBI full-length 16S reference sequences headers."
     params: organism=ORGANISM
     shell: "./all_required_asv_tree_building_files_FL16SDB/reduce_ncbi_fasta_headers.sh {params.organism}"


rule classify_fl_and_v4_reads:
     input: expand("{orgn}_ncbi_v4.fasta", orgn=ORGANISM)
     output: expand("{orgn}_ncbi_v4.utax", orgn=ORGANISM)
     threads: 16
     message: "Taxonomically classifying FL and V4 reads."
     params: organism=ORGANISM
     shell: "./all_required_asv_tree_building_files_FL16SDB/V4/classify_all_reads.sh {params.organism}"

rule filter_combine_reads:
     input: expand("{orgn}_ncbi_v4.fasta", orgn=ORGANISM)
     output: "all_filtered_seqs_w_id_v4.fasta", "all_reads_info_v4.tsv"
     threads: 1
     message: "Filtering and combining NCBI and ASV sequences, creating all_reads_info.tsv file"
     shell: "python3 all_required_asv_tree_building_files_FL16SDB/V4/filter_combine_asv_reps_ncbi_refs_and_make_tsv_FL16SDB_v4.py"

rule make_annotation_tab:
     input: "all_reads_info_v4.tsv"
     output: "annotation_tab_v4.tsv"
     threads: 1
     message: "Creating annotation table from all_reads_info.tsv"
     shell: "python3 all_required_asv_tree_building_files_FL16SDB/V4/make_annotation_tab_V4.py"
     
rule align_reads:
     input: "all_filtered_seqs_w_id_v4.fasta"
     output: expand("{orgn}_v4.aln", orgn=ORGANISM)
     params: organism=ORGANISM
     threads: 16
     message: "Aligning multifasta file with mafft einsi algorithm"
     shell: "einsi all_filtered_seqs_w_id_v4.fasta > {params.organism}_v4.aln"

rule edit_aln_header_info:
     input: expand("{orgn}_v4.aln", orgn=ORGANISM)
     output: "asv_alignment_w_edited_header_v4.aln"
     threads: 1
     message: "Removing extraneous header info from alignment for constructing newick string"
     shell: "python3 all_required_asv_tree_building_files_FL16SDB/V4/edit_fasta_header_info_for_newick_V4.py"

rule build_tree:
     input: "asv_alignment_w_edited_header_v4.aln"
     output: expand("{orgn}_v4.raxml.bestTree", orgn=ORGANISM)
     threads: 16
     message: "Building tree from alignment. NOTE: if this step fails, you might have duplicate headers. Try running command in terminal."
     params: organism=ORGANISM
     shell: "raxml-ng --all --threads 16 --msa {input} --prefix {params.organism}_v4 --model GTR+G --tree pars{{10}} --bs-trees 100 --msa-format fasta"

#if you get, "list index out of range" you might have duplicate ASV names.
#You can use FastTree in the build_tree rule, too, or whatever you want to use.

rule vizualize_with_R:
     input: tree = expand("{orgn}_v4.raxml.bestTree", orgn=ORGANISM), annot_tab = "annotation_tab_v4.tsv" #tree = expand("{orgn}.tree", orgn=ORGANISM), annot_tab = "annotation_tab.tsv"
     output: expand("{orgn}_v4_tree.png", orgn=ORGANISM)
     threads: 1
     message: "Importing data to R, midpoint rooting tree, vizualizing with ggtree"
     shell: "Rscript --vanilla all_required_asv_tree_building_files_FL16SDB/create_tree_image.R {input.tree} {input.annot_tab}"

