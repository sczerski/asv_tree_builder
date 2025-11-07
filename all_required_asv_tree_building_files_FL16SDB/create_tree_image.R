#!/usr/bin/env Rscript

#This Rscript is used to create a ggtree image file 
#usage: Rscript --vanilla create_tree_image.R tree_file annotation_rdata tip_color_variable
#e.g. Rscript --vanilla create_tree_image.R all_Azotobacter_FL16S.tree joined_Azotobacter.txt Species

args = commandArgs(trailingOnly=TRUE)

#This is a bit of an incantation to install pacman, which confusingly helps you install libraries if you need them
#if (!require("pacman")) install.packages("pacman", "CRAN")
#library(pacman)

#source("http://bioconductor.org/biocLite.R")
#biocLite("ggtree")
#biocLite("phyloseq")
#We will use pacman to install and load packages:
#pacman::p_load(ape, ggtree, phangorn, doParallel, seqinr, ggrepel, svglite, colorspace)
pacman::p_load(ape, ggtree, ggplot2, phangorn, svglite, colorspace)

#Lets add some paralellization up in here
#registerDoParallel(cores = 8)
ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw() }


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

#Initialize variables
tree_file = args[1]
annotation = args[2]
#tips = args[3]

tree_annotation <- readRDS(paste0(annotation))

tree <- read.tree(tree_file)
tree_mid <- midpoint(tree)


ggtree(tree_mid) %<+% tree_annotation +
  geom_tippoint(aes(color=type)) +
  theme(legend.position = "right") +
  ggtitle(tree_file)+ 
  geom_treescale(x=0,y=200, offset=2) +geom_tiplab(aes(label=label_w_strain_name))

out_name <- gsub(".tree", "", tree_file)

ggsave(file=paste0("final_utax_species_tree_images/",out_name, "_utax_species.png"))

