#Author: Sam Czerski
#Date: 8/26/2021
#Description: This script takes input as a multiple sequence
#alignment file (eg. otu_45_anaerococcus.aln) and remove extra
#header info for constructing a tree in newick format.

#imports
import glob

#Takes an ASV fasta file as input and returns a fasta file of filtered records based on a chosen "size" parameter
def remove_extra_header_info(alignment_file):
    handles = open(alignment_file, "r")
    handles_list = []
    for line in handles:
        if line.startswith(">"):
            handles_list.append(line)

    handles.close()
    edit_headers(handles_list)



def edit_headers(handles_list):
    new_headers = []
    ncbi_headers = []
    asv_headers = []

    for element in handles_list:
        if element.startswith(">NR"):
            handle_id_items = element.split(sep=';')
            ncbi_headers.append(handle_id_items[0])

        elif element.startswith(">"):
            handle_id_items = element.split(sep=';')
            asv_headers.append(handle_id_items[0])

    handles = open(aln_file[0], "r")
    i = 0
    j = 0
    with open("asv_alignment_w_edited_header.aln", "w") as fr:
        for line in handles:
            if line.startswith(">NR") == True:
                fr.write(line.replace(str(line), ncbi_headers[i]) + "\n")
                i+=1
            elif line.startswith(">") == True:
                fr.write(line.replace(str(line), asv_headers[j]) + "\n")
                j+=1
            else:
                fr.write(line)
#RIGHT NOW, THE ID I WANT AS THE HEADER IS BEING ADDED TO THE FIRST LINE OF THE ALIGNMENT AND THE OLD HEADER IS NOT REPLACED...
    #NEED TO FIGURE OUT HOW TO FIX THIS.

    handles.close()

    return new_headers

#main function
if __name__ == '__main__':
    #Get Input
    aln_file = glob.glob("*.aln")
    remove_extra_header_info(aln_file[0])
    print("\nFasta headers have been edited and wrote to asv_alignment_w_edited_header.aln")
    print("Fin")
