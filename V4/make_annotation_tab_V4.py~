#imports
import glob
import pandas as pd


def create_annotation_tab(all_reads_info_tsv):
    with open(all_reads_info_tsv, 'r') as a_tsv:
        labels_w_strain_names = []

        #for row in a_tsv:
            #print(row)
            #values = row.split('\t')
            #name = values[0]
            #genus_species = values[3]
            #label = values[6]
            #new_value = genus_species + "_" + label
            #labels_w_strain_names.append(new_value)

        #Create Headers
        header_list = ["name", "type", "genus_species", "passes", "num_identical", "label_w_strain_name"]
        #header_list ["name", "type", "genus_species", "passes", "num_identical", "note", "label"]

        #Add Headers to a new tsv file. I want to leave the original data untouched
        a_old_tsv = pd.read_csv(all_reads_info_tsv, delimiter = "\t", header = None)
        a_old_tsv.to_csv("annotation_tab.tsv", sep = "\t", header = header_list, index = False)

        #write new data
        #new_header = ["label_w_strain_name"]
        #df = pd.DataFrame(labels_w_strain_names)
        #df.to_csv("labels_w_strain_names.tsv", sep = "\t", header = new_header, index = False)

        #combine new column to data frame tsv with headers
        #pd.concat([pd.read_csv("pre_annotation_tab.tsv", delimiter="\t"), pd.read_csv("labels_w_strain_names.tsv", delimiter="\t")], axis = 1).to_csv("annotation_tab.tsv", sep = "\t", header = True, index = False)

#main function
if __name__ == '__main__':
    all_reads_info = glob.glob("all_reads_info.tsv")
    create_annotation_tab(all_reads_info[0])
