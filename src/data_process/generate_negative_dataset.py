""" module: generate the negative dataset
    the RNA fragments of 88 base long are retrieved
    the psuedo  miRNA-mRNA pairs are from random selection avoiding 
    the combinations reported in miRTarBase 
"""
import os
import random
import pandas as pd
import numpy as np
from xlrd import open_workbook

# function: get the cdna location from whole genome cdna file
def parse_cdna_file(file_path,species_name):
    cdna_info_list = []  # store all the cdna information
    fd = open(file_path)
    # num_cdna = 0
    for line in fd:
        if line.startswith(">"):
            cdna_info = []   # species,transcript id, location, gene name
            cdna_info.append(species_name)
            # num_cdna += 1
            line = line.split()
    
            transcript_id = line[3].split(":")[1]  # ensamble transcript name
            cdna_info.append(transcript_id)
        
            cdna_info.append(line[2])       #location
        
            gene_name = line[6].split(":")[1]     # gene name
            cdna_info.append(gene_name)
        
            cdna_info_list.append(cdna_info)
    fd.close()
    return cdna_info_list


compliment_dict = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
# function: get the reverse compliment sequence of input sequence
def reverse_compliment(m_str):
    m_str = m_str[::-1]
    m_str = m_str.upper()
    result = ""
    for base in m_str:
        result += compliment_dict[base]
    return result


# function:get the sequence fragment from 3UTR,5UTR and ORF with proportion of 7:2:1
# the input seq should be longer than 400, the retrieved fragment is 88 nt
def get_random_fragment(seq):
    random_choice =random.choice([1,1,1,1,1,1,1,2,2,3])
    fragment_length = 88       # pseudo target sequence length
    random_start_point = random.randrange(1,50)
    if random_choice ==1:
        random_fragment = seq[-(fragment_length+random_start_point):-random_start_point] #3UTR
    elif random_choice ==2:
        random_fragment = seq[random_start_point:random_start_point+fragment_length] # 5UTR
    else:
        mid_point = int(len(seq)/2)
        mid_point_offset = mid_point + random_start_point
        random_fragment = seq[mid_point_offset:mid_point_offset+fragment_length]
    return random_fragment
    

# open cdna_info csv file and retrieve RNA fragment of 88 nt
# return csv file (columns:species,gene name,sequence fragment of transcript 
# (88 nt from 3UTR,5UTR and ORF with proportion of 7:2:1)
# make sure that you have installed samtools and indexed the genome. 
# You should modify the paths when run this mofunction on your machine
def retrieve_seq(input_csv_file,genome_file_path,max_num_seq):
    cdna_df = pd.read_csv(input_csv_file)

    cdna_np_array = np.array(cdna_df)   # dataframe to numpy array
    cdna_list = cdna_np_array.tolist()  # numpy array to list
    
    new_cdna_list = []    # store the new list with 3' sequence
    
    random.shuffle(new_cdna_list)   # shuffle the order
    for each_cdna_info in cdna_list:   # order,species,transcript id, location, gene name
        if len(new_cdna_list) > max_num_seq:
            break
        else:
            location_info = each_cdna_info[3].split(":")
            chromesome = location_info[2]
            transcript_start = int(location_info[3])
            transcript_end =  int(location_info[4])
            transcript_length = transcript_end - transcript_start
     
            # filter the length of transcript and chromesome
            if transcript_length < 400 or \
                   not str(chromesome).startswith\
                            (("1","2","3","4","5","7","8","9","M","X","Y","I","II","III","IV","V")):
                 continue
    
            # get the sequence from reference genome with samtools    
            location_for_samtools = location_info[2] + ":" + location_info[3] + "-" + location_info[4]
            # pass the parameters
            seq_each = os.popen("samtools faidx %s %s" % (genome_file_path,location_for_samtools)).read( ) 

            seq_each = seq_each.split()[1:]  # remove the first line
            seq_each = ''.join(seq_each)
            if seq_each.find("N") >=0 or seq_each.find("n") >=0:
            # print(seq_each)
                continue                     # remove sequence containing "N" or "n"
        
            # coding strand is the complimentary strand
            if location_info[5]=="-1":
                seq_each = reverse_compliment(seq_each)
        
            seq_fragment = get_random_fragment(seq_each)

            # store in the new list
            new_each_cdna_info = each_cdna_info
            new_each_cdna_info.append(seq_fragment)
        
            new_cdna_list.append(new_each_cdna_info[1:])

        
    return new_cdna_list


# function: parse the miRTarBase_MTI.xlsx file and output the positive_combination dictionary
def generate_positive_combination():
    workbook = open_workbook(r"../../data/miRTarBase_MTI.xlsx")  # open xls file
    sheet_name= workbook.sheet_names()  # print the names of sheet
    interation_data_sheet = workbook.sheet_by_name('miRTarBase')
    # print(sheet_name)
    # workbook = open_workbook(r"../../data/test_miRTarBase_MTI.xlsx")  # open xls file
    # interation_data_sheet = workbook.sheet_by_name('Sheet1')

    # store the target gene, miRNAs in a dictionary
    positive_combination = {}
    for row in range(1,interation_data_sheet.nrows):
        # print(row)
        gene_name = str(interation_data_sheet.cell(row,3).value)
        miRNA_name =  str(interation_data_sheet.cell(row,1).value)
        if gene_name not in positive_combination.keys():
            miRNA_name_list = []
            miRNA_name_list.append(miRNA_name)
            positive_combination[gene_name] = miRNA_name_list
        else:
            positive_combination[gene_name].append(miRNA_name)
            
    return positive_combination


# function:read mature_miRNA.fa and store in a dictionary
def generate_miRNA_dic():
    miRNA_name_list = []
    miRNA_seq_list = []
    miRNA_name_seq_dict = {}
    with open("../../data/mature_miRNA.fa") as fd:
        for line in fd:
            if line.startswith(">"):
                line = line.split()
                line = line[0][1:]       # get the miRNA name
                miRNA_name_list.append(line)
            else:
                miRNA_seq_list.append(line.replace("\n",""))
    # print(len(miRNA_name_list))
    # print(len(miRNA_seq_list))

    for i in range(len(miRNA_name_list)):
        miRNA_name_seq_dict[miRNA_name_list[i]] = miRNA_seq_list[i]
    return miRNA_name_seq_dict


# function:read the human_cdna_fragment.csv, mouse_cdna_fragment.csv.
# select the miRNA not in positive_combination dictionary[specific_gene_name]
# generate the negative data set,writing to csv file
# column = [species,miRNA_name, miRNA_seq, taget_gene,target_seq,miRNA_target_seq,classification]
def generate_pseudo_miRNA_target(file_path,miRNA_name_seq_dict,positive_combination,\
                                 start_word,species_name,target_file_path):
    cdna_fragment_df = pd.read_csv(file_path)
    cdna_fragment_array = np.array(cdna_fragment_df)    # convert to numpy array
    cdna_fragment_list = cdna_fragment_array.tolist()   # convert to list

    # find the pseudo miRNA which target the cdna of human and store in a new list 
    # generate pseudo miRNA-target combination
    pseudo_miRNA_target_list = []
    miRNA_name_list = list(miRNA_name_seq_dict.keys())
    # select three pseodo miRNAs for each gene
    for i in range(3):
        for each_fragment in cdna_fragment_list:
            pseudo_miRNA_target = []
            gene_fragmen_name = each_fragment[4] # gene name
         
            # some gene names may not in the positive_combination dict
            if gene_fragmen_name not in positive_combination.keys(): 
                miRNA_choice = random.choice(miRNA_name_list)
            else:
                miRNA_list_target = positive_combination[gene_fragmen_name]        # miRNAs list
                while True:
                    miRNA_choice = random.choice(miRNA_name_list)
                    if miRNA_choice.startswith(start_word) and miRNA_choice not in miRNA_list_target:
                        break
            # write pseudo miRNA_target combination to a list       
            pseudo_miRNA_target.append(species_name)                          # species
            pseudo_miRNA_target.append(miRNA_choice)                          # miRNA_name,
            pseudo_miRNA_target.append(miRNA_name_seq_dict[miRNA_choice])     # miRNA_seq
            pseudo_miRNA_target.append(gene_fragmen_name)                     # pseudo_taget_gene
            pseudo_miRNA_target.append(each_fragment[5])                      # pseudo_target_seq
            # pseudo_miRNA_target_seq
            pseudo_miRNA_target.append(miRNA_name_seq_dict[miRNA_choice]+each_fragment[5])
            pseudo_miRNA_target.append(0)                                        # classification
    
            pseudo_miRNA_target_list.append(pseudo_miRNA_target)

    print(pseudo_miRNA_target_list)
    print("write to file start:")
    pseudo_miRNA_target_df = pd.DataFrame(pseudo_miRNA_target_list,columns=\
                                            ["species","miRNA_name", "miRNA_seq", "pseudo_taget_gene",\
                                             "pseudo_target_seq", "pseudo_miRNA_target_seq","classification"])

    pseudo_miRNA_target_df.to_csv(target_file_path)   # writer to file
    print("done")




if __name__ == "__main__":
    """
    # the file path in my computer, you should download the file from emsamble and give the right path
    human_cdna_file = "/home/zheng/reference_data/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa"
    human_cdna_info_list = parse_cdna_file(human_cdna_file,"Homo sapiens")
    # convert the human_cdna_info_new into dataframe: species,transcript id, location and gene name
    human_cdna_info_df = pd.DataFrame(human_cdna_info_list,\
                                 columns=["species","transcript_id", "location", "gene_name"])
    human_cdna_info_df.to_csv("../../data/human_cdna_location.csv")   # writer to file

    # the file path in my computer, you should download the file from emsamble and give the right path
    mouse_cdna_file = "/home/zheng/reference_data/Mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa"
    mouse_cdna_info_list = parse_cdna_file(mouse_cdna_file,"Mus musculus")
    # convert the human_cdna_info_new into dataframe: species,transcript id, location and gene name
    mouse_cdna_info_df = pd.DataFrame(mouse_cdna_info_list,\
                                 columns=["species","transcript_id", "location", "gene_name"])
    mouse_cdna_info_df.to_csv("../../data/mouse_cdna_location.csv")   # writer to file

    # generate 40000 human cdna fragment and write to csv file
    print("retrieve human transcripts fragment start:")
    new_human_cdna_list = retrieve_seq("../../data/human_cdna_location.csv",\
                    "/home/zheng/reference_data/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa"\
                    ,40000)
    print("human data handled!")
    new_human_cdna_df = pd.DataFrame(new_human_cdna_list,columns=["species","transcript_id", "location", \
                                                                  "gene_name","transcript_fragment"])
    new_human_cdna_df.to_csv("../../data/human_cdna_fragment.csv")   # writer to file
    print("human_cdna_fragment.csv file saved")
    new_human_cdna_list = []   # free the memory


    # generate 10000 mouse cdna fragment and write to csv file 
    new_mouse_cdna_list = retrieve_seq("../../data/mouse_cdna_location.csv",\
                      "/home/zheng/reference_data/Mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.1.fa"\
                       ,10000)
    print("mouse data handled!")
    new_mouse_cdna_df = pd.DataFrame(new_mouse_cdna_list,columns=["species","transcript_id", "location", \
                                                                  "gene_name","transcript_fragment"])
    new_mouse_cdna_df.to_csv("../../data/mouse_cdna_fragment.csv")   # writer to file
    print("mouse_cdna_fragment.csv file saved")
    """

    positive_combination_dict = generate_positive_combination()
    print("positive_combination_dict done")

    miRNA_name_seq_dict = generate_miRNA_dic()
    print("miRNA_name_seq_dict generated!")

    # human_pseudo_miRNA_target csv file
    generate_pseudo_miRNA_target("../../data/human_cdna_fragment.csv",miRNA_name_seq_dict,\
                                 positive_combination_dict,"hsa","homo_sapiens",\
                                 "../../data/human_pseudo_miRNA_target.csv")
    print("human_pseudo_miRNA_target csv file generated!")
    
    # mouse_pseudo_miRNA_target csv file
    generate_pseudo_miRNA_target("../../data/mouse_cdna_fragment.csv",miRNA_name_seq_dict,\
                                 positive_combination_dict,"mmu","Mus_musculus",\
                                "../../data/mouse_pseudo_miRNA_target.csv")
    print("mouse_pseudo_miRNA_target csv file generated!")
