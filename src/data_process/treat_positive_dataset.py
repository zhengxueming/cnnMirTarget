""" module: intergrate the positive datasets.
    read all the miRNA-target datasets of different format,merge into one and store in a csv file
"""

import pandas as pd
from xlrd import open_workbook

# function:read the txt file and convert into dataframe
def read_text_to_pd(input_file_path):
    fd = open(input_file_path,"r")
    # store a list of data(miNRA_species,miRNA_name,miRNA_seq,target_gene, target_seq, miRNA_target _seq)
    data_to_scv = []  
    i =1
    for data in fd: 
        i+=1
        if i <33:
            continue    # remove the head
        else:
            new_data = []
            data = data.split()
        
            new_data.append("Homo sapiens")   # species
            new_data.append(data[1])   # miRNA name
        
            miRNA_seq = data[4]
            new_data.append(miRNA_seq)           # miRNA sequence 
        
            new_data.append(data[5])  # target gene name
        
            target = data[8]
            new_data.append(target)  # target sequence
        
            miRNA_target_seq = miRNA_seq + target
            new_data.append(miRNA_target_seq)    # miRNA_target_seq
            new_data.append(1)               # indicate positive data
        
            data_to_scv.append(new_data)   
        
    fd.close()

    df = pd.DataFrame(data_to_scv, columns=["miNRA_species","miRNA_name","miRNA_seq",\
                                        "target_gene", "target_seq", "miRNA_target_seq","classification"])
    return df


# function:read the xls file and convert into dataframe
def read_xls_to_pd(file_path):
    workbook = open_workbook(file_path)  # open xls file
    sheet_name= workbook.sheet_names()  # print the names of sheet
    # print(sheet_name)

    worm_data = workbook.sheet_by_name('Sheet2')  # C. elegans miRNA_target data
    human_data = workbook.sheet_by_name('Sheet3')  # human miRNA_target data
    viras_data = workbook.sheet_by_name('Sheet4')  # viral miRNA_target data
    mouse_data = workbook.sheet_by_name('Sheet5')  # mouse miRNA_target data

    species_sheet = (worm_data,human_data,viras_data,mouse_data)
    species_name = ("Caenorhabditis elegans","Homo sapiens","virus","Mus musculus")   #species name of each sheet

    # store a list of data(miNRA_species,miRNA_name,miRNA_seq,target_gene, target_seq, miRNA_target_seq)
    data_to_scv = []  
    i = 0
    for sheet in species_sheet:    
        for row in range(sheet.nrows-4):   # for skipping the first 4 rows
            new_data = []
        
            new_data.append(species_name[i])   # species
            new_data.append(sheet.cell(row+4,6).value)   # miRNA name
        
            miRNA_seq = sheet.cell(row+4,7).value  
            new_data.append(miRNA_seq)           # miRNA sequence 
        
            target_location = sheet.cell(row+4,0).value + ":" + str(sheet.cell(row+4,1).value)\
                               + "-" + str(sheet.cell(row+4,2).value)
            new_data.append(target_location)  # target gene name,target gene location
            target = sheet.cell(row+4,8).value
            new_data.append(target)  # target sequence
        
            miRNA_target_seq = miRNA_seq + target
            new_data.append(miRNA_target_seq)    # miRNA_target_seq
            new_data.append(1)               # indicate positive data
        
            data_to_scv.append(new_data)
        i += 1     
        
    df = pd.DataFrame(data_to_scv, columns=["miNRA_species","miRNA_name","miRNA_seq",\
                                        "target_gene", "target_seq", "miRNA_target_seq","classification"])
    return df


# function:generate dictionary for miRNA name and sequence
def txt_to_dictionary(file_path):
    # read the data
    name_seq_list = []
    fd = open(file_path,'r')
    for line in fd:
        line = line.split()
        line[0] = line[0].replace(">","")
        name_seq_list.append(line[0])
    # print(name_seq_list)
    fd.close()

    # generate miRNA_name:miRNA_seq dictionary
    name_seq_dict = {}

    for i in range(len(name_seq_list)-1):
        name_seq_dict[name_seq_list[i]] = name_seq_list[i+1]
    # print(len(name_seq_dict))
    # print(name_seq_dict["hsa-miR-122-5p"]) 
    return name_seq_dict

# function:read the xlsx file and convert into dataframe
def read_xlsx_to_pd(file_path,name_seq_dict):
    workbook = open_workbook(file_path)  # 打开xls文件
    sheet_name= workbook.sheet_names()  # 打印所有sheet名称，是个列表
    # print(sheet_name)
    sheet = workbook.sheet_by_index(0)  # 根据sheet索引读取sheet中的所有内容

    # print(sheet.name, sheet.nrows, sheet.ncols)  # sheet的名称、行数、列数
    # print("")

    # store a list of data(miNRA_species,miRNA_name,miRNA_seq,target_gene, target_seq, miRNA_target_seq)
    data_to_scv = []  
    num_no_reference = 0
    for row in range(sheet.nrows-1):
        if sheet.cell(row+1,1).value not in name_seq_dict.keys():
            num_no_reference += 1
        else:
            new_data = []
            new_data.append(sheet.cell(row+1,2).value)   # species
            new_data.append(sheet.cell(row+1,1).value)   # name
        
            miRNA_seq = name_seq_dict[sheet.cell(row+1,1).value]  
            miRNA_seq = miRNA_seq.replace("U","T")  # replace U with T
            new_data.append(miRNA_seq)           # miRNA sequence 
        
            new_data.append(sheet.cell(row+1,3).value)  # target gene name
            target = sheet.cell(row+1,6).value
            new_data.append(target)  # target sequence
        
            miRNA_target_seq = miRNA_seq + target
            new_data.append(miRNA_target_seq)    # miRNA_target_seq
            new_data.append(1)               # indicate positive data
        
            data_to_scv.append(new_data)
        
    # print(data_to_scv)
    # print(num_no_reference)
    df = pd.DataFrame(data_to_scv, columns=["miNRA_species","miRNA_name","miRNA_seq",\
                                        "target_gene", "target_seq", "miRNA_target_seq","classification"])
    return df


# function:merge the dataset
def merge_csv(file_path1,file_path2,file_path3):
    csv_list = [file_path1,file_path2,file_path3]

    frames = []
    for file_path in csv_list:
        df_each = pd.read_csv(file_path)
        frames.append(df_each)

    df = pd.concat(frames,ignore_index=True)
    # print (len(df))
    # df['target_seq'].str.len().describe()

    # statistics on miRNA_target_seq
    # print(df['miRNA_target_seq'].str.len().describe())
    # df['miRNA_target_seq'].describe()

    # remove miRNA_target_seq >110 sequences
    # df = df[df['miRNA_target_seq'].str.len() <=110]
    # print(len(df))
    return df

if __name__ == "__main__":
    # read txt file
    txt_file_path = "../../data/mmc1.txt"
    df_from_txt = read_text_to_pd(txt_file_path)
    df_from_txt.to_csv("../../data/new_data_from_Aleksandra.csv")   # writer to csv file

    # read the xls files
    xls_file_path = "../../data/NIHMS600705-supplement-03.xls"
    df_from_xls = read_xls_to_pd(xls_file_path)
    df_from_xls.to_csv("../../data/new_data_from_NIHMS600705.csv")   # writer to file

    # load the dictionary
    miRNA_seq_name_dict = txt_to_dictionary("../../data/mature_miRNA.fa")
    # read the xls files
    xlsx_file_path = "../../data/MicroRNA_Target_Sites.xlsx"
    df_from_xlsx = read_xlsx_to_pd(xlsx_file_path,miRNA_seq_name_dict)
    df_from_xlsx.to_csv("../../data/new_data_from_miRTarBase.csv")   # writer to file
    
    # merge the dataset
    merged_df = merge_csv("../../data/new_data_from_Aleksandra.csv",\
                          "../../data/new_data_from_NIHMS600705.csv",\
                          "../../data/new_data_from_miRTarBase.csv")
    merged_df = merged_df[merged_df["target_seq"].str.len() >= 4]   # remove target_seq with len <4 

    merged_df = merged_df.drop_duplicates(["miRNA_target_seq"])     # remove duplicates

    merged_df.to_csv("../../data/all_positive_data.csv")

    print("Write file completed！")
