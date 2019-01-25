"""module: vectorize the miRNA_target_seq (X) for deep learning
"""
import pandas as pd
import numpy as np
import random

x_cast = {"A":[1,0,0,0],"U":[0,1,0,0],\
          "T":[0,1,0,0],"G":[0,0,1,0],\
          "C":[0,0,0,1],"N":[0,0,0,0]}


# function: get the reverse_complement sequence
def reverse_complement(seq):
    complement_dict = {"A":"T","T":"A","U":"A","G":"C","C":"G"}
    reverse_complement = ""
    seq = seq[::-1]   # reverse 
    for base in seq:
        reverse_complement += complement_dict[base]
    return reverse_complement
        

# function: randomly select bases for padding without 4 continuous nt pairing with miRNA sequence 
def selected_padding(merged_seq):
    merged_seq_len = len(merged_seq)
    SEQ_LEN = 110
    head_seq = merged_seq[0:10]
    head_reverse_complement = reverse_complement(head_seq)
    
    head_reverse_complement = head_reverse_complement.replace("U","T")
    while True:
#         print("outloop")
        flag = 0
        padding_seq = ""
        for i in range(SEQ_LEN-merged_seq_len):
#             print("innerloop")
            temp_base = random.choice("ATGC")
            padding_seq += temp_base
#             print(padding_seq)
        if len(padding_seq) < 4:
            break
        for j in range(len(padding_seq)-4):
            if head_reverse_complement.find(padding_seq[j:j+4]) >= 0:
                flag = 1
                break
        if flag ==0:
            break
    return padding_seq


# function: padding all the miRNA_target_seq to 110 nt and vectorize with one-hot encoding
def transform_xdata(df_column):
    x_dataset = []
    for line in df_column:
        line = line.strip()
        line = line.replace("X","")   # remove "X" in the sequence
        line = line.upper()   # upper case
        # padding sequence to SEQ_LEN
        padding_seq = selected_padding(line)   # padding with no-4-pairing sequence
        # padding_seq = "N"*(110-len(line))    # padding with "N"
        line = line + padding_seq
        # vectorization
        temp_list = []
        for base in line:
            temp_list.append(x_cast[base])
        x_dataset.append(temp_list)
    return x_dataset


if __name__ == "__main__":
    # read the file
    df = pd.read_csv("../../data/test_dataset.csv")
    print(df["miRNA_target_seq"].str.len().describe())
    # transformation
    X_test_list = transform_xdata(df["miRNA_target_seq"])
    #print(X_test_list)
    X_test_array = np.array(X_test_list)    
    print(X_test_array.shape)
