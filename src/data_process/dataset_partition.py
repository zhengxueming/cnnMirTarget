"""module: read the positive and negative dataset and merge into the whole training dataset
   partition the dataset into training and test dataset
"""

import pandas as pd
import numpy as np
from sklearn.utils import shuffle
import random

# function: merge all the dataset,shuffle the order and return the taining dataset and testdataset
def train_test_partition(state_num,test_num):
    # generate test and train dataset
    df_positive = pd.read_csv("../../data/all_positive_data.csv")
    df_positive = df_positive.iloc[:,2:]    # remove the first two column
    print(len(df_positive))
    df_positive = df_positive[df_positive["miRNA_target_seq"].str.len() <= 110]    # remove some items
    print(len(df_positive))
#    df_positive.to_csv("positive.csv")

    df_negative1 = pd.read_csv("../../data/human_pseudo_miRNA_target.csv")
    df_negative1 = df_negative1.iloc[:,1:]
    # replace the column names to merge the datasets
    df_negative1.rename(columns={'species':'miNRA_species', \
                                 'pseudo_taget_gene':'target_gene', \
                                 'pseudo_target_seq':'target_seq',\
                                 'pseudo_miRNA_target_seq':'miRNA_target_seq'}, inplace = True)
    print(len(df_negative1))
    df_negative1 = df_negative1[df_negative1["miRNA_target_seq"].str.len() <= 110]    # remove some items
    print(len(df_negative1))
 #   df_negative1.to_csv("negative_human.csv")

    df_negative2 = pd.read_csv("../../data/mouse_pseudo_miRNA_target.csv")
    df_negative2 = df_negative2.iloc[:,1:]
    # replace the column names to merge the datasets
    df_negative2.rename(columns={'species':'miNRA_species', \
                                'pseudo_taget_gene':'target_gene', \
                                'pseudo_target_seq':'target_seq',\
                                'pseudo_miRNA_target_seq':'miRNA_target_seq'}, inplace = True)

    print(len(df_negative2))
    df_negative2 = df_negative2[df_negative2["miRNA_target_seq"].str.len() <= 110]    # remove some items
    print(len(df_negative2))
  #  df_negative2.to_csv("negative_mouse.csv")
    # merger all the datasets
    df_merged = pd.concat([df_positive,df_negative1,df_negative2])
    print(df_merged.head())
    print (len(df_merged))
    
    df_merged = df_merged[df_merged["miRNA_target_seq"].str.len() <= 110]    # remove some items
    print (len(df_merged))
    # print (df_merged.columns.values)
    # print (df_merged.iloc[0:3,:])

    df_merged = shuffle(df_merged,random_state = state_num)
    
    test_df = df_merged.iloc[0:test_num,:]
    train_df = df_merged.iloc[test_num:,:]

    return train_df, test_df



if __name__ == "__main__":
    training_df,test_df = train_test_partition(1,5000)
    # write to csv files
    # training_df.to_csv("../../data/training_dataset.csv")
    # test_df.to_csv("../../data/test_dataset.csv")
    print("done!")
