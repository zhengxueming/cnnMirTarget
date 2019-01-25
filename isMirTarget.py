""" miRNA target prediction giving a miRNA sequence and target sequence
"""
#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

import sys, getopt
import os
import pandas as pd
import numpy as np
import sys
sys.path.append("./src/data_process")
from dataset_vectorization import transform_xdata

# function: handle the input sequence
def seq_process(seq):
    # remove line break   
    seq = seq.strip('\n')
    # remove sapce
    seq =seq.replace(' ', '')
    # onvert the string to all uppercase
    seq = seq.upper()
    # print(seq)
    # check correctness of the RNA sequence
    for char in seq:
        # print (char)
        if char not in ["A","T","U","G","C"]:
            print("Please input the right RNA sequence")
            exit(1)
    return seq

def usage():
    print("""
    USAGE: 
    python isMirTarget.py --miRNA sequence --mRNA mRNA_file_path 
    for example: python3 isMirTarget.py -s  TCAGTGCATCACAGAACTTTGT -t TCCGTCTTCGTCTTCCGTCTTCTTCTTCCG 
             or  python3 isMirTarget.py --miRNA  UCAGUGCAUCACAGAACUUUGU  --mRNA examples/test_mRNA.fasta

     """)


# get the parameters
def parse_opt(argv):
    miRNA_seq = ""
    target_seq = ""
    mRNA_file = ""
    try:
        opts, args = getopt.getopt(argv,"hs:t:",["miRNA=","target=","mRNA="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    if len(opts) < 1:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(0)
        elif opt in ("-s","--miRNA"):
            miRNA_seq = arg
            miRNA_seq = seq_process(miRNA_seq)
        elif opt in ("-t","--target"):
            target_seq = arg
            target_seq = seq_process(target_seq)
        
        elif opt == "--mRNA":
            mRNA_file = arg
            

    return miRNA_seq, target_seq, mRNA_file

# function: retrieve sequence fragments for prediction
#def get_merged_fragment(miRNA_seq,target_seq,step):
#   merged_fragment_list = []
#   miRNA_seq_len = len(miRNA_seq)
#    target_seq_len = len(target_seq)
#   window_size = 110 - miRNA_seq_len
#   for i in range(0,target_seq_len-window_size,step):
#       end =  min(i+window_size,target_seq_len)
#        target_fragment = target_seq[i:end]
#        merged_fragment_list.append(miRNA_seq + target_fragment)
#    return merged_fragment_list
    
    
# read the mRNA file (fasta sequence) and return mRNA sequence list ([[name1,seq1],[name2,seq2],...])
def read_mRNA_file(file_path):
    mRNA_name_seq_list = []
    name = ""
    seq = ""
    with open(file_path) as fd:
        for line in fd:
            line = line.strip()
            line = line.replace("\n","")
            if line.startswith(">"):
                mRNA_name_seq_list.append([name,seq])
                name = line
                seq = ""
            else:
                seq += line
        mRNA_name_seq_list.append([name,seq])   # add the last one
    if len(mRNA_name_seq_list) == 1:
        print ("read mRNA sequence file error!")
        sys.exit(1)
    #print("break3")
    #print (mRNA_name_seq_list)
    return mRNA_name_seq_list[1:]


# calculate MFE using RNAhybrid and get the candidate sites. return start positions of sites
def get_candidate_pos(m_seq,mi_seq):
    # write mi_seq to temp file
    with open("temp_miRNA.fasta","w") as fd:
        fd.write(">miRNA\n")
        fd.write(mi_seq)

    # using RNAhybrid to get the candidate sites
    
#    result = os.popen\
 #            ("./bin/RNAhybrid -d -t %s -q temp_miRNA.fasta -e -10 -m 80 -c|awk -F ':' '{ print $0,'\n'}'"\
  #            %(m_seq)).readlines()
    result = os.popen\
             ('./bin/RNAhybrid -d -t %s -q temp_miRNA.fasta -e -20 -m 25 -c'\
              %(m_seq)).readlines()
    #print(m_seq)
    #print("break4")
    #print (result)
    os.system("rm temp_miRNA.fasta")    # delete the temp file
    pos_list = []
    for line in result:
        line = line.strip()
        if line != "":
     #       print(line)
      #      print("===============")
            line = line.split(":")    
            each_pos = int(line[6])   # position
       #     print("break2")
        #    print(each_pos) 
            pos_list.append(each_pos)
    return pos_list


# print the prediction output
def print_output(prediction,threshold):
    # print the results based on the prediction
    probability = 0
    for i in range(len(prediction)):
        # get the highest probability
        if prediction[i] > probability:
            probability = prediction[i]
    
    if probability >= threshold:
        print("True,  ",end = "")
        print("probability: {}".format(probability))
    else:
        print("False,no predicted target sequence!")


# get the target sequence, return a list of target
def get_target_list(mRNA_seq,pos_list):
    target_seq_list = []
    # check and process the sequence
    mRNA_seq = seq_process(mRNA_seq)
    
    mRNA_len = len(mRNA_seq)
    for start in pos_list:
#        print(type(start))
        target_seq = mRNA_seq[max(0,start-1):min(mRNA_len,start+90)]
        target_seq_list.append(target_seq)
#        target_seq = mRNA_seq[max(0,start-6):min(mRNA_len,start+76)]   # move 5 nt backward 
#       target_seq_list.append(target_seq)
#       target_seq = mRNA_seq[max(0,start-11):min(mRNA_len,start+71)]  # move 10 nt backward
#       target_seq_list.append(target_seq)
#       target_seq = mRNA_seq[max(0,start-16):min(mRNA_len,start+65)]  # move 15 nt backward
#       target_seq_list.append(target_seq)
#       target_seq = mRNA_seq[max(0,start-21):min(mRNA_len,start+59)]  # move 20 nt backward
#       target_seq_list.append(target_seq)
        
       
    return target_seq_list


def main(argv):
    #parse_opt(argv)
    miRNA_seq, target_seq, mRNA_file = parse_opt(argv)
    print("")
    # analysis of the args
    if target_seq and mRNA_file:      # both provided 
        #print("break1")
        usage()
        sys.exit(2)

    if len(miRNA_seq) < 16:
        print("Please input the right miRNA sequence")
        sys.exit(1)

    if mRNA_file:    # only mRNA_file is provided
        # check file exist
        if not os.path.exists(mRNA_file):
            print("file doesn't exist!")
            sys.exit(1)
        # open the file and read the sequence
        mRNA_name_seq_list = read_mRNA_file(mRNA_file)    # get the list of [name,seq]
        num_mRNA = len(mRNA_name_seq_list)
        for i in range(num_mRNA):
            mRNA_name = mRNA_name_seq_list[i][0]
            mRNA_seq = mRNA_name_seq_list[i][1]
            mRNA_len = len(mRNA_seq)
            candidate_pos_list = get_candidate_pos(mRNA_seq,miRNA_seq)
            target_seq_list = get_target_list(mRNA_seq,candidate_pos_list)
            if len(target_seq_list) == 0:
                print (mRNA_name)
                print("False,no site accessibility!")
            # prediction using the trained model
            else:
                prediction = predict_results(target_seq_list,miRNA_seq)
                   
                # output the result
                print (mRNA_name)    # mRNA name
                print_output(prediction, THRESHOLD)       # prediction result
                print("")
    else: 
        if len(target_seq) < 4:
            print ("length of target sequence is at least 4 nt")
            sys.exit(1)

        else:       
            if len(target_seq) <= 110 -len(miRNA_seq):
                target_seq_list = [target_seq,target_seq]    # make a list
            else:
                candidate_pos_list = get_candidate_pos(target_seq,miRNA_seq)
                target_seq_list = get_target_list(target_seq,candidate_pos_list)
            
            # prediction using the trained model
            prediction = predict_results(target_seq_list,miRNA_seq)
                   
            # output the result
            print_output(prediction,THRESHOLD)       # prediction result
            
        print("")


# function: calculate the result based on the model
def predict_results(target_seq_list,miRNA_seq):
    num_seq_list = len(target_seq_list)
    merged_seq_list = []
    for i in range(num_seq_list):
        merged_seq = miRNA_seq + target_seq_list[i]
        if len(merged_seq) < 110:
            merged_seq += "N"*(110-len(merged_seq))
        if len(merged_seq) >= 110:
            merged_seq = merged_seq[0:110]
        merged_seq_list.append(merged_seq)
   # print("break 5")
   # print(merged_seq_list)
    merged_seq_df = pd.DataFrame(merged_seq_list,columns = ["miRNA_target_seq"])
    # print(merged_seq_df)
    # vectorization
    merged_seq_vector = transform_xdata(merged_seq_df["miRNA_target_seq"]) 
    
    merged_seq_vector_array = np.array(merged_seq_vector)
    # reload the model
    from keras.models import load_model
    try:
        model = load_model('data/cnn_model_preTrained.h5')
    except IOEEor:
        print("model file doesn't exit")
        
    # prediction
    result =  model.predict(merged_seq_vector_array)

    return result

if __name__ == "__main__":
    THRESHOLD = 0.5
    main(sys.argv[1:])
    exit(0)

"""
while True:
    print("miNRA sequence:(upper or lower case); @ to quit")
 print("....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8")
    miRNA_seq = input()
    if miRNA_seq.find("@") >= 0:
        break
    print("target sequence(upper or lower case); @ to quit")
    print("....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8")
    target_seq = input()
    if target_seq.find("@") >= 0:
        break
    else:
        print("in process...")
    print("")
print ("end")

"""
