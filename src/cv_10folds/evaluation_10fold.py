"""Train and evaluate the performance for 10-fold cross evaluation 
"""

#from keras.models import load_model
#from keras.models import Sequential
import numpy as np
import pandas as pd
#import math
import sys
sys.path.append("../data_process")
sys.path.append("../model_construction")
sys.path.append("../model_train")
sys.path.append("../model_evaluation")
import dataset_partition
import dataset_vectorization 
import cnn_model
import model_performance 

def cnn_train_10folds(X_data,y_data):
    model = cnn_model.cnn_model()
    model.fit(X_data,y_data,batch_size = 200, epochs = 20)
    return model


training_df,test_df = dataset_partition.train_test_partition(1,0)    # all in train_df
print(len(training_df))
each_num = 15938    # number of each fold in 10-fold cross valuation
max_num = 159380
for i in range(10):
    #transformation
    X_list = dataset_vectorization.transform_xdata(training_df["miRNA_target_seq"])
    X_array = np.array(X_list)
    print(X_array.shape)
    y_array = np.array(training_df["classification"])
    print(y_array.shape)
    print("fold {}".format(i+1))
    
    X_train_array =  np.concatenate((X_array[0:i*each_num],X_array[(i+1)*each_num:]),axis = 0)
    y_train_array =  np.concatenate((y_array[0:i*each_num],y_array[(i+1)*each_num:]),axis = 0)
    
    X_test_array = X_array[i*each_num:min(max_num,(i+1)*each_num)]
    y_test_array = y_array[i*each_num:min(max_num,(i+1)*each_num)]


    model = cnn_train_10folds(X_train_array,y_train_array)

    print("model train over")
    model_path = "cnn_model_preTrained" +str(i+1) + ".h5"
    model.save(model_path)
    print("The model is saved as","in the current directory:",model_path)
    sensitivity,specifity,f1_score,mcc,accuracy =\
           model_performance.test_evaluation(model_path,X_test_array,y_test_array)

    with open("performance_results","a+") as fd:
         fd.write("fold {}".format(i+1))
         fd.write("\n")
         fd.write("Sensitivity/recall on the test data is :{}".format(sensitivity))
         fd.write("\n")
         fd.write("specifity on the test data is :{}".format(specifity))
         fd.write("\n")
         fd.write("accuracy on the test data is :{}".format(accuracy))
         fd.write("\n")
         fd.write("f1_score on the test data is :{}".format(f1_score))
         fd.write("\n")
         fd.write("mcc on the test data is :{}".format(mcc))
         fd.write("\n")
         fd.write("\n")
print ("finished!")
