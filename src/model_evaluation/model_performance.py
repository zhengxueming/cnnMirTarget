""" module: evaluate the performance of trained model on test dataset
"""

"""Evaluate the performance of the trained model using the test dataset
"""
from keras.models import load_model
from keras.models import Sequential
import numpy as np
import pandas as pd
import math
import sys
sys.path.append("../data_process")
import dataset_vectorization


# calculate TP,TN, FP and FN
def predict_comparision(y_predict,y_test,cutoff):
    tp,tn,fp,fn = 0,0,0,0
    m = len(y_predict)
    for i in range(m):
        if y_predict[i] >= cutoff:
            if y_test[i] == 1:
                tp +=1    
            else:
                fp +=1
        else:
            if y_test[i] == 0:
                tn +=1
            else:
                fn += 1 
    return tp,tn,fp,fn


# calculate SENS, SPEC, ACC and Matthews Correlation Coefficient (MCC)
def test_evaluation(model_path,X_test_dataset, y_test_dataset):
    print("load the model")
    try:
        model = load_model(model_path)
    except Exception:
        print("The model file doesn't exist!")
        exit(1)
    predict_result = model.predict(X_test_dataset)
#     calculate tp,tn,fp,fn
    tp,tn,fp,fn = predict_comparision(predict_result,y_test_dataset,0.5)
    
#    calculate sens, spec, f1,mcc and acc based on tp,tn,fp,fn   
    try:
        sensitivity = float(tp)/(float(tp) + float(fn))
        specifity = float(tn)/(float(tn) + float(fp))
        precision = float(tp)/(float(tp) + float(fp))
        accuracy = (float(tp) + float(tn))/(float(tp) + float(fp) + \
                                            float(fn) + float(tn))
        f1_score = (2 * (precision * sensitivity)) / (precision + sensitivity)
        mcc = ((float(tp) * float(tn)) - (float(fp) * float(fn))) /\
                math.sqrt((float(tp) + float(fp)) * (float(tp) + float(fn))*\
                (float(tn) + float(fp)) * (float(tn) + float(fn)))
    except ZeroDivisionError as err:
        print("Exception:",err)
        exit(1)
    print("Sensitivity/recall on the test data is :{}".format(sensitivity)) 
    print("specifity on the test data is :{}".format(specifity)) 
    print("precision on the test data is :{}".format(precision))
    print("accuracy on the test data is :{}".format(accuracy))
    print("f1_score on the test data is :{}".format(f1_score))
    print("mcc on the test data is :{}".format(mcc))

    return sensitivity,specifity,f1_score,mcc,accuracy


if __name__ == "__main__":
    # read the training dataset
    test_df = pd.read_csv("../../data/test_dataset.csv")
    X_test_vector_list = dataset_vectorization.transform_xdata(test_df["miRNA_target_seq"])
    # tansform to numpy array
    X_test_vector_array = np.array(X_test_vector_list)
    y_test_vector_array = np.array(test_df["classification"])

    model_path = "../../data/cnn_model_preTrained.h5"
    sensitivity,specifity,f1_score,mcc,accuracy =\
        test_evaluation(model_path,X_test_vector_array,y_test_vector_array)


