""" Plotting the ROC curve of our trained model on the test dataset
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc  ###计算roc和auc
#from sklearn import cross_validation
import sys
sys.path.append("../data_process")
from keras.models import load_model
from dataset_vectorization import transform_xdata

def ROC_curve():
    df = pd.read_csv("../../data/test_dataset.csv")
    X_test_list = transform_xdata(df["miRNA_target_seq"])
    X_test_array = np.array(X_test_list)
    y_test_list = np.array(df["classification"])
    CNN_model_path = "../../data/cnn_model_preTrained.h5"

    print("load the model")
    try:
        CNN_model = load_model(CNN_model_path)
    except Exception:
        print("The model file doesn't exist!")
        exit(1)
    cnn_predict_result = CNN_model.predict(X_test_array)
#    print(predict_result)

    # Compute ROC curve and ROC area for each class
    cnn_fpr,cnn_tpr,cnn_threshold = roc_curve(y_test_list,cnn_predict_result) 
     ## calculate the AUC value
    cnn_roc_auc = auc(cnn_fpr,cnn_tpr) 
    # plotting
    plt.figure(figsize=(10,10))
    plt.plot(cnn_fpr, cnn_tpr, '-',\
         linewidth=2, label='CNN model-AUC:%0.4f)' %cnn_roc_auc) 
   # plt.plot([0, 1], [0, 1], color='navy', linewidth=2, linestyle='--')
    plt.xlim([0.0, 1.0])
   # plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
   # plt.title('Receiver operating characteristic')
    plt.legend(loc = "center right")
    plt.savefig("ROC_curve.png",dpi=600)
    plt.show()
if __name__  == "__main__":
    ROC_curve()

