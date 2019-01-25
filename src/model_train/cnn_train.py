"""module: train the model with training dataset
"""
import pandas as pd
import numpy as np
import sys
sys.path.append("../model_construction")
sys.path.append("../data_process")
import os
import cnn_model
import dataset_vectorization


# function: train the model with training dataset
def cnn_train(X_data,y_data):
    model = cnn_model.cnn_model()
    if os.path.exists("../../data/cnn_model_preTrained.h5"):
        print("load the weights")
        model.load_weights("../../data/cnn_model_preTrained.h5")
    model.fit(X_data,y_data,batch_size = 200, epochs = 2,\
          validation_split = 0.032)
    return model


# read the training dataset
training_df = pd.read_csv("../../data/training_dataset.csv")
X_data_vector_list = dataset_vectorization.transform_xdata(training_df["miRNA_target_seq"])
# tansform to numpy array
X_training_data_array = np.array(X_data_vector_list)
y_training_data_array = np.array(training_df["classification"])
print(X_training_data_array.shape)
print(y_training_data_array.shape)



model = cnn_train(X_training_data_array,y_training_data_array)
print("model train over")
model_path = "../../data/cnn_model_preTrained.h5"
model.save(model_path)
print("The model is saved as","in the current directory:",model_path)
