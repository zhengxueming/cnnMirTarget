"""
   module: construct the CNN deep learning model
"""

import keras
from keras import regularizers
from keras.models import Sequential
from keras.layers import Dense,Activation,Dropout,Flatten
from keras.layers import Conv1D,MaxPooling1D
from keras.optimizers import Adam

#function: design the model structure
def cnn_model():
    model = Sequential()
    #first layer of convolution and max-pooling
    model.add(Conv1D(16,2,activation = 'relu',padding = 'same',\
              input_shape = (110,4)))
    model.add(MaxPooling1D(pool_size = 2,padding = 'same'))
    #second layer of convolution and max-pooling
    model.add(Conv1D(32,3,activation = 'relu',padding = 'same'))
    model.add(MaxPooling1D(pool_size = 2,padding = 'same'))
    #third layer of convolution and max-pooling
    model.add(Conv1D(64,4,activation = 'relu',padding = 'same'))
    model.add(MaxPooling1D(pool_size = 2,padding = 'same'))
    #fourth layer of convolution and max-pooling
    model.add(Conv1D(128,5,activation = 'relu',padding = 'same'))
    model.add(MaxPooling1D(pool_size = 2,padding = 'same'))
   
    model.add(Flatten())
    model.add(Dropout(0.5))
    model.add(Dense(128,activation = 'relu',kernel_regularizer = regularizers.l2(0.1)))
    model.add(Dropout(0.5))
    model.add(Dense(1,activation = 'sigmoid'))
    adam = Adam()
    model.compile(loss = 'binary_crossentropy',optimizer = adam,\
              metrics = ['accuracy'])
    # print(model.summary())
    return model


if __name__ == "__main__":
    model = cnn_model()
    print(model.summary())
