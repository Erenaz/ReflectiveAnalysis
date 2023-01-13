#Current training sets are ~1k events
#May need much larger training sets depending upon difficulty of problem
#train on 75%, validate on 25%, see how good the result is


import os
import glob
from matplotlib import pyplot as plt
import numpy as np
from tensorflow import keras
from keras import optimizers
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, Reshape, GlobalAveragePooling1D, Activation, GlobalAveragePooling2D
from keras.layers import Conv2D, MaxPooling2D, Conv1D, MaxPooling1D
from keras.utils import np_utils
import random

round = '3rdpass'
path = f'DeepLearning/data/{round}/'
#Take first 750 events, only first 4 channels are being used
#RCR = np.load(os.path.join(path, "ReflCR_7436events_part0.npy"))[0:7000,0:4] #input a subset of the data here so that you can validate on the other set
RCR = np.load(os.path.join(path, "ReflCR_67950events_part0.npy"))[0:10000,0:4] #input a subset of the data here so that you can validate on the other set
TrainCut = len(RCR)
#Nu = np.load(os.path.join(path, "Nu_283events.npy"))[0:100,0:4] #make sure the signal and noise subset of data are the same size
Noise = np.load(os.path.join(path, "Station13_Data_500000events_part0.npy"))[0:10000,0:4] #make sure the signal and noise subset of data are the same size
print(Noise.shape)
index = np.arange(0, len(Noise), 1)
new_index = random.shuffle(index)
Noise = Noise[new_index]
print(Noise.shape)
Noise = Noise[0:TrainCut,0:4]
print(Noise.shape)


#make signal the same shape as the noise data, if needed
#Reuse one set multiple times to match larger dataset of the other
# signal = np.vstack((signal,signal,signal,signal))
# signal = signal[0:noise.shape[0]]

print(RCR.shape)
#print(Nu.shape)
print(Noise.shape)

x = np.vstack((RCR, Noise))  # shape is (200000, 1, 240)

n_samples = x.shape[2]
n_channels = x.shape[1]
x = np.expand_dims(x, axis=-1)
#Zeros are noise, 1 signal
#y is output array
y = np.vstack((np.zeros((RCR.shape[0], 1)), np.ones((Noise.shape[0], 1))))
s = np.arange(x.shape[0])
np.random.shuffle(s)
x = x[s]
y = y[s]
print(x.shape)

BATCH_SIZE = 32
#Iterate over many epochs to see which has lowest loss
#Then change epochs to be at lowest for final result
EPOCHS = 100

#This automatically saves when loss increases over a number of patience cycles
callbacks_list = [keras.callbacks.EarlyStopping(monitor='val_loss', patience=2)]
def training(j):
  model = Sequential()
#Convolutions have worked better for our data so far than connected NN
#First num = number of kernel/filters
#Second pair of numbers is kernel/filter size
#4 channels, 10 samples wide
#can widen the samples to get different set of features, could go up to 50 samples
  model.add(Conv2D(10, (4, 10), activation='relu', input_shape=(n_channels, n_samples, 1)))
  model.add(Dropout(0.5))
  # model.add(MaxPooling2D(pool_size=(1, 10)))
  model.add(Flatten())
#If doing 3 options, change activation to 'softmax'
  model.add(Dense(1, activation='sigmoid'))
  model.compile(optimizer='Adam',
                loss='binary_crossentropy',
                metrics=['accuracy'])

  model.fit(x, y, validation_split=0.2, epochs=EPOCHS, batch_size=BATCH_SIZE, verbose=1,callbacks=callbacks_list)
  model.summary()

  #input the path and file you'd like to save the model as (in h5 format)
  model.save(f'DeepLearning/h5_models/{round}_trained_CNN_1l-10-8-10_do0.5_fltn_sigm_valloss_p4_measNoise0-20k_0-5ksigNU-Scaled_shuff_monitortraining_{j}.h5')
  
#can increase the loop for more trainings is you want to see variation
for j in range(1):
  training(j)
