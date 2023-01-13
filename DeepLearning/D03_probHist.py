import os
import numpy as np
from numpy import save, load
import keras
import time
#can do tensorflow.keras if they have tensorflow
from keras import optimizers
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, Reshape, GlobalAveragePooling1D, Activation, GlobalAveragePooling2D
from keras.layers import Conv2D, MaxPooling2D, Conv1D, MaxPooling1D
from keras.utils import np_utils
from matplotlib.lines import Line2D
import random

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def plotTrace(traces, title, saveLoc, show=False):
    f, ax = plt.subplots(4,1)
    for chID, trace in enumerate(traces):
        ax[chID].plot(trace)
    plt.suptitle(title)
    if show:
        plt.show()
    else:
        plt.savefig(saveLoc, format='png')
    return


#######input path and file names#######
round = '3rdpass'
path = f'DeepLearning/data/{round}/'
RCR = np.load(os.path.join(path, "ReflCR_67950events_part0.npy"))[10000:,0:4] #input a subset of the data here so that you can validate on the other set
#Noise = np.load(os.path.join(path, "Station13_Data_500000events_part0.npy"))[67950:,0:4] #make sure the signal and noise subset of data are the same size
Noise = np.load(os.path.join(path, "Station13_Data_500000events_part0.npy"))[10000:,0:4] #make sure the signal and noise subset of data are the same size
#RCR = np.load(os.path.join(path, "ReflCR_499events.npy"))[100:,0:4] #validate on the other set of data
#Nu = np.load(os.path.join(path, "Nu_283events.npy"))[100:,0:4]

RCR = np.reshape(RCR, (RCR.shape[0], RCR.shape[1],RCR.shape[2],1))
Noise = np.reshape(Noise, (Noise.shape[0], Noise.shape[1],Noise.shape[2],1))


#######input path to trained h5 model#######
model = keras.models.load_model(f'DeepLearning/h5_models/{round}_trained_CNN_1l-10-8-10_do0.5_fltn_sigm_valloss_p4_measNoise0-20k_0-5ksigNU-Scaled_shuff_monitortraining_0.h5')

prob_RCR = model.predict(RCR)
prob_Noise = model.predict(Noise)

plotRandom = True
if plotRandom == True:
    for iE, Noi in enumerate(Noise):
        if iE % 10000 == 0:
            output = prob_Noise[iE][0]
            print(f'output {output}')
            plotTrace(Noi, f"Noise {iE}, Output {output:.2f}",f"DeepLearning/data/3rdpass/Noise_{iE}_Output_{output:.2f}.png")
    for iE, rcr in enumerate(RCR):
        if iE % 1000 == 0:
            output = prob_RCR[iE][0]
            plotTrace(rcr, f"RCR {iE}, Output {output:.2f}",f"DeepLearning/data/3rdpass/RCR_{iE}_Output_{output:.2f}.png")

quit()

"""
mask = np.zeros_like(prob_Nu).astype(bool)
max_amp_Nu = np.zeros_like(prob_Nu)
for iE, modeled_Nu in enumerate(prob_Nu):
    for chId, trace in enumerate(Nu[iE]):
        mask[iE] = mask[iE] or np.any(Nu[iE])
        max_amp_Nu[iE] = np.max([max_amp_Nu[iE], np.max(trace)])
    if not mask[iE]:
        continue
    continue
#    if modeled_Nu < 0.8:
    if modeled_Nu > 0.97:
        lowNu = Nu[iE]
        print(f'shape lowNu {np.shape(lowNu)}')

        f, ax = plt.subplots(4,1)
        print(f'modeled nu')
#        print(f'lowNu is {lowNu}')
        for chID, trace in enumerate(lowNu):
            ax[chID].plot(trace)
        print(modeled_Nu)
        plt.suptitle(f'Nu output {modeled_Nu}')
        plt.show()

max_amp_RCR = np.zeros_like(prob_RCR)
for iE, modeled_RCR in enumerate(prob_RCR):
    for chId, trace in enumerate(RCR[iE]):
        max_amp_RCR[iE] = np.max([max_amp_RCR[iE], np.max(trace)])

    f, ax = plt.subplots(4,1)
    print(f'modeled nu')
#        print(f'lowNu is {lowNu}')
    for chID, trace in enumerate(RCR[iE]):
        ax[chID].plot(trace)
#    print(modeled_Nu)
    plt.suptitle(f'RCR output {modeled_RCR}')
    plt.show()


mask = mask.astype(bool)

plt.scatter(max_amp_Nu[mask], prob_Nu[mask], label='Nu')
plt.scatter(max_amp_Nu[~mask], prob_Nu[~mask], label='Zero Trace')
plt.scatter(max_amp_RCR, prob_RCR, label='RCR')
plt.legend()
plt.title(f'Nu')
#plt.xscale('log')
plt.ylabel('Nu Output')
plt.xlabel('all traces max amp')
plt.show()

print(f'shape mask {np.shape(mask)} and nu {np.shape(prob_Nu)}')
"""

max_amps = np.zeros(len(prob_Noise))
for iC, trace in enumerate(Noise):
    max_amps[iC] = np.max(trace)
max_amps_RCR = np.zeros(len(prob_RCR))
for iC, trace in enumerate(RCR):
    max_amps_RCR[iC] = np.max(trace)

prob_Noise = 1 - prob_Noise
prob_RCR = 1 - prob_RCR
plt.scatter(max_amps, prob_Noise, label='Trig Events')
plt.scatter(max_amps_RCR, prob_RCR, label='RCR')
plt.ylabel('Network Output - 1=RCR')
plt.xlabel('Max amp of channels')
plt.legend()
plt.title('Station 13 Training')
plt.savefig(path+f'MaxAmps_Output.png', format='png')
plt.clf()
#plt.show()



fig = plt.figure()
ax = fig.add_subplot(111)

dense_val = False
ax.hist(prob_Noise, bins=20, range=(0, 1), histtype='step', color='red', linestyle='solid', label='Triggered Events', density=dense_val)
ax.hist(prob_RCR, bins=20, range=(0, 1), histtype='step',color='blue', linestyle='solid',label='RCR',density = dense_val) 


plt.xlabel('network output', fontsize=18)
plt.ylabel('events', fontsize=18)
plt.xticks(size=18)
plt.yticks(size=18)
plt.yscale('log')
plt.title('Station 13')
handles, labels = ax.get_legend_handles_labels()
new_handles = [Line2D([], [], c=h.get_edgecolor(), linestyle=h.get_linestyle()) for h in handles]
plt.legend(loc='upper center', handles=new_handles, labels=labels, fontsize=18)
#plt.show()
plt.savefig(path+f'AnalysisOutput.png', format='png')
