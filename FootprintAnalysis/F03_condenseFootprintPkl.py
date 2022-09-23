import sys
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt


def saveFile(condensed_output, filename, suffix):
    save_file_as = filename + suffix + '.pkl'
    with open(save_file_as, 'wb') as fout:
        print(f'saved file as {save_file_as}')
        pickle.dump(condensed_output, fout, protocol=pickle.HIGHEST_PROTOCOL)
    fout.close()


low = 0
#high = 49
high = 99		#SP Footprints
#low = 0
#high = 499		#MB Footprints


condensed_output = {}

CoREAS_mode = 'direct'
CoREAS_mode = 'refracted'
type = 'MB'
config = 'MB_old'
noise = True
amp = True

depthLayer = 300.0
dB = 0.0

cores = 1600
spacing = 7

amp_type = '100'
#amp_type = None


if type == 'SP':
    low = 0
    high = 49
#    low = 1000
#    high = 1049
    limit = 2199
elif type == 'MB':
    depthLayer = 576.0
    dB = 0.0

    low = 0
    high = 4
#    high = 24
#    high = 49
#    high = 99
    limit = 1000
#    limit = 4000
elif type == 'IceTop':
    low = 16.0
    high = 16.1
    limit = 18.7    
    numIceTop = 20
    iceTopSin = 0
#-1 makes a file with all sin bins
#    iceTopSin = -1


iter = (high - low) + 1
if type == 'IceTop':
    iter = high - low



if False:
    lpda_coinc = 2
    channels = 4
    config += f'_{lpda_coinc}of{channels}'

if noise:
    config += '_wNoise'
if amp:
    config += f'_wAmp'
config += f'{amp_type}s'

if type == 'IceTop':
    config += '_InfAir'

config += '_FinalStat'

save_file_as = f'FootprintAnalysis/data/CoREAS_{CoREAS_mode}_{config}_Layer{depthLayer}m_{dB}dB_Area{spacing:.2f}_{cores}cores'
suffix = '_part'
split = True
split_iD = 0

new_runid = 0
#while high < 4000:	#MB
#while high < 2099:	#SP
while high < limit:
    if type == 'IceTop':
        low = round(low, 1)
        high = round(high, 1)

###Do I need this anymore?
#   if low == 1000:
#        low += iter
#        high += iter
#        continue


    to_open = f'output/CRFootprintRates/CoREAS_{CoREAS_mode}_{config}'
    if type == 'IceTop':
#        to_open += f'_REALInfAir'
        to_open += f'_{iceTopSin:.1f}Sin_{numIceTop}PerBin'
    to_open += f'_Layer{depthLayer}m_{dB}dB_Area{spacing:.2f}_{cores}cores_id{low}_{high}.pkl'

#    with open(f'output/CRFootprintRates/CoREAS_{CoREAS_mode}_{config}_Layer{depthLayer}m_{dB}dB_Area{spacing:.2f}_{cores}cores_id{low}_{high}.pkl', 'rb') as fin:
    if not os.path.exists(to_open):
        print(f'File does not exist, {to_open}, skipping')
        if type != 'IceTop':
            low += iter
            high += iter
        else:
            if iceTopSin == -1:
                low += iter
                high += iter
            else:
                if iceTopSin >= 1:
                    iceTopSin = 0
                    low += iter
                    high += iter
                else:
                    iceTopSin += 0.1
        continue
    with open(to_open, 'rb') as fin:
        output = pickle.load(fin)

    print(f'runid of {new_runid}')
    for runid in output:
        print(runid)
        
        condensed_output[new_runid] = {}
        condensed_output[new_runid]['n'] = output[runid]['n']
        condensed_output[new_runid]['n_dip_dir'] = output[runid]['n_dip_dir']
        condensed_output[new_runid]['n_dip_refl'] = output[runid]['n_dip_refl']
        condensed_output[new_runid]['n_lpda_refl'] = output[runid]['n_lpda_refl']
        condensed_output[new_runid]['n_lpda_dir'] = output[runid]['n_lpda_dir']
        condensed_output[new_runid]['energy'] = output[runid]['energy']
        condensed_output[new_runid]['zenith'] = output[runid]['zenith']
        condensed_output[new_runid]['azimuth'] = output[runid]['azimuth']
        condensed_output[new_runid]['x_dir_lpda'] = output[runid]['x']
        condensed_output[new_runid]['y_dir_lpda'] = output[runid]['y']
        condensed_output[new_runid]['dip_dir_mask'] = output[runid]['dip_dir_mask']
        condensed_output[new_runid]['dip_refl_mask'] = output[runid]['dip_refl_mask']
        condensed_output[new_runid]['lpda_refl_mask'] = output[runid]['lpda_refl_mask']
        condensed_output[new_runid]['lpda_dir_mask'] = output[runid]['lpda_dir_mask']
        condensed_output[new_runid]['ant_zen'] = output[runid]['ant_zen']
        condensed_output[new_runid]['dip_dir_SNR'] = output[runid]['dip_dir_SNR']
        condensed_output[new_runid]['dip_refl_SNR'] = output[runid]['dip_refl_SNR']
        condensed_output[new_runid]['lpda_dir_SNR'] = output[runid]['lpda_dir_SNR']
        condensed_output[new_runid]['lpda_refl_SNR'] = output[runid]['lpda_refl_SNR']
        
#        for run in runid:
        """
        n.append(output[runid]['n'])
        n_trig_up.append(output[runid]['n_triggered_upward'])
        n_trig.append(output[runid]['n_triggered'])
        n_trig_att.append(output[runid]['n_att_triggered'])
        energy.append(output[runid]['energy'])
        ss_energy.append(output[runid]['ss_energy'])
        zenith.append(output[runid]['zenith'])
        azimuth.append(output[runid]['azimuth'])
        x.append(output[runid]['x'])
        y.append(output[runid]['y'])
        up_mask.append(output[runid]['up_mask'])
        refl_mask.append(output[runid]['refl_mask'])
        """
        new_runid += 1

    if type != 'IceTop':
        low += iter
        high += iter
    else:
        if iceTopSin == -1:
            low += iter
            high += iter
        else:
            if iceTopSin >= 1:
                iceTopSin = 0
                low += iter
                high += iter
            else:
                iceTopSin += 0.1

    if sys.getsizeof(condensed_output) > 10000:
        saveFile(condensed_output, save_file_as, suffix+f'{split_iD}')
        split_iD += 1
        condensed_output = {}


print(f'ended and saving')
saveFile(condensed_output, save_file_as, suffix+f'{split_iD}')
quit()

print(f'size of dict is {sys.getsizeof(condensed_output)}')

save_file_as = f'FootprintAnalysis/data/CoREAS_{CoREAS_mode}_{config}_Layer{depthLayer}m_{dB}dB_Area{spacing:.2f}_{cores}cores.pkl'

with open(save_file_as, 'wb') as fout:
    print(f'saved file as {save_file_as}')
    pickle.dump(condensed_output, fout, protocol=pickle.HIGHEST_PROTOCOL)

