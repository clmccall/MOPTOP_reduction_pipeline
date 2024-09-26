from calendar import c
import moptop_settings as mopset_1
import moptop_dicts as mopinf
import glob
import pandas as pd
import numpy as np
import shutil
import ast
import scipy.stats as stats
from astropy.stats import sigma_clip
import csv 
import os
from pylab import *
import warnings
from IPython.display import clear_output
from pathlib import Path
from astropy.io import fits
from scipy.optimize import curve_fit
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
pd.options.mode.chained_assignment = None  # default='warn'

pol_sources = ['HD251204','BD+64 106','VICyg12','BD+57 2615','HD155197', 'HILT960']
zpol_sources = ['HD14069','BD+32 3739','GD319','HD212311','G191B2B','BD+28 4211']

all_zpol_data = []
for source in zpol_sources:
    try:
        data = pd.read_csv(mopset_1.dir1+source+'/all_data.csv')
        data = data[data['photometric']=="PHOTOMETRIC"]
        data['Source'] = source
        all_zpol_data.append(data)
    except:
        continue
combined_zpol_data = pd.concat(all_zpol_data, ignore_index=True)
combined_zpol_data['epoch'] = np.nan

for i in range(len(combined_zpol_data['mjd'])):
    if combined_zpol_data['mjd'][i] < mopset_1.change_pol_constants[0]:
        combined_zpol_data['epoch'][i] = 1
    elif combined_zpol_data['mjd'][i] >= mopset_1.change_pol_constants[-1]:
        combined_zpol_data['epoch'][i] = len(mopset_1.change_pol_constants)+1
    else:
        for j in range(len(mopset_1.change_pol_constants) - 1):
            if mopset_1.change_pol_constants[j] <= combined_zpol_data['mjd'][i] < mopset_1.change_pol_constants[j+1]:
                combined_zpol_data['epoch'][i] = j + 2
                break

qu_zero = pd.DataFrame()
for epoch in unique(combined_zpol_data['epoch']):
    for filt in unique(combined_zpol_data['wave']):
        data = combined_zpol_data[(combined_zpol_data['epoch'] == epoch) & (combined_zpol_data['wave'] == filt)]
        data = data.reset_index(drop=True)

        q_zero = median(data['q_avg'])
        q_zero_err = np.std(data['q_avg'])/np.sqrt(len(data))
        u_zero = median(data['u_avg'])
        u_zero_err = np.std(data['u_avg'])/np.sqrt(len(data))

        qu_zero = qu_zero._append({'Epoch': int(epoch), 'Filter': filt, 'q_zero': f"{q_zero:.4f} ± {q_zero_err:.4f}", 'u_zero': f"{u_zero:.4f} ± {u_zero_err:.4f}"}, ignore_index=True)

qu_zero['Filter'] = pd.Categorical(qu_zero['Filter'], categories=['B','V','R','I','L'], ordered=True)
qu_zero = qu_zero.sort_values('Filter').reset_index(drop=True)
print(qu_zero)

all_pol_data = []
for source in pol_sources:
    try:
        data = pd.read_csv(mopset_1.dir1+source+'/all_data.csv')
        data = data[data['photometric']=="PHOTOMETRIC"]
        data['Source'] = source
        all_pol_data.append(data)
    except:
        continue
combined_pol_data = pd.concat(all_pol_data, ignore_index=True)
combined_pol_data['epoch'] = np.nan

for i in range(len(combined_pol_data['mjd'])):
    if combined_pol_data['mjd'][i] < mopset_1.change_pol_constants[0]:
        combined_pol_data['epoch'][i] = 1
    elif combined_pol_data['mjd'][i] >= mopset_1.change_pol_constants[-1]:
        combined_pol_data['epoch'][i] = len(mopset_1.change_pol_constants)+1
    else:
        for j in range(len(mopset_1.change_pol_constants) - 1):
            if mopset_1.change_pol_constants[j] <= combined_pol_data['mjd'][i] < mopset_1.change_pol_constants[j+1]:
                combined_pol_data['epoch'][i] = j + 2
                break

k_values = pd.DataFrame()
for epoch in unique(combined_pol_data['epoch']):
    for filt in unique(combined_pol_data['wave']):
        data = combined_pol_data[(combined_pol_data['epoch'] == epoch) & (combined_pol_data['wave'] == filt)]
        data = data.reset_index(drop=True)

        q = data['q_avg'] - float(qu_zero[(qu_zero['Epoch'] == epoch) & (qu_zero['Filter'] == filt)]['q_zero'].values[0].split(' ± ')[0])
        u = data['u_avg'] - float(qu_zero[(qu_zero['Epoch'] == epoch) & (qu_zero['Filter'] == filt)]['u_zero'].values[0].split(' ± ')[0])

        theta = np.degrees(0.5*np.arctan2(u,q)) % 180
        theta_cor = (data['rotskypa'] + theta) % 180
        theta_cor_act = data.apply(lambda row: mopinf.source_info[row['Source']]['theta_cor_act_' + filt], axis=1)
        k_value = (theta_cor_act - theta_cor) % 180
        k_value = k_value.apply(lambda x: x - 180 if x > 90 else x) % 180
        k = median(k_value)
        k_err = np.std(k_value)/np.sqrt(len(data))
        k_values = k_values._append({'Epoch': int(epoch), 'Filter': filt, 'k_value': f"{k:.2f} ± {k_err:.2f}"}, ignore_index=True)

k_values['Filter'] = pd.Categorical(k_values['Filter'], categories=['B','V','R','I','L'], ordered=True)
k_values = k_values.sort_values('Filter').reset_index(drop=True)
print(k_values)

depol_values = pd.DataFrame()
for epoch in unique(combined_pol_data['epoch']):
    for filt in unique(combined_pol_data['wave']):
        data = combined_pol_data[(combined_pol_data['epoch'] == epoch) & (combined_pol_data['wave'] == filt)]
        data = data.reset_index(drop=True)

        q = data['q_avg'] - float(qu_zero[(qu_zero['Epoch'] == epoch) & (qu_zero['Filter'] == filt)]['q_zero'].values[0].split(' ± ')[0])
        u = data['u_avg'] - float(qu_zero[(qu_zero['Epoch'] == epoch) & (qu_zero['Filter'] == filt)]['u_zero'].values[0].split(' ± ')[0])

        p = np.sqrt(q**2+u**2)*100
        p_act = data.apply(lambda row: mopinf.source_info[row['Source']]['per_p_act_' + filt], axis=1)
        
        def linear_func(x, a):
            return a * x

        m, c = curve_fit(linear_func, p_act, p)
        depol = m[0]
        depol_err = np.sqrt(np.diag(c))[0]
        depol_values = depol_values._append({'Epoch': int(epoch), 'Filter': filt, 'depol': f"{depol:.2f} ± {depol_err:.2f}"}, ignore_index=True)

depol_values['Filter'] = pd.Categorical(depol_values['Filter'], categories=['B','V','R','I','L'], ordered=True)
depol_values = depol_values.sort_values('Filter').reset_index(drop=True)
print(depol_values)
