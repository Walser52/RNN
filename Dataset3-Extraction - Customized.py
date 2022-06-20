# -*- coding: utf-8 -*-
"""
Created on Fri May 27 08:14:02 2022

@author: fes33
"""

import pandas as pd
import os
import re
import numpy as np
from scipy.stats import linregress


df_res = pd.DataFrame(columns=['cycle', 'Voltages', 'C_rate', 'D_rate', 'Tem', 'Capacity'])

files = os.listdir('../Dataset_3_NCM_NCA_battery/')
for file in range(len(files)):
    #General data
    Tem = int(files[file][2:4]) #Temperautre from file name
    
    data_r = pd.read_csv(os.path.join('../Dataset_3_NCM_NCA_battery/', files[file])) #Complete data
    k = np.min(data_r['cycle number'].values) #First cycle
    data_k = data_r[data_r['cycle number'] == k] #Data for first cycle
    Q_p = np.max(np.array(data_k['Q discharge/mA.h'])) #Discharge during first cycle
    delta = 1
    
    #Relaxation voltage data for each cycle 
    for i in range(int(np.min(data_r['cycle number'].values)), int(np.max(data_r['cycle number'].values))+1):
        filename = files[file]
        data_i = data_r[data_r['cycle number'] == i]    #Get ith cycle out of complete set
        Ecell = np.array(data_i['Ecell/V'])             #Get All the cell voltages
        Q_dis = np.array(data_i['Q discharge/mA.h'])    #All discharge capacities
        Current = np.array(data_i['<I>/mA'])            #All currents
        control = np.array(data_i['control/V/mA'])      #All control values
        cr = np.array(data_i['control/mA'])[1] / 2500   #C-Rate (Charging)
        cr_d = int(files[file][8])                      #C-Rate (Discharging)
        
        if np.max(Q_dis) < 1650 or np.max(Q_dis) > 2510:#Q_dis too high/too low skip.
            delta = delta + 1
            continue
        # Remove points where capacity changes too quickly
        if np.abs(np.max(Q_dis) - Q_p) > delta * 10:
            delta = delta + 1
            continue                                    #Move to next cycle.
        delta = 1
        Q_p = np.max(Q_dis)
        index = np.where(np.abs(control) == 0)          #Find points of relaxation (control == 0)
        if index[0][0] > 0:                             #index[0]=array, index[1]=type     
            start = index[0][0]                         #1st Relaxation Phase
        else:                                           #Is this ever used?
            start = index[0][14]
            print(i)
            
#Find remaining data     
    #Extract from CC
        fil_charge_cc = (data_i['control/mA']>0)
        len_cc = len(data_i[fil_charge_cc])
        vstart_cc = data_i['Ecell/V'].iloc[0]
        tstart_cc = data_i['time/s'].iloc[0]
        tend_cc = data_i[fil_charge_cc]['time/s'].iloc[len_cc-1]
        
        #Features CC
        vgap_cc = np.max(data_i[fil_charge_cc]['Ecell/V'])-vstart_cc
        v63perup_cc = vstart_cc + vgap_cc * 0.63
        t_v63perup_cc = data_i[data_i['Ecell/V'] > v63perup_cc]['time/s'].iloc[0]-tstart_cc
        totaltime_cc = tend_cc - tstart_cc

    #Extract from CV
        fil_charge_cv = (data_i['control/V']>0)
        len_cv = len(data_i[fil_charge_cv])
        istart_cv = data_i[fil_charge_cv]['<I>/mA'].iloc[0]
        iend_cv= data_i[fil_charge_cv]['<I>/mA'].iloc[len_cv-1]
        tstart_cv = data_i[fil_charge_cv]['time/s'].iloc[0]    
        tend_cv = data_i[fil_charge_cv]['time/s'].iloc[len_cv-1]
        
        
        #Features CV
        igap_cv         = istart_cv - iend_cv
        i63perdown_cv   = istart_cv - igap_cv * 0.63
        
        fil_currLT63per = (data_i['<I>/mA'] < i63perdown_cv)
        
        t_igap63perdown_cv  = data_i[fil_charge_cv & fil_currLT63per]['time/s'].iloc[0]-tstart_cv
        totaltime_cv    = tend_cv - tstart_cc
    
    #Extract from 1st Relaxaxtion phase
        if control[start + 19] == 0:   
            a = linregress(np.linspace(0,59,59, endpoint = False), y=Ecell[start:start + 59])
        #linear regression parameters + p-value + std error
        
        if control[start + 19] == 0:                    #If this is not the 2nd phase then write the data.
            df_res = df_res.append(
                {#General
                 'file': filename, 
                 'C_rate': cr, 
                 'D_rate': cr_d, 
                 'Tem': Tem,
                 #cycle
                 'Capacity': np.max(Q_dis),
                 'cycle': i, 
                 #CC
                 'V_gap(V)_cc':vgap_cc,
                 'V_63perup(V)_cc':v63perup_cc,
                 't_gap63perup_cc':t_v63perup_cc,
                 'time_cc':totaltime_cc, 
                 
                 #CV
                 'I_gap(mA)_cv': igap_cv,
                 'I_63perdown(mA)_cv': i63perdown_cv,
                 't_gap63perdown_cv': t_igap63perdown_cv,
                 'time_cv': totaltime_cv,
                 #Relaxation
                 'Voltages': Ecell[start:start + 59],
                 'slope_rv': a.slope, 
                 'intercept_rv': a.intercept, 
                 'rval_rv': a.rvalue, 
                 'pval_rv': a.pvalue, 
                 'stderr_rv': a.stderr
                 
                 }, ignore_index=True)

# Save to excel file
#df_res.to_excel('Dataset_3_NCM_NCA_battery.xlsx', index=False)
# Or save to csv file
df_res.to_csv('../Dataset_3_NCM_NCA_battery-customized.csv', index=False)
print('Features extraction is done')