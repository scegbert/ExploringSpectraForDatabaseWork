#%% -------------------------------------- load some libraries -------------------------------------- 

import numpy as np
import pickle 
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import os
from sys import path
path.append(os.path.abspath('..')+'\\modules')

import pldspectrapy as pld
import td_support as td # time doamain support
import linelist_conversions as db

import clipboard_and_style_sheet
clipboard_and_style_sheet.style_sheet()


#%% -------------------------------------- setup model parameters -------------------------------------- 

file_which = 'pure CO' # 'pure CO', 'diluted CO', 'pure CH4', 'diluted CH4'
wvn_step = 0.0033 # 0.03 = 1 GHz, 0.0066 = 200 MHz, 0.0033 = 100 MHz

molecule = file_which.split()[-1]

if molecule == 'CO':

    molecule_id = 5 # id's from HITRAN
    wvn_2 = [1850, 2325]

    PL = 0.1 # cm
    
    T = [300, 500, 700, 900, 1100, 1300] # temperature in K

    if file_which.split()[0] == 'pure': 
        P = [8,4,2,1,0.5] # pressure in Torr
        y = 1

    elif file_which.split()[0] == 'diluted': 
        P = [640,320,160,80,40] # pressure in Torr
        y = 0.05
    
if molecule == 'CH4':

    molecule_id = 6 # id's from HITRAN
    wvn_2 = [2700, 3250]
    
    PL = 0.1 # cm
    
    T = [300, 500, 700, 900, 1100, 1300] # temperature in K

    if file_which.split()[0] == 'pure': 
        P = [16,8,4,2,1] # pressure in Torr
        y = 1

    elif file_which.split()[0] == 'diluted': 
        P = [640,320,160,80,40] # pressure in Torr
        y = 0.05



#%% -------------------------------------- run the model -------------------------------------- 

wvn = np.arange(wvn_2[0],wvn_2[1],wvn_step) # wavenumber range
pld.db_begin('linelists')  # load the linelists into Python
 
mod, pars = td.spectra_single_lmfit() # this makes a single-path H2O cell path to model using lmfit, which is available online

pars['pathlength'].value = PL # pathlength in cm
pars['molefraction'].value = y # mole fraction

pars['mol_id' ].value = molecule_id # water = 1, ch4 = 6 (hitran molecular code)

absorbance = np.zeros((len(T), len(P),len(wvn)))

for i, T_i in enumerate(T): 
    pars['temperature'].value = T_i # temperature in K
    
    for j, P_j in enumerate(P): 
        pars['pressure'].value = P_j/760 # pressure in atm 
                
        absorbance[i,j,:] = np.fft.rfft(mod.eval(xx=wvn, params=pars, name=molecule)) 
        print(str(T_i) +'     '+ str(P_j) +'     '+ str(y))

f = open('MID-IR database test - '+file_which+'.pckl', 'wb')
pickle.dump([absorbance, T, P, y, PL, wvn, molecule], f)
f.close() 


#%% -------------------------------------- read in the data -------------------------------------- 


# f = open('MID-IR database test - '+file_which+'.pckl', 'rb')
# [absorbance, T, P, y, PL, wvn, molecule] = pickle.load(f)
# f.close()



#%% -------------------------------------- plot cu's stuff -------------------------------------- 

colors = ['#0028ff','#0080af','#117d11','#be961e', '#ff4000','#ff0000',     '#e6ab02', '#fee9ac']

transmission = np.exp(-absorbance)

wvl = 10000000 / wvn

plt.figure(figsize=(6, 4), dpi=200, facecolor='w', edgecolor='k')  
plt.rcParams.update({'figure.autolayout':True,'lines.linewidth':0.8})
plt.title(file_which)

P_j = max(P)
j = P.index(P_j)
for i, T_i in enumerate(T):  
    plt.plot(wvl, transmission[i,j,:], color=colors[i], label='P={}Torr, T={}K'.format(P_j, T_i))

plt.xlabel('wavelength (nm)')
plt.ylabel('transmission')
plt.legend(loc='lower right')





plt.figure(figsize=(6, 4), dpi=200, facecolor='w', edgecolor='k')  
plt.rcParams.update({'figure.autolayout':True,'lines.linewidth':0.8})
plt.title(file_which)

T_i = min(T)
i = T.index(T_i)
for j, P_j in enumerate(P):  

    plt.plot(wvl, transmission[i,j,:], color=colors[j], label='P={}Torr, T={}K'.format(P_j, T_i))
    
plt.xlabel('wavelength (nm)')
plt.ylabel('transmission')
plt.legend(loc='lower right')



