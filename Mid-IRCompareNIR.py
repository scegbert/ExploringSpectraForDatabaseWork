#%% -------------------------------------- load some libraries -------------------------------------- 

import numpy as np
import pickle 
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from sys import path
path.append(r'C:\Users\scott\Documents\1-WorkStuff\code\scottcode\modules')

from packfind import find_package
find_package('pldspectrapy')
import pldspectrapy as pld
import td_support as td # time doamain support
import linelist_conversions as db

import clipboard_and_style_sheet
clipboard_and_style_sheet.style_sheet()

#%% -------------------------------------- setup model parameters -------------------------------------- 

file_which = 'NIR'

nu2_NIR = [6600, 7600]
nu2_cu = [1500, 3500]

nu_step = 0.005 # 0.03 = 1 GHz, Ron at 0.3 = 10 GHz

T = [300] # [2000, 2500, 3000] # temperature in K
P = [0.8] # [3, 2.5, 2] # pressure in atm
y = [1e-4] # [0.2, 0.1, 0.05, 0.02] # see above

PL = 100 # 1.27 # cm shock tube is 0.5" ID

molecule = 'H2O'
molecule_id = 1 # id's from HITRAN

#%% -------------------------------------- download the database -------------------------------------- 

r'''
# get the database from online, as needed
sim = pld.SpectraSim()

sim.def_environment(300, 1, 1, nu_start, nu_stop) # T[K], P[atm], pathlength [m], wvn start, wvn stop

for molecule in molecules: 
    sim.def_molecules([molecule], [1]) # molecule and mole fraction

sim.print_Hitran()

sim.simulate()
please = stophere
r'''

#%% -------------------------------------- run the model -------------------------------------- 

r'''

if file_which == 'NIR': 
    nu = np.arange(nu2_NIR[0],nu2_NIR[1],nu_step) # wavenumber range
    pld.db_begin('data - HITRAN 2020')  # load the linelists into Python
else: 
    nu = np.arange(nu2_cu[0],nu2_cu[1],nu_step) # wavenumber range
    pld.db_begin('linelists')  # load the linelists into Python
 
mod, pars = td.spectra_single_lmfit() # this makes a single-path H2O cell path to model using lmfit, which is available online

pars['pathlength'].value = PL # pathlength in cm
pars['mol_id' ].value = molecule_id # water = 1, ch4 = 6 (hitran molecular code)

absorbance = np.zeros((len(T), len(P), len(y),len(nu)))

for i in range(len(T)): 
    pars['temperature'].value = T[i] # temperature in K
    
    for j in range(len(P)): 
        pars['pressure'].value = P[j] # pressure in atm 
        
        for k in range(len(y)): 
            pars['molefraction'].value = y[k] # mole fraction
            
            absorbance[i,j,k,:] = np.fft.rfft(mod.eval(xx=nu, params=pars, name=molecule)) # 'H2O' or 'CH4' needs to be in db_begin directory   
            print(str(T[i]) +'     '+ str(P[j]) +'     '+ str(y[k]))

f = open('MID-IRCombustionProducts - '+file_which+' background.pckl', 'wb')
pickle.dump([absorbance, T, P, y, PL, nu, molecule], f)
f.close() 

r'''

#%% -------------------------------------- read in the data -------------------------------------- 


f = open('MID-IRCombustionProducts - '+file_which+' fine.pckl', 'rb')
[absorbance_cu, T, P, y, PL_cu, nu_cu, molecules] = pickle.load(f)
f.close()

f = open('MID-IRCombustionProducts - '+file_which+' 1 GHz.pckl', 'rb')
[absorbance_cu_res, _, _, _, _, nu_cu_res, _] = pickle.load(f)
f.close()

f = open('MID-IRCombustionProducts - '+file_which+' background.pckl', 'rb')
[absorbance_cu_back, _, _, _, _, _, _] = pickle.load(f)
f.close()

colors = ['navy', 'g', 'firebrick', 'orangered']	


#%% -------------------------------------- plot cu's stuff -------------------------------------- 

T_which = [2000, 2500, 3000] # [2000, 2500, 3000] # temperature in K
P_which = 3 # [3, 2.5, 2] # pressure in atm
y_which = 0.02 # [0.2, 0.1, 0.05, 0.02] # concentration


if type(T_which) is list: 
    transmission_cu = np.exp(-absorbance_cu[:, P.index(P_which), y.index(y_which), :])
    transmission_cu_res = np.exp(-absorbance_cu_res[:, P.index(P_which), y.index(y_which), :])
    prop_which = 'T'
    prop_units = 'K'
    prop_vals = T_which
    plot_title = 'H2O absorption at P = '+str(P_which)+' and y = '+str(y_which)
    
elif type(P_which) is list: 
    transmission_cu = np.exp(-absorbance_cu[T.index(T_which), :, y.index(y_which), :])
    transmission_cu_res = np.exp(-absorbance_cu_res[T.index(T_which), :, y.index(y_which), :])
    prop_which = 'P'
    prop_units = 'atm'
    prop_vals = P_which
    plot_title = 'H2O absorption at T = '+str(T_which)+' and y = '+str(y_which)
    
elif type(y_which) is list: 
    transmission_cu = np.exp(-absorbance_cu[T.index(T_which), P.index(P_which), :, :])
    transmission_cu_res = np.exp(-absorbance_cu_res[T.index(T_which), P.index(P_which), :, :])
    prop_which = 'y'
    prop_units = '%'
    prop_vals = y_which
    plot_title = 'H2O absorption at T = '+str(T_which)+' and P = '+str(P_which)

transmission_cu_back = np.exp(-absorbance_cu_back[0,0,0,:])

lambda_cu = 10000000 / nu_cu
lambda_cu_res = 10000000 / nu_cu_res

plt.figure(figsize=(6, 4), dpi=200, facecolor='w', edgecolor='k')  
plt.rcParams.update({'figure.autolayout':True,'lines.linewidth':0.8})
plt.title(plot_title)

for i in range(len(prop_vals)):  

    plt.plot(lambda_cu, transmission_cu[i,:], color=colors[i], label=str(prop_vals[i]) + ' ' + prop_units)
    plt.plot(lambda_cu_res, transmission_cu_res[i,:], '-.', color=colors[i])
    
plt.plot(lambda_cu, transmission_cu_back, '--', color='k', label='purged bkgd (1 m, 1% RH)')
plt.plot(lambda_cu, 1-(transmission_cu[i,:]-transmission_cu_back))

plt.xlabel('wavelength (nm)')
plt.ylabel('transmission')
plt.legend(loc='lower right')


# plt.xlim([3200, 3900])
# plt.ylim([0.993,1.0001])

# plt.xlim([3338.2, 3340.8])
# plt.ylim([0.9962,1.0001])

# plt.xlim([4000, 5000])
# plt.ylim([0.94,1.001])
# plt.legend(loc='lower left')

# plt.xlim([4782.2, 4785.8])
# plt.ylim([0.972,1.001])




