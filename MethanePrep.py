## -------------------------------------- load some libraries -------------------------------------- 

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

##% -------------------------------------- download the database -------------------------------------- 

r'''
nu_stop = 7800
nu_start = 6500

# get the database from online, as needed
sim = pld.SpectraSim()

sim.def_environment(300, 1, 1, nu_start, nu_stop) # T[K], P[atm], pathlength [m], wvn start, wvn stop
sim.def_molecules(["CH4"], [1]) # molecule and mole fraction
sim.print_Hitran()

sim.simulate()
please = stophere

r'''

##% -------------------------------------- setup model parameters -------------------------------------- 

r'''
nu = np.arange(2200,3500,0.03) # wavenumber range

T = 300 # temperature in K
P = 500 # pressure in torr
PL = 15 # cm length of cell in furnace (double pass)
ych4 = 1.0 # fraction ch4

##% -------------------------------------- run the model -------------------------------------- 


pld.db_begin('linelists')  # load the linelists into Python
mod, pars = td.spectra_single_lmfit() # this makes a single-path H2O cell path to model using lmfit, which is available online

pars['pathlength'].value = PL # pathlength in cm
pars['molefraction'].value = ych4 # mole fraction
pars['mol_id' ].value = 6 # water = 1, ch4 = 6 (hitran molecular code)

pars['temperature'].value = T # temperature in K
pars['pressure'].value = P / 760 # pressure in atm (converted from Torr)

absorbance = np.fft.rfft(mod.eval(xx=nu, params=pars, name='CH4')) # 'H2O' or 'CH4' needs to be in db_begin directory
        

f = open('PureMethaneNate.pckl', 'wb')
pickle.dump([absorbance, T, P, PL, nu], f)
f.close() 

please = stop
r'''

##% -------------------------------------- functions to visualize the data -------------------------------------- 

def plotspectra(nu, absorbance, T_all, P_all, PL, T_plot, P_plot, lgnd='upper left'):

    if type(T_plot) is list: 
     
        j = np.where(P_all == P_plot)[0][0]
        
        for Ti in T_plot:
        
            i = np.where(T_all == Ti)[0][0]
            
            plt.plot(nu, absorbance[i,j,:], label = str(T_all[i])+'K, '+str(P_all[j])+'T')
            
        if np.max(absorbance[1,1,:]) > 1:
            plt.title('Pure Water Absorbance at '+str(P_all[j])+' torr and '+str(PL)+' cm')
            plt.ylabel('Absorbance')
        else:
            plt.title('Pure Water Absorption at '+str(P_all[j])+' torr and '+str(PL)+' cm')
            plt.ylabel('Absorption')
            
    else:
        
        i = np.where(T_all == T_plot)[0][0]
        
        for Pj in P_plot:
            
           j = np.where(P_all == Pj)[0][0]
            
           plt.plot(nu, absorbance[i,j,:], label='P = '+str(P_all[j])+' torr')
            
        if np.max(absorbance[1,1,:]) > 1.1:
            plt.title('Pure Water Absorbance at '+str(T_all[i])+' Kelvin and '+str(PL)+' cm')
            plt.ylabel('Absorbance')
        else:
            plt.title('Pure Water Transmission at '+str(T_all[i])+' Kelvin and '+str(PL)+' cm')
            plt.ylabel('Transmission')
            
    plt.rcParams.update({'figure.autolayout':True,'lines.linewidth':1.5})
    
    if type(lgnd) is str:
        plt.legend(loc=lgnd)
    
    plt.xlabel('Wavenumber ($cm^{-1}$)')
    
    return

##% -------------------------------------- execute functions to visualize the data -------------------------------------- 

f = open('PureMethaneNate.pckl', 'rb')
[absorbance, T, P, PL, nu] = pickle.load(f)
f.close()


lamb = 10000000 / nu

transmission = np.exp(-absorbance)

plt.figure(figsize=(6, 4), dpi=200, facecolor='w', edgecolor='k')  
plt.rcParams.update({'figure.autolayout':True,'lines.linewidth':0.8})

plt.plot(lamb, transmission)













