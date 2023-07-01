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

## -------------------------------------- download the database -------------------------------------- 
r'''
# get the database from online, as needed
sim = pld.SpectraSim()

sim.def_environment(300, 1, 10e-2, 10e6/nu_stop, 10e6/nu_start) # Temp_Kelvin, Press_atm, pathlength_meter, nm_start, nm_stop
sim.def_molecules(["H2O"], [1])
sim.print_Hitran()

sim.simulate()
r'''
## -------------------------------------- setup model parameters -------------------------------------- 

nu = np.arange(6300,8000,0.005) # wavenumber range

T = np.array([300,500,700,900,1100]) # temperature in K
P = np.array([20,40,60,80,120,160,320,600]) # pressure in torr
PL = 91.4 # cm length of cell in furnace (double pass)
yh2o = 0.02 # water fraction

nu_scott = [6500,7700]

nu_combustor = [6700, 7000]
nu_isolator = [6875, 7250]
nu_silmaril = [6385,7985]
nu_paul = [6801, 7188]

## -------------------------------------- run the model -------------------------------------- 
r'''
pld.db_begin('data')  # load the linelists into Python
mod, pars = td.spectra_single_lmfit() # this makes a single-path H2O cell path to model using lmfit, which is available online

pars['pathlength'].value = PL # pathlength in cm
pars['molefraction'].value = yh2o # mole fraction
pars['mol_id' ].value = 1 # water = 1, ch4 = 6 (hitran molecular code)

absorbance = np.zeros((len(T),len(P),len(nu)))

for i in range(0,len(T)):

    pars['temperature'].value = T[i] # temperature in K
    
    for j in range(0,len(P)):
        
        pars['pressure'].value = P[j] / 760 # pressure in atm (converted from Torr)
            
        print(pars['pathlength'].value, pars['molefraction'].value, 
              pars['temperature'].value, pars['pressure'].value * 760)
        
        absorbance[i,j,:] = np.fft.rfft(mod.eval(xx=nu, params=pars, name='H2O')) # 'H2O' needs to be in db_begin directory
        

f = open('WaterAir2.pckl', 'wb')
pickle.dump([absorbance, T, P, PL, nu], f)
f.close() 
r'''
## -------------------------------------- functions to visualize the data -------------------------------------- 

def plotspectra(nu, absorbance, T_all, P_all, PL, T_plot, P_plot, lgnd='upper left', 
                xlim=[6500, 7700], ylim=None, line50MHz=False, lineMeas=False, lineYmax=False):
    
    plt.figure(figsize=(6, 4), dpi=200, facecolor='w', edgecolor='k')  
    plt.rcParams.update({'figure.autolayout':True,'lines.linewidth':0.8})

    if type(T_plot) is list or type(T_plot) is np.ndarray: 
     
        j = np.where(P_all == P_plot)[0][0]
        
        for Ti in T_plot:
        
            i = np.where(T_all == Ti)[0][0]
            
            plt.plot(nu, absorbance[i,j,:], label='T = '+str(T_all[i])+' Kelvin')
            
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
            
        if np.max(absorbance[1,1,:]) > 1:
            plt.title('Pure Water Absorbance at '+str(T_all[i])+' Kelvin and '+str(PL)+' cm')
            plt.ylabel('Absorbance')
        else:
            plt.title('Pure Water Absorption at '+str(T_all[i])+' Kelvin and '+str(PL)+' cm')
            plt.ylabel('Absorption')
            
    plt.rcParams.update({'figure.autolayout':True,'lines.linewidth':1.5})
    
    if lineMeas:
        plt.axvline(x=nu_combustor[0],color='k', linestyle='--', label='AFRL Combustor')
        plt.axvline(x=nu_combustor[1],color='k', linestyle='--')
        
        plt.axvline(x=nu_isolator[0],color='k', linestyle='-.', label='AFRL Isolator')
        plt.axvline(x=nu_isolator[1],color='k', linestyle='-.')
        
        plt.axvline(x=nu_paul[0],color='k', linestyle=':', label='Paul Database')
        plt.axvline(x=nu_paul[1],color='k', linestyle=':')
    
    if line50MHz:
        plt.axvline(x=nu_silmaril[0],color='r', linestyle='-', label='Sub 50 MHz')
        plt.axvline(x=nu_silmaril[1],color='r', linestyle='-')
        
    if lineYmax:
        plt.axhline(y=lineYmax,color='k', linestyle='-')
    
    if type(lgnd) is str:
        plt.legend(loc=lgnd)
    
    if type(ylim) is list:
        plt.ylim(ylim)
    
    plt.xlim(xlim)
    plt.xlabel('Wavenumber ($cm^{-1}$)')
    
    return

def plotAbsorbance(nu, absorbance, T_all, P_all, nu_plot, T_plot, lgnd='upper left', xlim=None, ylim=None):
    # this function might not resolve/distinguish all peaks (depending on resolution of the model) but it's a good estimate
    
    plt.figure(figsize=(6, 4), dpi=200, facecolor='w', edgecolor='k')  
    plt.rcParams.update({'figure.autolayout':True,'lines.linewidth':0.8})

    delta = int(0.1 / (nu[1] - nu[0]))
    
    if type(T_plot) is list: 
     
        k = np.argmin(abs(nu-nu_plot))
        
        a = 0
        abs_max = np.zeros([len(T_plot),len(P_all)])
                
        for Ti in T_plot:
        
            i = np.where(T == Ti)[0][0]
                        
            for Pj in P_all:
                    
                j = np.where(P_all == Pj)[0][0]
                
                abs_max[a,j] = max(absorbance[i,j,k-delta:k+delta]) # max value of absorption for given feature at T and P
                 
            plt.plot(P_all, abs_max[a,:], label='T = '+str(Ti)+' Kelvin')
            
            a+=1
            
        plt.title('pure water absorbance at '+str(nu_plot)+' $cm^{-1}$ and '+str(PL)+' cm')

    else:
     
        i = np.where(T == T_plot)[0][0]   
       
        a = 0      
        abs_max = np.zeros([len(nu_plot),len(P_all)])
        
        for nuk in nu_plot:
            
            k = np.argmin(abs(nu-nuk))        

            for Pj in P_all:
            
                j = np.where(P_all == Pj)[0][0]
                
                abs_max[a,j] = max(absorbance[i,j,k-delta:k+delta]) # max value of absorption for given feature at T and P
                 
            plt.plot(P_all, abs_max[a,:], label='\u03BD = '+str(nuk)+' $cm^{-1}$')
            
            a+=1
                
        plt.title('pure water absorbance at '+str(T_plot)+' Kelvin and '+str(PL)+' cm')
        
    plt.axhline(y=1.5, color='k', linestyle='--')
    plt.axhline(y=2.6, color='k', linestyle='-.')
    
    if type(lgnd) is str:
        plt.legend(loc=lgnd)
    
    if type(ylim) is list:
        plt.ylim(ylim)
    
    if type(xlim) is list:
        plt.xlim(xlim)
    
    plt.xlabel('Pressure (torr)')
    plt.ylabel('Peak Absorbance')
        
    return

def saturatedPeaks(nu, absorbance, T_all, P_all, T_test, P_test, nu_range=[7000,7400], cutoff_val=1.5):
    
    
    plt.figure(figsize=(6, 4), dpi=200, facecolor='w', edgecolor='k')  
    plt.rcParams.update({'figure.autolayout':True,'lines.linewidth':0.8}) 
    
    hitran = db.par_to_df(r'C:\Users\scott\Documents\1-WorkStuff\code\pldspectrapy-nam\ScottFiles\data\H2O.data')
    hitran_large = hitran[hitran.sw > 1e-21]
    
    imin = np.argmin(abs(nu-nu_range[0]))
    imax = np.argmin(abs(nu-nu_range[1]))
    
    test = absorbance > cutoff_val
    
    nu_filter = nu[imin:imax]
    test_filter = 1 * test[:,:,imin:imax]
    test_filter_shift =  1 * test[:,:,imin+1:imax+1]
    
    change = test_filter - test_filter_shift # value of 1 for each time absoprtion passes cuttoff_val (up and down)
    
    nu_saturated = []
    count_saturated = np.zeros([len(T_test),len(P_test)])
    
    i_test = -1
    
    for Ti in T_test:
                
        i = np.where(T_all == Ti)[0][0]
        i_test += 1
    
        j_test = -1
        
        for Pj in P_test:
            
            j = np.where(P_all == Pj)[0][0]        
            j_test += 1
            
            changeTP = change[i,j,:]
            nu_TP = nu_filter
                 
            c = 0
            
            while max(abs(changeTP)) > 0.5: # while there are values left to snag
                
                k_low = np.argmin(changeTP) # find next -1 (crossing cutoff_val going up)
                k_high = np.argmax(changeTP) # find next 1 (crossing cutoff_val going down)
                
                k_center = int(np.mean([k_low,k_high]))
                nu_center = nu_TP[k_center]
                
                k_ht = np.argmin(abs(hitran_large.nu - nu_center))
                
                saturated = np.array(hitran_large.nu)[k_ht]
               
                if saturated not in nu_saturated:
                    
                    nu_saturated.append(saturated)
                           
                changeTP = changeTP[k_high+1:-1] # remove this feature so you can find the next one
                nu_TP = nu_TP[k_high+1:-1]
    
                c+=1 
                
            count_saturated[i_test,j_test] = c
                
        plt.plot(P_test, count_saturated[i_test,:], 'x-', label='T = '+str(T_all[i])+' Kelvin')   
            
    plt.legend()
    
    plt.xlabel('Pressure (torr)')
    plt.ylabel('Number of Saturated Features (\u03B1 > '+str(cutoff_val)+')')    
        
    return count_saturated, nu_saturated

## -------------------------------------- execute functions to visualize the data -------------------------------------- 

f = open('WaterAir2.pckl', 'rb')
[absorbance, T, P, PL, nu] = pickle.load(f)
f.close()

T_plot = T[0]
P_plot = [40,600]
nu_plot = [6600, 7600]
plotspectra(nu, np.exp(-absorbance), T, P, PL, T_plot, P_plot, xlim=nu_plot, lgnd='upper right', line50MHz=False, lineMeas=False)



r'''
nu_plot = [6300,9650]
plotspectra(nu, absorbance, T, P, PL, T_plot, P_plot, xlim=nu_plot, lgnd='upper right', line50MHz=True, lineMeas=True)


T_plot = T[0]
P_plot = np.flipud(P)
plotspectra(nu, absorbance, T, P, PL, T_plot, P_plot, lineMeas=True)

ysat = 1.5
ylim = [0,2.6]
nu_saturated = [7000,7450]
plotspectra(nu, absorbance, T, P, PL, T_plot, P_plot, xlim=nu_saturated, lgnd=None, ylim=ylim, lineYmax=ysat, lineMeas=True)


nu_sat = [7047.684236, 7047.687745, 7071.48447, 7071.49757, 7094.68301, 7095.85518, 
          7104.61877, 7105.85466, 7117.42026, 7117.75474, 7139.08913, 7139.61009,
          7142.56137, 7145.08447, 7145.21725, 7160.23543, 7161.4101, 7165.82054, 
          7168.43701, 7181.15578, 7182.94962] # Paul's saturated features


T_plot = 300
plotAbsorbance(nu, absorbance, T, P, nu_sat, T_plot, ylim=[0,3], lgnd = None)

T_plot = T
nu_sat = 7047.684236
#plotAbsorbance(nu, absorbance, T, P, nu_sat, T_plot, ylim=[0,3])

       


nu_test = nu_scott
T_plot = T
P_plot = P[-1]
plotspectra(nu, absorbance, T, P, PL, T_plot, P_plot, xlim=nu_test, ylim=[0,2.6], lgnd=None)
T_test = [T[0],T[0]]
P_test = [P[5],P[5]]
cutoff_val = 1.5
count_saturated, nu_saturated = saturatedPeaks(nu, absorbance, T, P, T_test, P_test, nu_range=nu_test, cutoff_val=cutoff_val)

r'''














