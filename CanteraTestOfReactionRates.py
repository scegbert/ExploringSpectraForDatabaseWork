r'''
file to test reaction rates for CO in H2 and CH4 in H2 
before blowing myself up trying to measure Jovian conditions
r'''

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

gas = ct.Solution('gri30.yaml')
gas.basis = 'molar'

X = 0.02

target = 'CO' # 'CH4', 'CO'
bath = 'air' # 'air' 'Ar', 'H2'

if bath == 'air': 
    gas.X = {target:X, 'O2':(1-X)*0.21, 'N2':(1-X)*0.78}    
else: 
    gas.X = {target:X, bath:1-X}


T_all = np.linspace(300, 1500) # K
P_all = [10, 30, 90, 180, 400, 600] # Torr (converted to Pa later)

co_left = np.zeros((len(P_all), len(T_all)))

for i_P, P in enumerate(P_all): 

    for i_T, T in enumerate(T_all): 
    
        gas.TP = T, P*133.3
        
        d_destruct = gas.destruction_rates[14] # kmol / m3 / s
        v_CO = gas.v * X # specifiv volume of CO in m3 / kmol 
        
        d_destruct = d_destruct*3600 * v_CO # destruction rate of CO per hour
        
        co_left_i = 1 - d_destruct
        
        co_left[i_P, i_T] = co_left_i
    
    plt.plot(T_all, co_left[i_P, :], label='{} Torr'.format(int(P)))

# plt.yscale('log')
plt.legend()





















