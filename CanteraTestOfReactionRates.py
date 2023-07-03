r'''
file to test reaction rates for CO in H2 and CH4 in H2 
before blowing myself up trying to measure Jovian conditions
r'''

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

gas = ct.Solution('gri30.yaml')
gas.basis = 'molar'

X = 0.05

target = 'CH4' # 'CH4', 'CO'
bath = 'air' # 'air' 'Ar', 'H2'

if bath == 'air': 
    gas.X = {target:X, 'O2':(1-X)*0.21, 'N2':(1-X)*0.78}    
else: 
    gas.X = {target:X, bath:1-X}


T_all = np.linspace(300, 1300) # K
P_all = [640] # [40, 80, 160, 320, 640] # Torr (converted to Pa later)

target_left = np.zeros((len(P_all), len(T_all)))
bath_left = np.zeros((len(P_all), len(T_all)))


plt.figure(figsize=(8,4))

for i_P, P in enumerate(P_all): 

    for i_T, T in enumerate(T_all): 
    
        gas.TP = T, P*133.3
        
        # d_destruct_H2 = gas.destruction_rates[0] # kmol / m3 / s of H2
        d_destruct_H2 = gas.destruction_rates[3] # kmol / m3 / s of H2
        v_bath = gas.v * (1-X) # specific volume of CO in m3 / kmol 
        d_destruct_bath = d_destruct_H2*3600 * v_bath # destruction rate of CO per hour

        
        # d_destruct_target = gas.destruction_rates[14] # kmol / m3 / s of CO
        d_destruct_target = gas.destruction_rates[13] # kmol / m3 / s of CH4

        v_sample = gas.v * X # specific volume of CO in m3 / kmol 
        d_destruct_target = d_destruct_target*3600 * v_sample # destruction rate of CO per hour

        
        bath_left[i_P, i_T] = 1 - d_destruct_bath
        target_left[i_P, i_T] = 1 - d_destruct_target

    
    plt.plot(T_all, target_left[i_P, :], label='{} Torr, {}'.format(int(P), target))
    plt.plot(T_all, bath_left[i_P, :], '--', label='{} Torr, {}'.format(int(P), bath))

# plt.yscale('log')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('Rate of Destruction (% remaining/ hour)')

plt.ylim((0.999, 1.0001))





















