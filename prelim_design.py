import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# starting parameters
tit = 1100 # K
T_inlet = 288 # K
C_a = 150 # m/s
rho_infty = 1.22500 # kg/m3
m_desired = 20 # desired mass flow, kg/sec

# need to pin down: hub to tip ratio
hub_tip_ratio = 0.5 # about average, idk
r_tip = np.sqrt(m_desired * (np.pi*rho_infty*C_a*(1-hub_tip_ratio**2))**(-1))
r_hub = r_tip*hub_tip_ratio
print(f'radius of the tip: {r_tip}'
      f'r hub: {r_hub}')

A = lambda r_major, r_minor: np.pi*(r_major**2 - r_minor**2) # annulus dimensions

print(f'speed of sound: {np.sqrt(1.4*288.0*T_inlet)}')
desired_tip_speed = 340
rpm = 60*desired_tip_speed/(2*np.pi*r_tip) # m/s
print(f'rpm : {rpm}')

print(f'tip mach number: {desired_tip_speed/np.sqrt(1.4*288.0*T_inlet)}')

fig, ax = plt.subplots()
ax.add_patch(patches.Rectangle((1.0,1.0),1.0,1.0))
plt.show()