""" 
Plot the tangential and azimuthal loads along the wind turbine blade, for
the DTU10MW reference wind turbine, from HAWC2S and PyWakeEllipSys.
"""
from dtu10mw_des_bladeloads import dtu10mw_bladeloading
import numpy as np

import matplotlib.pyplot as plt


rho = 1.225
UH = 8

# %% Load HAWC2 bladeloads
hawc2_ind_path = 'D:\AE EWEM MSc\T H E S I S\\6_code\code repository\my_dtu_10mw\DTU_10MW_rigid_hawc2s_flattened_u8000.ind'

data = np.loadtxt(hawc2_ind_path)
s = data[:, 0]      # Blade span [m]
S = s[-1] + 0.001   # Maximum blade span (w/ tip fix) [m]
s = np.append(s, S)  # Tip fix
Ft = data[:, 6]      # Thrust force per unit length [N/m]
Ft = np.append(Ft, 0)  # Tip fix
Fn = data[:, 7]      # Azumuthal force per unit length [N/m]
Fn = np.append(Fn, 0)  # Tip fix

# Normalisation of forces, non-dimensionalisation of radius
s_nd_h2 = s / S           # Blade span, non-dimensioned [-]
ft_h2 = Ft / (rho * S * UH**2)  # Normalised tangential force [-]
fn_h2 = Fn / (rho * S * UH**2)  # Normalised azimuthal force [-]

# %% Loads PWES bladeloads (ugly af coding, 14 days left to give in lol)

jou_path = r'D:\AE EWEM MSc\T H E S I S\\6_code\code repository\run_data\thrust\joukowsky_0.35\adbladeloads.dat'
dtu10mw_path = r'D:\AE EWEM MSc\T H E S I S\\6_code\code repository\run_data\thrust\airfoil_dtu10mw\adbladeloads.dat'
wdr10mw_path = r'D:\AE EWEM MSc\T H E S I S\\6_code\code repository\run_data\thrust\airfoil_wdr10mw\adbladeloads.dat'

# Joukowsky
loads_data_jou = np.loadtxt(jou_path)

generator_loss = 0.064  # [%], already accounted for!
# WHY is it used for thrust curve as well?
correction = 1 / (1 - generator_loss)
# anyway, if it is used for the thrust curve, then it should only be applied to the thrust force, not azimuthal
# Manipulate data
r = loads_data_jou[:, 0]  # Blade radius [m]
Ft = loads_data_jou[:, 3]  # Tangential force per unit length [N/m] (assumed)
Fn = loads_data_jou[:, 4]  # Azimuthal force per unit length [N/m] (assumed)
# Normalisation of forces, non-dimensionalisation of radius
R = r[-1]  # Rotor max radius [m]
r_nd = r/R  # Non-dimensionalised blade radius r/R [-]
Ft_jou = Ft / (rho * R * UH**2)  # * correction
Fn_jou = Fn / (rho * R * UH**2)  # * correction
print(f"AD joukowsky R={R}m")

# DTU10MW
loads_data_dtu = np.loadtxt(dtu10mw_path)

# Manipulate data
r = loads_data_dtu[:, 0]  # Blade radius [m]
Ft = loads_data_dtu[:, 3]  # Tangential force per unit length [N/m] (assumed)
Fn = loads_data_dtu[:, 4]  # Azimuthal force per unit length [N/m] (assumed)
# Normalisation of forces, non-dimensionalisation of radius
R = r[-1]  # Rotor max radius [m]
r_nd = r/R  # Non-dimensionalised blade radius r/R [-]
Ft_dtu = Ft / (rho * R * UH**2)
Fn_dtu = Fn / (rho * R * UH**2) * correction
print(f"AD air R={R}m")
# WDR10MW
loads_data_wdr = np.loadtxt(wdr10mw_path)

# Manipulate data
r = loads_data_wdr[:, 0]  # Blade radius [m]
Ft = loads_data_wdr[:, 3]  # Tangential force per unit length [N/m] (assumed)
Fn = loads_data_wdr[:, 4]  # Azimuthal force per unit length [N/m] (assumed)
# Normalisation of forces, non-dimensionalisation of radius
R = r[-1]  # Rotor max radius [m]
r_nd = r/R  # Non-dimensionalised blade radius r/R [-]
Ft_wdr = Ft / (rho * R * UH**2)
Fn_wdr = Fn / (rho * R * UH**2)

# DES
# Non-dimensionalised blade radius r/R [-]
r_nd_val = dtu10mw_bladeloading[:, 0]
# Non-dimensionalised tangential force [-] (assumed axis system
ft_val = -dtu10mw_bladeloading[:, 1]
# Non-dimensionalised azimuthal force [-] (assumed)
fn_val = dtu10mw_bladeloading[:, 3]

fig, ax = plt.subplots(1, 1, sharex=True, figsize=[13, 5])

#ax.plot(r_nd_val, fn_val, label="DES DTU10MW")

ax.plot(s_nd_h2, fn_h2, label="HAWC2S DTU10MW")

#ax.plot(r_nd, Fn_jou, label="Joukowsky AD $\overline{\delta}=0.35$")

Fn_dtu[-2] = 0.02
Fn_wdr[-2] = 0.02
ax.plot(r_nd, Fn_dtu, label="AD Airfoil DTU10MW")
ax.plot(r_nd, Fn_wdr, label="AD Airfoil WDR10MW")


###
# ax[0].plot(r_nd_val, ft_val, label="DES DTU10MW")
# ax[1].plot(r_nd_val, fn_val, label="DES DTU10MW")
# ax[0].plot(s_nd_h2, ft_h2, label="HAWC2S DTU10MW")
# ax[1].plot(s_nd_h2, fn_h2, label="HAWC2S DTU10MW")
# ax[0].plot(r_nd, Ft_jou, label="AD Joukowsky $\overline{\delta}=0.35$")
# ax[1].plot(r_nd, Fn_jou, label="AD Joukowsky $\overline{\delta}=0.35$")
# ax[0].plot(r_nd, Ft_dtu, label="AD Airfoil DTU10MW")
# ax[1].plot(r_nd, Fn_dtu, label="AD Airfoil DTU10MW")
#ax[0].plot(r_nd, Ft_wdr, label="AD Airfoil WDR10MW")
#ax[1].plot(r_nd, Fn_wdr, label="AD Airfoil WDR10MW")

lfsize = 28

ax.set_xlabel("$r/R$ $[-]$", fontsize=lfsize)
ax.set_ylabel("$f_T$ $[-]$", fontsize=lfsize)
###
# ax[1].set_xlabel("$r/R$ $[-]$", fontsize=lfsize)
# ax[0].set_ylabel("$f_T$ $[-]$", fontsize=lfsize)
# ax[1].set_ylabel(r"$f_\theta$ $[-]$", fontsize=lfsize)

ax.set_xlim([0, 1])
ax.set_ylim([0, 0.8])
# ax[0].set_xlim([0, 1])
plt.legend(loc='best')

fig.savefig('plots/ad_bladeloads_thrust_wdr.pdf', bbox_inches='tight')
