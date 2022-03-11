"""
This file was made for increased ease of automating PyWakeEllipSys simulations.
"""
import numpy as np
import matlotlib as mpl
import matplotlib.pyplot as plt

#%% Plot style, generously donated by Mads Christian Baungaard
mpl.style.use('classic')
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams["legend.scatterpoints"] = 1
plt.rcParams["legend.numpoints"] = 1
plt.rcParams['grid.linestyle'] = ':' # Dotted gridlines
mpl.rcParams['lines.linewidth'] = 2
plt.rcParams['axes.grid']=True
yd = dict(rotation=0,ha='right') # I couldn't find a way to customize these, so use a dict everytime..
plt.close('all')


#%% Wind turbine definition
from py_wake_ellipsys.wind_farm_models.ellipsys_lib.ellipsys_wind_turbines import EllipSysOneTypeWT
from py_wake_ellipsys.examples.data.turbines.ADairfoil import ADairfoil_path


# DTU10MW for Joukowski `2xxx` AD forcing
from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_ct_curve
from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_power_curve
from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_rpm_curve
from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_pitch_curve
from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_bladeloading

class DTU10MW(EllipSysOneTypeWT):
    def __init__(self):
        EllipSysOneTypeWT.__init__(self,
        name='DTU10MW',
        diameter=178.3,
        hub_height=119.0,
        cutin=4.0,
        cutout=25.0,
        dws=1.0, # what is this?
        tilt=0.0,
        rated_power=10000.0, # 10e6 [W], but this is in [kW] due to our choice of `power_unit`
        rotdir='cw',
        power_unit='kW',
        ct_func=self.dtu10mw_ct,
        power_func=self.dtu10mw_power,
        rpm_func=self.dtu10mw_rpm_curve,
        pitch_func=self.dtu10mw_pitch,
        bladeloading_func=self.dtu10mw_bladeloading,
        airfoildata_file='dtu10mw.pc',
        airfoildata_path=ADairfoil_path,
        airfoildata_file_type=2,  # 1=flex5, 2=hawc2
        bladegeo_file='dtu10mw.geo',
        bladegeo_path=ADairfoil_path
        ) # Below are re-definitions for ease of access
        self.D=178.3
        self.zRef=119.0
    
    def dtu10mw_ct(self, u):
        return np.interp(u, dtu10mw_ct_curve[:, 0], dtu10mw_ct_curve[:, 1])

    def dtu10mw_power(self, u):
        return np.interp(u, dtu10mw_power_curve[:, 0], dtu10mw_power_curve[:, 1])

    def dtu10mw_rpm_curve(self, u):
        return np.interp(u, dtu10mw_rpm_curve[:, 0], dtu10mw_rpm_curve[:, 1])
    
    def dtu10mw_pitch(self, u):
        return np.interp(u, dtu10mw_pitch_curve[:, 0], dtu10mw_pitch_curve[:, 1])
    
    def dtu10mw_bladeloading(self, r):
        Qx = np.interp(r, dtu10mw_bladeloading[:, 0], dtu10mw_bladeloading[:, 1])
        Qy = np.interp(r, dtu10mw_bladeloading[:, 0], dtu10mw_bladeloading[:, 2])
        Qz = np.interp(r, dtu10mw_bladeloading[:, 0], dtu10mw_bladeloading[:, 3])
        return Qx, Qy, Qz

from py_wake_ellipsys.examples.data.turbines.nrel5mw import nrel5mw_ct_curve
from py_wake_ellipsys.examples.data.turbines.nrel5mw import nrel5mw_power_curve
from py_wake_ellipsys.examples.data.turbines.nrel5mw import nrel5mw_rpm_curve

class NREL5MW(EllipSysOneTypeWT):
    def __init__(self):
        EllipSysOneTypeWT.__init__(self, name='NREL5MW',
        diameter=126.0,
        hub_height=90.0,
        cutin=4.0,
        cutout=25.0,
        dws=1.0,
        tilt=0.0,
        rated_power=5000.0,
        rotdir='cw',
        power_unit='kW',
        ct_func=self.nrel5mw_ct,
        power_func=self.nrel5mw_power,
        rpm_func=self.nrel5mw_rpm,
        # input on pitch function etc is not needed for non-3xxx forcings, but they are required positional arguments
        pitch_func='',
        bladeloading_func='',
        airfoildata_file='',
        airfoildata_path='',
        airfoildata_file_type='',  # 1=flex5, 2=hawc2
        bladegeo_file='',
        bladegeo_path=''
        ), # Below are re-definitions for ease of access
        self.D=126.0
        self.zRef=90.0

    def nrel5mw_ct(self, u):
        return np.interp(u, nrel5mw_ct_curve[:, 0], nrel5mw_ct_curve[:, 1])

    def nrel5mw_power(self, u):
        return np.interp(u, nrel5mw_power_curve[:, 0], nrel5mw_power_curve[:, 1])

    def nrel5mw_rpm(self, u):
        return np.interp(u, nrel5mw_rpm_curve[:, 0], nrel5mw_rpm_curve[:, 1])

# Grid resolution settings
# l: low, m: medium, f: fine (production), uf: ultra-fine (for grid study)

l_res_grid = {
    'grid_cells1_D': 2, # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 8, # Number of nodes in radial direction
    'adgrid_ntheta': 8, # Number of nodes in polar direction
    'ad_out': True, # Required for post_ad_bladeloads() method used later
    'ad_out_reduce_r': 1, # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_a': 1, # Reduction of azimuthal grid points for AD output (integer)
    }

m_res_grid = {
    'grid_cells1_D': 8, # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 64, # Number of nodes in radial direction
    'adgrid_ntheta': 64, # Number of nodes in polar direction
    'ad_out': True, # Required for post_ad_bladeloads() method used later
    'ad_out_reduce_r': 1, # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_a': 8, # Reduction of azimuthal grid points for AD output (integer)
    }

f_res_grid = {
    'grid_cells1_D': 16, # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 64, # Number of nodes in radial direction
    'adgrid_ntheta': 64, # Number of nodes in polar direction
    'ad_out': True, # Required for post_ad_bladeloads() method used later
    'ad_out_reduce_r': 1, # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_a': 8, # Reduction of azimuthal grid points for AD output (integer)
    }

uf_res_grid = {
    'grid_cells1_D': 20, # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 64, # Number of nodes in radial direction
    'adgrid_ntheta': 64, # Number of nodes in polar direction
    'ad_out': True, # Required for post_ad_bladeloads() method used later
    'ad_out_reduce_r': 1, # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_a': 8, # Reduction of azimuthal grid points for AD output (integer)
    }

# Wind farm settings
from py_wake.examples.data.hornsrev1 import Hornsrev1Site
from py_wake.site._site import UniformWeibullSite
   

class Single_Hornsrev1_wt(Hornsrev1Site):
    def __init__(self):
        Hornsrev1Site.__init__(self)
        self.Ti=0.06
        self.wt_x=[0]
        self.wt_y=[0]
        self.type_i=np.array([0])

class TwobyTwo_Hornsrev1_wf(Hornsrev1Site):
    def __init__(self, wt):
        Hornsrev1Site.__init__(self)
        self.Ti=0.06
        self.wt_x=[0, 5*wt.D, 0, 5*wt.D]
        self.wt_y=[0, 0, 5*wt.D, 5*wt.D]
        self.type_i=np.array([0, 0, 0, 0])

# Flow cases, [wind directions], [wind speeds]
flow_case_1 = [270.0], [8.0] # Normal scenario
flow_case_2 = [270.0], [11.0] # High thrust

# Code basement (currently unused)

# dtu10mw_data = { # Might come in handy for expansion to more WTs
#     'name': "DTU10MW",
#     'diameter': 178.3, # Rotor diameter [m]
#     'hub_height': 119.0, # [m]
#     'zRef': 119.0, # Reference hub height [m]
#     'cutin': 4.0, # Cut-in wind speed [m]
#     'cutout': 25.0, # Cut-out wind speed [m]
#     'rated_power': 10e6, # [W]'
#     }

# # Ti calculation - check later, reference - given by Paul in Teams chat on 3/2/22.
# zRef = hub_height # Reference height at which the turbulence is defined [m]
# z0 = 2e-04 # Roughness length, open sea [m]
# kappa = 0.4 # von Karman constant
# cmu = 0.03 # ???
# Ti = kappa*np.sqrt(2.0/3.0) / (cmu**0.25 * np.log((zRef + z0) / z0)) # Turbulence intensity [-]
# print(f"For this simulation, turbulence intensity Ti = {Ti}.") # ~0.059 for DTU10MW in offshore conditions

# # Check all EllipSys input variables
# variables = vars(e3d_model)
# print("\n EllipSys variables \n")
# [print(var, variables[var]) for var in variables]
# # OR: check only grid_* variables
# print("\n Grid variables \n")
# [print(var, variables[var]) for var in variables if var[0:4] == 'grid']