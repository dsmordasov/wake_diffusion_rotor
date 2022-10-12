"""
This file was made for increased ease of automating PyWakeEllipSys simulations,
HAWC2S work, post-processing and more. Keep in folder with post-processing 
code locally, or if running simulations on gbar, drop in the same folder as the
job script. 
"""
from sys import platform  # Check if running locally, or on gbar
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler

import os

# %% Plot style, generously donated by Mads Christian Baungaard
mpl.style.use('classic')

# Use below rcParams if you have many lines to plot
# mpl.rcParams['axes.prop_cycle'] = cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', 'k'])

plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams["legend.scatterpoints"] = 1
plt.rcParams["legend.numpoints"] = 1
plt.rcParams['grid.linestyle'] = ':'  # Dotted gridlines
mpl.rcParams['lines.linewidth'] = 2
plt.rcParams['axes.grid'] = True
# This is the best way to customize the reversed velocity deficit plots, trust me
yd = dict(rotation=0, ha='right')
plt.close('all')

# %% To enable post-processing locally, not on gbar, we can not import PWES files
# This block imports PWES if running on gbar, if ran locally on a Windows
# distribution for post-processing, then it does not import PWES

if not platform == 'win32':
    # Wind turbine imports
    from py_wake_ellipsys.wind_farm_models.ellipsys_lib.ellipsys_wind_turbines import EllipSysOneTypeWT
    from py_wake_ellipsys.examples.data.turbines.ADairfoil import ADairfoil_path

    # DTU10MW
    from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_ct_curve
    from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_power_curve
    from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_rpm_curve
    from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_pitch_curve
    from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_bladeloading

    # NREL5MW
    from py_wake_ellipsys.examples.data.turbines.nrel5mw import nrel5mw_ct_curve
    from py_wake_ellipsys.examples.data.turbines.nrel5mw import nrel5mw_power_curve
    from py_wake_ellipsys.examples.data.turbines.nrel5mw import nrel5mw_rpm_curve

    # Sites
    from py_wake.examples.data.hornsrev1 import Hornsrev1Site

    # %% Wind turbine definition
    class DTU10MW(EllipSysOneTypeWT):
        def __init__(self):
            EllipSysOneTypeWT.__init__(self,
                                       name='DTU10MW',
                                       diameter=178.3,
                                       hub_height=119.0,
                                       cutin=4.0,
                                       cutout=25.0,
                                       dws=1.0,
                                       tilt=0.0,
                                       # 10e6 [W], but this is in [kW] due to our choice of `power_unit`
                                       rated_power=10000.0,
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
                                       bladegeo_path=ADairfoil_path,
                                       )  # Below are re-definitions for ease of access
            self.D = 178.3
            self.zRef = 119.0

        def dtu10mw_ct(self, u):
            return np.interp(u, dtu10mw_ct_curve[:, 0], dtu10mw_ct_curve[:, 1])

        def dtu10mw_power(self, u):
            return np.interp(u, dtu10mw_power_curve[:, 0], dtu10mw_power_curve[:, 1])

        def dtu10mw_rpm_curve(self, u):
            return np.interp(u, dtu10mw_rpm_curve[:, 0], dtu10mw_rpm_curve[:, 1])

        def dtu10mw_pitch(self, u):
            return np.interp(u, dtu10mw_pitch_curve[:, 0], dtu10mw_pitch_curve[:, 1])

        def dtu10mw_bladeloading(self, r):
            Qx = np.interp(
                r, dtu10mw_bladeloading[:, 0], dtu10mw_bladeloading[:, 1])
            Qy = np.interp(
                r, dtu10mw_bladeloading[:, 0], dtu10mw_bladeloading[:, 2])
            Qz = np.interp(
                r, dtu10mw_bladeloading[:, 0], dtu10mw_bladeloading[:, 3])
            return Qx, Qy, Qz

    # Wind farm definions
    class Single_Hornsrev1_wt(Hornsrev1Site):  # Single wind turbine
        def __init__(self):
            Hornsrev1Site.__init__(self)
            self.Ti = 0.06
            self.wt_x = [0]
            self.wt_y = [0]
            self.type_i = np.array([0])

    class Row3_Hornsrev1_wf(Hornsrev1Site):  # Three wind turbine row
        def __init__(self, wt):
            Hornsrev1Site.__init__(self)
            self.Ti = 0.06
            self.wt_x = np.arange(0, 15, 5)*wt.D
            self.wt_y = np.zeros(3)
            self.type_i = np.zeros(3)

    class Row10_Hornsrev1_wf(Hornsrev1Site):  # Ten wind turbine row
        def __init__(self, wt):
            Hornsrev1Site.__init__(self)
            self.Ti = 0.06
            self.wt_x = np.arange(0, 50, 5)*wt.D
            self.wt_y = np.zeros(10)
            self.type_i = np.zeros(10)


# Grid resolution settings
# l: low, m: medium, f: fine (production), uf: ultra-fine (for grid study)

test_res_grid = {
    'grid_cells1_D': 2,  # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 4,  # Number of nodes in radial direction
    'adgrid_ntheta': 4,  # Number of nodes in polar direction
    'ad_out': True,  # Required for post_ad_bladeloads() method used later
    # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_r': 1,
    # Reduction of azimuthal grid points for AD output (integer)
    'ad_out_reduce_a': 1,
}


l_res_grid = {
    'grid_cells1_D': 2,  # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 8,  # Number of nodes in radial direction
    'adgrid_ntheta': 8,  # Number of nodes in polar direction
    'ad_out': True,  # Required for post_ad_bladeloads() method used later
    # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_r': 1,
    # Reduction of azimuthal grid points for AD output (integer)
    'ad_out_reduce_a': 1,
}

m_res_grid = {
    'grid_cells1_D': 8,  # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 64,  # Number of nodes in radial direction
    'adgrid_ntheta': 64,  # Number of nodes in polar direction
    'ad_out': True,  # Required for post_ad_bladeloads() method used later
    # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_r': 1,
    # Reduction of azimuthal grid points for AD output (integer)
    'ad_out_reduce_a': 8,
}

# --- PRODUCTION GRID ---
f_res_grid = {
    'grid_cells1_D': 16,  # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 64,  # Number of nodes in radial direction
    'adgrid_ntheta': 64,  # Number of nodes in polar direction
    'ad_out': True,  # Required for post_ad_bladeloads() method used later
    # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_r': 1,
    # Reduction of azimuthal grid points for AD output (integer)
    'ad_out_reduce_a': 8,
}

uf_res_grid = {
    'grid_cells1_D': 20,  # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 64,  # Number of nodes in radial direction
    'adgrid_ntheta': 64,  # Number of nodes in polar direction
    'ad_out': True,  # Required for post_ad_bladeloads() method used later
    # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_r': 1,
    # Reduction of azimuthal grid points for AD output (integer)
    'ad_out_reduce_a': 8,
}

# --- GRID STUDY GRIDS ---
gs_5_res_grid = {
    'grid_cells1_D': 16,  # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 128,  # Number of nodes in radial direction
    'adgrid_ntheta': 128,  # Number of nodes in polar direction
    'ad_out': True,  # Required for post_ad_bladeloads() method used later
    # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_r': 1,
    # Reduction of azimuthal grid points for AD output (integer)
    'ad_out_reduce_a': 8,
}

gs_4_res_grid = {
    'grid_cells1_D': 16,  # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 64,  # Number of nodes in radial direction
    'adgrid_ntheta': 64,  # Number of nodes in polar direction
    'ad_out': True,  # Required for post_ad_bladeloads() method used later
    # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_r': 1,
    # Reduction of azimuthal grid points for AD output (integer)
    'ad_out_reduce_a': 8,
}

gs_3_res_grid = {
    'grid_cells1_D': 16,  # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 32,  # Number of nodes in radial direction
    'adgrid_ntheta': 32,  # Number of nodes in polar direction
    'ad_out': True,  # Required for post_ad_bladeloads() method used later
    # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_r': 1,
    # Reduction of azimuthal grid points for AD output (integer)
    'ad_out_reduce_a': 8,
}

gs_2_res_grid = {
    'grid_cells1_D': 16,  # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 16,  # Number of nodes in radial direction
    'adgrid_ntheta': 16,  # Number of nodes in polar direction
    'ad_out': True,  # Required for post_ad_bladeloads() method used later
    # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_r': 1,
    # Reduction of azimuthal grid points for AD output (integer)
    'ad_out_reduce_a': 8,
}

gs_1_res_grid = {
    'grid_cells1_D': 16,  # Cells per rotor diameter in the inner domain [D]
    'adgrid_nr': 8,  # Number of nodes in radial direction
    'adgrid_ntheta': 8,  # Number of nodes in polar direction
    'ad_out': True,  # Required for post_ad_bladeloads() method used later
    # Reduction of radial grid points for AD output (integer)
    'ad_out_reduce_r': 1,
    # Reduction of azimuthal grid points for AD output (integer)
    'ad_out_reduce_a': 8,
}


# Flow cases: [wind directions], [wind speeds]
flow_case_1 = [270.0], [8.0]  # High thrust scenario
flow_case_2 = [270.0], [11.0]  # Near rated scenario

# Code basement (uncomment to use)

# # Turbulence intensity calculation
# zRef = hub_height # Reference height at which the turbulence is defined [m]
# z0 = 2e-04 # Roughness length, open sea [m]
# kappa = 0.4 # von Karman constant [-]
# cmu = 0.03 # turbulence model constant, calibrated for neutral atmospheric conditiosn [-]
# TI = kappa*np.sqrt(2.0/3.0) / (cmu**0.25 * np.log((zRef + z0) / z0)) # Turbulence intensity [-]
# print(f"For this simulation, turbulence intensity TI = {TI}.") # ~0.059 for DTU10MW in offshore conditions

# # Check all EllipSys input variables
# variables = vars(e3d_model)
# print("\n EllipSys variables \n")
# [print(var, variables[var]) for var in variables]
# # OR: check only grid_* variables
# print("\n Grid variables \n")
# [print(var, variables[var]) for var in variables if var[0:4] == 'grid']
