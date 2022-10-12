"""
This file was made in order to run simulations a ten wind turbine row 
simulation using the Joukowsky `2xxx` forcing model of PyWakeEllipSys, with a
user-chosen \overline{\delta} parameter value.
"""

# Imports
from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_bladeloading
from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_pitch_curve
from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_rpm_curve
from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_power_curve
from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_ct_curve
from py_wake_ellipsys.wind_farm_models.ellipsys_lib.ellipsys_wind_turbines import EllipSysOneTypeWT
import numpy as np
import time
import os
from py_wake_ellipsys.wind_farm_models.ellipsys import EllipSys

import config

# Force method 2112 does not need force calibration simulations.
create_ad_grid = True
create_cal_grid = True
create_wf_grid = True
run_calibration = True

design_name = 'jou10_0.35_1'

baseline_power = 7.293024  # MW, flattened 3-wt row PWES at TI=0.06

# DTU10MW
class DTU10MW(EllipSysOneTypeWT):
    def __init__(self):
        EllipSysOneTypeWT.__init__(self,
                                   name='DTU10MW',
                                   diameter=178.3,
                                   hub_height=119.0,
                                   cutin=4.0,
                                   cutout=25.0,
                                   dws=1.0,  # what is this?
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
                                   airfoildata_path=config.ADairfoil_path,
                                   airfoildata_file_type=2,  # 1=flex5, 2=hawc2
                                   bladegeo_file='baseline_blade.geo',
                                   bladegeo_path=f'{os.getcwd()}/',
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


def run_pwes_simulation(
        ad_jou_rootdelta=0.35,
        ad_jou_roota=2.335,
        ad_jou_rootb=4.0,
        grid=config.f_res_grid,  # Grid resolution [l, m, f], see config.py
        e3d_reslim=1e-5,  # Convergence criteria
        wf=config.Row10_Hornsrev1_wf(DTU10MW()),  # Wind farm
):

    start_time = time.time()

    wt = DTU10MW()
    wd = [270]  # [deg]
    ws = [8]  # [m/s]

    e3d_model = EllipSys(wf, wt, wf.Ti, wt.zRef,  # Required positional arguments
                         subcasename=f"_{design_name}",
                         ad_force='2111',  # Forcing method
                         # Cells per rotor diameter in the inner domain [D], from 2 to 8
                         grid_cells1_D=grid['grid_cells1_D'],
                         # Height of the first cell [D]
                         grid_zFirstCell_D=0.5/wt.D,
                         # Number of nodes in radial direction
                         adgrid_nr=grid['adgrid_nr'],
                         # Number of nodes in polar direction
                         adgrid_ntheta=grid['adgrid_ntheta'],
                         # Root correction constant \overline{delta} for AD Joukowsky
                         ad_jou_rootdelta=ad_jou_rootdelta,
                         ad_jou_roota=ad_jou_roota,  # Root correction constant `a` for AD Joukowsky
                         ad_jou_rootb=ad_jou_rootb,  # Root correction constant `b` for AD Joukowsky
                         casename='dtu10mw',
                         # 3 inner margin, west, after rotation [D]
                         grid_m1_w_D=3.0,
                         # 12 inner margin, east, after rotation [D]
                         grid_m1_e_D=12.0,
                         grid_m1_s_D=1.5,  # 1.5
                         grid_m1_n_D=1.5,  # 1.5
                         # Convergence criteria of EllipSys3D flow solver, 1e-5 is normal, 1e-4 is ok, 1e-2 is ridiculous
                         e3d_reslim=e3d_reslim,
                         grid_radius_D=100,  # Outer distance in the grid [D]
                         run_machine='gbar',
                         maxnodes=4,
                         corespernode=24,
                         queue='hpc',
                         e3d_turbmodel='kefP',
                         run_write_restart=True,  # Write restart file per flow case
                         walltime_wf_grid='1:00:00',
                         walltime_wf_run='4:00:00',
                         walltime_cal_run='4:00:00',  # Important! Default 10 mins too small
                         ad_out=True,  # Required for post_ad_bladeloads() method used later
                         # Reduction of radial grid points for AD output (integer)
                         ad_out_reduce_r=grid['ad_out_reduce_r'],
                         # Reduction of azimuthal grid points for AD output (integer)
                         ad_out_reduce_a=grid['ad_out_reduce_a'],
                         )

    # Create wind farm grid
    if create_ad_grid:
        e3d_model.create_adgrid(wf.type_i[0])
        print("create_adgrid() started")

    if create_cal_grid:
        e3d_model.create_calibration_grid(wf.type_i[0])
        print("create_calibration_grid started")

    if create_wf_grid:
        e3d_model.create_windfarm_grid(wf.wt_x, wf.wt_y, wf.type_i)
        print("create_windfarm_grid() started")

    if run_calibration:
        e3d_model.run_calibration(wf.type_i[0])
        print("run_calibration() started")

    e3d_model.run_windfarm(wf.wt_x, wf.wt_y, wd, ws, wf.type_i)
    print("run_windfarm() started")
    WS_eff_ilk, power_ilk, ct_ilk = e3d_model.post_windfarm(
        wf.type_i, wd, ws)  # what does ilk stand for?
    print("post_windfarm() started")

    fname = f"{design_name}_row3_power.csv"
    np.savetxt(fname, power_ilk, delimiter="\t")

    new_power = np.sum(power_ilk)/1e6  # [MW]
    percentage_difference = (new_power - baseline_power) / baseline_power * 100
    print(f"{'Baseline P: ':<20}" + f"{str(np.round(baseline_power, 4))} MW.")
    print(f"{'{design_name} P: ':<20}" + f"{str(np.round(new_power, 4))} MW.")
    print(f"{'Delta P: ':<19}" + f"{str(np.round(percentage_difference, 2))} %.")

    end_time = time.time()
    runtime = round(end_time - start_time, 2)
    print(f"Simulation ran in {runtime}s.")

    # Store 3D flow data as a netCDF file for the one flow case
    # Clip the netCDF data
    # [-x < xyz[0] < x, -y < xyz[1] < y, xyz[2] < z]
    xyz_3 = [55 * wt.D, 2 * wt.D, wt.zRef + wt.D]
    # On gbar/hpc, post_windfarm_flow() with more than one nodes sometimes fails
    e3d_model.maxnodes = 1
    e3d_model.post_windfarm_flow(
        wd, ws, outputformat='netCDF', clip_netCDF=xyz_3)
    print("post_windfarm_flow() started")
    # Calculate the averaged loads of the first AD, output results to <>adbladeloads.dat
    e3d_model.post_ad_bladeloads(wd, ws, np.array(wf.type_i))
    print("post_ad_bladeloads() started")

    # Post-process the ADVD.

    # pps.pp_advd(design_name=design_name)

    return  # np.array(power_ilk, runtime)


if __name__ == '__main__':
    run_pwes_simulation()
