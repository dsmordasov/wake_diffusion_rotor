"""
This file was made in order to run simulations for Joukowsky and RANS 
DTU10MW reference wind turbine, in order to validate these methods against 
present DES data: U = 8.0 m/s, rho = 1.225 kg/m^3, RPM = 6.4259, TI = 0.00
"""

# Imports
import numpy as np
import time
import os
from py_wake_ellipsys.wind_farm_models.ellipsys import EllipSys

import config

run_grid = False
design_name = 'new_design'


def run_pwes_simulation(
        ad_jou_rootdelta=0.35,
        ad_jou_roota=2.335,
        ad_jou_rootb=4.0,
        grid=config.f_res_grid,  # Grid resolution [l, m, f], see config.py
        e3d_reslim=1e-5,  # Convergence criteria
        wt=config.DTU10MW(),  # Wind turbine [DTU10MW()]
        wf=config.Single_Hornsrev1_wt(),  # Wind farm
):

    start_time = time.time()

    wf.Ti = 1e-2  # Set the turbulence intensity as close to 0 as in the DES simulation
    wd = [270]  # [deg]
    ws = [8]  # [m/s]

    if design_name:
        wt.bladegeo_file = f"{design_name}_blade.geo"
        wt.bladegeo_path = os.getcwd()
        print(f"bladegeo_path: {wt.bladegeo_path}")

    e3d_model = EllipSys(wf, wt, wf.Ti, wt.zRef,  # Required positional arguments
                         ad_force='3110',  # Forcing method
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
                         walltime_wf_run='1:00:00',
                         ad_out=True,  # Required for post_ad_bladeloads() method used later
                         # Reduction of radial grid points for AD output (integer)
                         ad_out_reduce_r=grid['ad_out_reduce_r'],
                         # Reduction of azimuthal grid points for AD output (integer)
                         ad_out_reduce_a=grid['ad_out_reduce_a'],
                         )

    # Create wind farm grid
    if run_grid:
        e3d_model.create_adgrid(wf.type_i[0])
        print("create_adgrid() started")
        e3d_model.create_windfarm_grid(wf.wt_x, wf.wt_y, wf.type_i[0])
        print("create_windfarm_grid() started")

    e3d_model.run_windfarm(wf.wt_x, wf.wt_y, wd, ws, wf.type_i)
    print("run_windfarm() started")
    WS_eff_ilk, power_ilk, ct_ilk = e3d_model.post_windfarm(
        wf.type_i, wd, ws)  # what does ilk stand for?
    print("post_windfarm() started")

    print(f"\n Power: {float(power_ilk)/1e6} MW \n")
    end_time = time.time()
    runtime = round(end_time - start_time, 2)
    print(f"Simulation ran in {runtime}s.")

    # Store 3D flow data as a netCDF file for the one flow case
    # Clip the netCDF data
    # [-x < xyz[0] < x, -y < xyz[1] < y, xyz[2] < z]
    xyz = [10 * wt.D, 2 * wt.D, wt.zRef + wt.D]
    # On gbar/hpc, post_windfarm_flow() with more than one nodes sometimes fails
    e3d_model.maxnodes = 1
    e3d_model.post_windfarm_flow(
        wd, ws, outputformat='netCDF', clip_netCDF=xyz)
    print("post_windfarm_flow() started")
    # Calculate the averaged loads of the first AD, output results to <>adbladeloads.dat
    e3d_model.post_ad_bladeloads(wd, ws, np.array(wf.type_i))
    print("post_ad_bladeloads() started")

    return np.array([float(power_ilk), runtime])


if __name__ == '__main__':
    run_pwes_simulation()
