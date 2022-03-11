"""
This code was made for the purposes of studying the effect of varying 
\overline{\delta} in the Joukowsky `2xxx` forcing model of PyWakeEllipSys.
"""
# Imports
import numpy as np
import time 
from py_wake_ellipsys.wind_farm_models.ellipsys import EllipSys 

import config

new_grid = True
run = True
delta_variants = np.arange(0.15, 0.55, 0.05)

def run_delta_study(ad_jou_rootdelta=0.25,
                        ad_jou_roota=2.335,
                        ad_jou_rootb=4.0,
                        grid=config.f_res_grid, # Grid resolution [l, m, f, uf], see config.py
                        e3d_reslim=1e-5, # Convergence criteria
                        wt=config.DTU10MW(), # Wind turbine [DTU10MW()]
                        wf=config.Single_Hornsrev1_wt(), # Wind farm 
                        fc=config.flow_case_1,): # Flow case
    start_time = time.time()
    
    # Define the wind farm simulation
    e3d_model = EllipSys(wf, wt, wf.Ti, wt.zRef, # Required positional arguments
        ad_force='2110', # Forcing method
        grid_cells1_D=grid['grid_cells1_D'], # Cells per rotor diameter in the inner domain [D], from 2 to 8
        grid_zFirstCell_D=0.5/wt.D, # Height of the first cell [D]
        adgrid_nr=grid['adgrid_nr'], # Number of nodes in radial direction
        adgrid_ntheta=grid['adgrid_ntheta'], # Number of nodes in polar direction
        ad_jou_rootdelta=ad_jou_rootdelta, # Root correction constant \overline{delta} for AD Joukowsky
        ad_jou_roota=ad_jou_roota, # Root correction constant `a` for AD Joukowsky
        ad_jou_rootb=ad_jou_rootb, # Root correction constant `b` for AD Joukowsky
        casename='dtu10mw', #
        grid_m1_w_D=3.0, # Inner margin, west, after rotation [D]
        grid_m1_e_D=12.0, # Inner margin, east, after rotation [D]
        grid_m1_s_D=1.5, # Inner margin, south, after rotation [D]
        grid_m1_n_D=1.5, # Inner margin, north, after rotation [D]
        e3d_reslim=e3d_reslim, # Convergence criteria of EllipSys3D flow solver
        grid_radius_D=100, # Outer distance in the grid [D]
        run_machine='gbar',
        maxnodes=4,
        corespernode=24,
        queue='hpc',
        e3d_turbmodel='kefP',
        run_write_restart=True, # Write restart file per flow case
        walltime_wf_run='1:30:00',
        ad_out=True, # Required for post_ad_bladeloads() method used later
        ad_out_reduce_r=grid['ad_out_reduce_r'], # Reduction of radial grid points for AD output (integer)
        ad_out_reduce_a=grid['ad_out_reduce_a'], # Reduction of azimuthal grid points for AD output (integer)
        )

    # ---Run step by step---

    if new_grid:
        # Create AD polar grid
        e3d_model.create_adgrid(wf.type_i[0])
        # Create wind farm grid
        e3d_model.create_windfarm_grid(wf.wt_x, wf.wt_y, wf.type_i[0])
        print("create_windfarm_grid() started")
        
    if run:
        # Run flow cases
        for delta in delta_variants:
            e3d_model.maxnodes = 4
            e3d_model.ad_jou_rootdelta = delta
            subcasename = f"_d{round(delta, 2)}"
            e3d_model.subcasename=subcasename
            wd, ws = fc # fc is a tuple in this script atm, setting wind direction & speed
            e3d_model.run_windfarm(wf.wt_x, wf.wt_y, wd, ws, wf.type_i)
            
            # Post process flow cases
            WS_eff_ilk, power_ilk, ct_ilk = e3d_model.post_windfarm(wf.type_i, wd, ws) # what does ilk stand for?
            print(f"\n Power: {float(power_ilk)/1e6} MW \n")
            end_time = time.time()
            runtime = round(end_time - start_time, 2)
            print(f"Simulation ran in {runtime}s.")
            # Store 3D flow data as a netCDF file for the one flow case
            # Clip the netCDF data 
            xyz = [10 * wt.D, 2 * wt.D, wt.zRef + wt.D] # [-x < xyz[0] < x, -y < xyz[1] < y, xyz[2] < z]
            e3d_model.maxnodes = 1 # On gbar/hpc, post_windfarm_flow() with more than one nodes sometimes fails
            e3d_model.post_windfarm_flow(wd, ws, outputformat='netCDF', clip_netCDF=xyz)
            # Calculate the averaged loads of the first AD, output results to <>adbladeloads.dat
            e3d_model.post_ad_bladeloads(wd, ws, np.array(wf.type_i))

if __name__ == '__main__':
    run_delta_study()
    