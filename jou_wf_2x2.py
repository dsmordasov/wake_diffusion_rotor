"""
This code was made for the purposes of studying the effect of varying 
\overline{\delta} in the Joukowsky `2xxx` forcing model of PyWakeEllipSys.
"""
# Imports
import numpy as np
import time 
from py_wake_ellipsys.wind_farm_models.ellipsys import EllipSys 

import config



wt = config.NREL5MW()

#%% Wind farm definition
from py_wake.site._site import UniformWeibullSite
# Below from https://gitlab.windenergy.dtu.dk/TOPFARM/cuttingedge/pywake/pywake_ellipsys/-/blob/master/py_wake_ellipsys/examples/Multirotor/mr.py#L185
class UniformSite(UniformWeibullSite):
    def __init__(self, A, k, nsector=12.0, shear=None): # TODO:  Add shear!
        UniformWeibullSite.__init__(self,
                                    np.ones(int(nsector)) / float(nsector), # Sector frequencies [-]
                                    [A] * int(nsector), # Weibull scale parameter [m/s]
                                    [k] * int(nsector), # Weibull shape parameter
                                    0.06, # Turbulence intensity
                                    shear=shear)
wf = UniformSite(A=9.0, k=2.2) # Horns Rev 2 values, apparently
wf.Ti = 0.06 # General offshore
wf.wt_x=[0, 5*wt.D, 0, 5*wt.D]
wf.wt_y=[0, 0, 5*wt.D, 5*wt.D]
wf.type_i=np.array([0, 0, 0, 0])

#%% Run parameters
create_grid = True
run = True
aep = True

delta_variants = np.array([0.10, 0.25, 0.4])

wd = np.arange(0.0, 360.0, 15.0)
ws = np.arange(4.0, 25.1, 4.0)
wdsym=[0.0, 90.0, 'rsym']
wscutoff = 17 # No wake effects are assumed for U = 17 m/s and above to reduce number of flow cases

# %%
def run_jou_wf_2x2(ad_jou_rootdelta=0.25,
                        ad_jou_roota=2.335,
                        ad_jou_rootb=4.0,
                        grid=config.l_res_grid, # Grid resolution [l, m, f, uf], see config.py
                        e3d_reslim=1e-2, # Convergence criteria
                        wt=wt, # Wind turbine
                        wf=wf, # Wind farm
                        ):
    start_time = time.time()
    
    # Define the wind farm simulation
    e3d_model = EllipSys(wf, wt, wf.Ti, wt.zRef, # Required positional arguments
        ad_force='2111', # Forcing method
        grid_cells1_D=grid['grid_cells1_D'], # Cells per rotor diameter in the inner domain [D], from 2 to 8
        grid_zFirstCell_D=0.5/wt.D, # Height of the first cell [D]
        adgrid_nr=grid['adgrid_nr'], # Number of nodes in radial direction
        adgrid_ntheta=grid['adgrid_ntheta'], # Number of nodes in polar direction
        ad_jou_rootdelta=ad_jou_rootdelta, # Root correction constant \overline{delta} for AD Joukowsky
        ad_jou_roota=ad_jou_roota, # Root correction constant `a` for AD Joukowsky
        ad_jou_rootb=ad_jou_rootb, # Root correction constant `b` for AD Joukowsky
        casename='nrel5mw', #
        #grid_m1_w_D=3.0, # Inner margin, west, after rotation [D]
        #grid_m1_e_D=12.0, # Inner margin, east, after rotation [D]
        #grid_m1_s_D=1.5, # Inner margin, south, after rotation [D]
        #grid_m1_n_D=1.5, # Inner margin, north, after rotation [D]
        e3d_reslim=e3d_reslim, # Convergence criteria of EllipSys3D flow solver
        grid_radius_D=5e4/wt.D, # Outer distance in the grid [D]
        e3d_runwf_wdsym='auto', # Wind farm symmetry, set to 'auto'?
        e3d_runwf_wscutoff='auto', 
        e3d_runwf_march_dir='ws', # Ask if faster than default
        run_machine='gbar',
        maxnodes=4,
        corespernode=24,
        queue='hpc',
        e3d_turbmodel='kefP',
        run_write_restart=True, # Write restart file per flow case
        walltime_cal_grid='0:30:00',
        walltime_wf_grid='0:30:00',
        walltime_cal_run='0:30:00',
        walltime_wf_run='2:00:00',
        #ad_out=True, # Required for post_ad_bladeloads() method used later
        #ad_out_reduce_r=grid['ad_out_reduce_r'], # Reduction of radial grid points for AD output (integer)
        #ad_out_reduce_a=grid['ad_out_reduce_a'], # Reduction of azimuthal grid points for AD output (integer)
        )

    # ---Run step by step---

    if create_grid:
        # Create AD polar grid
        e3d_model.create_adgrid(wf.type_i[0])
        # Create calibration grid
        e3d_model.create_calibration_grid(wf.type_i[0])
        # Create wind farm grid
        e3d_model.create_windfarm_grid(wf.wt_x, wf.wt_y, wf.type_i)
        print("create_windfarm_grid() started")
    
    if run:
        # Change delta
        for delta in delta_variants:
            e3d_model.ad_jou_rootdelta = delta
            subcasename = f"_d{round(delta, 2)}"
            e3d_model.subcasename=subcasename
            
            # Run force control calibration
            e3d_model.run_calibration(wf.type_i[0])
            
            # Run flow case
            e3d_model.run_windfarm(wf.wt_x, wf.wt_y, wd, ws, wf.type_i)
            
            # Post process flow case
            WS_eff_ilk, power_ilk, ct_ilk = e3d_model.post_windfarm(wf.type_i, wd, ws) # what does ilk stand for?
            print(f"\n Power:{(np.sum(power_ilk))/1e6} MW \n")
            print(f"\n Power:{power_ilk/1e6} MW \n")
            
    end_time = time.time()
    print(f"Simulation ran in {round((end_time-start_time)/60, 2)} minutes")

if __name__ == '__main__':
    run_jou_wf_2x2()