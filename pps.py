""" Post-processing suite/toolbox, single wind turbine.
"""

import glob
import time
import sys

import numpy as np
import matplotlib.pyplot as plt
import xarray

# Run parameters, DTU10MW
rho = 1.225 # Air density [kg/m^3]
UH = 8 # Simulation velocity [m/s]
D = 178.3 # Rotor diameter [m]
zh = 119 # Reference hub height [m]

advd_x_probe_D = 5 # Distance x at which we probe for velocity deficit [D]

# Baseline design parameters

def calc_UAD(data, x, y, z, D, n=128, var='U'):
    '''
    Calculate <U_AD> at each position.
    A n by n rectangular grid is masked by sqrt(y^2 + z^2) < R.
    I recommend using n=128 (can increase to minimize geometric error).
    Same resolution and diameter used for all disks.
    Important that x, y and z have the same length!
    Input:
        - data: netcdf data       (xarray Dataset)
        - x: streamwise positions (array)
        - y: lateral positions    (array)
        - z: vertical positions   (array)
        - n: resolution           (scalar)
        - D: diameter             (scalar)
        
    Written by Mads Christian Baungaard, taken from plotnetCDFad.py
    '''
    R = D / 2
    results = np.zeros(len(x))

    for i in range(len(x)):
        # Extract rectangular data
        y_rec = np.linspace(-R, R, n) + y[i]
        z_rec = np.linspace(-R, R, n) + z[i]
        A = (y_rec[1] - y_rec[0]) * (z_rec[1] - z_rec[0])    # Rectangular shape element
        rec_data = data.interp(x=x[i], y=y_rec, z=z_rec)  # This is fast

        # Mask operation. If mask condition not fulfilled, set value to 0.0.
        rec_data_masked = rec_data.where(np.sqrt((rec_data.y - y[i])**2 + (rec_data.z - z[i])**2) < R,
                                         other=0.0)

        # Each velocity contribution is weighted by the shape element area
        results[i] = np.sum(rec_data_masked[var].values * A) / (np.pi * R**2)

    return results

def pp_advd(design_name="new_design"):
    """
    Post process actuator disk velocity deficit - data extraction and plotting.
    Returns the velocity deficit at distance `advd_x_probe_D` defined earlier
    in the script, x=5D is a common value in research literature.
    Also from plotnetCDFad.py by Mads Christian Baungaard.
    """
    baseline_path = "/work3/s203350/pwes_tests/baseline_data/flowdata.nc"
    
    glob_expression = f"*{design_name}*/*/flowdata.nc"
    flowdata_path = glob.glob(glob_expression, recursive=True)
    print("Waiting for post_wf to finish:")
    while not flowdata_path:
        flowdata_path = glob.glob(glob_expression, recursive=True)
        time.sleep(5)
        sys.stdout.write('.', )
        sys.stdout.flush()
    flowdata_path = flowdata_path[0] # Otherwise it's a list, and doesn't open
    print(f"Flowdata path:\n{flowdata_path}")
    
    x = np.arange(-7, 7.01, 0.25) * D
    y = np.zeros(x.shape)
    z = np.ones(x.shape) * zh
    
    baseline_data = xarray.open_dataset(baseline_path)
    baseline_ADvar = calc_UAD(baseline_data, x, y, z, D)
    baseline_advd = np.interp(advd_x_probe_D, x / D, baseline_ADvar / UH)

    # Calculate disk-averaged quantities
    data = xarray.open_dataset(flowdata_path)  
    ADvar = calc_UAD(data, x, y, z, D)
    probed_advd = np.interp(advd_x_probe_D, x / D, ADvar / UH)
    
    fig, ax = plt.subplots(figsize=[10, 6])
    
    percentage_difference = (probed_advd - baseline_advd) / baseline_advd * 100
    print(f"{'Baseline ADVD: ':<20}" + f"{str(np.round(baseline_advd, 4))} [-].")
    print(f"{'{design_name} ADVD: ':<20}" + f"{str(np.round(probed_advd, 4))} [-].")
    print(f"{'Delta ADVD: ':<19}" + f"{str(np.round(percentage_difference, 2))} %.")
    
    ax.plot(x / D, baseline_ADvar / UH, label='Baseline')
    ax.plot(x / D, ADvar / UH, label=f'{design_name}')
    
    yd = dict(rotation=0, ha='right')     # Shortcut to flip y-label, credit to Mads
    ax.set_ylabel(r'$\dfrac{ \langle U_{AD} \rangle}{U_{H}}$ [-]', yd)
    ax.set_xlabel('$x/D$')
    ax.grid()
    ax.legend(loc='upper right')
    
    fig.savefig('AD_deficits.pdf', bbox_inches='tight')
    
    return design_name, probed_advd

# for ad blade loads
# 