"""
This code was made in order to post process the results of `delta_study.py`,
focusing on the effects of varying \overline{\delta} and thus the force 
distribution over a wind turbine blade on the wake of a single wind turbine.
"""

import glob
import re

import numpy as np
import matplotlib.pyplot as plt
import xarray

# Run options
parameter_name = "delta" # Set to either of: ["grid", "delta"]
check_against_DES = False # Set to True to check against DTU10MW loading
pp_adbladeloads_option = True # Set to True to post process AD blade loads
pp_power_option = True # Set to True to post process power
pp_advd_option = True # Set to True to post process velocity deficits (takes a while)
advd_x_probe_D = 5 # Distance x at which we probe for velocity deficit [D]

# Run parameters, DTU10MW
rho = 1.225 # Air density [kg/m^3]
UH = 11 # Simulation velocity [m/s]
D = 178.3 # Rotor diameter [m]
zh = 119 # Reference hub height [m]

def numerical_sort_key(parsed_string): 
    """
    Helper key function to sort a list of pathnames based on studied 
    parameter's value. 
    """
    if parameter_name == "grid":
        regex_expression = r'\d\d?'
    elif parameter_name == "delta":
        regex_expression = r'0.\d\d?'
    return list(map(float, re.findall(regex_expression, parsed_string)))

def determine_parameter_value(filename_path):
    """
    Helper function to grab parameter value using regex.
    Change expression based on required parameter.
    """
    if parameter_name == "grid":
        regex_magic = re.search('(\d\d?cD){1}', filename_path)
        if regex_magic:
            grid_fineness = re.search('(\d\d?){1}', regex_magic.group(1))
            parameter_value = np.int(grid_fineness[0])
    elif parameter_name == "delta":
        regex_magic = re.search('(0.\d\d?){1}', filename_path)
        if regex_magic:
            parameter_value = round(np.float(regex_magic.group(1)), 2)
    if not regex_magic:
        print("Cannot discern the parameter value from filename.\n")
        print("Check your `parameter_name` input or regex expression.")
        return
    return parameter_value

def plot_adbladeloads(filename_path):
    """
    Function to plot <>adbladeloads.dat normalised tangential and azimuthal loads
    """
    # Read file
    loads_data = np.loadtxt(filename_path)
    parameter_value = determine_parameter_value(filename_path)
    
    # Manipulate data
    r  = loads_data[:, 0] # Blade radius [m]
    ft = loads_data[:, 1] # Tangential force per unit area [N/m^2] (assumed)
    fn = loads_data[:, 2] # Azimuthal force per unit area [N/m^2] (assumed)
    Ft = loads_data[:, 3] # Tangential force per unit length [N/m] (assumed)
    Fn = loads_data[:, 4] # Azimuthal force per unit length [N/m] (assumed)
    ut = loads_data[:, 5] # Tangential velocity [m/s] (assumed)
    un = loads_data[:, 6] # Azumithal velocity [m/s] (assumed)
    
    # Normalisation of forces, non-dimensionalisation of radius
    R = r[-1] # Rotor max radius [m]
    r_nd = r/R # Non-dimensionalised blade radius r/R [-]
    Ft = Ft / (rho * R * UH**2)
    Fn = Fn / (rho * R * UH**2)

    ax[0].plot(r_nd, Ft, label=fr"${parameter_name}$={parameter_value}")
    ax[1].plot(r_nd, Fn, label=fr"${parameter_name}$={parameter_value}")
    return

def calc_UAD(data, x, y, z, D, n=32, var='U'):
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

def pp_advd(filename_path):
    """
    Post process actuator disk velocity deficit - data extraction and plotting.
    Returns the velocity deficit at distance `advd_x_probe_D` defined earlier
    in the script, x=5D is a common value in research literature.
    Also from plotnetCDFad.py by Mads Christian Baungaard.
    """
    
    x = np.arange(-7, 7.01, 0.25) * D
    y = np.zeros(x.shape)
    z = np.ones(x.shape) * zh
    data = xarray.open_dataset(filename_path)
    parameter_value = determine_parameter_value(filename_path)    
    
    # Calculate disk-averaged quantities
    ADvar = calc_UAD(data, x, y, z, D)
    probed_advd = np.interp(advd_x_probe_D, x / D, ADvar / UH)
    ax.plot(x / D, ADvar / UH, label=f'{parameter_name}={parameter_value}')
    return parameter_value, probed_advd

def pp_power(filename_path):
    """
    Post process power data of a single wind turbine.
    """
    # Read file
    power_data = np.loadtxt(filename_path, usecols=4)
    parameter_value = determine_parameter_value(filename_path)
    return parameter_value, power_data 
        
# Create lists of found pathnames
adbladeloads_paths = glob.glob("**/*_adbladeloads.dat", recursive=True)
flowdata_paths = glob.glob("**/flowdata.nc", recursive=True)
power_paths = glob.glob("**/Power*.dat", recursive=True)

if not adbladeloads_paths:
    print("Found no adbladeloads.dat files!")    
if not flowdata_paths:
    print("Found no flowdata.nc files!")
if not power_paths:
    print("Found no power files!")

# Sort pathnames, print for checking    
adbladeloads_paths.sort(key=numerical_sort_key)
flowdata_paths.sort(key=numerical_sort_key)
power_paths.sort(key=numerical_sort_key)
[print(path) for path in adbladeloads_paths]
[print(path) for path in flowdata_paths]
[print(path) for path in power_paths]

#%% Plot blade loads
if pp_adbladeloads_option:
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=[10,6])
    
    for filename_path in adbladeloads_paths:
        plot_adbladeloads(filename_path)
        
    # Validation (optional)    
    if check_against_DES:
        # Loads DTU10MW DES blade loads and plots them to validate against
        # The run parameters of the DES simulation are in the file
        from py_wake_ellipsys.examples.data.turbines.dtu10mw import dtu10mw_bladeloading
        r_nd_val = dtu10mw_bladeloading[:, 0] # Non-dimensionalised blade radius r/R [-]
        ft_val = -dtu10mw_bladeloading[:, 1] # Non-dimensionalised tangential force [-] (assumed axis system
        fn_val = dtu10mw_bladeloading[:, 3] # Non-dimensionalised azimuthal force [-] (assumed)
        ax[0].plot(r_nd_val, ft_val, label="DES")
        ax[1].plot(r_nd_val, fn_val, label="DES")
        
    
    ax[0].set_ylabel('ft $[-]$')
    ax[0].legend(loc=1)
    ax[1].set_ylabel('fn $[-]$')
    ax[1].set_xlabel('Blade radius r/R [-]')
    plot_title = "$\overline{\delta}$ parameter study\n"
    fig.suptitle(plot_title, y=1.01)
    plt.tight_layout()
    
    for axis in ax:
        axis.grid()
        
    if check_against_DES:
        fig.savefig('AD_b_loads_val.pdf', bbox_inches='tight')
    else:
        fig.savefig('AD_b_loads.pdf', bbox_inches='tight')
    print("AD blade loads plotted.")

#%% Power post-processing

if pp_power_option:
    
    power_results = np.empty([len(power_paths), 2])
    
    power_iterator = 0
    for filename_path in power_paths:
        power_result = pp_power(filename_path)
        power_results[power_iterator, :] = power_result
        power_iterator +=1 
    
    #np.savetxt('power_results.txt', power_results)

#%% Velocity deficit post-processing

if pp_advd_option:
    fig, ax = plt.subplots(figsize=[10, 6])
    
    # Initialise empty results array
    probed_advd_results = np.empty([len(flowdata_paths), 2])
    
    advd_iterator = 0 # This takes a while, so it's nice to see progress in the console
    for filename_path in flowdata_paths:
        print(f"Processing velocity deficits: {advd_iterator+1}/{len(flowdata_paths)}")
        probed_advd_result = pp_advd(filename_path)
        print(probed_advd_result)
        probed_advd_results[advd_iterator, :] = probed_advd_result
        advd_iterator +=1
    
    yd = dict(rotation=0, ha='right')     # Shortcut to flip y-label
    ax.set_ylabel(r'$\dfrac{ \langle U_{AD} \rangle}{U_{H}}$ [-]', yd)
    ax.set_xlabel('$x/D$')
    ax.grid()
    ax.legend(loc='upper right')
    
    fig.savefig('AD_b_deficits.pdf', bbox_inches='tight')
    #np.savetxt('advd_results.txt', probed_advd_results)
#%% 

if pp_power_option and pp_advd_option:
    results_matrix = np.empty([len(power_paths), 3])
    results_matrix[:, [0,1]] = power_results
    results_matrix[:, 2] = probed_advd_results[:,1]
    
    header = f"{parameter_name}                            Power [W]   Velocity deficit [-]"
    np.savetxt('results.txt', results_matrix, header=header)





