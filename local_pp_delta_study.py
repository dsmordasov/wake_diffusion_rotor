"""
This code was made in order to post process the results of `delta_study.py`,
focusing on the effects of varying \overline{\delta} and thus the force 
distribution over a wind turbine blade on the wake of a single wind turbine.
"""

import glob
import re
import os

import config

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray

# Run parameters, DTU10MW
rho = 1.225 # Air density [kg/m^3]
UH = 8 # Simulation velocity [m/s]
D = 178.3 # Rotor diameter [m]
R = D / 2 # Rotor radius [m]
#R = 86.36599 # Rotor radius, airfoil [m]
zh = 119 # Reference hub height [m]

# Let the hardcoding begin, change path if required
os.chdir(r"D:\AE EWEM MSc\T H E S I S\6_code\code repository\run_data\joukowski")

# Run options
parameter_name = "delta" # Set to either of: ["grid", "delta", design]

pp_adbladeloads_option = False # Set to True to post process AD blade loads
check_against_DES = True # Set to True to check against DTU10MW loading

pp_advd_option = False # Set to True to post process velocity deficits (takes a while)
advd_x_probe_D = 5 # Distance x at which we probe for velocity deficit [D]

# Hardcoding for Gaussian hats (GH) velocity profile analysis
analysed_spacing = 2.5 # Distance between probed Gaussian hat distances [D]
n_hats = 6 # How many Gaussian hats would you like?
first_probed_distance = -analysed_spacing
last_probed_distance = (n_hats - 1) * analysed_spacing 
analysed_flowdata_paths = ['d0.15/flowdata.nc', 'd0.35/flowdata.nc', 'd0.60/flowdata.nc']


gh_option = True # Set to True to get a gaussian hat graph of velocity deficits
if gh_option:
    analysed_deltas = np.array([0.15, 0.35, 0.60])
    colors = ['r', 'g', 'b'] # Set as many as there are analysed deltas
    analysed_downstream_xs = np.arange(first_probed_distance, last_probed_distance, analysed_spacing)

gh_single_option = False
if gh_single_option:
    analysed_deltas = np.array([0.15])
    analysed_downstream_xs = 0

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
    return format(parameter_value, '.2f')

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

    ax[0].plot(r_nd, Ft, label="$\overline{\delta}$="+parameter_value)
    ax[1].plot(r_nd, Fn, label="$\overline{\delta}$="+parameter_value)
    return

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
    ax.plot(x / D, ADvar / UH, label="$\overline{\delta}$="+parameter_value)
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

def pp_single_gh(filename_path, n=128):
    """
    Post process netcdf file to produce gaussian hat plots.
    """
    U_profiles = np.zeros([1, n])
    
    z_line = 1.2 * np.linspace(-R, R, n) + zh # Vertical line of the velocity profile

    data = xarray.open_dataset(filename_path)
    pos = analysed_downstream_xs
    
    print(f"x/D = {pos}")
    U_profiles = (data.U.interp(x=pos*D, y=0, z=z_line) / UH) # Velocity deficit [-] (0 none, 1 fully stopped flow)
    print(U_profiles)
    U_hat =U_profiles #+ pos #- analysed_spacing
    print(np.mean(U_hat))
    #print((U_hat - pos ) / analysed_spacing)
    ax.plot(U_hat + pos, (z_line - 29.85) / D) # 29.85m  is the distance to ground
    ax.axvline(pos, color='k', linestyle='dotted')

def pp_gh(filename_path, n=128, color='r'):
    """
    Post process netcdf file to produce gaussian hat plots.
    """
    U_profiles = np.zeros([len(analysed_downstream_xs), n])
    
    z_line = 1.5 * np.linspace(-R, R, n) + zh # Vertical line of the velocity profile

    data = xarray.open_dataset(filename_path)
    
    for i, pos in enumerate(analysed_downstream_xs):
        print(f"x/D = {pos}")
        U_profiles[i] = analysed_spacing * (data.U.interp(x=pos*D, y=0, z=z_line) / UH) # Velocity deficit [-] (0 none, 1 fully stopped flow)
        U_hat = U_profiles[i] 
        
        # Plot the gaussian hat, along with text denoting the position probed
        ax.plot(U_hat + pos - analysed_spacing, (z_line - 29.85) / D, color=color) # 29.85m  is the distance to ground
        ax.axvline(pos, color='k', linestyle='dotted')
        text_x_offset = -1e-1 * analysed_spacing
        text_y_offset = -2e-2 * analysed_spacing
        ax.text(pos + text_x_offset, text_y_offset, f"x/D = {pos}", rotation=90)
        

# Sort pathnames, print for checking    
adbladeloads_paths.sort(key=numerical_sort_key)
flowdata_paths.sort(key=numerical_sort_key)

[print(path) for path in adbladeloads_paths]
[print(path) for path in flowdata_paths]

#%% Plot blade loads
if pp_adbladeloads_option:
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=[10,6])
    
    for filename_path in adbladeloads_paths:
        plot_adbladeloads(filename_path)
        
    # Validation (optional)    
    if check_against_DES:
        # Plots loaded DTU10MW RWT DES values to validate against
        # The run parameters of the DES simulation are in the file
        from dtu10mw_des_bladeloads import dtu10mw_bladeloading
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


#%% Velocity deficit post-processing

if pp_advd_option:
    fig, ax = plt.subplots(figsize=[10, 6])
    
    # Initialise empty results array
    probed_advd_results = np.empty([len(flowdata_paths), 2])
    
    # This takes a while, so it is nice to see progress in the console
    for advd_iterator, filename_path in enumerate(flowdata_paths):
        print(f"Processing velocity deficits: {advd_iterator+1}/{len(flowdata_paths)}")
        probed_advd_result = pp_advd(filename_path)
        print(probed_advd_result)
        probed_advd_results[advd_iterator, :] = probed_advd_result
    
    yd = dict(rotation=0, ha='right')     # Shortcut to flip y-label, credit to Mads
    ax.set_ylabel(r'$\dfrac{ \langle U_{AD} \rangle}{U_{H}}$ [-]', yd)
    ax.set_xlabel('$x/D$')
    ax.grid()
    ax.legend(loc='upper right')
    
    fig.savefig('AD_b_deficits.pdf', bbox_inches='tight')
    #np.savetxt('advd_results.txt', probed_advd_results)

#%% Gaussian hats

if gh_single_option:
    fix, ax = plt.subplots(figsize=[6, 6])
    colors = ['r', 'g', 'b']
    
    for gh_iterator, filename_path in enumerate(analysed_flowdata_paths):
        color = colors[gh_iterator]
        print(color)
        pp_single_gh(filename_path) #color=color)
    
    # TODO: Set ylim to [0, 1]
    #plt.xticks(analysed_downstream_xs)
    plt.ylabel("z/D [-]")
    plt.xlabel("$U/U_{inf}$ [-]")


if gh_option:
    fix, ax = plt.subplots(figsize=[12, 4])
    
    for gh_iterator, filename_path in enumerate(analysed_flowdata_paths):
        color = colors[gh_iterator]
        print(f"Drawing gaussian hats for delta={analysed_deltas[gh_iterator]}...")
        print(color)
        pp_gh(filename_path, color=color)
    
    # Cut the graph to show the wake beyond the wind turbine width
    wake_margin = 0.15
    plt.ylim([0 - wake_margin, 1 + wake_margin])
    
    # The following block makes the 1|0, 0.5, 0|1... ticks on the x axis
    n_ticks = len(analysed_downstream_xs)/2 + 2 # Number of ticks we will have
    x_ticks = np.arange(first_probed_distance, last_probed_distance, analysed_spacing / 2)[:-1]
    internal_ticks = int(n_ticks) * ["1|0", "0.5"]
    x_labels = [*internal_ticks, "1"]
    latex_fontsize = 18
    plt.xticks(ticks=x_ticks, labels = x_labels)
    plt.ylabel("$z/D \ [-]$", fontsize=latex_fontsize)
    plt.xlabel("$U/U_{inf} \ [-]$", fontsize=latex_fontsize)
    
    # Now let us create a hardcoded legend for this graph
    hardcode_legend_elements = [mpl.lines.Line2D([0], [0], color=colors[0], lw=2, label=r"$\hat{\delta}$ = " + str(analysed_deltas[0])),
                                mpl.lines.Line2D([0], [0], color=colors[1], lw=2, label=r"$\hat{\delta}$ = " + str(analysed_deltas[1])),
                                mpl.lines.Line2D([0], [0], color=colors[2], lw=2, label=r"$\hat{\delta}$ = " + str(analysed_deltas[2]))]
    ax.legend(handles=hardcode_legend_elements, loc="upper left")