"""
This code was made in order to post process the results of `delta_study.py`,
focusing on the effects of varying \overline{\delta} and thus the force 
distribution over a wind turbine blade on the wake of a single wind turbine.
"""

import glob
import re
import os

import config # import for plotting style

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.pyplot import cm
import xarray

# Run parameters, DTU10MW
rho = 1.225 # Air density [kg/m^3]
UH = 8 # Simulation velocity [m/s]
D = 178.3 # Rotor diameter [m]
R = D / 2 # Rotor radius [m]
#R = 86.36599 # Rotor radius, airfoil [m]
zh = 119 # Reference hub height [m]
z_norm = zh / D # Normalising hub height [D]

# Let the hardcoding begin, change path if required
os.chdir(r"D:\AE EWEM MSc\T H E S I S\6_code\code repository\run_data\joukowski")
analysed_flowdata_paths = ['d0.15/flowdata.nc', 'd0.35/flowdata.nc', 'd0.60/flowdata.nc']

# Run options
parameter_name = "delta" # Set to either of: ["grid", "delta", design]

pp_adbladeloads_option = False # Set to True to post process AD blade loads
check_against_DES = False # Set to True to check against DTU10MW loading

pp_advd_option = False # Set to True to post process velocity deficits (takes a while)
advd_x_probe_D = 5 # Distance x at which we probe for velocity deficit [D]

#%% Hardcoding for Gaussian hats (gh) velocity profile analysis
analysed_spacing = 2.5 # Distance between probed Gaussian hat distances [D]
n_hats = 5 # How many Gaussian hats would you like in your graph?
first_probed_distance = -analysed_spacing
last_probed_distance = (n_hats - 1) * analysed_spacing 
n_deltas = len(analysed_flowdata_paths) # Number of deltas investigated

analysed_deltas = np.array([0.15, 0.35, 0.60])

gh_option = False # Set to True to get a gaussian hat graph of velocity deficits
if gh_option:
    colors = ['r', 'g', 'b'] # Set as many as there are analysed deltas
    analysed_downstream_xs = np.arange(first_probed_distance, last_probed_distance, analysed_spacing)

tke_option = False
tke_cmap = cm.jet

# Momentum transport terms 
# 0: d(p)/dx
# 1: d\overline(u'u')/dx
# 2: d\overline(u'v')/dy
# 3: d\overline(u'w')/dw
mtt_option = False # range(4) for mtt terms, range(6) to plot u'v' and u'w' too
mtt_cmap = cm.RdBu

vel_option = True

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
        #ax.plot(U_hat + pos - analysed_spacing, (z_line - 29.85) / D, color=color) # 29.85m  is the distance to ground
        ax.plot(U_hat + pos - analysed_spacing, (z_line - (29.85 + D/2)) / D, color=color) # 29.85m  is the distance to ground
        ax.axvline(pos, color='k', linestyle='dotted')
        text_x_offset = -1e-1 * analysed_spacing
        text_y_offset = -2.5e-1 * analysed_spacing
        ax.text(pos + text_x_offset, text_y_offset, f"x/D = {pos}", rotation=90)

def calc_sigma(data):
    '''
    Input: xarray object
    Calculate flow field variables:
        tau = turbulent time scale = k/eps
        s_ij = normalized symmetric velocity gradient = 0.5*k/eps*(dU_i/dxj + dU_j/dxi)
        omega_ij = normalized anti-symmetric velocity gradient = 0.5*k/eps*(dU_i/dxj - dU_j/dxi)
        eta1 = 1st invariant = tr(s^2)
        eta2 = 2nd invariant = tr(omega^2)
        sigma = local shear parameter = k/eps*sqrt((dU_i/dx_j)^2) = sqrt(eta1 - eta2)
        
        courtesy of Mads Christian Baungaard
    '''
   
    # Store variables in this dictionary
    data_vars = {}
   
    # Also, include the original data:
    data_vars['U'] = data['U']
    data_vars['V'] = data['V']
    data_vars['W'] = data['W']
    data_vars['P'] = data['P']
    data_vars['muT'] = data['muT'] # Dynamic viscosity [kg m^-1 s^-1]
    data_vars['tke'] = data['tke'] # Turbulent kinetic energy [m^2 s^-2]
    data_vars['epsilon'] = data['epsilon'] 
   
    # 1st derivatives
    diffx = data.differentiate(coord='x')
    diffy = data.differentiate(coord='y')
    diffz = data.differentiate(coord='z')

    # Calculate turbulent time scale
    data_vars['tau'] = data['tke']/data['epsilon']

    # Calculate Reynolds stresses (mu / rho = nu === kinematic viscosity)
    # Using Boussinesq eddy viscosity hypothesis 
    uv = -data['muT'] / rho * (diffy['U'] + diffx['V'])
    vu = uv
    uw = -data['muT'] / rho * (diffz['U'] + diffx['W'])
    wu = uw
    vw = -data['muT'] / rho * (diffz['V'] + diffy['W'])
    wv = vw
    uu = 2. / 3. * data['tke'] - data['muT'] / rho * 2 * diffx['U']
    vv = 2. / 3. * data['tke'] - data['muT'] / rho * 2 * diffy['V']
    ww = 2. / 3. * data['tke'] - data['muT'] / rho * 2 * diffz['W']
    
    # Calculate momentum transport terms in RANS equation
    data_vars['term0'] = data_vars['P'].differentiate(coord='x')
    data_vars['term1'] = uu.differentiate(coord='x')
    data_vars['term2'] = uv.differentiate(coord='y')
    data_vars['term3'] = uw.differentiate(coord='z')
    data_vars['term4'] = uv
    data_vars['term5'] = uw
   
    # Calculate TKE production to dissipation ratio
    data_vars['Peps'] = (- (uu * diffx['U'] + uv * diffy['U'] + uw * diffz['U']) \
                     - (vu * diffx['V'] + vv * diffy['V'] + vw * diffz['V']) \
                     - (wu * diffx['W'] + wv * diffy['W'] + ww * diffz['W']))/data['epsilon']

   
    # Calculate eta1
    s11 = 0.5*data_vars['tau']*2*diffx['U']
    s22 = 0.5*data_vars['tau']*2*diffx['V']
    s33 = -(s11+s22) # s_ij is traceless, i.e. tr(s)=0
    s12 = 0.5*data_vars['tau']*(diffy['U'] + diffx['V'])
    s13 = 0.5*data_vars['tau']*(diffz['U'] + diffx['W'])
    s23 = 0.5*data_vars['tau']*(diffz['V'] + diffy['W'])
    ss11=s11**2+s12**2+s13**2  # Helper variable 1
    ss22=s12**2+s22**2+s23**2  # Helper variable 2
    ss33=s13**2+s23**2+s33**2  # Helper variable 3
    data_vars['eta1'] = ss11 + ss22 + ss33  # This is the EllipSys way!
   
    # Calculate eta2
    o12 = 0.5*data_vars['tau']*(diffy['U'] - diffx['V'])
    o13 = 0.5*data_vars['tau']*(diffz['U'] - diffx['W'])
    o23 = 0.5*data_vars['tau']*(diffz['V'] - diffy['W'])
    oo11=-o12**2 - o13**2  # Helper vairable 1
    oo22=-o12**2 - o23**2  # Helper vairable 2
    oo33=-o13**2 - o23**2  # Helper vairable 3
    data_vars['eta2'] = oo11 + oo22 + oo33   # Again, this is the EllipSys way.
   
    # Calculate local shear parameter
    data_vars['sigma'] = np.sqrt(data_vars['eta1'] - data_vars['eta2'])
   
    # Return new netcdf dataset
    sigma_data = xarray.Dataset(data_vars=data_vars)
    return sigma_data

def pp_tke(sigma_data, axis_counter):
    ''' Post-process turbulent kinetic energy.
    
    Produces an [M x 2] sized graph, with left side containing M y-plane cuts
    of normalised TKE, and the right side containing M z-plane cuts of 
    normalised TKE.
    
    Input:
        sigma_data = xarray object
        axis_counter = 
    '''
    current_label = r"$\overline{\delta}$ = " + str(analysed_deltas[axis_counter])
    
    x_limit = (-0.5, 5) # [D]
    y_limit = (-0.7, 0.7) # [D]
    y_limit_zplane = (-0.65, 0.9)

    # Left graphs: y plane
    zplane_data = sigma_data.interp(z=zh) # Slice at hub height
    X, Y = np.meshgrid(zplane_data.x, zplane_data.y)
    tim_z = np.sqrt(zplane_data['tke'].T * (2/3)) / UH # Turbulence intensity measure
    left_plot = axes[axis_counter, 0].contourf(X / D, Y / D, tim_z,
                                               cmap=tke_cmap,
                                               levels=np.linspace(0.04, 0.20, 11))
    axes[axis_counter, 0].text(x_limit[0] + 1e-1,
                               y_limit[1] + 0.5e-1,
                               current_label)
    
    
    # Right graphs: z plane
    yplane_data = sigma_data.interp(y=0) # Slice vertically
    tim_y = np.sqrt(yplane_data['tke'].T * (2/3)) / UH # Turbulence intensity measure
    X, Z = np.meshgrid(yplane_data.x, yplane_data.z)
    right_plot = axes[axis_counter, 1].contourf(X / D, Z / D - z_norm, tim_y,
                                                cmap=tke_cmap,
                                                levels=np.linspace(0.04, 0.20, 11))
    axes[axis_counter, 1].text(x_limit[0] + 1e-1,
                               y_limit_zplane[1] + 1e-1,
                               current_label)

    plt.xlim(x_limit)
    axes[axis_counter, 0].set_ylim(y_limit)
  
    return left_plot, right_plot

def pp_vel(sigma_data, axis_counter):
    ''' Post-process turbulent kinetic energy.
    '''
    current_label = r"$\overline{\delta}$ = " + str(analysed_deltas[axis_counter])
    
    x_limit = (-0.5, 5) # [D]
    y_limit = (-0.7, 0.7) # [D]
    y_limit_zplane = (-0.65, 0.9)

    # Left graphs: y plane
    zplane_data = sigma_data.interp(z=zh) # Slice at hub height
    X, Y = np.meshgrid(zplane_data.x, zplane_data.y)
    vel_z = zplane_data['U'].T/UH # Streamwise elocity
    left_plot = axes[axis_counter, 0].contourf(X / D, Y / D, vel_z,
                                               cmap=tke_cmap,)
                                               #levels=np.linspace(0.04, 0.20, 11))
    axes[axis_counter, 0].text(x_limit[0] + 1e-1,
                               y_limit[1] + 0.5e-1,
                               current_label)
    
    
    # Right graphs: z plane
    yplane_data = sigma_data.interp(y=0) # Slice vertically
    vel_y = yplane_data['U'].T/UH # Streamwise elocity
    X, Z = np.meshgrid(yplane_data.x, yplane_data.z)
    right_plot = axes[axis_counter, 1].contourf(X / D, Z / D - z_norm, vel_y,
                                                cmap=tke_cmap,)
                                                #levels=np.linspace(0.04, 0.20, 11))
    axes[axis_counter, 1].text(x_limit[0] + 1e-1,
                               y_limit_zplane[1] + 1e-1,
                               current_label)

    plt.xlim(x_limit)
    axes[axis_counter, 0].set_ylim(y_limit)
  
    return left_plot, right_plot
    
def pp_mtt(sigma_data, mtt_term_index, axis_counter):
    ''' Post-process momentum transport terms.
            Momentum transport terms:
        0: d(p)/dx
        1: d\overline(u'u')/dx
        2: d\overline(u'v')/dy
        3: d\overline(u'w')/dw
    '''
    # Left graphs: y plane
    zplane_data = sigma_data.interp(z=zh) # Slice at hub height
    X, Y = np.meshgrid(zplane_data.x, zplane_data.y)
    
    # Right graphs: z plane
    yplane_data = sigma_data.interp(y=0) # Slice vertically
    X_z, Z = np.meshgrid(yplane_data.x, yplane_data.z)

    # Depending on the currently plotted term
    if mtt_term_index == 0:
        current_xlabel = r"$\dfrac{\partial p}{\partial x}$"
        plot_data_l = zplane_data['term0'].T
        plot_data_r = yplane_data['term0'].T
        current_levels_l = np.linspace(-1.5, 1.5, 10) # -1.5 -> 0.5
        current_levels_r = np.linspace(-1.8, 1.8, 10) # -1.8 -> 0.6
    elif mtt_term_index == 1:
        current_xlabel = r"$\dfrac{\partial \overline{u' u'}}{\partial x}$"
        plot_data_l = zplane_data['term1'].T
        plot_data_r = yplane_data['term1'].T
        current_levels_l = np.linspace(-0.02, 0.02, 10) # -0.02 -> 0.02
        current_levels_r = np.linspace(-0.02, 0.02, 10) # -0.02 -> 0.02
    elif mtt_term_index == 2:
        current_xlabel = r"$\dfrac{\partial \overline{u' v'}}{\partial y}$"
        plot_data_l = zplane_data['term2'].T
        plot_data_r = yplane_data['term2'].T
        current_levels_l = np.linspace(-0.06, 0.06, 10) # -0.06 -> 0.06
        current_levels_r = np.linspace(-0.048, 0.048, 10) # -0.016 -> 0.048
    elif mtt_term_index == 3:
        current_xlabel = r"$\dfrac{\partial \overline{u' w'}}{\partial z}$"
        plot_data_l = zplane_data['term3'].T
        plot_data_r = yplane_data['term3'].T
        current_levels_l = np.linspace(-0.048, 0.048, 10) # -0.016 -> 0.048
        current_levels_r = np.linspace(-0.15, 0.15, 10) # -0.2 -> 0.08
    elif mtt_term_index == 4:
        current_xlabel = r"$\overline{u' v'}$"
        plot_data_l = zplane_data['term4'].T
        plot_data_r = yplane_data['term4'].T
        current_levels_l = np.linspace(-1.0, 1.0, 10) # -1.00 -> 1.00
        current_levels_r = np.linspace(-0.12, 0.12, 10) # -0.12 -> 0.09
    elif mtt_term_index == 5:
        current_xlabel = r"$\overline{u' w'}$"
        plot_data_l = zplane_data['term5'].T
        plot_data_r = yplane_data['term5'].T
        current_levels_l = np.linspace(-0.36, 0.36, 10) # -0.36 -> 0.12
        current_levels_r = np.linspace(-1.2, 1.2, 10) # -1.2 -> 1.2
        
    current_label = r"$\overline{\delta}$ = " + str(analysed_deltas[axis_counter])
    
    x_limit = (-0.5, 5) # [D]
    y_limit = (-0.7, 0.7) # [D]
    y_limit_zplane = (-0.65, 0.9)
    
    left_plot = axes[axis_counter, 0].contourf(X / D, Y / D, plot_data_l,
                                               cmap=mtt_cmap,
                                               levels=current_levels_l)#, vmin=vmin_left, vmax=vmax_left)
    axes[axis_counter, 0].text(x_limit[0] + 1e-1,
                               y_limit[1] + 0.5e-1,
                               current_label)
    
    right_plot = axes[axis_counter, 1].contourf(X_z / D, Z / D - z_norm, plot_data_r,
                                                cmap=mtt_cmap,
                                                levels=current_levels_r)#, vmin=vmin_right, vmax=vmax_right)
    axes[axis_counter, 1].text(x_limit[0] + 1e-1,
                               y_limit_zplane[1] + 1e-1,
                               current_label)

    plt.xlim(x_limit)
    axes[axis_counter, 0].set_ylim(y_limit)
    #axes[axis_counter, 1].set_ylim(y_limit)
  
    return left_plot, right_plot, current_xlabel
#%% Find, sort and print pathnames for various data files
# Create lists of found pathnames
adbladeloads_paths = glob.glob("**/*_adbladeloads.dat", recursive=True)
flowdata_paths = glob.glob("**/flowdata.nc", recursive=True)
power_paths = glob.glob("**/Power*.dat", recursive=True)

if not adbladeloads_paths:
    print("Found no adbladeloads.dat files!")    
if not flowdata_paths:
    print("Found no flowdata.nc files!")

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


#%% Post-process velocity deficits

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

#%% Plot velocity deficit Gaussian hats

if gh_option:
    fig, ax = plt.subplots(figsize=[12, 4])
    
    for gh_iterator, filename_path in enumerate(analysed_flowdata_paths):
        color = colors[gh_iterator]
        print(f"Drawing gaussian hats for delta={analysed_deltas[gh_iterator]}...")
        print(color)
        pp_gh(filename_path, color=color)
    
    # Cut the graph to show the wake beyond the wind turbine width
    wake_margin = 0.
    plt.ylim([-0.6674 - wake_margin, 0.6674 + wake_margin]) # (29.85)/D + 0.5 = 0.6674, height of hub [D]
    
    # The following block makes the 1|0, 0.5, 0|1... ticks on the x axis
    n_ticks = len(analysed_downstream_xs)/2 + 2 # Number of ticks we will have
    x_ticks = np.arange(first_probed_distance, last_probed_distance, analysed_spacing / 2)[:-1]
    internal_ticks = int(n_ticks) * ["1|0", "0.5"]
    x_labels = [*internal_ticks, "1"]
    latex_fontsize = 18
    plt.xticks(ticks=x_ticks, labels = x_labels)
    plt.ylabel("$z/D \ [-]$", fontsize=latex_fontsize)
    plt.xlabel("$U/U_{H} \ [-]$", fontsize=latex_fontsize)
    
    # Now let us create a hardcoded legend for this graph
    hardcode_legend_elements = [mpl.lines.Line2D([0], [0], color=colors[0], lw=2, label=r"$\hat{\delta}$ = " + str(analysed_deltas[0])),
                                mpl.lines.Line2D([0], [0], color=colors[1], lw=2, label=r"$\hat{\delta}$ = " + str(analysed_deltas[1])),
                                mpl.lines.Line2D([0], [0], color=colors[2], lw=2, label=r"$\hat{\delta}$ = " + str(analysed_deltas[2]))]
    ax.legend(handles=hardcode_legend_elements, loc="upper left")
    
    ax.set_xlim([-analysed_spacing*1.3, (n_hats-2)*analysed_spacing])
    fig.savefig('../../plots/jou_gaussian_hats.pdf', bbox_inches='tight')


#%% Plot turbulent kinetic energy

if tke_option:
    fig, axes = plt.subplots(3, 2, sharex=True, figsize=(12, 6))
    
    for tke_iterator, filename_path in enumerate(analysed_flowdata_paths):
        print(f"Ploting turbulence intensity measure {tke_iterator + 1}/{n_deltas}")
        
        data = xarray.open_dataset(filename_path)
        pped_data = calc_sigma(data)
        
        left_plot, right_plot = pp_tke(pped_data, tke_iterator)
        
        current_delta = analysed_deltas[tke_iterator]
    
    colorbar_left = fig.colorbar(left_plot, ax=axes[:, 0])
    colorbar_right = fig.colorbar(right_plot, ax=axes[:, 1])
    
    
    
    colorbar_left.ax.set_xlabel("    TKE $[-]$")
    colorbar_right.ax.set_xlabel("    TKE $[-]$")
    
    axes[n_deltas-1, 0].set_xlabel("x/D $[-]$")
    axes[n_deltas-1, 1].set_xlabel("x/D $[-]$")
    
    axes[1, 0].set_ylabel("y/D $[-]$")
    axes[1, 1].set_ylabel("z/D $[-]$")
    fig.savefig('../../plots/jou_tke.pdf', bbox_inches='tight')

#%% Plot RANS equation terms

if mtt_option:    
    for mtt_index in mtt_option:
        fig, axes = plt.subplots(3, 2, sharex=True, figsize=(12, 6))
        
        for mtt_iterator, filename_path in enumerate(analysed_flowdata_paths):
            print(f"Ploting momentum transport term {mtt_index}:  {mtt_iterator + 1}/{n_deltas}")
            
            data = xarray.open_dataset(filename_path)
            pped_data = calc_sigma(data)
            
            left_plot, right_plot, current_xlabel = pp_mtt(pped_data, mtt_index, mtt_iterator)
            
            current_delta = analysed_deltas[mtt_iterator]
        
        
        colorbar_left = fig.colorbar(left_plot, ax=axes[:, 0])
        colorbar_right = fig.colorbar(right_plot, ax=axes[:, 1])
        
        colorbar_left.ax.set_xlabel(current_xlabel)
        colorbar_right.ax.set_xlabel(current_xlabel)
        
        axes[n_deltas-1, 0].set_xlabel("x/D $[-]$")
        axes[n_deltas-1, 1].set_xlabel("x/D $[-]$")
        
        axes[1, 0].set_ylabel("y/D $[-]$")
        axes[1, 1].set_ylabel("z/D $[-]$")
        
        fig.savefig(f'../../plots/jou_mtt_{mtt_index}.pdf', bbox_inches='tight')

#%% Plot velocity

if vel_option:
    fig, axes = plt.subplots(3, 2, sharex=True, figsize=(12, 6))
    
    for vel_iterator, filename_path in enumerate(analysed_flowdata_paths):
        print(f"Ploting velocity {vel_iterator + 1}/{n_deltas}")
        
        data = xarray.open_dataset(filename_path)
        pped_data = calc_sigma(data)
        
        left_plot, right_plot = pp_vel(pped_data, vel_iterator)
        
        current_delta = analysed_deltas[vel_iterator]
    
    colorbar_left = fig.colorbar(left_plot, ax=axes[:, 0])
    colorbar_right = fig.colorbar(right_plot, ax=axes[:, 1])
    
    colorbar_left.ax.set_xlabel(r"    $\dfrac{U}{U_H}$ $[-]$")
    colorbar_right.ax.set_xlabel(r"    $\dfrac{U}{U_H}$ $[-]$")
    
    axes[n_deltas-1, 0].set_xlabel("x/D $[-]$")
    axes[n_deltas-1, 1].set_xlabel("x/D $[-]$")
    
    axes[1, 0].set_ylabel("y/D $[-]$")
    axes[1, 1].set_ylabel("z/D $[-]$")
    fig.savefig('../../plots/jou_vel.pdf', bbox_inches='tight')

# 30/5/22 - 1 hour work