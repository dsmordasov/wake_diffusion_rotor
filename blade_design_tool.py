import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import shutil

import config # Import no matter what, it loads up the figure style

import h2s

U = 8 # Tested velocity [m/s]

design_name = "flattened"

#%% Load original blade design
# Read .ae and .htc files
ae_file_pathname = r"D:\AE EWEM MSc\T H E S I S\6_code\code repository\my_dtu_10mw\data\DTU_10MW_RWT_ae.dat"
htc_file_pathname = r"D:\AE EWEM MSc\T H E S I S\6_code\code repository\my_dtu_10mw\DTU_10MW_rigid_hawc2s_flattened.htc"

# Get blade aerodynamic data from the ae.dat file
ae_data = np.loadtxt(ae_file_pathname, skiprows=2, usecols=(0, 1, 2))
radius = ae_data[:, 0]  # Span [m]
chord = ae_data[:, 1]  # Chord length [m]
rel_thickness = ae_data[:, 2]  # Relative thickness (%)

N = np.size(radius) # Number of sections

# Where the twist data can be found in the .htc file
start_line = 110
n_sec = 27
end_line = start_line + n_sec

# Get blade twist data from the input .htc file
twist = np.genfromtxt(htc_file_pathname, skip_header=start_line, max_rows=27, usecols=(4, 5), comments=";",)
twist_radii = twist[:, 0] # Positions at which twist is defined in the .htc file
    
# Create a .geo blade equivalent format matrix
geo_mat = np.zeros([N, 4])
geo_mat[:, 0] = radius
geo_mat[:, 1] = np.interp(geo_mat[:, 0], twist[:, 0], twist[:, 1])
geo_mat[:, 2] = chord
geo_mat[:, 3] = rel_thickness

# Plot chord and twist
R = geo_mat[:, 0][-1] # Maximum radius [m]
radius_nd = geo_mat[:, 0] / R # Non-dimensionalised radius r/R [-]

fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].plot(radius_nd, geo_mat[:, 2], label="Original")
ax[0].set_ylabel("Chord [m]")

ax[1].plot(radius_nd, geo_mat[:, 1], label="Original")
ax[1].set_ylabel("Twist [deg]")
ax[1].set_xlabel("Blade radius r/R [-]")

#%% Changing the design
working_mat = geo_mat.copy()

cos_coords = np.arange(0, np.pi, 0.01)
cos_multiplier = np.cos(cos_coords)
# Linear chord addition
linear_chord_addition = False
if linear_chord_addition:
    linear_multiplier = radius_nd - 1.0
    chord_addition = -1.0 # [m]
    working_mat[:, 2] = working_mat[:, 2] + chord_addition * linear_multiplier
    
untwist_the_root = False
if untwist_the_root:
    linear_multiplier = radius_nd - 1.0
    twist_addition = 14.5 # [deg]
    working_mat[:, 1] = working_mat[:, 2] + twist_addition * linear_multiplier

# Plot new chord and twist on top of old
ax[0].plot(radius_nd, working_mat[:, 2], label="New")
ax[1].plot(radius_nd, working_mat[:, 1], label="New")
plt.legend(loc=2)

#%% Creating the new design 
new_htc_filename = rf"D:\AE EWEM MSc\T H E S I S\6_code\code repository\my_dtu_10mw\{design_name}.htc"
new_ae_filename = rf"D:\AE EWEM MSc\T H E S I S\6_code\code repository\my_dtu_10mw\data\{design_name}_ae.dat" 

try:
    shutil.copy(htc_file_pathname, new_htc_filename)
    shutil.copy(ae_file_pathname, new_ae_filename)
    print("Copied original .htc and ae.dat files successfully.")
except:
    print("Copying original .htc and ae.dat files FAILED. Script shutdown.")
    quit()

with open(new_htc_filename, 'r') as new_htc_file:
    content = new_htc_file.readlines()

with open(new_htc_filename, 'r+') as new_htc_file:
    for i, line in enumerate(new_htc_file):
        if "ae_filename" in line:
            content[i] = f"\tae_filename ./data/{design_name}_ae.dat;\n"
        if i >= start_line and i < end_line:
            line = line.strip().split(" ")
            if "" in line: line.remove("") # In single-digit sections there can be a whitespace
            line[4] = str(working_mat[i-start_line, 1]) + ";\n"
            line = "\t\t" + ' '.join(line)
            content[i] = line
            
with open(new_htc_filename, 'w') as new_htc_file:
    new_htc_file.writelines(content)
#%% Run HAWC2S

h2s.run_hawc2s(design_name)

h2s.pp_hawc2s_ind(design_name, U=U)

h2s.pp_hawc2s_pwr(design_name)

h2s.hawc2s_files_to_geo(design_name)

    


