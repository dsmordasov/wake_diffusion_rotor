import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import shutil

import config # Import no matter what, it loads up the figure style

import h2s

U = 8 # Tested velocity [m/s]

design_name = "new_design"
save_design = True

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

# Where the twist data can be found in the .htc file, and data in ae_dat file
htc_start_line = 110
n_sec = 27
htc_end_line = htc_start_line + n_sec
ae_start_line = 2

# Get blade twist data from the input .htc file
twist = np.genfromtxt(htc_file_pathname, skip_header=htc_start_line, max_rows=27, usecols=(4, 5), comments=";",)
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
ax[0].plot(radius_nd, geo_mat[:, 2], label="Baseline")
ax[0].set_ylabel("Chord [m]")

ax[1].plot(radius_nd, geo_mat[:, 1], label="Baseline")
ax[1].set_ylabel("Twist [deg]")
ax[1].set_xlabel("Blade radius r/R [-]")

#%% Changing the design
working_mat = geo_mat.copy()


# Linear chord addition
root_linear_chord_addition = False
if root_linear_chord_addition:
    linear_multiplier = radius_nd - 1.0
    chord_addition = -2.0 # [m]
    changes_until = 0.4 # r/R until which we change the chord
    i = 0
    for radius in radius_nd:
        if radius <= changes_until:
            working_mat[i, 2] = working_mat[i, 2] - chord_addition * linear_multiplier[i]
        i += 1

tip_cos_chord_addition = True
if tip_cos_chord_addition:
    ncos_multiplier = - np.cos(radius_nd * np.pi) # Negative cos multiplier, from [-1, 1]
    chord_addition = - 0.5 # [m]
    changes_from = 0.6 # r/R, counting down from 1.0, until which we change the chord
    i = 0
    for radius in radius_nd:
        if radius  >= changes_from:
            working_mat[i, 2] = working_mat[i, 2] + chord_addition * ncos_multiplier[i]
        i += 1

cos_chord_addition = False
if cos_chord_addition:
    ncos_multiplier = - np.cos(radius_nd * np.pi) # Negative cos multiplier, from [-1, 1]
    print(ncos_multiplier)
    chord_addition = 2.0 # [m]
    working_mat[:, 2] = working_mat[:, 2] + chord_addition * ncos_multiplier
    
twist_the_root = False
if twist_the_root:
    linear_multiplier = - (radius_nd - 1.0)
    twist_addition = -10.0 # Positive value adds twist [deg]
    changes_until = 0.25 # r/R until which we change the chord
    i = 0
    for radius in radius_nd:
        if radius <= changes_until:
            working_mat[i, 1] = working_mat[i, 1] + twist_addition * linear_multiplier[i]
        i += 1

make_root_cylindrical = False
if make_root_cylindrical:
    changes_until = 0.25 # r/R until which we use cylindrical airfoil
    i = 0
    for radius in radius_nd:
        if radius <= changes_until:
            working_mat[i, 3] = 100.0
        i += 1
            


# Plot new chord and twist on top of old
ax[0].plot(radius_nd, working_mat[:, 2], label=f"{design_name}")
ax[1].plot(radius_nd, working_mat[:, 1], label=f"{design_name}")
plt.legend(loc=2)

# Format the working_mat twists into the .htc version
htc_twist = np.interp(twist_radii, working_mat[:, 0], working_mat[:, 1])


#%% Creating the new design .htc and ae.dat files
new_htc_filename = rf"D:\AE EWEM MSc\T H E S I S\6_code\code repository\my_dtu_10mw\{design_name}.htc"
new_ae_filename = rf"D:\AE EWEM MSc\T H E S I S\6_code\code repository\my_dtu_10mw\data\{design_name}_ae.dat" 

try:
    shutil.copy(htc_file_pathname, new_htc_filename)
    shutil.copy(ae_file_pathname, new_ae_filename)
    print("Copied original .htc and ae.dat files successfully.")
except:
    print("Copying original .htc and ae.dat files FAILED. Script shutdown.")
    quit()
    
# .htc file creation
# Read the old file
with open(new_htc_filename, 'r') as new_htc_file:
    content = new_htc_file.readlines()

# Create what to write in the new file
with open(new_htc_filename, 'r+') as new_htc_file:
    for i, line in enumerate(new_htc_file):
        if "ae_filename" in line:
            content[i] = f"\tae_filename ./data/{design_name}_ae.dat;\n"
        if i >= htc_start_line and i < htc_end_line:
            line = line.strip().split(" ")
            if "" in line: line.remove("") # In single-digit sections there can be a whitespace
            line[4] = str(htc_twist[i-htc_start_line]) + ";\n"
            line = "\t\t" + ' '.join(line)
            content[i] = line
            
# Write the new file
with open(new_htc_filename, 'w') as new_htc_file:
    new_htc_file.writelines(content)
    
# ae.dat file creation
# Read the old file
with open(new_ae_filename, 'r') as new_ae_file:
    content = new_ae_file.readlines()

# Create what to write in the new file
with open(new_ae_filename, 'r+') as new_ae_file:
    for i, line in enumerate(new_ae_file):
        if i >= ae_start_line:
            line = line.strip().split("\t")
            line[0] = str(working_mat[i - ae_start_line, 0])
            line[1] = str(working_mat[i - ae_start_line, 2])
            line[2] = str(working_mat[i - ae_start_line, 3])
            line = '\t'.join(line) + "\n"
            content[i] = line

# Write the new file            
with open(new_ae_filename, 'w') as new_ae_file:
    new_ae_file.writelines(content)

#%% Run HAWC2S

h2s.run_hawc2s(design_name)

h2s.pp_hawc2s_ind(design_name, U=U)

h2s.pp_hawc2s_pwr(design_name)

h2s.pp_hawc2s_bladepower(design_name)

if save_design:
    h2s.hawc2s_files_to_geo(design_name)

    


