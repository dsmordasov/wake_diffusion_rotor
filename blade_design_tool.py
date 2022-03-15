import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import shutil

design_name = "new_design"

# Read .ae and .htc files
ae_file_pathname = r"D:\AE EWEM MSc\T H E S I S\6_code\code repository\my_dtu_10mw\data\DTU_10MW_RWT_ae.dat"
htc_file_pathname = r"D:\AE EWEM MSc\T H E S I S\6_code\code repository\my_dtu_10mw\DTU_10MW_rigid_hawc2s.htc"


# Get blade aerodynamic data from the ae.dat file
ae_data = np.loadtxt(ae_file_pathname, skiprows=2, usecols=(0, 1, 2))
radius = ae_data[:, 0]  # Span [m]
chord = ae_data[:, 1]  # Chord length [m]
rel_thickness = ae_data[:, 2]  # Relative thickness (%)

N = np.size(radius)
    
# Get blade twist data from the input .htc file
twist = np.genfromtxt(
    htc_file_pathname,
    # --------------------------------------------------
    skip_header=110,  # This line may need to be changed.
    # --------------------------------------------------
    max_rows=27,
    usecols=(4, 5),
    comments=";",
)
twist_radii = twist[:, 0] # Positions at which twist is defined in the .htc file
    
# Create a .geo blade equivalent format matrix
geo_mat = np.zeros([N, 4])
geo_mat[:, 0] = radius
geo_mat[1:, 1] = np.interp(geo_mat[1:, 0], twist[:, 0], twist[:, 1])
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

# Change chord and twist (?)

working_mat = geo_mat.copy()

cos_coords = np.arange(0, np.pi, 0.01)
cos_multiplier = np.cos(cos_coords)

# Linear chord addition
linear_chord_addition = True
if linear_chord_addition:
    linear_multiplier = radius_nd - 1.0
    chord_addition = -1. # [m]
    working_mat[:, 2] = working_mat[:, 2] + chord_addition * linear_multiplier

# Plot new chord and twist on top of old
ax[0].plot(radius_nd, working_mat[:, 2], label="New")
ax[1].plot(radius_nd, working_mat[:, 1], label="New")
plt.legend(loc=2)

# Output a an .htc and a ae.dat file
# .htc file
new_htc_filename = rf"D:\AE EWEM MSc\T H E S I S\6_code\code repository\my_dtu_10mw\{design_name}.htc"
new_ae_filename = rf"D:\AE EWEM MSc\T H E S I S\6_code\code repository\my_dtu_10mw\data\{design_name}_ae.dat" 

try:
    shutil.copy(htc_file_pathname, new_htc_filename)
    shutil.copy(ae_file_pathname, new_ae_filename)
    print("Copied original .htc and ae.dat files successfully.")
except:
    print("")