import glob

import numpy as np

# Make sure that there is only ONE ae.dat and .htc file in the given folder
ae_file_pathname = glob.glob("**/*ae.dat", recursive=True)[0]
htc_file_pathname = glob.glob("**/*.htc", recursive=True)[0]

ae_data = np.loadtxt(ae_file_pathname, skiprows=2, usecols=(0, 1, 2))
radius = ae_data[:, 0]  # Span [m]
chord = ae_data[:, 1]  # Chord length [m]
rel_thickness = ae_data[:, 2]  # Relative thickness (%)

twist = np.genfromtxt(
    htc_file_pathname,
    # --------------------------------------------------
    skip_header=110,  # This line may need to be changed.
    # --------------------------------------------------
    max_rows=27,
    usecols=(4, 5),
    comments=";",
)

