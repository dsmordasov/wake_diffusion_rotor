"""Short script to turn HAWC2S blade files into a PyWakeEllipSys blade file."""

import glob
import numpy as np

if __name__ == "__main__":
    # Make sure that there is only ONE ae.dat and .htc file in the given folder
    ae_file_pathname = glob.glob("**/*ae.dat", recursive=True)[0]
    htc_file_pathname = glob.glob("**/*.htc", recursive=True)[0]
    
    ae_data = np.loadtxt(ae_file_pathname, skiprows=2, usecols=(0, 1, 2))
    radius = ae_data[:, 0]  # Span [m]
    chord = ae_data[:, 1]  # Chord length [m]
    rel_thickness = ae_data[:, 2]  # Relative thickness (%)
    
    N = np.size(radius)
    
    twist = np.genfromtxt(
        htc_file_pathname,
        # --------------------------------------------------
        skip_header=110,  # This line may need to be changed.
        # --------------------------------------------------
        max_rows=27,
        usecols=(4, 5),
        comments=";",
    )
    
    geo_mat = np.zeros([N, 4])
    geo_mat[:, 0] = radius
    geo_mat[1:, 1] = np.interp(geo_mat[1:, 0], twist[:, 0], twist[:, 1])
    geo_mat[:, 2] = chord
    geo_mat[:, 3] = rel_thickness
    
    header = " # Number of positions, cols: r [m]  twist [g]  chord [m] rel_thickness [%]"
    
    np.savetxt("blade.geo", geo_mat, comments=str(N), header=header)