import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import os
import glob

import config

def hawc2s_blade_to_geo(save=True):
    """Turn HAWC2S blade files into a PyWakeEllipSys blade file."""

    # Make sure that there is only ONE ae.dat and .htc file in the given folder
    ae_file_pathname = glob.glob("**/*ae.dat", recursive=True)[0]
    htc_file_pathname = glob.glob("**/*.htc", recursive=True)[0]
    
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
    
    # Create and save the final .geo blade format matrix
    geo_mat = np.zeros([N, 4])
    geo_mat[:, 0] = radius
    geo_mat[1:, 1] = np.interp(geo_mat[1:, 0], twist[:, 0], twist[:, 1])
    geo_mat[:, 2] = chord
    geo_mat[:, 3] = rel_thickness
    
    header = " # Number of positions, cols: r [m]  twist [g]  chord [m] rel_thickness [%]"
    
    if save:
        np.savetxt("blade.geo", geo_mat, comments=str(N), header=header)
        print(".geo blade file created!")
    
    return geo_mat


def plot_c_and_theta(geo_mat):
    """Plot chord and twist distributions of the given .geo matrix"""
    R = geo_mat[:, 0][-1] # Maximum radius [m]
    radius_nd = geo_mat[:, 0] / R # Non-dimensionalised radius r/R [-]
    
    fig, ax = plt.subplots(2, 1, sharex=True)
    ax[0].plot(radius_nd, geo_mat[:, 2], label="Chord")
    ax[0].set_ylabel("Chord [m]")
    
    ax[1].plot(radius_nd, geo_mat[:, 1], label ="Twist")
    ax[1].set_ylabel("Twist [deg]")
    ax[1].set_xlabel("Blade radius r/R [-]")
    
    
def run_hawc2s(design_name):
    """Runs a HAWC2S simulation for an .htc input file with design_name in its filename"""
    design_name = 'flattened'

    # Reminder to self: combining the backlashes and quotation marks and everything
    # in order to make runnable CMD commands is such pain
    hawc2s_path = '\"D:\AE EWEM MSc\T H E S I S\\6_code\code repository\hawc2s_hawcstab2_2.15_x64\HAWC2S_x64.exe\"'
    
    cwd = os.getcwd()
    if 'my_dtu_10mw' not in cwd:
        os.chdir('my_dtu_10mw')
    
    
    htc_path = glob.glob(f"*{design_name}*.htc", recursive=True)[0]
    if not htc_path:
        ".htc file not found!"
    htc_path = f'"{htc_path}\"'
    
    
    stream = os.popen(fr'cmd /k "{hawc2s_path} {htc_path}"')
    stream = stream.read()
    print(f"CMD OUTPUT FOR THE '{design_name}' DESIGN\n-------------------------------------\n{stream}")
    