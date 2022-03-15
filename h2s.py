"""HAWC2S script suite/toolbox."""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import os
import glob

import config

def hawc2s_blade_to_geo(design_name, save=True):
    """Turn HAWC2S blade files into a PyWakeEllipSys blade file."""

    # Make sure that there is only ONE ae.dat file in the given folder
    ae_file_pathname = glob.glob("**/*ae.dat", recursive=True)[0] # TODO: Take `design_name` into consideration
    htc_file_pathname = glob.glob(f"**/*{design_name}*.htc", recursive=True)[0]
    
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
        os.chdir('my_dtu_10mw') # HAWC2S must be ran from the folder with the .htc file
    
    
    htc_path = glob.glob(f"*{design_name}*.htc", recursive=True)[0]
    if not htc_path:
        ".htc file not found!"
    htc_path = f'"{htc_path}\"'
    
    
    stream = os.popen(fr'cmd /k "{hawc2s_path} {htc_path}"')
    stream = stream.read()
    print(f"CMD OUTPUT FOR THE '{design_name}' DESIGN\n-------------------------------------\n{stream}")


def pp_hawc2s_ind(design_name, U=8, rho=1.225):
    """Post process HAWC2S .ind  file.
    
    HAWC2S goes bonkers if its input files include the tip of the blade,
    so for plotting purposes, this was written to add a zero tangential and
    azimuthal force at the blade tip (denoted `tip fix` in comments).
    """
    hawc2s_ind_path = glob.glob(f"**/*{design_name}_u*.ind", recursive=True)[0]
    if not hawc2s_ind_path:
        print("Found no HAWC2S .ind files!")
    
    data = np.loadtxt(hawc2s_ind_path)
    s       = data[:, 0]      # Blade span [m]
    S       = s[-1] + 0.001   # Maximum blade span (w/ tip fix) [m]
    s       = np.append(s, S) # Tip fix
    Ft      = data[:, 6]      # Tangential force per unit length [N/m]
    Ft      = np.append(Ft, 0)# Tip fix
    Fn      = data[:, 7]      # Azumuthal force per unit length [N/m]
    Fn      = np.append(Fn, 0)# Tip fix
    
    # Normalisation of forces, non-dimensionalisation of radius
    s_nd    = s / S           # Blade span, non-dimensioned [-]
    ft = Ft / (rho * S * U**2) # Normalised tangential force [-]
    fn = Fn / (rho * S * U**2) # Normalised azimuthal force [-]
    
    # Some other useful values you can grab fron the .ind file
    # C_L     = data[:, 16]     # Lift coefficient [-]
    # C_D     = data[:, 17]     # Drag coefficient [-]
    # C_L_C_D = C_L/C_D         # C_L/C_D [-]
    # C_P     = data[:, 33]     # Power coefficient [-]
    # C_T     = data[:, 32]     # Thrust coefficient [-]
    # aoa     = np.rad2deg(data[:, 4])     # Angle of attack [deg]
    
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=[10,6])
    
    ax[0].plot(s_nd, ft, label="HAWC2S")
    ax[1].plot(s_nd, fn, label="HAWC2S")
    
    ax[0].set_ylabel('ft $[-]$')
    ax[0].legend(loc=1)
    ax[1].set_ylabel('fn $[-]$')
    ax[1].set_xlabel('Blade radius r/R [-]')
    plt.tight_layout()
    
    
def pp_hawc2s_pwr(design_name):
    """Post process HAWC2S .pwr file.
    """
    
    hawc2s_pwr_path = glob.glob(f"**/*{design_name}.pwr", recursive=True)[0]
    if not hawc2s_pwr_path:
        print("Found no HAWC2S .pwr files!")
        
    data = np.loadtxt(hawc2s_pwr_path)
    power = data[1] / 1000 # [MW]
    print(f"Power: {np.round(power, 4)} MW.")
    
    