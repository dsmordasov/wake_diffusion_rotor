"""
This code was made for the purposes of studying the effect of varying 
the grid resolution, along with the residuals limit for PyWakeEllipSys.
"""
import os

import numpy as np
import config

from run_pwes_simulation import run_pwes_simulation

parent_working_directory = os.getcwd()

# Choose whether you are running grid resolution or residuals simulations
#parameter_name = "grid"
#parameter_variants = [config.gs_1_res_grid, config.gs_2_res_grid, config.gs_3_res_grid, config.gs_4_res_grid, config.gs_5_res_grid, ]
parameter_name = "res"
parameter_variants = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
n_of_runs = len(parameter_variants)

results = np.empty([n_of_runs, 4])

if __name__ == '__main__':
    run_number = 0
    master_folder_name = f"{parameter_name}_study_results"
    os.mkdir(master_folder_name, 0o777)
    os.chdir(master_folder_name)
    try:
        for parameter in parameter_variants:
            design_name = parameter_name + str(run_number)
            print(
                f"Running for {parameter_name}={run_number + 1}/{n_of_runs}:")
            parameter_folder_name = f"{parameter_name}_{run_number}"
            os.mkdir(parameter_folder_name, 0o777)
            os.chdir(parameter_folder_name)
            # ---RUN---
            #aep, advd, runtime = run_pwes_simulation(grid=parameter, design_name=design_name)
            aep, advd, runtime = run_pwes_simulation(
                e3d_reslim=parameter, design_name=design_name)
            print(
                f"{parameter_name} fineness: {run_number}, AEP: {aep}, runtime: {runtime}")
            results[run_number, :] = run_number, aep, advd, runtime
            run_number += 1
            # ---RUN---
            os.chdir(os.path.join(parent_working_directory, master_folder_name))
    except Exception as e:
        print(f"Exception occured:\n{e}")
        os.chdir(parent_working_directory)
    os.chdir(parent_working_directory)

print(results)

# CLEANUP:
# shutil.rmtree(master_folder_name)
