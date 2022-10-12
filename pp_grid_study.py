"""
This code was made in order to post process the results of the grid study.
Excuse the heavy hardcoding, as this was performed very close to the thesis 
deadline.
"""
import numpy as np
import config

import matplotlib.pyplot as plt

graph_marker = '-.D'
save_figs = True

# cells per diameter [-], power [MW], AD velocity deficit [-]
wake_grid_results = np.array([[2, 4.674398, 0.691722747],
                              [4, 4.380171, 0.693197909],
                              [8, 4.154834, 0.713697409],
                              [16, 4.074896, 0.726278678],
                              [20, 4.064092, 0.728962293]])

# x cells * x cells [-], power [MW], AD velocity deficit [-]
ad_grid_results = np.array([[8, 3.791462, 0.741725320],
                            [16, 4.019477, 0.729535184],
                            [32, 4.060197, 0.726818417],
                            [64, 4.073563, 0.726308684],
                            [128, 4.074896, 0.726278678]])

# resolution limit 1e-x [-], power [MW], AD velocity deficit [-]
res_results = np.array([[2, 4.059669, 0.972886074],
                        [3, 9.977986, 0.693827555],
                        [4, 4.066159, 0.722025149],
                        [5, 4.073563, 0.726308684],
                        [6, 4.073977, 0.726419348]])

res_x_ticks = res_results[:,0]

# Compute relative errors
wake_grid_power = (wake_grid_results[:, 1] - wake_grid_results[-1, 1]) / wake_grid_results[-1, 1]
wake_grid_advd = (wake_grid_results[:, 2] - wake_grid_results[-1, 2]) / wake_grid_results[-1, 2]
ad_grid_power = (ad_grid_results[:, 1] - ad_grid_results[-1, 1]) / ad_grid_results[-1, 1]
ad_grid_advd = (ad_grid_results[:, 2] - ad_grid_results[-1, 2]) / ad_grid_results[-1, 2]
res_power = (res_results[:, 1] - res_results[-1, 1]) / res_results[-1, 1]
res_advd = (res_results[:, 2] - res_results[-1, 2]) / res_results[-1, 2]


fig, ax = plt.subplots(2, 1, sharex=True, figsize=[10,6])
ax[0].plot(wake_grid_results[:, 0], wake_grid_power, graph_marker)
ax[1].plot(wake_grid_results[:, 0], wake_grid_advd, graph_marker)

ax[0].set_ylim([-0.02, 0.16])
ax[1].set_ylim([-0.06, 0.01])
plt.xticks(wake_grid_results[:, 0])

if save_figs:
    fig.savefig('plots/wake_grid_relel.pdf', bbox_inches='tight')

fig, ax = plt.subplots(2, 1, sharex=True, figsize=[10,6])
ax[0].plot(ad_grid_results[:, 0], ad_grid_power, graph_marker)
ax[1].plot(ad_grid_results[:, 0], ad_grid_advd, graph_marker)

ax[0].set_ylim([-0.07, 0.01])
ax[1].set_ylim([-0.005, 0.025])
plt.xticks(ad_grid_results[:, 0])

if save_figs:
    fig.savefig('plots/ad_grid_relel.pdf', bbox_inches='tight')

fig, ax = plt.subplots(2, 1, sharex=True, figsize=[10,6])
ax[0].plot(res_x_ticks, res_power, graph_marker)
ax[1].plot(res_x_ticks, res_advd, graph_marker)

plt.xticks(res_x_ticks)

if save_figs:
    fig.savefig('plots/res_relel.pdf', bbox_inches='tight')