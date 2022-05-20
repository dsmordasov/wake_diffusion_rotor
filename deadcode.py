#%% local_pp_delta_study.py
# momentum transfer graphs - velocity differentials

# Hardcoding for momentum transport analysis
mt_option = True
differentials_cmap = cm.RdBu

def pp_mt(sigma_data, axis_counter):
    ''' Post-process momentum transfer. 
    '''
    current_label = r"$\hat{\delta}$ = " + str(analysed_deltas[axis_counter])
    
    x_limit = (-0.5, 3)
    y_limit = (-0.7, 0.7)
    #y_limit_dWdz = (-0.6, 0.6)
    #z_limit = () # TODO: Normalize z with zh
    
    # Hardcoded colorbar minima and maxima
    vmin_left = -0.06
    vmax_left = 0.08
    vmin_right = -0.11
    vmax_right = 0.09
    
    # # Set up limits for the colorbar
    # global left_min
    # global left_max
    # global right_min
    # global right_max
    
    # Left: dV/dy
    zplane_data = sigma_data.interp(z=zh)
    X, Y = np.meshgrid(zplane_data.x, zplane_data.y)
    # left_current_min = np.min(zplane_data['dV/dy'])
    # left_current_max = np.max(zplane_data['dV/dy'])
    left_plot = axes[axis_counter, 0].contourf(X / D, Y / D, zplane_data['dV/dy'].T,
                                               cmap=differentials_cmap, vmin=vmin_left, vmax=vmax_left)
    axes[axis_counter, 0].text(x_limit[0] + 1e-1, y_limit[1] + 0.5e-1, current_label)
    
    
    # Right: dW/dz
    yplane_data = sigma_data.interp(y=0)
    X, Z = np.meshgrid(yplane_data.x, yplane_data.z)
    # right_current_min = np.min(yplane_data['dW/dz'])
    # right_current_max = np.max(yplane_data['dW/dz'])
    right_plot = axes[axis_counter, 1].contourf(X / D, Z / D - z_norm, yplane_data['dW/dz'].T,
                                                cmap=differentials_cmap, vmin=vmin_right, vmax=vmax_right)
    axes[axis_counter, 1].text(x_limit[0] + 1e-1, y_limit[1] + 0.5e-1, current_label)

    plt.xlim(x_limit)
    axes[axis_counter, 0].set_ylim(y_limit)
    axes[axis_counter, 1].set_ylim(y_limit)

# =============================================================================
#     # Set colorbar limits for left dV/dy and right dW/dz plots
#     if axis_counter == 0:
#         left_min = left_current_min
#         left_max = left_current_max
#         right_min = right_current_min
#         right_max = right_current_max
#     
#     # Update the minimum/maximum if current value is lower/higher
#     if left_current_min < left_min:
#         left_min = left_current_min
#     if left_current_max > left_max:
#         left_max = left_current_max
#     if right_current_min < right_min:
#         right_min = right_current_min
#     if right_current_max > right_max:
#         right_max = right_current_max
# =============================================================================
    
    return left_plot, right_plot
    #fig.colorbar(right_plot, ax=axes.ravel().tolist())
    #if (axis_counter + 1) == len(analysed_flowdata_paths):

if mt_option:
    fig, axes = plt.subplots(3, 2, sharex=True, figsize=(12, 5))
    
    for mt_iterator, filename_path in enumerate(analysed_flowdata_paths):
        print(f"Ploting velocity differentials {mt_iterator + 1}/{n_deltas}")
        
        data = xarray.open_dataset(filename_path)
        pped_data = calc_sigma(data)
        
        left_plot, right_plot = pp_mt(pped_data, mt_iterator)
        
        current_delta = analysed_deltas[mt_iterator]
        #current_label = r"$\hat{\delta}$ = " + str(current_delta)
        #plt.legend(current_label, loc="upper left", fontsize=5)
    
    
    colorbar_left = fig.colorbar(left_plot, ax=axes[:, 0])
    colorbar_right = fig.colorbar(right_plot, ax=axes[:, 1])
    
    colorbar_left.ax.set_xlabel("$\dfrac{dV}{dy}$")
    colorbar_right.ax.set_xlabel("$\dfrac{dW}{dz}$")
    
    axes[n_deltas-1, 0].set_xlabel("x [D]")
    axes[n_deltas-1, 1].set_xlabel("x [D]")
    
    axes[1, 0].set_ylabel("y [D]")
    axes[1, 1].set_ylabel("z [D]")