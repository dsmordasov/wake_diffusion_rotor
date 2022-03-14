import h2s

blade_matrix = h2s.hawc2s_blade_to_geo(save=False)

h2s.plot_c_and_theta(blade_matrix)

h2s.run_hawc2s("flattened")

# Run HAWC2S

# Plot the thrust curve, save the power 

