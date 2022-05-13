import h2s

design_name = "flattened"
U = 8 # Tested velocity [m/s]

blade_matrix = h2s.hawc2s_blade_to_geo(design_name, save=False)

h2s.plot_c_and_theta(blade_matrix)

h2s.run_hawc2s(design_name)

h2s.pp_hawc2s_ind(design_name, U=U)

h2s.pp_hawc2s_pwr(design_name)

