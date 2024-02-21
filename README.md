Index of options for using MD_drivers.py:

# Flags for running (required)
Versions = 
'run_v1' : False 
'run_cpy' : True
'run_c' : True
'run_f' : True

# Required
dynamics_args = 
'static' : False, 
'sin' : True, 
'from_file' : False, 
 # if sin is true
'period' : 5, 
'A' : 1, 
'axis' : 4 # 0 -> x, 1 -> y, 2 -> z, 3 -> rx , 4 -> ry, 5 -> rz

# Required
run_args = 
'simulate' : True,
'plot' : True,
'del_logs' : False,
'rootname' : 'USFLOWT_updated',
'extension' : '.dat', 
'path' : 'MooringTest/', 
'tMax' : 2,  # simulation duration (s)
'dof' : 6, # DOF of coupled objects: 3 DOF for lines, points, connections, 6 DOF for bodies and rods (for no coupled objects, set to 3). TODO: make this work for coupled bodies and rods at the same time
'dt_scaling' : 1, # Scaling factor from dtM to dtC
'printing' : 1 # verbosity for driver script printing. Recommended is 1

# Required if plot is true in run_args
plot_args = 
'display': True, 
'save': False,    
'plot_channels': True,
'animate_all': False,
'animate_start_end': False,
'plot_individual_start_end': True,
'plot_all_start_end': True,
'plot3d': True,
'plot2d': True,
'from_saved_runs': False,
'outputs_dir' : None,
'v1' : False,
'one_dataset' : False,
'plot_tRange': 'All' # needs to be [double, double] or 'All'