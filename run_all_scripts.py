import MD_drivers

#------------------- Run All Scripts -----------------------------

# Flags for running (required)
versions = {'run_v1' : False, 'run_cpy' : False, 'run_c' : True, 'run_f' : False}

# Required
dynamics_args = {'static' : True, 
                    'sin' : False, 
                    'from_file' : False, 
                    # if sin is true
                    'period' : 5, 
                    'A' : 1, 
                    'axis' : 4 # 0 -> x, 1 -> y, 2 -> z, 3 -> rx , 4 -> ry, 5 -> rz
                    }

# Required
run_args = {'simulate' : True,
            'plot' : True,
            'del_logs' : False,
            'rootname' : 'viv_cable',
            'extension' : '.txt', 
            'path' : 'MooringTest/', 
            'tMax' : 20,  # simulation duration (s)
            'dof' : 3, # DOF of coupled objects: 3 DOF for lines, points, connections, 6 DOF for bodies and rods (for no coupled objects, set to 3). TODO: make this work for coupled bodies and rods at the same time
            'dt_scaling' : 1, # Scaling factor from dtM to dtC
            'printing' : 1 # verbosity for driver script printing. Recommended is 1
            } 

# Required if plot is true in run_args
plot_args = {}
if run_args['plot']: # for plotting moordyn outputs
    plot_args = {'display': True, 
                    'save': False,    
                    'plot_channels': True,
                    'animate_all': False,
                    'animate_start_end': False,
                    'plot_individual_start_end': False,
                    'plot_all_start_end': True,
                    'plot3d': True,
                    'plot2d': True,
                    'from_saved_runs': False,
                    'outputs_dir' : None,
                    'v1' : False,
                    'one_dataset' : False,
                    'plot_tRange': 'All' # needs to be [double, double] or 'All'
                    } 

instance = MD_drivers.run_infile(plot_args, dynamics_args, versions)
instance.run(run_args = run_args)