import run_all_scripts
import os

if __name__ == "__main__":

    # Flags for running
    versions = {'run_v1' : False, 'run_dev2' : False, 'run_cpy' : False, 'run_c' : True, 'run_f' : True}
    
    # Required
    dynamics_args = {'static' : True, 
                     'sin' : False, 
                     'from_file' : False, 
                     # if sin is true
                     'period' : 10, 
                     'A' : 6, 
                     'axis' : 5 # 0 -> x, 1 -> y, 2 -> z, 3 -> rx , 4 -> ry, 5 -> rz
                     }

    # Required
    run_args = {'debug' : False, 
                'simulate' : True,
                'plot' : True,
                'rootname' : '' ,
                'extension' : '', 
                'path' : 'MooringTest/', 
                'tMax' : 20,  # simulation duration (s)
                'dof' : 3} # DOF of coupled objects: 3 DOF for lines, points, connections, 6 DOF for bodies and rods (for no coupled objects, set to 3). TODO: make this work for coupled bodies and rods at the same time

    # Required if plot is true in run_args
    plot_args = {}
    if run_args['plot']: # for plotting moordyn outputs
        plot_args = {'display': False, 
                     'save': True,    
                     'line_rmse': False, 
                     'ten_rmse': False, 
                     'plot_channels': True,
                     'lines_and_tens': False, 
                     'animate_all': False,
                     'animate_start_end': False,
                     'plot_individual_start_end': False,
                     'plot_all_start_end': True,
                     'plot3d': True,
                     'plot2d': False,
                     'from_saved_runs': False,
                     'outputs_dir' : None,
                     'v1' : False,
                     'one_dataset' : False,
                     'plot_tRange': 'All'} # needs to be double

    # Read in a list of files in MooringTest
    os.system('OS_scripts/clean_outputs')
    file_list = os.listdir(run_args['path'])
    if '.DS_Store' in file_list:
        file_list.remove('.DS_Store')

    for in_file in file_list:
        if ("dev2" in in_file) or ("Cpy" in in_file) or ("C" in in_file) or ("_mod" in in_file) or ("motions" in in_file):
            file_list.remove(in_file)
        # if ("cable" in in_file):
        #     file_list.remove(in_file)
    # Notes about inputs files:
    #   case4.dat and lines.txt are 3d systems
    #   v1 calculates incorrect tensions for case4.dat () 
    
    print("Running for files: ", file_list)

    for in_file in file_list:
        fname = in_file.split('.')
        run_args['rootname'] = fname[0]
        run_args['extension'] = '.'+fname[1]
        # if '2D' in in_file:
        #     plot_args['plot2d'] = True
        #     plot_args['plot3d'] = False
        # else:
        #     plot_args['plot2d'] = False
        #     plot_args['plot3d'] = True

        instance = run_all_scripts.run_infile(plot_args, dynamics_args, versions)
        instance.run(run_args = run_args)
        print("Sucessfully run for ", fname[0]+'.'+fname[1])
        print("----------------------------------------------")