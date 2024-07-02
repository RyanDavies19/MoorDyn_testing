import MD_drivers
import os

if __name__ == "__main__":

    # Flags for running
    versions = {'run_v1' : False, 'run_cpy' : False, 'run_c' : True, 'run_f' : True}
    
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
                'path' : 'reg_tests/', 
                'tMax' : 60,  # simulation duration (s)
                'dof' : 6} # DOF of coupled objects: 3 DOF for lines, points, connections, 6 DOF for bodies and rods (for no coupled objects, set to 3). TODO: make this work for coupled bodies and rods at the same time

    # Required if plot is true in run_args
    plot_args = {}
    if run_args['plot']: # for plotting moordyn outputs
        plot_args = {'display': True, 
                     'save': True,    
                     'plot_channels': True, 
                     'animate_all': False,
                     'animate_start_end': False,
                     'plot_individual_start_end': False,
                     'plot_all_start_end': True,
                     'plot3d': True,
                     'plot2d': True,
                     'from_saved_runs': False,
                     'outputs_dir' : 'reg_tests/',
                     'v1' : False,
                     'one_dataset' : False,
                     'plot_tRange': 'All'} # needs to be double

    # Read in a list of files in MooringTest
    os.system('OS_scripts/clean_outputs > OS_output.txt')
    dir_list = os.listdir(run_args['path'])
    if '.DS_Store' in dir_list:
        dir_list.remove('.DS_Store')

    file_list = []
    for in_file in dir_list:
        if not (("Cpy" in in_file) or ("C" in in_file) or ("_mod" in in_file) or ("F" in in_file) or ("water" in in_file) or ("current" in in_file) or ("wave" in in_file) or (".py" in in_file)):
            file_list.append(in_file)
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

        # os.system("cp MooringTest/ptfm_motions_{}.dat MooringTest/ptfm_motions.dat".format(fname[0]))

        if 'point' in in_file:
            run_args['dof'] = 3
        else:
            run_args['dof'] = 6

        instance = MD_drivers.run_infile(plot_args, dynamics_args, versions)
        instance.run(run_args = run_args)
        print("Sucessfully run for ", fname[0]+'.'+fname[1])
        # os.system("rm MooringTest/ptfm_motions.dat")
        print("----------------------------------------------")