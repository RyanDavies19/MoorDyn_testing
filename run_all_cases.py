import run_all_scripts
import os

if __name__ == "__main__":

    # Flags for running
    versions = {'run_v1' : False, 'run_dev2' : False, 'run_cpy' : True, 'run_c' : False, 'run_f' : True}
    
    dynamics_args = {'static' : True, 
                     'sin' : False, 
                     'from_file' : False, 
                     'period' : 10, 
                     'A' : 10, # Amplitude of driving funtion
                     'axis' : 0} # Axis of oscillation x: 0, y: 1, z: 2


    run_args = {'debug' : False, 
                'simulate' : True,
                'plot' : True,
                'rootname' : '' ,
                'extension' : '', 
                'path' : 'MooringTest/', 
                'tMax' : 10.0,  # simulation duration (s)
                'dof' : 3} # Size of X and XD vectors: 3 DOF for lines, points, connections, 6 DOF for bodies and rods. Ex for three points, size should be 9. 

    plot_args = {}
    if run_args['plot']:
        plot_args = {'display': True, 
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
                     'plot_tRange': 'All'}

    # Read in a list of files in MooringTest
    os.system('OS_scripts/clean_outputs')
    file_list = os.listdir(run_args['path'])
    if '.DS_Store' in file_list:
        file_list.remove('.DS_Store')

    for in_file in file_list:
        if ("dev2" in in_file) or ("Cpy" in in_file) or ("C" in in_file) or ("_mod" in in_file) or ("motions" in in_file):
            file_list.remove(in_file)
        if ("cable" in in_file):
            file_list.remove(in_file)
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