import run_all_scripts
import os

if __name__ == "__main__":

    # Flags for running
    versions = {'run_v1' : True, 'run_v2n' : True, 'run_v2o' : True}
    
    dynamics_args = {'static' : False, 
                     'x_sin' : True, 
                     'from_file' : False, 
                     'period' : 10, 
                     'A' : 10, # Amplitude of driving funtion
                     'axis' : 0} # Axis of oscillation x: 0, y: 1, z: 2

    run_args = {'debug' : False, 
                'simulate' : True,
                'plot' : True,
                'rootname' : '' ,
                'extension' : '', 
                'path' : "MooringTest/", 
                'tMax' : 300.0,  # simulation duration (s)
                'dof' : 3} # Size of X and XD vectors: 3 DOF for lines, points, connections, 6 DOF for bodies and rods. Ex for three points, size should be 9. 

    plot_args = {}
    if run_args['plot']:
        plot_args = {'display': False, 
                     'save': True,    
                     'line_rmse': False, 
                     'ten_rmse': False, 
                     'plot_ten': False, 
                     'lines_and_tens': True,
                     'animate_all': False,
                     'animate_start_end': False,
                     'plot_individual_start_end': False,
                     'plot_all_start_end': True,
                     'plot3d': False,
                     'plot2d': True,
                     'from_saved_runs': False,
                     'outputs_dir' : 'saved_runs/',
                     'plot_tRange': (200,250)}

    # Read in a list of files in MooringTest
    os.system('OS_scripts/clean_outputs')
    file_list = os.listdir(run_args['path'])
    file_list.remove('.DS_Store')
    file_list.remove('Outputs')

    for in_file in file_list:
        if ("dev2" in in_file) or ("v2new" in in_file) or ("v2old" in in_file) or ("_mod" in in_file):
            file_list.remove(in_file)

    # Notes about inputs files:
    #   case4.dat and lines.txt are 3d systems
    #   v1 calculates incorrect tensions for case4.dat () 
    
    print("Running for files: ", file_list)

    for in_file in file_list:
        fname = in_file.split('.')
        run_args['rootname'] = fname[0]
        run_args['extension'] = '.'+fname[1]
        if 'case4' in in_file or 'lines' in in_file:
            plot_args['plot2d'] = False
            plot_args['plot3d'] = True
        else:
            plot_args['plot2d'] = True
            plot_args['plot3d'] = False

        instance = run_all_scripts.run_infile(plot_args, dynamics_args, versions)
        instance.run(run_args = run_args)
        print("Sucessfully run for ", fname[0]+fname[1])
        print("----------------------------------------------")