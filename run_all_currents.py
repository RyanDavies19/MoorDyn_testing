import MD_drivers
import numpy as np
import os

if __name__ == "__main__":

    # Flags for running
    versions = {'run_v1' : False, 'run_cpy' : False, 'run_c' : True, 'run_f' : False}
    
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
    run_args = {'simulate' : True,
                'plot' : True,
                'del_logs' : False,
                'rootname' : 'catenary_riser_small',
                'extension' : '.txt', 
                'path' : 'MooringTest/', 
                'tMax' : 60,  # simulation duration (s)
                'dof' : 3, # DOF of coupled objects: 3 DOF for lines, points, connections, 6 DOF for bodies and rods (for no coupled objects, set to 3). TODO: make this work for coupled bodies and rods at the same time
                'dt_scaling' : 1, # Scaling factor from dtM to dtC
                'printing' : 0 # verbosity for driver script printing. Recommended is 1
                } 
    # Required if plot is true in run_args
    plot_args = {}
    if run_args['plot']: # for plotting moordyn outputs
        plot_args = {'display': False, 
                     'save': True,    
                     'plot_channels': True, 
                     'animate_all': False,
                     'animate_start_end': False,
                     'plot_individual_start_end': False,
                     'plot_all_start_end': False,
                     'plot3d': False,
                     'plot2d': False,
                     'from_saved_runs': False,
                     'outputs_dir' : 'MooringTest/',
                     'v1' : False,
                     'one_dataset' : False,
                     'plot_tRange': 'All'} # needs to be double
        

    currents = np.arange(0.1, 5.0, 0.1) # 0.1 m/s increments. Current must be greater than 0

    wtr_depth = 10
    z = np.linspace(-(wtr_depth+1), 0.0, num=50)

    root = run_args['rootname']

    for current in currents:
        os.system("rm MooringTest/current_profile.txt")
        os.system("cp "+run_args['path']+root+run_args['extension']+" "+run_args['path']+root+"_"+f"{current:.1f}"+run_args['extension'])
        run_args['rootname'] = root+"_"+f"{current:.1f}"
        with open('MooringTest/current_profile.txt', 'w') as f:
            f.write('--------------------- MoorDyn steady currents File ----------------------------------\n')
            f.write('Tabulated file with the water currents components\n')
            f.write('z (m), ux (m/s), uy (m/s), uz (m/s)\n')
            for zz in z:
                f.write(str(zz) + " " + str(current) + " 0 0\n")

        instance = MD_drivers.run_infile(plot_args, dynamics_args, versions)
        instance.run(run_args = run_args)
        print("Sucessfully run for",current,"m/s current")
        print("----------------------------------------------")
        os.system("rm "+run_args['path']+root+"_"+f"{current:.1f}"+run_args['extension'])
