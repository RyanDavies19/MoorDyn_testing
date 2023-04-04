import run_all_scripts
import os

if __name__ == "__main__":

    run_v1 = True
    run_v2n = True
    run_v2o = True
    static = False

    debug = False
    simulate = True
    plot = True

    fig_options = {}
    if plot:
    # Note: This only apply if plot boolean is true
    # Flags for plotting: can only plot indidvual or all, not both. If both are true, only individual runs
        fig_options = {'animate_all' : False, 'animate_start_end' : False, 'plot2d' : False, 'plot3d' : False, 'plot_individual_start_end' : False, 'plot_all_start_end' : False, 'display' : False, 'save_fig' : True, 'show_rmse' : True}
    
    path = "MooringTest/"

    tMax = 600.0  # simulation duration (s)
    
    dof = 3 # Size of X and XD vectors: 3 DOF for lines, points, connections, 6 DOF for bodies and rods. Ex for three points, size should be 9. 

    # Read in a list of files in MooringTest
    file_list = os.listdir(path)
    file_list.remove('.DS_Store')
    file_list.remove('Outputs')

    for in_file in file_list:
        if ("dev2" in in_file) or ("v2new" in in_file) or ("v2old" in in_file):
            file_list.remove(in_file)

    # Notes about inputs files:
    #   case4.dat segfaults on v1 and cannot be plotted 2d
    #   lines.txt cannot be plotted 2d


    for in_file in file_list:
        fname = in_file.split('.')
        rootname = fname[0]
        extension = '.'+fname[1]
        instance = run_all_scripts.run_infile()
        if fig_options['plot2d']:
            if rootname == 'case4':
                fig_options['plot2d'] = False
                fig_options['plot3d'] = True
                instance.run(rootname = rootname, extension = extension, path = path, tMax = tMax, dof = dof, debug = debug, run_v1 = run_v1, run_v2n = run_v2n, run_v2o = run_v2o, simulate = simulate, plot = plot, fig_options = fig_options, static = static)                
                fig_options['plot2d'] = True 
                fig_options['plot3d'] = False       
            elif rootname == 'lines':
                fig_options['plot2d'] = False
                fig_options['plot3d'] = True
                instance.run(rootname = rootname, extension = extension, path = path, tMax = tMax, dof = dof, debug = debug, run_v1 = run_v1, run_v2n = run_v2n, run_v2o = run_v2o, simulate = simulate, plot = plot, fig_options = fig_options, static = static)                
                fig_options['plot2d'] = True 
                fig_options['plot3d'] = False
            else: 
                instance.run(rootname = rootname, extension = extension, path = path, tMax = tMax, dof = dof, debug = debug, run_v1 = run_v1, run_v2n = run_v2n, run_v2o = run_v2o, simulate = simulate, plot = plot, fig_options = fig_options, static = static)                   
        else: 
            instance.run(rootname = rootname, extension = extension, path = path, tMax = tMax, dof = dof, debug = debug, run_v1 = run_v1, run_v2n = run_v2n, run_v2o = run_v2o, simulate = simulate, plot = plot, fig_options = fig_options, static = static)                

        print("Sucessfully run for ", rootname+extension)
        print("----------------------------------------------")