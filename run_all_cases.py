import run_all_scripts
import os

if __name__ == "__main__":

    # Flags for running
    debug = False 
    run_v1 = True
    simulate = True
    plot = True
    
    # Note: This only apply if plot boolean is true
    # Flags for plotting: can only plot indidvual or all, not both. If both are true, only individual runs
    fig_options = {'animate_all' : False, 'animate_start_end' : False, 'plot2d' : True, 'plot3d' : True, 'plot_individual_start_end' : True, 'plot_all_start_end' : True, 'display' : False, 'save_fig' : True}
    
    path = "MooringTest/" 

    tMax = 600.0  # simulation duration (s)
    
    dof = 3 # Size of X and XD vectors: 3 DOF for lines, points, connections, 6 DOF for bodies and rods. Ex for three points, size should be 9. 
    
    # Read in a list of files in MooringTest
    file_list = os.listdir(path)
    file_list.remove('.DS_Store')
    file_list.remove('Outputs')

    # Notes about inputs files:
    #   case4.dat segfaults on v1 and cannot be plotted 2d
    #   lines.txt cannot be plotted 2d


    for in_file in file_list:
        fname = in_file.split('.')
        rootname = fname[0]
        extension = '.'+fname[1]
        instance = run_all_scripts.run_infile()
        if rootname == 'case4':
            fig_options['plot2d'] = False
            instance.run(rootname = rootname, extension = extension, path = path, tMax = tMax, dof = dof, debug = debug, run_v1 = False, simulate = simulate, plot = plot, fig_options = fig_options)
            fig_options['plot2d'] = True        
        elif rootname == 'lines':
            fig_options['plot2d'] = False
            instance.run(rootname = rootname, extension = extension, path = path, tMax = tMax, dof = dof, debug = debug, run_v1 = run_v1, simulate = simulate, plot = plot, fig_options = fig_options)
            fig_options['plot2d'] = True
        else: 
            instance.run(rootname = rootname, extension = extension, path = path, tMax = tMax, dof = dof, debug = debug, run_v1 = run_v1, simulate = simulate, plot = plot, fig_options = fig_options)

        print("Sucessfully run for ", rootname+extension)
        print("----------------------------------------------")