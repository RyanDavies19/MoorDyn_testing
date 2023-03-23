"""
This is working without error, however this is pre vtk updates to MoorDyn. Python library no longer compiles locally 
"""


import moordyn
import numpy as np

def run (x_in, xd_in, tMax, dtC_py, vector_size, infilename, dylib = None): #only works for constant x and xd inputs

    print("==================================================")
    print("This runs the python wrapper of MoorDynV2, it does not reference the local MoorDyn copy that is being edited in ../MoorDyn")
    system = moordyn.Create(infilename)
    x = list(x_in)
    xd = list(xd_in)
    xold = np.zeros(vector_size)
    dx = [0]*vector_size
    ts = np.arange(0,tMax,dtC_py)
    moordyn.Init(system, x, dx)
    # loop through coupling time steps
    print("MoorDyn initialized - now performing calls to MoorDynStep...")
    for i in range(len(ts)):
        t = ts[i]
        dt = dtC_py
        
        # # # update position vector here (keeping at zero for now)
        # for j in range(vector_size):
        #     x[j] = 0.0
        #     if j == (vector_size-1):
        #         del j
            
        # # calculate velocities using finite difference
        # if i==0:
        #     for j in range(6):
        #         xd  [j] = 0.0
        #         xold[j] = 0.0
        # else:
        #     for j in range(6):
        #         xd  [j] = (x[j] - xold[j])/dtC
        #         xold[j] =  x[j]
                    
        # call the MoorDyn step function
        moordyn.Step(system, x, xd, t, dt)    #force value returned here in array

        
    print("Succesffuly simulated for {} seconds - now closing MoorDyn...".format(tMax))  
        

    # close MoorDyn simulation (clean up the internal memory, hopefully) when finished
    moordyn.Close(system)   

    print("New API v2 script executed successfully")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++")
    del system

if __name__ == "__main__":

    # Creating MoorDyn Instance with v2 API.
    infilename = "MooringTest/lines.txt" #Debugging
    #Size of X and XD vectors: 3 DOF for lines, points, connections, 6 DOF for bodies and rods. Ex for three points, size should be 3. 
    vector_size = int(9)
    print("==================================================")
    print("This runs the python wrapper of MoorDynV2, it does not reference the local MoorDyn copy that is being edited in ../MoorDyn")

    system = moordyn.Create(infilename)
    #Initializing MoorDyn

    dx = [0]*vector_size
    x = [0]*vector_size
    xd = [0]*vector_size
    xold = np.zeros(vector_size)  # for storing previous positions

    # lines.txt fairlead locations
    x[0] = 5.2
    x[1] = 0.0
    x[2] = -70.0
    x[3] = -2.6
    x[4] = 4.5
    x[5] = -70.0
    x[6] = -2.6
    x[7] = -4.5
    x[8] = -70.0

    for i in range(0,len(xd)):
        xd[i] = 0

    # parameters
    dtC = 0.02   # coupling time step size (s)
    tMax = 60.0  # simulation duration (s)
    ts = np.arange(0,tMax,dtC) #np.array of 30 seconds, 0.02 s timesteps

    moordyn.Init(system, x, dx)

    # loop through coupling time steps
    print("MoorDyn initialized - now performing calls to MoorDynStep...")
    for i in range(len(ts)):
        t = ts[i]
        dt = dtC
        
        # # # update position vector here (keeping at zero for now)
        # for j in range(vector_size):
        #     x[j] = 0.0
        #     if j == (vector_size-1):
        #         del j
            
        # # calculate velocities using finite difference
        # if i==0:
        #     for j in range(6):
        #         xd  [j] = 0.0
        #         xold[j] = 0.0
        # else:
        #     for j in range(6):
        #         xd  [j] = (x[j] - xold[j])/dtC
        #         xold[j] =  x[j]
                    
        # call the MoorDyn step function
        f = moordyn.Step(system, x, xd, t, dt)    #force value returned here in array

        
    print("Succesffuly simulated for {} seconds - now closing MoorDyn...".format(tMax))  
        

    # close MoorDyn simulation (clean up the internal memory, hopefully) when finished
    moordyn.Close(system)   

    print("New API v2 script executed successfully")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++")
