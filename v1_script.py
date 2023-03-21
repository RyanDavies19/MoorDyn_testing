"""
Workin process. Need to translate input file Mooring/lines.txt from v2 version to v1 version. Also need to check data types into functions 
"""



import numpy as np
import math as m
from ctypes import *

def run (xp, xdp, tMax, dtC_py, vector_size, dylib = None): #only works for constant x and xd inputs
    #Double vector pointer data type
    double_p = POINTER(c_double)
    # -------------------- load the MoorDyn DLL ---------------------

    # Make MoorDyn function prototypes and parameter lists (remember, first entry is return type, rest are args)
    MDInitProto = CFUNCTYPE(c_int, c_double*vector_size, c_double*vector_size) #need to add filename option here, maybe this c_char works? #need to determine char size 
    MDStepProto = CFUNCTYPE(c_int, c_double*vector_size, c_double*vector_size, c_double*vector_size, double_p, double_p)
    MDClosProto = CFUNCTYPE(c_int)

    MDInitParams = (1, "x"), (1, "xd")
    MDStepParams = (1, "x"), (1, "xd"), (2, "f"), (1, "t"), (1, "dtC") 

    if dylib != None:
        dylib_path = dylib
    else:
        dylib_path = "compileDYLIB/MoorDyn.dylib"
    MDdylib = CDLL(dylib_path) #load moordyn dylib

    MDInit = MDInitProto(("LinesInit", MDdylib), MDInitParams)
    MDStep = MDStepProto(("LinesCalc", MDdylib), MDStepParams)
    MDClose= MDClosProto(("LinesClose", MDdylib))
    # ------------------------ run MoorDyn ---------------------------
    print("==================================================")    
    # initialize some arrays for communicating with MoorDyn
    t  = double_p()    # pointer to t
    dt = double_p()     # pointer to dt
    xold = np.zeros(vector_size)  # for storing previous positions

    ts_py = np.arange(0,tMax,dtC_py) 

    # Converting to ctypes
    dtC = pointer(c_double(dtC_py))

    # initialize MoorDyn at origin
    MDInit(xp,xdp)

    # loop through coupling time steps
    print("MoorDyn initialized - now performing calls to MoorDynStep...")
    for i in range(len(ts_py)):
        
        t = pointer(c_double(ts_py[i]))
        dt = dtC
        MDStep(xp, xdp, t, dt)    #force value returned here in array    
    print("Succesffuly simulated for {} seconds - now closing MoorDyn...".format(tMax))  
        

    # close MoorDyn simulation (clean up the internal memory, hopefully) when finished
    MDClose()   
    print("v1 script executed successfully")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++")
    del MDdylib


if __name__ == "__main__":
    #Double vector pointer data type
    double_p = POINTER(c_double)

    #Size of X and XD vectors: 3 DOF for lines, points, connections, 6 DOF for bodies and rods. 
    vector_size = int(6)

    # -------------------- load the MoorDyn DLL ---------------------

    # Make MoorDyn function prototypes and parameter lists (remember, first entry is return type, rest are args)
    MDInitProto = CFUNCTYPE(c_int, c_double*vector_size, c_double*vector_size) #need to add filename option here, maybe this c_char works? #need to determine char size 
    MDStepProto = CFUNCTYPE(c_int, c_double*vector_size, c_double*vector_size, c_double*vector_size, double_p, double_p)
    MDClosProto = CFUNCTYPE(c_int)

    MDInitParams = (1, "x"), (1, "xd")
    MDStepParams = (1, "x"), (1, "xd"), (2, "f"), (1, "t"), (1, "dtC") 

    MDdylib = CDLL("v1_DYLIB/MoorDyn.dylib") #load moordyn dylib


    MDInit = MDInitProto(("LinesInit", MDdylib), MDInitParams)
    MDStep = MDStepProto(("LinesCalc", MDdylib), MDStepParams)
    MDClose= MDClosProto(("LinesClose", MDdylib))

    # ------------------------ run MoorDyn ---------------------------
    print("==================================================")    
    # initialize some arrays for communicating with MoorDyn
    t  = double_p()    # pointer to t
    dt = double_p()     # pointer to dt
    x  = (c_double*vector_size)()
    xd = (c_double*vector_size)()
    xold = np.zeros(vector_size)  # for storing previous positions

    # # For circular platform motion, not working
    # r=30.0
    # w=(2*3.14159)/250

    # lines.txt fairlead locations
    x[0] = 0.0
    x[1] = 0.0
    x[2] = 0.0
    x[3] = 0.0
    x[4] = 0.0
    x[5] = 0.0

    for i in range(0,len(xd)):
        xd[i] = 0

    # parameters
    dtC_py = 0.02   # coupling time step size (s)
    tMax = 60.0  # simulation duration (s)
    ts_py = np.arange(0,tMax,dtC_py) #np.array of 30 seconds, 0.02 s timesteps

    # Converting to ctypes
    dtC = pointer(c_double(dtC_py))

    # initialize MoorDyn at origin
    MDInit(x,xd)

    # loop through coupling time steps
    print("MoorDyn initialized - now performing calls to MoorDynStep...")
    for i in range(len(ts_py)):
        
        t = pointer(c_double(ts_py[i]))
        dt = dtC
        
        # # update position vector here (keeping at zero for now)
        # x[0] = c_double(r*m.cos(w*ts_py[i]))
        # x[1] = c_double(r*m.sin(w*ts_py[i]))



            
        # # calculate velocities using finite difference
        # if i==0:
        #     for j in range(vector_size):
        #         xd  [j] = 0.0
        #         xold[j] = 0.0
        #     del j    
        # else:
        #     for j in range(vector_size):
        #         xd  [j] = (x[j] - xold[j])/dtC_py
        #         xold[j] =  x[j]
        #     del j
                    
        # call the MoorDyn step function
        f = MDStep(x, xd, t, dt)    #force value returned here in array    
    print("Succesffuly simulated for {} seconds - now closing MoorDyn...".format(tMax))  
        

    # close MoorDyn simulation (clean up the internal memory, hopefully) when finished
    MDClose()   
    print("v1 script executed successfully")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++")
