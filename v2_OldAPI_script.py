"""
This is working without error, however output file does not contain fairlead tensions. Likely input file 
"""


import numpy as np
from ctypes import *

def run (x, xd, tMax, dtC_py, vector_size, filename, dylib = None): #only works for constant x and xd inputs
    print("==================================================")
    # -------------------- load the MoorDyn DLL ---------------------

    #Double vector pointer data type
    double_p = POINTER(c_double)

    # Make MoorDyn function prototypes and parameter lists (remember, first entry is return type, rest are args)
    MDInitProto = CFUNCTYPE(c_int, c_double*vector_size, c_double*vector_size, c_char_p) #need to add filename option here, maybe this c_char works? #need to determine char size 
    MDStepProto = CFUNCTYPE(c_int, c_double*vector_size, c_double*vector_size, c_double*vector_size, double_p, double_p)
    MDClosProto = CFUNCTYPE(c_int)

    MDInitParams = (1, "x"), (1, "xd"), (1, "infilename") #guessing the 1 flag here means input?
    MDStepParams = (1, "x"), (1, "xd"), (2, "f"), (1, "t"), (1, "dtC") 

    if dylib != None:
        dylib_path = dylib
    else:
        dylib_path = "../MoorDyn/compile/DYLIB/libmoordyn2.dylib"
    print("dylib path is ", dylib_path)
    MDdylib = CDLL(dylib_path) #load moordyn dylib


    MDInit = MDInitProto(("MoorDynInit", MDdylib), MDInitParams)
    MDStep = MDStepProto(("MoorDynStep", MDdylib), MDStepParams)
    MDClose= MDClosProto(("MoorDynClose", MDdylib))  
  # ------------------------ run MoorDyn ---------------------------
    # initialize some arrays for communicating with MoorDyn
    t  = double_p()    # pointer to t
    dt = double_p()     # pointer to dt
    xold = np.zeros(vector_size)  # for storing previous positions

    # parameters
    ts_py = np.arange(0,tMax,dtC_py) #np.array of 30 seconds, 0.02 s timesteps
    dtC = pointer(c_double(dtC_py))

    infile = c_char_p(bytes(filename, encoding='utf8'))

    # initialize MoorDyn at origin
    MDInit(x,xd,infile)

    # loop through coupling time steps
    print("MoorDyn initialized - now performing calls to MoorDynStep...")

    for i in range(len(ts_py)):
        t = pointer(c_double(ts_py[i]))
        dt = dtC
        MDStep(x, xd, t, dt)    
    print("Succesffuly simulated for {} seconds - now closing MoorDyn...".format(tMax))  

    # close MoorDyn simulation (clean up the internal memory, hopefully) when finished
    MDClose()   
    print("Old API v2 script executed successfully")
    print("++++++++++++++++++++++++++++++++++++++++++++++++++")
    del MDdylib
    

if __name__ == "__main__":
        
    #Double vector pointer data type
    double_p = POINTER(c_double)

    #Size of X and XD vectors: 3 DOF for lines, points, connections, 6 DOF for bodies and rods. Ex for three points, size should be 3. 
    vector_size = int(9)

    # initialize some arrays for communicating with MoorDyn
    t  = double_p()    # pointer to t
    dt = double_p()     # pointer to dt
    x  = (c_double*vector_size)()
    xd = (c_double*vector_size)()
    xold = np.zeros(vector_size)  # for storing previous positions

    # parameters
    dtC_py = 0.02   # coupling time step size (s)
    tMax = 30.0  # simulation duration (s)
    ts_py = np.arange(0,tMax,dtC_py) #np.array of 30 seconds, 0.02 s timesteps

    # Converting to ctypes
    dtC = pointer(c_double(dtC_py))

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

    infileroot = "MooringTest/lines.txt"


#     # # For circular platform motion, not currently working
#     # r=100
#     # w=(2*3.14159)/10

#     for i in range(len(ts_py)):
        
#         t = pointer(c_double(ts_py[i]))
#         dt = dtC
        
#         # # update position vector here (keeping at zero for now)
#         # for j in range(0,len(x)):
#         #     if j%3 == 0 :
#         #         x[j] = c_double(x[j] + r*m.cos(w*ts_py[i]))
#         #         x[j+1] = c_double(x[j+1] + r*m.sin(w*ts_py[i]))
#         # del j

            
#         # # calculate velocities using finite difference
#         # if i==0:
#         #     for j in range(vector_size):
#         #         xd  [j] = 0.0
#         #         xold[j] = 0.0
#         # else:
#         #     for j in range(vector_size):
#         #         xd  [j] = (x[j] - xold[j])/dtC_py
#         #         xold[j] =  x[j]
#         # del j

    run (x, xd, tMax, dtC_py, vector_size, infileroot, dylib = None)