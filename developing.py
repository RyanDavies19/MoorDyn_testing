import moordyn
import numpy as np

rootname = 'BodiesAndRods_mod2'
extension = '.dat'
path = 'MooringTest/'
tMax = 25.0
dtM = 0.001
time = np.arange(0, tMax, dtM)
size = (len(time), 6)

x = np.zeros(size)
xd = np.zeros(size)

print("==================================================")
print("This runs the python wrapper of MoorDynV2, it does not reference the local MoorDyn copy that is being edited in ../MoorDyn")
system = moordyn.Create(path+rootname+extension)
moordyn.Init(system, x[0,:], xd[0,:])
# loop through coupling time steps
print("MoorDyn initialized - now performing calls to MoorDynStep...")
for i in range(len(time)):
    # call the MoorDyn step function
    moordyn.Step(system, x[i,:], xd[i,:], time[i], dtM)    #force value returned here in array

print("Successfuly simulated for {} seconds - now closing MoorDyn...".format(tMax))  

# close MoorDyn simulation (clean up the internal memory, hopefully) when finished
moordyn.Close(system)   

print("New API v2 script executed successfully")
print("++++++++++++++++++++++++++++++++++++++++++++++++++")