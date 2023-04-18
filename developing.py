import numpy as np
from moorpy.helpers import read_mooring_file

dtC = 1
tMax = 600
time = np.arange(0, tMax, dtC)

data, ch, channels, units = read_mooring_file('saved_runs/cable_outputs/', 'cabledev2.out')

data1 = np.zeros((len(time), len(channels)))

# interpolation
ts = 0
for its in range(0,len(time)):
    t = its*dtC

    while ts < (len(data[:,0])-1):  # search through platform motion data time steps (from .out file)	
        if (data[ts+1,0] > t):				
            frac = ( t - data[ts,0] )/( data[ts+1,0] - data[ts,0] )		# get interpolation fraction
            for j in range(0, len(channels)):
                data1[its][j] = data[ts][j] + frac*( data[ts+1][j] - data[ts][j] ) # interpolate for each platform DOF
            break
        ts += 1