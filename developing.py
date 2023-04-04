import numpy as np

dtC = 0.1
tMax = 10
time = np.arange(0, tMax, dtC)


scaler = [1., 1., 1., np.pi/180., np.pi/180., np.pi/180.]  # for scaling platform position inputs
outFileName = "PtfmMotions.dat"
i=0  # file line number
t_in = []
Xp_in = []
myfile2 = open(outFileName, 'r');     # open an input stream to the line data input file
if myfile2:
    print(outFileName, " opened.")
    if not myfile2.closed: 
          
        for line2 in myfile2:
            # skip data in first two lines (headers)
            if (i < 2):
                i+=1
                continue
            
            #  split line by tabs
            datarow = list(line2.split())
            
            if (len(datarow) < 7): 
                print("Seems like we've hit a bad line or end of file. ")
                break;                  # break if we're on a last empty line
            
            t_in.append(float(datarow[0]))
            scaled_data = [0,0,0,0,0,0]			
            for j in range(5):
                scaled_data[j] = float(datarow[j+1])*scaler[j] # add platform positions
            Xp_in.append(scaled_data)
            i += 1

    myfile2.close()
    print("Done reading PtfmMotions.dat. Last line read: ", i)

Xp_in = np.array(Xp_in)


# interpolator for platform positions: t_in is vector of time steps from position input file. xp_in is dof
ts = 0
xp = np.zeros((len(time),len(Xp_in[0])))
for its in range(0, len(time)):

    t = its*dtC
    
    # interpolate platform positions from .out file data, and approximate velocities

    while ts < (len(t_in)-1):  # search through platform motion data time steps (from .out file)	
        if (t_in[ts+1] > t):				
            frac = ( t - t_in[ts] )/( t_in[ts+1] - t_in[ts] )		# get interpolation fraction
            for j in range(0, len(Xp_in[0])):
                xp[its][j] = Xp_in[ts][j] + frac*( Xp_in[ts+1][j] - Xp_in[ts][j] ) # interpolate for each platform DOF
            break
        ts += 1

xd = np.zeros((len(time),len(Xp_in[0])))
xold = np.zeros(len(Xp_in[0]))
# calculate velocities using finite difference
for i in range(len(time)):
    xd [i] = (xp[i] - xold)/dtC
    xold =  xp[i]