import os
import sys
import numpy as np
import matplotlib.pyplot as plt

import moorpy as mp
import moorpy.line as Line
from moorpy.helpers import set_axes_equal, read_mooring_file

rootname = "case5"
dirname = "saved_runs/{}_outputs/".format(rootname)
versionname1 = "dev2"
versionname2 = "v2old"

def compare_line (line):

    data1, ch1, channels1, units1 = read_mooring_file(dirname, rootname+versionname1+line+'.out') # remember number starts on 1 rather than 0
    data2, ch2, channels2, units2 = read_mooring_file(dirname, rootname+versionname2+line+'.out') # remember number starts on 1 rather than 0

    channels = []
    for channel in channels1:
        if 'Seg' not in channel:
            if 'px' in channel:
                channel = (channel.split('p'))[0]
                channels.append(channel)

    # channels: list of strings of channel names (1st line of file)
    # units: list of strings of units (2nd line of file)
    # ch: disctionary where each channel name corresponds to the column number, easier indexing with data array
    # data: np array of all lines below first two

    # RMSE for single timestep and single node
    # find difference between both datasets for x, y, and z at each timestep for each node
    # square the 3 differences and sum squares
    # divide by 3 (number of datapoints)
    # sqrt whole thing

    if (len(data1[:,0]) != len(data2[:,0])) or (len(channels1) != len(channels2)):
        print("Error, make sure both datasets are same length")
        # return

    time = data1[:,0]
    size = (len(time), int(len(channels)-1))
    rmse_array = np.zeros(size)


    for j in range(0,len(channels),3):
        i = 0
        for i in range(0,len(time)):
            if i == 0:
                i = 2
                continue

            position1 = data1[i,j:j+3]
            position2 = data2[i,j:j+3]
            if i-1 % 4000 == 0:
                print('Timestep = ', i-1, ' Node = ', j/3)
                print("position1 (cont): ", position1)
                print("position2 (test): ", position2)

            diff = np.subtract(position1, position2)

            rmse = np.sqrt(np.sum(np.square(diff))/len(diff))
            rmse_array[i,j] = rmse

    rmse_totals = np.zeros((len(time),2))
    rmse_totals [:,0] = time
    i = 0
    for i in range(0,len(time)):
        rmse_totals[i,1] = np.sqrt(np.sum(np.square(rmse_array[i,:]))/len(rmse_array[i,:]))

    # np.set_printoptions(threshold=sys.maxsize)

    # rmse_totals = np.flip(rmse_totals)

    return rmse_totals

    # Working! Need to find a good way to smooth the data though, becasue the high time resolution means really messy plotting 

def compare_lines(lines):

    for i in range(0,len(lines)):
        print("Line ", i)
        rmse_line = compare_line(lines[i])
        if i == 0:
            time = rmse_line[:,0]
            rmse = np.zeros((len(time),len(lines)))
        rmse [:,i] = rmse_line[:,1]
    
    rmse_final = np.zeros(len(time))
    for j in range(0,len(time)):
        rmse_final[j] = np.sqrt(np.sum(np.square(rmse[j,:]))/len(rmse[j,:]))
    
    return time, rmse_final


if __name__ == "__main__":

    lines = ["_Line1", "_Line2", "_Line3"]

    time, rmse_out = compare_lines(lines)

    plt.plot(time, rmse_out)
    plt.xlabel('time')
    plt.ylabel('rmse')
    plt.title('RMSE between {} and {}'.format(rootname+versionname1, rootname+versionname2))
    # plt.ylim(0,0.005)
    plt.show()