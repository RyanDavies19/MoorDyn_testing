import numpy as np
import os
import matplotlib.pyplot as plt

# dtWave = 0.15 # needs to match dtWave from MDC input file
# T = 19.22 # wave period pulse
# T = 20 # wave period
# H = 12.0 # max wave height
# UX = 0 # currents
wtr_depth = 200
# Tmax = 100.0
# plot = True
x = np.linspace(-50, 50.0, num=200)
y = np.linspace(-50, 50, num=2)
z = np.linspace(-100, 15.0, num=50)
# t = np.linspace(0.0, Tmax, num=1001)
# theta = (2.0 * np.pi * (t-(Tmax/2)) / T)
# # zeta = H * np.exp(-(theta**2)) * np.sin(theta) # one wave pulse at t = 10
# zeta = (H/2) * np.sin(2.0 * np.pi * t / T) # sin wave

with open('water_grid.txt', 'w') as f:
    f.write('--------------------- MoorDyn Waves grid File ----------------------------------\n')
    f.write('List of grid points, in 3 blocks (x, y, z)\n')
    f.write('Each block starts with a 2 (i.e. list of coords), and then the coordinates (m)\n')
    for l in (x, y, z):
        txt = '1\n'
        for p in l:
            txt = txt + str(p) + ' '
        f.write(txt.strip() + '\n')

# with open('wave_elevation.txt', 'w') as f:
#     for i, tt in enumerate(t):
#         f.write(str(tt) + ' ' + str(zeta[i]) + '\n')

# with open('current_profile.txt', 'w') as f:
#     f.write('--------------------- MoorDyn steady currents File ----------------------------------\n')
#     f.write('Tabulated file with the water currents components\n')
#     f.write('z (m), ux (m/s), uy (m/s), uz (m/s)\n')
#     for zz in z:
#         ux = 0.0
#         if zz > -wtr_depth and zz <= -10.0:
#             ux = UX * (1.0 - (-zz - 10.0) / (wtr_depth*0.8))
#         elif zz > -10.0 and zz < -H:
#             ux = UX * (1.0 - (zz + 10.0) / (10.0 - H))
#         f.write(str(zz) + " " + str(ux) + " 0 0\n")

# MDF stuff:

# if UX > 0:
#     current = '1'
# else:
#     current = '0'

# with open('WaterKin.txt', 'w') as f:
#     f.write('MoorDyn Waves and Currents input file\n')
#     f.write('space for notes\n')
#     f.write('--------------------------- WAVES -------------------------------------\n')
#     f.write('3                    WaveKinMod  - type of wave input {0 no waves; set up grid of wave data based on time series} \n')
#     f.write('"WaveElev.txt"       WaveKinFile - file containing wave elevation time series at 0,0,0\n')
#     f.write('{}                  dtWave      - time step to use in setting up wave kinematics grid (s)\n'.format(dtWave))
#     f.write('0                    WaveDir     - wave heading (deg)\n')
    
#     # write wave grid MDF
#     coords = ['X','Y','Z']
#     i = 0
#     for l in (x, y, z):
#         f.write('2     - {} wave input type (0: not used; 1: list values in ascending order; 2: uniform specified by -{}lim, {}lim, num) \n'.format(coords[i],coords[i],coords[i]))
#         f.write(str(np.min(l)) + ', ' + str(np.max(l)) + ', ' + str(len(l)) + '\n')
#         i=+1

#     # for l in (x, y, z):
#     #     f.write('1     - {} wave input type (0: not used; 1: list values in ascending order; 2: uniform specified by -{}lim, {}lim, num) \n'.format(coords[i],coords[i],coords[i]))
#     #     for node in l:
#     #         f.write(str(node)+', ')
#     #     f.write('\n')
#     #     i=+1

#     # Write currents MDF (Only in X direction for now)
#     f.write('--------------------------- CURRENT -------------------------------------\n')
#     f.write(current+'                    CurrentMod  - type of current input {0 no current; 1 steady current profile described below} \n')
#     f.write('z-depth     x-current      y-current\n')
#     f.write('(m)           (m/s)         (m/s)\n')
#     for zz in z:
#         ux = 0.0
#         if zz > -50.0 and zz <= -10.0:
#             ux = UX * (1.0 - (-zz - 10.0) / 40.0)
#         elif zz > -10.0 and zz < -H:
#             ux = UX * (1.0 - (zz + 10.0) / (10.0 - H))
#         f.write(str(zz) + "     " + str(ux) + "     0\n")  # zero in string represents y direction current   
#     f.write('--------------------- need this line ------------------\n')

# os.system('cp wave_elevation.txt ../MD_fortran_input/WaveElev.txt') # Format is same for MDC and MDF
# os.system('mv WaterKin.txt ../MD_fortran_input/WaterKin.txt') # Format is same for MDC and MDF

# if plot:
#     plt.plot(t,zeta)
#     plt.xlabel("time (s)")
#     plt.ylabel("wave elevation (m)")
#     plt.show()