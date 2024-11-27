import numpy as np
import matplotlib.pyplot as plt

UX = 1.0 # current speed
sheared = True
depth = 200
z = np.linspace(-depth, 0.0, num=int(depth/2))
current_profile = np.zeros((2,len(z)))

with open('current_profile.txt', 'w') as f:
    f.write('--------------------- MoorDyn steady currents File ----------------------------------\n')
    f.write('Tabulated file with the water currents components\n')
    f.write('z (m), ux (m/s), uy (m/s), uz (m/s)\n')
    for i, zz in enumerate(z):
        if sheared:
            ux = 0.0
            if zz > -depth and zz <= -10.0:
                ux = UX * (1.0 - (-zz - 10.0) / 190.0)
            elif zz > -10.0 and zz < 0:
                ux = UX * (1.0 - (zz + 10.0) / (10.0))
        else:
            ux = UX
        f.write(str(zz) + " " + str(ux) + " 0 0\n")
        current_profile[0][i] = zz
        current_profile[1][i] = ux
    
plt.plot(current_profile[1],current_profile[0])
plt.show()
