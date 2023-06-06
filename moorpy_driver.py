# MoorPy Example Script:
# Example of manually setting up a mooring system in MoorPy and solving equilibrium.

import numpy as np
import matplotlib.pyplot as plt
import moorpy as mp


# ----- set up the mooring system and floating body -----

# Create new MoorPy System and set its depth
ms = mp.System(file='MooringTest/case4.dat')

# ----- run the model to demonstrate -----

ms.initialize()                                             # make sure everything's connected

# find equilibrium of just the mooring lines (the body is 'coupled' and will by default not be moved)
ms.solveEquilibrium()                                       # equilibrate
                                   

T = ms.getTensions()
print(T)

# fig, ax = ms.plot()    # plot the system in original configuration  
# plt.show()