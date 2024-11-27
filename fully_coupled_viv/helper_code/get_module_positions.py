import numpy as np

# ---------- Inputs ----------

# Parameters
can_l = 0.9 # can length (m)
N = 6 # number of Modules
rAnchor = [200,0,-200] # cable endpoint location (to find plane of mooring system)

# Module centers (m)
x = np.array([90.95760369,
82.83485989,
72.0160812,
63.99998884,
58.8866394,
55.27623366])
y = np.array([0,
0,
0,
0,
0,
0])
z = np.array([-149.488172,
-141.9186593,
-142.0586469,
-149.7505959,
-159.7269454,
-170.3557543]) 

# estimated module orientation relative to negative horizontal - I know this is strange but its what works, its becasue I have my A and B backwards here - (rad)
orientations = (180/np.pi) * np.array([45,
180,
180,
-45,
-45,
-45])

# ---------- Calculations ----------
"""
Coordinate conventions from David J. Griffiths - Introduction to Electrodynamics - Cambridge University Press (2017)
End A and B positions are calculated using spherical coordinates from the defined can center point and estimated orentation
"""

# check lengths
if not (N == len(x) == len(y) == len(z) == len(orientations)):
    exit("Error: array lengths must all be equal")

# Cable profile plane angle (rad)
phi = np.arctan2(rAnchor[1],rAnchor[0]) # offset from x axis

# calculate end A and B 
rA = [] # closer to hangoff
rB = [] # closer to anchor
rod_index = 3
for i in range(N):

    theta = (np.pi/2) - orientations[i] # offset from z axis (from vertical)

    # spherical coordinates stuff 
    offset_x = (can_l/2) * np.cos(phi) * np.sin(theta)
    offset_y = (can_l/2) * np.sin(phi) * np.sin(theta)
    offset_z = (can_l/2) * np.cos(theta)
    rA.append( np.array([x[i]-offset_x, y[i]-offset_y, z[i]-offset_z]) ) # closer to hangoff
    rB.append( np.array([x[i]+offset_x, y[i]+offset_y, z[i]+offset_z]) ) # closer to anchor
    # print("Can: ", i)
    # print("rA: ", rA[i])
    # print("rB: ", rB[i])
    # MoorDyn format, with second set of coords as unit vector of module
    # print("{}     connector  Free    {:.3f} {:.3f} {:.3f}  {:.3f} {:.3f} {:.3f}       0        p".format(rod_index, rA[i][0], rA[i][1], rA[i][2], rA[i][0]+(np.cos(phi) * np.sin(theta)), rA[i][1] + (np.sin(phi) * np.sin(theta)), rA[i][2] + (np.cos(theta))))
    # print("{}     connector  Free    {:.3f} {:.3f} {:.3f}  {:.3f} {:.3f} {:.3f}       0        p".format(rod_index+1, rB[i][0], rB[i][1], rB[i][2], rB[i][0]+(np.cos(phi) * np.sin(theta)) , rB[i][1] + (np.sin(phi) * np.sin(theta)), rB[i][2] + (np.cos(theta))))

    print("{}     buoy  Free    {:.3f} {:.3f} {:.3f}  {:.3f} {:.3f} {:.3f}       1        p".format(rod_index, rA[i][0], rA[i][1], rA[i][2], rB[i][0], rB[i][1], rB[i][2]))
    rod_index += 1