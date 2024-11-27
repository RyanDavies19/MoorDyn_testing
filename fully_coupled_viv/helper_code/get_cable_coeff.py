import numpy as np

C_dnc = 1.2
C_dnb = 0.7
C_dab1 = 1
C_dab2 = 0
C_dac = 0
C_anb = 1
C_anc = 1
C_aab = 0.5 
C_aac = 0
'''

C_dnc : Normal drag coeff for the cable. The default is 1.2.
C_dnb : Normal drag coeff for the buoyancy module. The default is 0.7.
C_dab1 : Drag coefficient for exposed ends of buoyancy module. The default is 1.
C_dab2 : Axial drag coefficient for buoyancy module (skin friction). The default is 0.
C_dac : Axial drag coefficient for cable (skin friction). The default is 0.
C_anb : Normal added mass coefficient for buoyancy module. The default is 1.
C_anc : Normal added mass coefficient for cable. The default is 1.
C_aab : Axial added mass coefficient for buoyancy module. The default is 0.5 (assumed sphere added mass coeff).
C_aac : Axial added mass coefficient for cable. The default is 0.

'''

# # Buoyancy module properties smoothed out along section
# deq = 0.3 # volume equiv diameter for buoy section
# dc = 0.1607 # diameter of cable
# db = 0.865 # diameter of buoy
# Lbs = 11.23 # length of segment for per length calcs (or spacing)
# if Lbs == 0:
#     ValueError('Buoyancy module spacing is zero')
# Lb = 0.9 # buoy length

# Hydro coefficients of modules not smoohed out
deq = 0.865 # volume equiv diameter for buoy section
dc = 0.1607 # diameter of cable
db = 0.865 # diameter of buoy
Lbs = 0.9 # length of segment for per length calcs (or spacing)
if Lbs == 0:
    ValueError('Buoyancy module spacing is zero')
Lb = 0.9 # buoy length

Cd = 1/(Lbs * deq)*(C_dnc * dc * (Lbs - Lb) + C_dnb * db * Lb)
CdAx = 1 / (Lbs * deq) *(C_dab1 * (db**2 - dc**2)/4 + C_dab2 * db * Lb + C_dac * dc * (Lbs - Lb))
Ca = 1 / (Lbs * deq**2) * (C_anb * db**2 * Lb + C_anc * dc**2 *(Lbs - Lb))
CaAx =  1 / (Lbs * deq**2) * (C_aab * db**2 * Lb + C_aac * dc**2 *(Lbs - Lb))

print("Cd: ", Cd)
print("CdAx: ", CdAx)
print("Ca: ", Ca)
print("CaAx: ", CaAx)