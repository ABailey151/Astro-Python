from math import cos, gamma
import numpy as np

u = 398600.4

rft = np.array([4.1852E7,6.2778E7,10.463E7])
vft = np.array([2.5936E4,5.1972E4,0])

r = 0.0003048 * rft
rmag = np.linalg.norm(r)

v = 0.0003048 * vft
vmag = np.linalg.norm(v)

h = np.cross(r,v)
hmag = np.linalg.norm(h)

spe = ((vmag**2)/2) - (u/rmag)

fpa = np.rad2deg(np.arccos(hmag/(rmag*vmag)))

print('The specific Mechanical Energy is',spe,'km^2/s^2\nThe Specific Angular Momentum is',hmag,'km^2/s\nThe Flight Path Angle is',fpa,'degrees')