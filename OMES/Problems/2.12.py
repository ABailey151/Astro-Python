from cmath import pi
import numpy as np

ur = np.array([0.26726,0.53452,0.80178])
v = np.array([-4,3,-5])

vr = (np.dot(v,ur))

vt = np.sqrt((np.linalg.norm(v)**2)-(np.linalg.norm(vr)**2))

fparad = (pi/2) - np.arccos((np.dot(v,ur))/(np.linalg.norm(v)*np.linalg.norm(ur)))

fpa = np.rad2deg(fparad)

print("Radial Velocity Component",vr,"km/s\nTangential Velocity Component",vt,"km/s\nFlight Path Angle",fpa,"Degrees")