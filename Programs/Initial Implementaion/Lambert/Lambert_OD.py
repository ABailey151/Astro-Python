from sys import path
path.append('./Tools')
import numpy as np
import Planetary_data as pd
from Kep_Anom import KU_Cz,KU_Sz
from Lambert_Functions import *
from sv_coe import coe_from_sv,true_anom

cb = pd.earth
mu = cb['mu']

dt = 60*60

r1 = np.array([5000, 10000, 2100])
r2 = np.array([-14600, 2500, 7000])
print('Postion vectors over',dt,'second period:','\nr1 =',r1,'\nr2 =',r2)

r1mag = np.linalg.norm(r1)
r2mag = np.linalg.norm(r2)

C12 = np.cross(r1,r2)

#//Orbit direction
#print('1 = Prograde, 2 = Retrograde')
dir = 1 #int(input())
#//

#//Determine change in true anomaly
#For prograde orbit
if dir == 1:
        if C12[2] >= 0:
            del_ta = np.arccos(np.dot(r1,r2)/(r1mag*r2mag))      
        else:
              del_ta = 2*np.pi - np.arccos(np.dot(r1,r2)/(r1mag*r2mag))
#For retrograde orbit
else:
        if C12[2] < 0:
            del_ta = np.arccos(np.dot(r1,r2)/(r1mag*r2mag))      
        else:
              del_ta = 2*np.pi - np.arccos(np.dot(r1,r2)/(r1mag*r2mag))
#//

#//Calculate the constant A
A = np.sin(del_ta) * np.sqrt((r1mag*r2mag)/(1-np.cos(del_ta)))
print('\nA (km) =',A,'\ndel_ta (Deg) =',np.rad2deg(del_ta))
#//

#//Find initial guess for z
z0 = 0
Fz0 = lambert_F_z(r1mag,r2mag,mu,dt,A,z0)

while Fz0 < 0:
    cz0 = KU_Cz(z0)
    sz0 = KU_Sz(z0)
    yz0 = lambert_y_z(r1mag,r2mag,A,z0)
    Fz0 = lambert_F_z(r1mag,r2mag,mu,dt,A,z0)
    z0 = z0 + 0.01
#//

#//Iterate z0 using Newton's method
tol = 1e-8
ratio = 1
z = z0
while abs(ratio) > tol:
      ratio = lambert_F_z(r1mag,r2mag,mu,dt,A,z)/lambert_dF_z(r1mag,r2mag,A,z)
      z = z - ratio
print('z =',z)
#//

#//Calculate Lagrange Functions
cz = KU_Cz(z)
sz = KU_Sz(z)
yz = lambert_y_z(r1mag,r2mag,A,z)

lag_f = 1 - yz/r1mag
lag_g = A*np.sqrt(yz/mu)
lag_gdot = 1 - yz/r2mag
#//

#//Calculate Velocities at both postions
v1 = (1/lag_g)*(r2 - lag_f*r1)
v2 = (1/lag_g)*(lag_gdot*r2 - r1)
#//

#//Determine the orbital elements
a,h,hmag,i,lon_an,e,emag,pearg,ta1 = coe_from_sv(mu,r1,v1)
print('\nOrbital Elements\nh (km^2/s) =',hmag,'\na (km) =',a,'\ni (Deg) =',np.rad2deg(i),'\nRAAN (Deg) =',np.rad2deg(lon_an),'\ne =',emag,'\nArg of Pe (Deg) =',np.rad2deg(pearg),'\nTa (Deg) =',np.rad2deg(ta1),'\n')
#//

#//Perigee altitude and time since passage
pe_h = ((hmag**2)/mu)*(1/(1+emag)) - cb['radius']

#Time since Perigee passage for the first sighting
E1 = 2*np.arctan(np.sqrt((1-emag)/(1+emag))*np.tan(ta1/2))
Me1 = E1 - emag*np.sin(E1)
tp1 = ((hmag**3)/(mu**2))*((1-(emag**2))**(-3/2))*Me1
#Time since Perigee passage for the first sighting
ta2 = true_anom(r2,v2,e,emag)
E2 = 2*np.arctan(np.sqrt((1-emag)/(1+emag))*np.tan(ta2/2))
Me2 = E2 - emag*np.sin(E2)
tp2 = ((hmag**3)/(mu**2))*((1-(emag**2))**(-3/2))*Me2

print('Perigee Altitude (km) =',pe_h,'\nTime since Perigee (s), r1 =',tp1,'r2 =',tp2,'\n')
#//
