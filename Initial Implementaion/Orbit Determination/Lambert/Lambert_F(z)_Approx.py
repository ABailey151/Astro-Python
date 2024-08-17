from sys import path
path.append('./Tools')
import numpy as np
import matplotlib.pyplot as plt
import Planetary_data as pd
from Kep_Anom import KU_Sz,KU_Cz

cb = pd.earth
mu = cb['mu']

dt = 60*60

r1 = np.array([5000, 10000, 2100])
r2 = np.array([-14600, 2500, 7000])

r1mag = np.linalg.norm(r1)
r2mag = np.linalg.norm(r2)

C12 = np.cross(r1,r2)

#Orbit direction
print('1 = Prograde, 2 = Retrograde')
dir = 1#int(input())


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
#//
print('A =',A,'',np.rad2deg(del_ta))

z = np.linspace(-np.pi,np.pi,200)
Fz = np.zeros_like(z)
for m in range(0,200):
    zm = z[m]
    czm = KU_Cz(zm)
    szm = KU_Sz(zm)
    yzm = (r1mag + r2mag + A*((zm*szm - 1)/np.sqrt(czm)))
    Fz[m] = szm*((yzm/czm)**(3/2)) + A*np.sqrt(yzm) - dt*np.sqrt(mu)

fig, ax = plt.subplots()
ax.plot(z,Fz,'b')
ax.grid()
ax.set_title('Lambert F(z) Plot')
plt.show()
