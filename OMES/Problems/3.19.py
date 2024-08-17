from sys import path
path.append('./Tools')
from U_anom import *

u = 398600.4
r0mag = 7200
vr0 = 1
alpha = 1/10000
dt = 60*60

X = kepler_U(u,dt,r0mag,vr0,alpha)

print(X)