from sys import path
path.append('./Tools')
import Planetary_data as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from sv_coe import *
from plot_func import *
from misc_maths import time_arr


#Orbit Initial conditions
cb = pd.moon
ra = cb['radius'] + 2000
rp = cb['radius'] + 2000
a = (ra+rp)/2
e = (ra-rp)/(ra+rp)
if e < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(e**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(e**2)))
i = np.deg2rad(10)
lon_an = np.deg2rad(0)
pearg = np.deg2rad(270)
ta0 = np.deg2rad(0)
T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))

r0,v0 = sv_from_coe(cb['mu'],e,hmag,i,lon_an,pearg,ta0)

#Spacecraft properties
m0 = 1000
isp = 10000
Th = 2500e-6

def f(t,y,Th,isp):

    rx,ry,rz,vx,vy,vz,m = y
    r = np.array([rx,ry,rz])
    v = np.array([vx,vy,vz])
    ag = -(cb['mu']/(np.linalg.norm(r)**3)) * r
    th = (Th/m)*(v/np.linalg.norm(v))
    a = ag + th
    mdot = -(Th/(isp*9.804e-3))
    dydt = np.array([v[0],v[1],v[2],a[0],a[1],a[2],mdot]) #np.concatenate((v,a), axis=0)
    return dydt

y0 = np.array([r0[0],r0[1],r0[2],v0[0],v0[1],v0[2],m0]) #np.concatenate((r0,v0), axis=0)
t0 = 0

#Set time period and step length/amount
ttotal = 100*(60**2)
tstep = 100
N = int(np.ceil(ttotal/tstep))
print(N,'Steps')
n = 1

#Initialise arrays and write initial conditions
ts = time_arr(N,tstep)
rs = np.zeros([N,3])
vs = np.zeros([N,3])
es = np.zeros([N,1])
smas = np.zeros([N,1])

rs[0] = r0
vs[0] = v0
es[0] = eccen(cb['mu'],r0,v0)[1]
smas[0] = coe_from_sv(cb['mu'],r0,v0)[0]


prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
while n < N:

    prop.set_initial_value(y0,0).set_f_params(Th,isp)
    prop.integrate(tstep)
 
    rs[n] = np.array([prop.y[0],prop.y[1],prop.y[2]])
    vs[n] = np.array([prop.y[3],prop.y[4],prop.y[5]])

    es[n] = eccen(cb['mu'],rs[n],vs[n])[1]
    smas[n] = coe_from_sv(cb['mu'],rs[n],vs[n])[0]
    
    y0 = prop.y
    n = n + 1


Elems = plt.figure(figsize=(8,8))
ax2d1 = Elems.add_subplot(211)
ax2d1.grid()
ax2d1.plot((ts/60**2), smas,'r')
plt.title('Semi-Major Axis')
plt.xlabel('Time (Hours)')
plt.ylabel('a (km)')

ax2d2 = Elems.add_subplot(212,sharex=ax2d1)
ax2d2.grid()
ax2d2.plot((ts/60**2), es,'b')
plt.title('Eccentricity')
plt.xlabel('Time (Hours)')
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.94,hspace=0.25)

orb3d = plt.figure(figsize=(8,8))
ax3d = plt.subplot(projection='3d')
ax3d.set_title(str(cb['name'])+' Centered Inertial Frame')
plt.subplots_adjust(left=0.05,right=0.95,bottom=0.1,top=0.94)
tr_plot = Orbit_Plot(ax3d,cb['radius'],rs,ttotal,n,tr_col='r',cb_col=cb['col'])
ax3d.legend()
plt.show()
