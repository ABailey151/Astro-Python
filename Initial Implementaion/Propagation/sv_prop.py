from sys import path
path.append('./Tools')
import Planetary_data as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from sv_coe import sv_from_coe,coe_from_sv
from plot_func import *


#Orbit Initial conditions
cb = pd.earth
'''
ra = cb['radius'] + 1500
rp = cb['radius'] + 800
a = (ra+rp)/2
e = (ra-rp)/(ra+rp)
if e < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(e**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(e**2)))
i = np.deg2rad(90)
lon_an = np.deg2rad(0)
pearg = np.deg2rad(270)
ta0 = np.deg2rad(180)

T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))
'''
#r0,v0 = sv_from_coe(cb['mu'],e,hmag,i,lon_an,pearg,ta0)
r0 = np.array([ 3830.67827273, -2216.46778571,  6605.09381616]) #([5000, 10000,  2100]) #np.array([0,42164.1712,0])   
v0 = np.array([ 1.50425319, -4.56246159, -0.29217639]) #([-5.99249498,  1.9253664,   3.24563791]) #np.array([-np.sqrt((u)/np.linalg.norm(r0)),0,0])

a,h,hmag,i,lon_an,e,emag,pearg,ta = coe_from_sv(cb['mu'],r0,v0)
T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))
def f(t,y):

    rx,ry,rz,vx,vy,vz = y
    r = np.array([rx,ry,rz])
    v = np.array([vx,vy,vz])
    am = -(cb['mu']/(np.linalg.norm(r)**3))
    a = am*r
    dydt = np.array([v[0],v[1],v[2],a[0],a[1],a[2]]) #np.concatenate((v,a), axis=0)
    return dydt

y0 = np.array([r0[0],r0[1],r0[2],v0[0],v0[1],v0[2]]) #np.concatenate((r0,v0), axis=0)
t0 = 0

#Set time period and step length/amount
ttotal = 10*T
tstep = 10
N = int(np.ceil(ttotal/tstep) + 1)
n = 1

rs = np.zeros([N,3])
vs = np.zeros([N,3])
rs[0] = r0
vs[0] = v0

prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
while n < N:

    prop.set_initial_value(y0,0)
    prop.integrate(tstep)
 
    rs[n] = prop.y[0:3] #np.array([prop.y[0],prop.y[1],prop.y[2]])
    vs[n] = prop.y[3:6] #np.array([prop.y[3],prop.y[4],prop.y[5]])
    y0 = prop.y
    n = n + 1
    #print(n,'/',N)

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(projection='3d')
tr_plot = Orbit_Plot(ax,cb['radius'],rs,ttotal,n,tr_col='r',cb_col=cb['col'])
plt.show()