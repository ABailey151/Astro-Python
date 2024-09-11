from sys import path
path.append('./Tools')
import Planetary_data as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from sv_coe import sv_from_coe
from plot_func import *

#Orbit Initial conditions
cb = pd.earth
#'''
ra = cb['radius'] + 150000
rp = cb['radius'] + 800
a = (ra+rp)/2
emag = (ra-rp)/(ra+rp)
if emag < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(emag**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(emag**2)))
i = np.deg2rad(90)
lon_an = np.deg2rad(0)
pearg = np.deg2rad(270)
ta0 = np.deg2rad(180)

T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))

r0,v0 = sv_from_coe(cb['mu'],emag,hmag,i,lon_an,pearg,ta0)

#'''
#r0 = np.array([5000, 10000,  2100]) # #np.array([0,42164.1712,0])   
#v0 = np.array([-5.99249498,  1.9253664,   3.24563791]) # #np.array([-np.sqrt((u)/np.linalg.norm(r0)),0,0])

#a,h,hmag,i,lon_an,e,emag,pearg,ta = coe_from_sv(cb['mu'],r0,v0)

def f(t,y):

    r = y[0:3]
    v = y[3:]
    am = -(cb['mu']/(np.linalg.norm(r)**3))
    a = am*r
    dydt = np.append(v,a)
    return dydt

y0 = np.append(r0,v0)

#Set time period and step length/amount
ttotal = T/2
tstep = 10*(1.65**(-((hmag/np.linalg.norm(r0)) + 9.1))) + 20 #Note, this only considers the tangential velocity.
n = 0
t = 0

rs = r0
vs = v0

prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
while t < ttotal:

    prop.set_initial_value(y0,0)
    prop.integrate(tstep)
 
    rs = np.vstack((rs,prop.y[0:3]))
    vs = np.vstack((vs,prop.y[3:6]))
    y0 = prop.y
    
    n = n + 1
    t = t + tstep
    tstep = 10*(1.65**(-(hmag/np.linalg.norm(prop.y[0:3])) + 9.1)) + 20

print('\nOrbital Elements\nh (km^2/s) =',hmag,'\na (km) =',a,'\ni (Deg) =',np.rad2deg(i),'\nRAAN (Deg) =',np.rad2deg(lon_an),'\ne =',emag,'\nArg of Pe (Deg) =',np.rad2deg(pearg),'\nTa (Deg) =',np.rad2deg(ta0),'\n')

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(projection='3d')
tr_plot = Orbit_Plot(ax,cb['radius'],rs,ttotal,n+1,tr_col='r',cb_col=cb['col'],line=False)
ax.set_title(('Keplerian Propagation with variable step size. ('+str(n)+' Steps)'))
plt.show()
