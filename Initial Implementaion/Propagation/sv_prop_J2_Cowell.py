from sys import path
path.append('./Tools')
import Planetary_data as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from sv_coe import sv_from_coe, coe_from_sv
from plot_func import *
from misc_maths import time_arr

#Orbit Initial conditions
cb = pd.earth
ra = cb['radius'] + 3062
rp = cb['radius'] + 300
a = (ra+rp)/2
emag = (ra-rp)/(ra+rp)
if emag < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(emag**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(emag**2)))
i = np.deg2rad(28)
lon_an = np.deg2rad(45)
pearg = np.deg2rad(30)
ta0 = np.deg2rad(40)

T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))

r0,v0 = sv_from_coe(cb['mu'],emag,hmag,i,lon_an,pearg,ta0)



#Set time period and step length/amount
ttotal = 100*T#2*24*(60**2)
tstep = 600
N = int(np.ceil(ttotal/tstep) + 1)
print(tstep,'Seconds per',N,'Steps')
ts = time_arr(N,tstep)

def J2_Cowell(cb,N,tstep,r0,v0):

    y0 = np.array([r0[0],r0[1],r0[2],v0[0],v0[1],v0[2]])

    def f(t,y):

        rx,ry,rz,vx,vy,vz = y
        r = np.array([rx,ry,rz])
        rmag = np.linalg.norm(r)
        v = np.array([vx,vy,vz])

        ag = -(cb['mu']/(rmag**3)) * r

        j2c = ((3*cb['j2']*cb['mu']*(cb['radius']**2))/(2*rmag**4))
        aj2x = j2c * (r[0]/rmag) * ((5*((r[2]**2)/(rmag**2))) - 1)
        aj2y = j2c * (r[1]/rmag) * ((5*((r[2]**2)/(rmag**2))) - 1)
        aj2z = j2c * (r[2]/rmag) * ((5*((r[2]**2)/(rmag**2))) - 3)
        aj2 = np.array([aj2x, aj2y, aj2z])

        a = ag + aj2

        dydt = np.array([v[0],v[1],v[2],a[0],a[1],a[2]]) #np.concatenate((v,a), axis=0)
        return dydt


    rs = np.zeros([N,3])
    vs = np.zeros([N,3])

    rs[0] = r0
    vs[0] = v0

    n = 1

    prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
    while n < N:

        prop.set_initial_value(y0,0)
        prop.integrate(tstep)
    
        rs[n] = prop.y[0:3] #np.array([prop.y[0],prop.y[1],prop.y[2]])
        vs[n] = prop.y[3:6] 

        y0 = prop.y
        n = n + 1
    
    return rs,vs
    
rs,vs = J2_Cowell(cb,N,tstep,r0,v0)

lon_ans = np.zeros([N,1])
peargs =np.zeros([N,1])
lon_ans[0] = coe_from_sv(cb['mu'],r0,v0)[4]
peargs[0] = coe_from_sv(cb['mu'],r0,v0)[7]
i = 0
while i < N:
    lon_ans[i] = coe_from_sv(cb['mu'],rs[i],vs[i])[4]
    peargs[i] = coe_from_sv(cb['mu'],rs[i],vs[i])[7]
    i = i + 1

lon_an_rate = np.rad2deg((lon_ans[-1]-lon_ans[0])/(ttotal/3600))
pearg_rate = np.rad2deg((peargs[-1]-peargs[0])/(ttotal/3600))
print('Lon_AN deg/hour =',lon_an_rate,'Cowell(5000 Steps) delta =',lon_an_rate - -0.17230583)
print('Pe_Arg deg/hour =',pearg_rate,'Cowell(5000 Steps) delta =',pearg_rate - 0.28216778)
#np.savetxt('rs_Cowell.csv',np.concatenate((ts,lon_ans),axis=1),delimiter=',')

#Element Plots
Elems = plt.figure(figsize=(8,8))
Elems.suptitle('J2 Perturbations using Cowell\'s Method', fontsize=16)
ax2d1 = Elems.add_subplot(211)
ax2d1.grid()
ax2d1.plot((ts/60**2), lon_ans,'r')
plt.title('Right ascension of the Acsending Node')
plt.xlabel('Time (Hours)')
plt.ylabel('(Degrees)')

ax2d2 = Elems.add_subplot(212,sharex=ax2d1)
ax2d2.grid()
ax2d2.plot((ts/60**2), peargs,'b')
plt.title('Argument of Periapsis')
plt.ylabel('(Degrees)')
plt.xlabel('Time (Hours)')
plt.subplots_adjust(left=0.08,right=0.95,bottom=0.07,top=0.91,hspace=0.25)

#3d Plot
fig = plt.figure(figsize=(8,8))
ax = plt.subplot(projection='3d')
tr_plot = Orbit_Plot(ax,cb['radius'],rs,ttotal,N,tr_col='r',cb_col=cb['col'])
plt.show()
