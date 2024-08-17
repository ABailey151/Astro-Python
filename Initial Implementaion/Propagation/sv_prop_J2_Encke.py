from sys import path
path.append('./Tools')
import Planetary_data as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from sv_coe import *
from plot_func import *
from Kep_Anom import rv_from_kepler_U
from misc_maths import time_arr,NeN

#Orbit Initial conditions
cb = pd.earth
ra = cb['radius'] + 3062
rp = cb['radius'] + 300
a = (ra+rp)/2
e = (ra-rp)/(ra+rp)
if e < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(e**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(e**2)))
i = np.deg2rad(28)
lon_an = np.deg2rad(45)
pearg = np.deg2rad(30)
ta0 = np.deg2rad(40)
T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))

r0,v0 = sv_from_coe(cb['mu'],e,hmag,i,lon_an,pearg,ta0)
y0 = np.array([r0[0],r0[1],r0[2],v0[0],v0[1],v0[2]])

#Set time period and step length/amount
ttotal = 480*(3600)
tstep = ttotal/10000
N = int(np.ceil(ttotal/tstep) + 1)
print(tstep,'Seconds per',N,'Steps')
ts = time_arr(N,tstep)

def J2_Encke(cb,N,tstep,r0,v0):

    rs = np.zeros([N,3])
    vs = np.zeros([N,3])
    zrs = np.zeros([N,3])
    zvs = np.zeros([N,3])
    rs[0] = r0
    vs[0] = v0
    zrs[0] = 0
    zvs[0] = 0
    n = 1

    del_y0 = np.zeros([6])

    def f(t,del_y0):
            
        del_r = np.zeros([3])
        del_v = np.zeros([3])

        rosc,vosc = rv_from_kepler_U(cb['mu'],tstep,r0,v0)

        rpp = rosc + del_r
        vpp = vosc + del_v
        rmag_osc = np.linalg.norm(rosc)
        rmag_pp = np.linalg.norm(rpp)

        j2c = ((3*cb['j2']*cb['mu']*(cb['radius']**2))/(2*rmag_pp**4))
        aj2x = j2c * (rpp[0]/rmag_pp) * ((5*((rpp[2]**2)/(rmag_pp**2))) - 1)
        aj2y = j2c * (rpp[1]/rmag_pp) * ((5*((rpp[2]**2)/(rmag_pp**2))) - 1)
        aj2z = j2c * (rpp[2]/rmag_pp) * ((5*((rpp[2]**2)/(rmag_pp**2))) - 3)
        aj2 = np.array([aj2x, aj2y, aj2z])

        f = NeN(del_r,rpp)#1 - (rmag_osc/rmag_pp)**3 
        del_a = ((-cb['mu'] * (del_r - (f*rpp)))/(rmag_osc**3)) + aj2

        dydt = np.concatenate([del_v,del_a], axis=0)

        return dydt

    del_prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
    
    while n < N:

        del_prop.set_initial_value(del_y0,0)
        #print(del_prop.y)
        del_prop.integrate(tstep)
        zr =  del_prop.y[0:3]  #np.array([del_prop.y[0],del_prop.y[1],del_prop.y[2]]).reshape(1,-1)
        zv =  del_prop.y[3:6]  #np.array([del_prop.y[3],del_prop.y[4],del_prop.y[5]]).reshape(1,-1)

        rosc, vosc = rv_from_kepler_U(cb['mu'],tstep,r0,v0)
        
        r0 = np.squeeze(rosc + zr)
        v0 = np.squeeze(vosc + zv)
        del_y0 = np.zeros([6])
        rs[n] = r0
        vs[n] = v0
        zrs[n-1] = zr
        zvs[n-1] = zv
        
        n = n + 1
    
    return rs,vs,zrs,zvs

rs,vs,zrs,zvs = J2_Encke(cb,N,tstep,r0,v0)

lon_ans = np.zeros([N,1])
peargs = np.zeros([N,1])
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

#np.savetxt('rs_Encke.csv',np.concatenate((ts,lon_ans),axis=1),delimiter=',')

#Element Plots
Elems = plt.figure(figsize=(8,8))
Elems.suptitle('J2 Perturbations using Encke\'s Method', fontsize=16)
ax2d1 = Elems.add_subplot(211)
ax2d1.grid()
ax2d1.plot((ts/60**2), lon_ans-lon_ans[0],'r')
plt.title('Right ascension of the Acsending Node')
plt.xlabel('Time (Hours)')
plt.ylabel('(Degrees)')

ax2d2 = Elems.add_subplot(212,sharex=ax2d1)
ax2d2.grid()
ax2d2.plot((ts/60**2), peargs-peargs[0],'b')
plt.title('Argument of Periapsis')
plt.ylabel('(Degrees)')
plt.xlabel('Time (Hours)')
plt.subplots_adjust(left=0.08,right=0.95,bottom=0.07,top=0.91,hspace=0.25)

#3d Plot
fig = plt.figure(figsize=(8,8))
ax = plt.subplot(projection='3d')
e_plot = Orbit_Plot(ax,cb['radius'],rs,ttotal,N,tr_col='r',cb_col=cb['col'])
#buh = Traj_simple(ax,zrs[N-1,1],zrs,ttotal,N)
plt.show()
