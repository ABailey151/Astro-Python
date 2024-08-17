from sys import path
path.append('./Tools')
import Planetary_data as pd
import numpy as np
import scipy as sp

#Import functions where needed
from sv_coe import *
#from geo_func import *
from Kep_Anom import rv_from_kepler_U
from plot_func import *
#from misc_maths import NeN
from spher_harm import j2_XYZ
from misc_maths import NeN
#from int_props import *

'''
def J2_Encke_variable_step(cb,dt,ddt,ttotal,r0,v0):
    rs = r0
    vs = v0
    zrs = np.zeros([3])
    
    
    rs[0] = r0
    vs[0] = v0
    zrs[0] = 0
    zvs[0] = 0
    

    def f(t,y0):
            
            del_r = y0[0:3]
            del_v = y0[3:6]

            #print(dt)
            #print(del_r,del_v)
            rosc,vosc = rv_from_kepler_U(cb['mu'],dt,r0,v0)

            rpp = rosc + del_r
            vpp = vosc + del_v
            rmag_osc = np.linalg.norm(rosc)
            rmag_pp = np.linalg.norm(rpp)

            aj2 = j2_XYZ(cb,rpp)

            fq = NeN(del_r,rpp)
            del_a = -((cb['mu'] * (del_r - (fq)))/(rmag_osc**3)) + aj2
            #del_a = aj2 + ((pd.earth['mu']/(rmag_osc**3)) * ((fq*rpp) - del_r))
            
            dydt = np.concatenate([del_v,del_a], axis=0)

            return dydt

    #Initialise integrator
    t0 = 0
    N = 0

    y0 = np.zeros([6])
    prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')

    while t0 < ttotal:

        rosc, vosc = rv_from_kepler_U(cb['mu'],t0,r0,v0)

        tol = 0
        #print(rosc,t0)

        while tol < 1e-3:
            prop.set_initial_value(y0,0)
            prop.integrate(dt)
            zr = prop.y[0:3]
            zv = prop.y[3:6]

            rpp = rosc + zr
            vpp = vosc + zv

            tol = np.linalg.norm(prop.y[0:3])/np.linalg.norm(rpp)
            
            y0 = prop.y
            dt = dt + ddt
            print(tol,dt)
        
        print('Rectified',t0,dt)
        rs = np.vstack([rs,rpp])
        vs = np.vstack([vs,vpp])
        zrs = np.vstack([zrs,prop.y[0:3]])
        t0 = t0 + dt
        dt = 0
        #y0 = np.zeros([6])
        N = N + 1
    
    return rs,vs,N
'''
def J2_Encke_variable_step(cb,dt,ddt,ttotal,r0,v0):
    rs = r0
    vs = v0
    zrs = np.zeros([3])
    '''
    rs[0] = r0
    vs[0] = v0
    zrs[0] = 0
    zvs[0] = 0
    '''

    def f(t,y0):

            del_r = y0[0:3]
            del_v = y0[3:6]

            #print(del_r,del_v)
            rosc,vosc = rv_from_kepler_U(cb['mu'],dt,r0,v0)

            rpp = rosc + del_r
            #vpp = vosc + del_v
            rmag_osc = np.linalg.norm(rosc)
            rmag_pp = np.linalg.norm(rpp)

            aj2 = j2_XYZ(cb,rpp)

            fq = NeN(del_r,rpp)
            #del_a = -((cb['mu'] * (del_r - (fq)))/(rmag_osc**3)) + aj2
            del_a = aj2 + ((pd.earth['mu']/(rmag_osc**3)) * ((fq) - del_r))
            
            dydt = np.concatenate([del_v,del_a], axis=0)

            return dydt

    #Initialise integrator
    t0 = 0
    N = 0

    y0 = np.zeros([6])
    prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')

    while t0 < ttotal:
        rosc, vosc = rv_from_kepler_U(cb['mu'],t0,r0,v0)

        tol = 0
        #print(rosc,t0)

        while tol < 1e-5:
            prop.set_initial_value(y0,0)
            prop.integrate(dt)

            rpp = rosc + prop.y[0:3]
            vpp = vosc + prop.y[3:6]

            tol = np.linalg.norm(prop.y[0:3])/np.linalg.norm(rpp)
            
            y0 = prop.y
            dt = dt + ddt
            print(tol,t0,dt)
        
        print('Rectified',t0,dt)
        rs = np.vstack([rs,rpp])
        vs = np.vstack([vs,vpp])
        zrs = np.vstack([zrs,prop.y[0:3]])
        t0 = t0 + dt
        dt = dt/2
        y0 = np.zeros([6])
        N = N + 1
    
    return rs,vs,N

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

#T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))

r0,v0 = sv_from_coe(cb['mu'],e,hmag,i,lon_an,pearg,ta0)
r0mag = np.linalg.norm(r0)

#Set time period and initial step size
ttotal = 2*(3600)
dt = 1
ddt = 1

rs,vs,N = J2_Encke_variable_step(cb,dt,ddt,ttotal,r0,v0)

elems_encke = sv_to_elems(cb,rs,vs,N)
lon_an_rate_encke = np.rad2deg((elems_encke[-1,3]-elems_encke[0,3])/(ttotal/3600))
pearg_rate_encke = np.rad2deg((elems_encke[-1,5]-elems_encke[0,5])/(ttotal/3600))
print('Encke Lon_AN deg/hour =',lon_an_rate_encke,)
print('Encke Pe_Arg deg/hour =',pearg_rate_encke,)

##Plots
fig = plt.figure(figsize=(10,8))
axes = plt.subplot(projection='3d')
p1 = Orbit_Plot(axes,cb['radius'],rs,ttotal,N,'r',cb['col'])
plt.subplots_adjust(left=0.05,right=0.95,bottom=0.09,top=0.95)
#zr_plot = Traj_simple(axes,zrs[N,1],zrs,ttotal,N)
#axes.quiver(0,0,0,h0[0],h0[1],h0[2],length=(r0mag/h0mag),color='y',arrow_length_ratio=0.1)
plt.show()
