import numpy as np
import scipy as sp
from Kep_Anom import  rv_from_kepler_U
import Planetary_data as pd
from misc_maths import NeN
from spher_harm import *
from Gauss_ve import *
from misc_maths import time_arr
from sv_coe import sv_from_coe



def J2_Encke(cb,N,tstep,r0,v0):

    rs = np.zeros([N,3])
    vs = np.zeros([N,3])
    rs[0] = r0
    vs[0] = v0
    n = 1

    del_y0 = np.zeros([6])

    def f(t,del_y0):

        #del_rx,del_ry,del_rz,del_vx,del_vy,del_vz = del_y0
        #del_r = np.array([del_rx,del_ry,del_rz])
        #del_v = np.array([del_vx,del_vy,del_vz])

        del_r = np.zeros([3])
        del_v = np.zeros([3])

        #print(del_r,del_v)

        rosc,vosc = rv_from_kepler_U(cb['mu'],tstep,r0,v0)

        rpp = rosc + del_r
        vpp = vosc + del_v
        rmag_osc = np.linalg.norm(rosc)
        rmag_pp = np.linalg.norm(rpp)

        aj2 = j2_XYZ(cb,rpp)

        fq = NeN(del_r,rpp)    #1 - (rmag_osc/rmag_pp)**3
        del_a = -((cb['mu'] * (del_r - (fq)))/(rmag_osc**3)) + aj2
        
        dydt = np.concatenate([del_v,del_a], axis=0) #np.array([del_v[0],del_v[1],del_v[2],del_a[0],del_a[1],del_a[2]])

        return dydt

    del_prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
    print('Propagating',N ,'Steps using the Encke Method...')
    while n < N:

        del_prop.set_initial_value(del_y0,0)
        del_prop.integrate(tstep)
        zr = np.array([del_prop.y[0],del_prop.y[1],del_prop.y[2]])#.reshape(1,-1)
        zv = np.array([del_prop.y[3],del_prop.y[4],del_prop.y[5]])#.reshape(1,-1)

        rosc, vosc = rv_from_kepler_U(cb['mu'],tstep,r0,v0)

        r0 = np.squeeze(rosc + zr)
        v0 = np.squeeze(vosc + zv)
        rs[n] = r0
        vs[n] = v0
        
        '''
        Use below to rectify only after a given tolerance is exceeded
        del_y0 = del_prop.y

        rtol = np.linalg.norm(zr)/np.linalg.norm(rosc)

        if rtol > 1e-7:

            del_y0 = np.zeros([6])
            #print('rec')
        '''
            
        #y0 = np.concatenate([r0,v0], axis=0)
        #print(del_y0)
        #print(np.linalg.norm(zr)/np.linalg.norm(r0))

        n = n + 1
        #print('Encke:',n,'/',N)
    
    print('Done')

    return rs,vs

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
        print('before',t0)
        rosc, vosc = rv_from_kepler_U(cb['mu'],t0,r0,v0)

        tol = 0
        #print(rosc,t0)

        while tol < 1e-3:
            prop.set_initial_value(y0,0)
            prop.integrate(dt)

            rpp = rosc + prop.y[0:3]
            vpp = vosc + prop.y[3:6]

            tol = np.linalg.norm(prop.y[0:3])/np.linalg.norm(rpp)
            
            y0 = prop.y
            dt = dt + ddt
            #print(tol,t0,dt)
        
        #print('Rectified',t0,dt)
        rs = np.vstack([rs,rpp])
        vs = np.vstack([vs,vpp])
        zrs = np.vstack([zrs,prop.y[0:3]])
        t0 = t0 + dt
        dt = dt/2
        y0 = np.zeros([6])
        N = N + 1
    
    return rs,vs,N

def J2_Cowell(cb,N,tstep,r0,v0):

    y0 = np.array([r0[0],r0[1],r0[2],v0[0],v0[1],v0[2]])

    def f(t,y):

        rx,ry,rz,vx,vy,vz = y
        r = np.array([rx,ry,rz])
        rmag = np.linalg.norm(r)
        v = np.array([vx,vy,vz])

        ag = -(cb['mu']/(rmag**3)) * r

        aj2 = j2_XYZ(cb,r)

        a = ag + aj2

        dydt = np.array([v[0],v[1],v[2],a[0],a[1],a[2]]) #np.concatenate((v,a), axis=0)

        return dydt


    rs = np.zeros([N,3])
    vs = np.zeros([N,3])

    rs[0] = r0
    vs[0] = v0

    n = 1

    prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
    print('Propagating',N ,'Steps using the Cowell Method...')
    while n < N:

        prop.set_initial_value(y0,0)
        prop.integrate(tstep)
    
        rs[n] = np.array([prop.y[0],prop.y[1],prop.y[2]])
        vs[n] = np.array([prop.y[3],prop.y[4],prop.y[5]])

        y0 = prop.y
        n = n + 1
    
    print('Done') 

    return rs,vs

def tb_prop(cb,N,tstep,r0,v0):

    y0 = np.array([r0[0],r0[1],r0[2],v0[0],v0[1],v0[2]])

    def f(t,y):

        rx,ry,rz,vx,vy,vz = y
        r = np.array([rx,ry,rz])
        rmag = np.linalg.norm(r)
        v = np.array([vx,vy,vz])

        a = -(cb['mu']/(rmag**3)) * r

        dydt = np.array([v[0],v[1],v[2],a[0],a[1],a[2]]) #np.concatenate((v,a), axis=0)
        return dydt

    rs = np.zeros([N,3])
    vs = np.zeros([N,3])

    rs[0] = r0
    vs[0] = v0

    n = 1

    prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
    #print('Propagating',N ,'Steps...')
    while n < N:

        prop.set_initial_value(y0,0)
        prop.integrate(tstep)
    
        rs[n] = np.array([prop.y[0],prop.y[1],prop.y[2]])
        vs[n] = np.array([prop.y[3],prop.y[4],prop.y[5]])

        y0 = prop.y
        n = n + 1
        #print(n,'/',N)
    
    #print('Done\n')

    return rs,vs

def gauss_ve_prop(cb,N,tstep,hmag,emag,ta,lon_an,i,pearg):

    y0 = np.array([hmag,emag,ta,lon_an,i,pearg])

    def f(t,y):
        hmag,emag,ta,lon_an,i,pearg = y
        rmag = (hmag**2) / (cb['mu'] * (1 + (emag * np.cos(ta))))
        
        df_dt = guass_dt(cb,rmag,hmag,emag,ta,i,pearg)
        return df_dt


    elems = np.zeros([N,6])
    elems[0] = y0

    n = 1

    prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
    print('Propagating',N ,'Steps using Gauss\' Variational Equations...')
    while n < N:

        prop.set_initial_value(y0,0)
        prop.integrate(tstep)

        elems[n] = prop.y

        y0 = prop.y
        n = n + 1
        #print('Gauss:',n,'/',N)

    print('Converting Elements to state vectors...')
    rs = np.zeros([N,3])
    vs = np.zeros([N,3])
    si = 0
    while si < N:
        hmag = np.squeeze(elems[si,0])
        emag = np.squeeze(elems[si,1])
        ta = np.squeeze(elems[si,2])
        lon_an = np.squeeze(elems[si,3])
        i = np.squeeze(elems[si,4])
        pearg = np.squeeze(elems[si,5])
        rs[si],vs[si] = sv_from_coe(cb['mu'],emag,hmag,i,lon_an,pearg,ta)
        si = si + 1
    
    print('Done')

    return elems,rs,vs