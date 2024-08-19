import numpy as np
import scipy as sp
from misc_maths import JulianD

def JD_Rm(JD):

    T0 = (JD - 2451545.0)/36525 #Number of Julian centuries since J2000

    ob = np.deg2rad(23.439 - (0.0130042*T0))

    LpCoefs = np.array(([0,     218.32,  481267.881, 0,     0,      0,         0.9508, 0,      0],
                        [6.29,  135.0,   477198.87,  5.13,  93.3,   483202.03, 0.0518, 135.0,  477198.87],
                        [-1.27, 259.3,  -413335.36,  0.28,  220.2,  960400.89, 0.0095, 259.3, -413335.38],
                        [0.66,  235.7,   890534.22, -0.28,  318.3,  6003.15,   0.0078, 253.7,  890534.22],
                        [0.21,  269.9,   954397.74, -0.17,  217.6, -407332.21, 0.0028, 269.9,  954397.70],
                        [-0.19, 357.5,   35999.05,  0,      0,     0,          0,      0,      0],
                        [-0.11, 106.5,   966404.03, 0,      0,     0,          0,      0,      0]))
    
    i = 0
    for m in range(0,7):
        i = i + (LpCoefs[m,0] * np.sin(np.deg2rad(LpCoefs[m,1] + (LpCoefs[m,2]*T0))))
    Lelon = np.deg2rad((LpCoefs[0,1] + (LpCoefs[0,2]*T0) + i)%360)
    
    j = 0
    for b in range(1,5):
        j = j + (LpCoefs[b,3] * np.sin(np.deg2rad(LpCoefs[b,4] + (LpCoefs[b,5]*T0))))
    Lelat = np.deg2rad(j%360)

    k = 0
    for f in range(1,5):
        k = k + (LpCoefs[f,6] * np.sin(np.deg2rad(LpCoefs[f,7] + (LpCoefs[f,8]*T0))))
    HP = np.deg2rad(LpCoefs[0,6] + k%360)
    
    '''
    Lelon = np.deg2rad((LpCoefs[0,1] + (LpCoefs[0,2]*T0) + (LpCoefs[0,0] * np.sin(LpCoefs[0,1] + (LpCoefs[0,2]*T0))) + (LpCoefs[1,0] * np.sin(LpCoefs[1,1] + (LpCoefs[1,2]*T0))) + (LpCoefs[2,0] * np.sin(LpCoefs[2,1] + (LpCoefs[2,2]*T0))) + (LpCoefs[3,0] * np.sin(LpCoefs[3,1] + (LpCoefs[3,2]*T0))) + (LpCoefs[4,0] * np.sin(LpCoefs[4,1] + (LpCoefs[4,2]*T0))) + (LpCoefs[5,0] * np.sin(LpCoefs[5,1] + (LpCoefs[5,2]*T0))) + (LpCoefs[6,0] * np.sin(LpCoefs[6,1] + (LpCoefs[6,2]*T0))))%360)
    np.deg2rad
    Lelat = np.deg2rad(((LpCoefs[1,3] * np.sin(LpCoefs[1,4] + (LpCoefs[1,5]*T0))) + (LpCoefs[2,3] * np.sin(LpCoefs[2,4] + (LpCoefs[2,5]*T0))) + (LpCoefs[3,3] * np.sin(LpCoefs[3,4] + (LpCoefs[3,5]*T0))) + (LpCoefs[4,3] * np.sin(LpCoefs[4,4] + (LpCoefs[4,5]*T0))))%360)

    HP = np.deg2rad((LpCoefs[0,6] + (LpCoefs[1,6] * np.sin(LpCoefs[1,7] + (LpCoefs[1,8]*T0))) + (LpCoefs[2,6] * np.sin(LpCoefs[2,7] + (LpCoefs[2,8]*T0))) + (LpCoefs[2,6] * np.sin(LpCoefs[2,7] + (LpCoefs[2,8]*T0))) + (LpCoefs[3,6] * np.sin(LpCoefs[3,7] + (LpCoefs[3,8]*T0))) + (LpCoefs[4,6] * np.sin(LpCoefs[4,7] + (LpCoefs[4,8]*T0))))%360)
    '''

    rm = 6378.14/np.sin(HP)

    Uxyz = np.array([np.cos(Lelat)*np.cos(Lelon),
                    (np.cos(ob)*np.cos(Lelat)*np.sin(Lelon)) - (np.sin(ob)*np.sin(Lelat)),
                    (np.sin(ob)*np.cos(Lelat)*np.sin(Lelon)) + (np.cos(ob)*np.sin(Lelat))])
    
    Rm = rm * Uxyz
    return Rm

'''
def prop_luna_3b(cb,N,tstep,r0,v0,epoch):
    def f(t,y,JD):

        rx,ry,rz,vx,vy,vz = y
        r = np.array([rx,ry,rz])
        v = np.array([vx,vy,vz])

        am = -(cb[0]['mu']/(np.linalg.norm(r)**3))*r

        Rm = JD_Rm(JD)
        Rms = Rm - r
        pm = cb[1]['mu'] * ((Rms/(np.linalg.norm(Rms)**3)) - (Rm/(np.linalg.norm(Rm)**3)))

        a = am + pm

        dydt = np.concatenate((v,a), axis=0)
        return dydt

    y0 = np.concatenate((r0,v0), axis=0)

    JD = JulianD(epoch)
    rs = np.zeros([N,3])
    vs = np.zeros([N,3])
    Rms = np.zeros([N,3])

    rs[0] = r0
    vs[0] = v0
    Rms[0] = JD_Rm(JD)
    n = 1
    

    prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
    while n < N:

        prop.set_initial_value(y0,0).set_f_params(JD)
        prop.integrate(tstep)
    
        rs[n] = prop.y[0:3] #np.array([prop.y[0],prop.y[1],prop.y[2]])
        vs[n] = prop.y[3:6] #np.array([prop.y[3],prop.y[4],prop.y[5]])
        Rms[n] = JD_Rm(JD)
        y0 = prop.y
        JD = JD + ((tstep)/86400)
        n = n + 1
        #print(n,'/',N)
    
    return rs,vs,Rms
'''

def prop_luna_3b(bods,tstep0,ttotal,r0,v0,epoch):
    def f(t,y,JD):

        rx,ry,rz,vx,vy,vz = y
        r = np.array([rx,ry,rz])
        v = np.array([vx,vy,vz])

        am = -(bods[0]['mu']/(np.linalg.norm(r)**3))*r

        Rm = JD_Rm(JD)
        RmS = Rm - r
        pm = bods[1]['mu'] * ((RmS/(np.linalg.norm(RmS)**3)) - (Rm/(np.linalg.norm(Rm)**3)))

        a = am + pm

        dydt = np.concatenate((v,a), axis=0)
        return dydt

    y0 = np.concatenate((r0,v0), axis=0)

    N = int(np.ceil(ttotal/tstep0) + 1)
    J0,UT = JulianD(epoch)
    JD = J0 + (UT/24)
    rs = np.zeros([N,3])
    vs = np.zeros([N,3])
    Rms = np.zeros([N,3])
    
    rs[0] = r0
    vs[0] = v0
    Rms[0] = JD_Rm(JD)
    tstep0 = 600
    t = 0
    n = 0
    tstep = tstep0

    prop = sp.integrate.ode(f,jac=None).set_integrator('dopri5')
    while t < ttotal:
        
        prop.set_initial_value(y0,0).set_f_params(JD)
        prop.integrate(tstep)

        rs[n] = prop.y[0:3]
        vs[n] = prop.y[3:6]

        Rms[n] = JD_Rm(JD)
        y0 = prop.y

        
        mba = -(bods[0]['mu']/(np.linalg.norm(rs[n])**3))*rs[n]
        RmS = Rms[n] - rs[n]
        sba = bods[1]['mu'] * ((RmS/(np.linalg.norm(RmS)**3)) - (Rms[n]/(np.linalg.norm(Rms[n])**3)))
        aratio = np.linalg.norm(sba)/np.linalg.norm(mba)
        
        tstep = (tstep0/2) + (tstep0 * aratio)

        if tstep > tstep0:
            tstep = tstep0
        
        print(aratio,tstep)

        JD = JD + ((tstep)/86400)
        n = n + 1
        t = t + tstep
        #print(n,'/',N)
    
    return rs,vs,Rms,n