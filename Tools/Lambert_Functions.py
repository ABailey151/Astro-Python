from sys import path
path.append('./Tools')
import numpy as np
from Kep_Anom import KU_Cz,KU_Sz

def lambert_z(r1mag,r2mag,mu,dt,A,z0,Fz0):

    while Fz0 < 0:
        cz0 = KU_Cz(z0)
        sz0 = KU_Sz(z0)
        yz0 = lambert_y_z(r1mag,r2mag,A,z0)
        Fz0 = lambert_F_z(r1mag,r2mag,mu,dt,A,z0)
        z0 = z0 + 0.01
    #//

    #//Iterate z0 using Newton's method
    tol = 1e-8
    ratio = 1
    z = z0
    while abs(ratio) > tol:
        ratio = lambert_F_z(r1mag,r2mag,mu,dt,A,z)/lambert_dF_z(r1mag,r2mag,A,z)
        z = z - ratio
    return z
    #//

def lambert_y_z(r1mag,r2mag,A,z):
    cz = KU_Cz(z)
    sz = KU_Sz(z)
    yz = (r1mag + r2mag + A*((z*sz - 1)/np.sqrt(cz)))
    return yz

def lambert_F_z(r1mag,r2mag,mu,dt,A,z):
    cz = KU_Cz(z)
    sz = KU_Sz(z)
    yz = lambert_y_z(r1mag,r2mag,A,z)
    Fz = sz*((yz/cz)**(3/2)) + A*np.sqrt(yz) - dt*np.sqrt(mu)
    return Fz
      
def lambert_dF_z(r1mag,r2mag,A,z):
    
    if z == 0:
            y0 = lambert_y_z(r1mag,r2mag,A,0)
            dFz = (np.sqrt(2)/40)*(y0**(3/2)) + (A/8)*(np.sqrt(y0) + A*np.sqrt(1/(2*y0)))
            return dFz
    else:
        cz = KU_Cz(z)
        sz = KU_Sz(z)
        yz = lambert_y_z(r1mag,r2mag,A,z)
        dFz = ((yz/cz)**(3/2))*((1/(2*z))*(cz - (3/2)*(sz/cz)) + (3/4)*((sz**2)/cz)) + (A/8)*(3*(sz/cz)*np.sqrt(yz) + A*np.sqrt(cz/yz))
        return dFz

def lambert_vs_from_rs(r1,r2,r1mag,r2mag,z,A,mu):
    #//Calculate Lagrange Functions
    yz = lambert_y_z(r1mag,r2mag,A,z)

    lag_f = 1 - yz/r1mag
    lag_g = A*np.sqrt(yz/mu)
    lag_gdot = 1 - yz/r2mag
    #//

    #//Calculate Velocities at both postions
    v1 = (1/lag_g)*(r2 - lag_f*r1)
    v2 = (1/lag_g)*(lag_gdot*r2 - r1)
    return v1,v2
    #//