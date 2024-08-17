from sys import path
path.append('./Tools')
import numpy as np
from Kep_Anom import KU_Cz,KU_Sz

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
