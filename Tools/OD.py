from sys import path
path.append('./Tools')
import numpy as np
from sv_coe import coe_from_sv

def Gibbs_OD(cb,r1,r2,r3):
    '''
    Return order, a,h,hmag,i,lon_an,e,emag,pearg,ta

    **Angles are in Radians**
    '''
    r1mag = np.linalg.norm(r1)
    r2mag = np.linalg.norm(r2)
    r3mag = np.linalg.norm(r3)

    C12 = np.cross(r1,r2)
    C23 = np.cross(r2,r3)
    C31 = np.cross(r3,r1)

    N = (r1mag*C23) + (r2mag*C31) + (r3mag*C12)
    D = C12 + C23 + C31
    S = r1*(r2mag-r3mag) + r2*(r3mag-r1mag) + r3*(r1mag-r2mag)

    v2 = np.sqrt(cb['mu']/(np.linalg.norm(N)*np.linalg.norm(D))) * ((np.cross(D,r2)/r2mag) + S)

    a,h,hmag,i,lon_an,e,emag,pearg,ta = coe_from_sv(cb['mu'],r2,v2)

    return a,h,hmag,i,lon_an,e,emag,pearg,ta