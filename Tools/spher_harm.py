import numpy as np

def j2_XYZ(cb,r):
    rmag = np.linalg.norm(r)
        
    j2c = (3/2) * ((cb['j2']*cb['mu']*(cb['radius']**2))/(rmag**4))
    aj2x = j2c * (r[0]/rmag) * ((5*((r[2]**2)/(rmag**2))) - 1)
    aj2y = j2c * (r[1]/rmag) * ((5*((r[2]**2)/(rmag**2))) - 1)
    aj2z = j2c * (r[2]/rmag) * ((5*((r[2]**2)/(rmag**2))) - 3)
    aj2 = np.array([aj2x, aj2y, aj2z])

    return aj2

def j2_rsw(cb,rmag,i,pearg,ta):
    u = pearg + ta

    j2c = -(3/2) * ((cb['j2']*cb['mu']*(cb['radius']**2))/(rmag**4))
    j2_r = j2c * (1 - (3 * (np.sin(i)**2) * (np.sin(u)**2)))
    j2_s = j2c * (np.sin(i)**2) * np.sin(2 * u)
    j2_w = j2c * np.sin(2 * i) * np.sin(u)
    j2_rsw = np.array([j2_r, j2_s, j2_w])

    return j2_rsw

def j3(cb,r):
    rmag = np.linalg.norm(r)

    j3c = ((cb['j3']*cb['mu']*(cb['radius']**3))/(2*rmag**5))

    aj3x = j3c * (r[0]/rmag) * 5 * ((7*((r[2]/rmag)**3)) - (3*(r[2]/rmag)))
    aj3y = j3c * (r[1]/rmag) * 5 * ((7*((r[2]/rmag)**3)) - (3*(r[2]/rmag)))
    aj3z = j3c * (r[2]/rmag) * 3 * (((35/3)*((r[2]/rmag)**4)) - (10*((r[2]/rmag)**2)) + 1)
    aj3 = np.array([aj3x, aj3y, aj3z])

    return aj3
    