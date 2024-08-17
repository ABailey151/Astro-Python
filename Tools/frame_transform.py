import numpy as np

def XYZ_to_rsw(Vec_XYZ,i,lon_an,pearg,ta):
    i = np.deg2rad(i)
    lon_an = np.deg2rad(lon_an)
    pearg = np.deg2rad(pearg)
    ta = np.deg2rad(ta)

    #Calculate the Transformation Matrix from the perifocal frame to the geo equ coords
    QXr = np.zeros((3,3))

    QXr[0,0] = -(np.sin(lon_an)*np.cos(i)*np.sin(pearg)) + (np.cos(lon_an)*np.cos(pearg))
    QXr[0,1] = -(np.sin(lon_an)*np.cos(i)*np.cos(pearg)) - (np.cos(lon_an)*np.sin(pearg))
    QXr[0,2] = np.sin(lon_an)*np.sin(i)

    QXr[1,0] = (np.cos(lon_an)*np.cos(i)*np.sin(pearg)) + (np.sin(lon_an)*np.cos(pearg))
    QXr[1,1] = (np.cos(lon_an)*np.cos(i)*np.cos(pearg)) - (np.sin(lon_an)*np.sin(pearg))
    QXr[1,2] = -np.cos(lon_an)*np.sin(i)

    QXr[2,0] = np.sin(i)*np.sin(pearg)
    QXr[2,1] = np.sin(i)*np.cos(pearg)
    QXr[2,2] = np.cos(i)

    #Transform the perifical vectors in the geocentric frame
    rx = QXr[0,0]*Vec_XYZ[0] + QXr[0,1]*Vec_XYZ[1] + QXr[0,2]*Vec_XYZ[2]
    ry = QXr[1,0]*Vec_XYZ[0] + QXr[1,1]*Vec_XYZ[1] + QXr[1,2]*Vec_XYZ[2]
    rz = QXr[2,0]*Vec_XYZ[0] + QXr[2,1]*Vec_XYZ[1] + QXr[2,2]*Vec_XYZ[2]

    Vec_rsw = np.array([rx, ry, rz])

    return Vec_rsw

def int_to_fix_rs(phi,dphi,rs,tstep,n):
    j = 0
    rsfixed = np.zeros((n,3))
    Rphi = np.zeros((3,3))

    while j < n:

        phi = phi + (tstep*dphi)

        #Define Rotation Matrix
        Rphi[0,0] = np.cos(phi)
        Rphi[0,1] = np.sin(phi)
        Rphi[0,2] = 0
        Rphi[1,0] = -np.sin(phi)
        Rphi[1,1] = np.cos(phi)
        Rphi[1,2] = 0
        Rphi[2,0] = 0
        Rphi[2,1] = 0
        Rphi[2,2] = 1

        rsfixed[j,0] = Rphi[0,0]*rs[j,0] + Rphi[0,1]*rs[j,1] + Rphi[0,2]*rs[j,2]
        rsfixed[j,1] = Rphi[1,0]*rs[j,0] + Rphi[1,1]*rs[j,1] + Rphi[1,2]*rs[j,2]
        rsfixed[j,2] = Rphi[2,0]*rs[j,0] + Rphi[2,1]*rs[j,1] + Rphi[2,2]*rs[j,2]

        j = j + 1
    
    return rsfixed

def int_to_fix(phi,r):

    Rphi = np.zeros((3,3))
    rfixed = np.zeros((3))

    #Define Rotation Matrix
    Rphi[0,0] = np.cos(phi)
    Rphi[0,1] = np.sin(phi)
    Rphi[0,2] = 0
    Rphi[1,0] = -np.sin(phi)
    Rphi[1,1] = np.cos(phi)
    Rphi[1,2] = 0
    Rphi[2,0] = 0
    Rphi[2,1] = 0
    Rphi[2,2] = 1

    rfixed[0] = Rphi[0,0]*r[0] + Rphi[0,1]*r[1] + Rphi[0,2]*r[2]
    rfixed[1] = Rphi[1,0]*r[0] + Rphi[1,1]*r[1] + Rphi[1,2]*r[2]
    rfixed[2] = Rphi[2,0]*r[0] + Rphi[2,1]*r[1] + Rphi[2,2]*r[2]

    return rfixed

