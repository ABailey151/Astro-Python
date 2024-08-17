import numpy as np

##State vector conversions

#Individual Elements
def ang_mom(r,v):
    #Calculate Specific Angular Momentum and its magnitude
    h = np.cross(r,v)
    hmag = np.linalg.norm(h)
    return h,hmag

def inclin(h,hmag):
    i = (np.arccos(h[2]/hmag))
    return i

def L_of_N(h):
    K = np.array([0, 0, 1]) #ECI Z axis

    N = np.cross(K,h)
    Nmag = np.linalg.norm(N)
    return N,Nmag

def LON_of_an(i,N,Nmag):
    if i == 0 or i == 180:
        lon_an = 0
    else:
        if np.squeeze(N[1]) >= 0:
            lon_an = (np.arccos(N[0]/Nmag))
        else:
            lon_an = (2*np.pi) - (np.arccos(N[0]/Nmag))
    return lon_an

def eccen(u,r,v):
    rmag = np.linalg.norm(r)
    vmag = np.linalg.norm(v)
    vr = (np.dot(r,v))/rmag

    l = (vmag**2)-(u/rmag)
    e = (1/u)*(l*r - (rmag*vr)*v)
    emag = np.linalg.norm(e)
    return e,emag

def arg_of_pe(i,e,emag,N,Nmag):
    I = np.array([1, 0, 0]) #ECI X axis

    if i == 0 or i == 180:
        pearg = 0
    else:
        if e[2] >= 0:
            pearg = (np.arccos((np.dot(N,e))/(Nmag*emag)))
        else:
            pearg = (2*np.pi) - (np.arccos((np.dot(N,e))/(Nmag*emag)))
    return pearg

def true_anom(r,v,e,emag):
    rmag = np.linalg.norm(r)
    vmag = np.linalg.norm(v)
    vr = (np.dot(r,v))/rmag

    if vr >= 0:
        ta = (np.arccos(np.dot((e/emag),(r/rmag))))
    else:
        ta = (2*np.pi) - (np.arccos(np.dot((e/emag),(r/rmag))))
    return ta

def sma(mu,hmag,emag):
    #Calulate apogee, perigee and Semi-major axis
    rp = ((hmag**2)/mu)*(1/(1+emag))
    ra = ((hmag**2)/mu)*(1/(1-emag))
    a = 0.5*(rp+ra)
    return a

def orb_period(cb,a):
    T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))
    return T

#Full conversions
def coe_from_sv(mu,r,v):
    '''
    Return order, a,h,hmag,i,lon_an,e,emag,pearg,ta

    **Angles are in Radians**
    '''
    #Calculate Specific Angular Momentum and its magnitude
    h, hmag = ang_mom(r,v)

    #Calculate Inclination
    i = inclin(h,hmag)

    #Calculate Line of Nodes Vector, Perp to h and k
    N,Nmag = L_of_N(h)

    #Calculate lon_an
    lon_an = LON_of_an(i,N,Nmag)

    #Calulate the eccentricity vetor and its magnitude
    e,emag = eccen(mu,r,v)

    #Calculate the argument of Perigee
    pearg = arg_of_pe(i,e,emag,N,Nmag)

    #Calculate the true anomaly
    ta = true_anom(r,v,e,emag)

    #Calulate apogee, perigee and Semi-major axis
    a = sma(mu,hmag,emag)

    return a,h,hmag,i,lon_an,e,emag,pearg,ta

def sv_to_elems(cb,rs,vs,N):

    elems = np.zeros([N,6])

    j = 0
    
    while j < N:
        _, _, elems[j,0], elems[j,4], elems[j,3], _, elems[j,1], elems[j,5], elems[j,2] = coe_from_sv(cb['mu'],rs[j],vs[j])
        j = j + 1
    
    return elems

def sv_from_coe(u,emag,hmag,i,lon_an,pearg,ta):
    
    #Calculate the Position vetor in the Perifocal Frame
    rpxcoef = ((hmag**2)/u) * (1/(1+(emag*np.cos(ta))))
    rpx = rpxcoef * np.array([np.cos(ta), np.sin(ta), 0])

    #Calculate the velocity vector in the Perifocal frame
    vpx = (u/hmag) * np.array([-np.sin(ta), (emag+np.cos(ta)), 0])

    #Calculate the Transformation Matrix from the perifocal frame to the geo equ coords
    QXx = np.zeros((3,3))

    QXx[0,0] = -(np.sin(lon_an)*np.cos(i)*np.sin(pearg)) + (np.cos(lon_an)*np.cos(pearg))
    QXx[0,1] = -(np.sin(lon_an)*np.cos(i)*np.cos(pearg)) - (np.cos(lon_an)*np.sin(pearg))
    QXx[0,2] = np.sin(lon_an)*np.sin(i)

    QXx[1,0] = (np.cos(lon_an)*np.cos(i)*np.sin(pearg)) + (np.sin(lon_an)*np.cos(pearg))
    QXx[1,1] = (np.cos(lon_an)*np.cos(i)*np.cos(pearg)) - (np.sin(lon_an)*np.sin(pearg))
    QXx[1,2] = -np.cos(lon_an)*np.sin(i)

    QXx[2,0] = np.sin(i)*np.sin(pearg)
    QXx[2,1] = np.sin(i)*np.cos(pearg)
    QXx[2,2] = np.cos(i)

    #Transform the perifical vectors in the geocentric frame
    rx = QXx[0,0]*rpx[0] + QXx[0,1]*rpx[1] + QXx[0,2]*rpx[2]
    ry = QXx[1,0]*rpx[0] + QXx[1,1]*rpx[1] + QXx[1,2]*rpx[2]
    rz = QXx[2,0]*rpx[0] + QXx[2,1]*rpx[1] + QXx[2,2]*rpx[2]

    r = np.array([rx, ry, rz])

    vx = QXx[0,0]*vpx[0] + QXx[0,1]*vpx[1] + QXx[0,2]*vpx[2]
    vy = QXx[1,0]*vpx[0] + QXx[1,1]*vpx[1] + QXx[1,2]*vpx[2]
    vz = QXx[2,0]*vpx[0] + QXx[2,1]*vpx[1] + QXx[2,2]*vpx[2]

    v = np.array([vx, vy, vz])

    return r,v
##