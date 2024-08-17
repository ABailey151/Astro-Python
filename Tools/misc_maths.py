import numpy as np

#Misc maths
def bisection(f,a,b,tol):
    i = 1
    Nmax = int(((np.log((b-a)/tol))/np.log(2)) -  1)

    if f(a)*f(b) >= 0:
        print("The Bisection Method fails for these values of a and b")

    while i < (Nmax+1) or f(c) > tol:
        c = (a+b)/2
        
        if f(a)*f(c) < 0:
            b = c
        
        else:
            a = c

        i = i+1
        #print(i,c)
    return c

def ra_and_dec_from_r(r):
    '''
    **OUTPUT IS IN RADIANS**
    '''

    rmag = np.linalg.norm(r)

    l = r[0]/rmag
    m = r[1]/rmag
    n = r[2]/rmag

    dec = np.arcsin(n)
    decdeg = np.rad2deg(dec)

    if m > 0:
        ra = np.arccos((l/np.cos(dec)))
        radeg = np.rad2deg(ra)

    elif m <= 0:
        ra = 2*np.pi - np.arccos((l/np.cos(dec)))
        radeg = np.rad2deg(ra)
    
    return ra, dec

def r_from_ra_and_dec(ra,dec,rmag):
    ra_rad = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)

    ur = np.array([(np.cos(dec_rad)*np.cos(ra_rad)), (np.cos(dec_rad)*np.sin(ra_rad)), (np.sin(dec_rad))])

    r = rmag * ur

    return r

def time_arr(n,tstep):
    ts = np.zeros((n,1))
    ti = 0

    while ti < n:

        ts[ti] = tstep*ti
        ti = ti + 1
    
    return ts

#'''
def NeN(a,b):

    q = (np.dot(a,((2*b)-a))) / (np.linalg.norm(b)**2)

    fq  = (((q**2)-(3*q)+3) / (1+((1-q)**(3/2))))*q

    return fq
'''


def NeN(a,b):

    q = (np.dot(a,((2*b)-a))) / (np.linalg.norm(b)**2)
    #q = (np.dot(a,a) - np.dot((2*a),b)) / (np.linalg.norm(b)**2)

    fq = (((q**2)+(3*q)+3) / (1+((1+q)**(1.5))))*q

    return fq
'''

def JulianD(gregdate):
    """
    Epoch format, array([D, M, Y, h, m, s])
    """

    UT = gregdate[3] + (gregdate[4]/60) + (gregdate[5]/3600)

    J0 = (367*gregdate[2]) - int((7*(gregdate[2]+int((gregdate[1]+9)/12)))/4) + int((275*gregdate[1])/9) + gregdate[0] + 1721013.5

    return J0,UT