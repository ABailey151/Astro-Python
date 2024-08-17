import numpy as np

#Solving Keplers equations for anomalies given time
def kepler_E(e,M):
    if M < np.pi:
        E = M+(e/2)

    else:
        E = M-(e/2)
    
    error = 1e-8

    ratio = float(1)

    while abs(ratio) > error:
        ratio = (E - e*np.sin(E)-M)/(1-e*np.cos(E))
        E = E - ratio

    return E

def kepler_F(e,M):
    error = 1e-8

    ratio = float(1)

    F = M/(e-1)

    while abs(ratio) > error:
        ratio = (e*np.sinh(F)-F-M)/(e*np.cosh(F)-1)
        F = F - ratio
    
    return F

def ta_from_F(e,M):
    F = kepler_F(e,M)

    a = np.sqrt((e+1)/(e-1))
    ta = (2 * np.arctan(a*np.tanh(F/2)))
    return F, ta

def ta_from_E(e,M):
    E = kepler_E(e,M)
    a = np.sqrt((1+e)/(1-e))
    if E > np.pi:
        ta = 2*np.pi - (2 * np.arctan(a*np.arctan(E/2)))
    else:
        ta = (2 * np.arctan(a*np.arctan(E/2)))
    return E, ta

def para_ta(M):
    a = np.sqrt(1+(3*M)**2)
    z = (3*M + a)**(1/3)

    ta = 2*np.arctan(z-(1/z))
    return ta

#Universal Anomaly
def KU_Cz(z):
    if z > 0:
        cz = (1 - np.cos(np.sqrt(z)))/z
        return cz

    elif z == 0:
        cz = 1/2
        return cz

    elif z < 0:
        cz = (np.cosh(np.sqrt(-z)) - 1)/(-z)
        return cz

def KU_Sz(z):
    if z > 0:
        sz = (np.sqrt(z) - np.sin(np.sqrt(z)))/np.sqrt(z)**3
        return sz
    
    elif z == 0:
        sz = 1/6
        return sz
    
    elif z < 0:
        sz = (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/np.sqrt(-z)**3
        return sz

def kepler_U(u,dt,rmag,vr,alpha):
    i = 1
    error = 1e-9
    ratio = float(1)

    X = np.sqrt(u)*abs(alpha)*dt

    while abs(ratio) > error:
        z = alpha*(X**2)

        cz = KU_Cz(z)
        sz = KU_Sz(z)

        a = (rmag*vr)/np.sqrt(u)
        b = 1 - alpha*rmag

        FX = a*(X**2)*cz + b*(X**3)*sz + rmag*X - np.sqrt(u)*dt
        FdX = a*X*(1 - alpha*(X**2)*sz) + b*(X**2)*cz + rmag

        ratio = FX/FdX
        X = X - ratio
        i = i + 1
    return X

def lag_f_g(u,X,dt,r,v):
    rmag = np.linalg.norm(r)
    vmag = np.linalg.norm(v)
    alpha = 2/rmag - (vmag**2)/u
    z = alpha*(X**2)

    f = 1 - ((X**2)/rmag)*KU_Cz(z)
    g = dt - (1/np.sqrt(u))*(X**3)*KU_Sz(z)
    R = f*r + g*v
    Rmag = np.linalg.norm(R)

    fdot = (np.sqrt(u)/(Rmag*rmag))*(alpha*(X**3)*KU_Sz(z) - X)
    gdot = 1 - ((X**2)/Rmag)*KU_Cz(z)

    return f,g,fdot,gdot

def rv_from_kepler_U(u,dt,r,v):
    rmag = np.linalg.norm(r)
    vmag = np.linalg.norm(v)
    vr = (np.dot(r,v))/rmag
    alpha = 2/rmag - (vmag**2)/u

    X = kepler_U(u,dt,rmag,vr,alpha)
    z = alpha*(X**2)

    f = 1 - ((X**2)/rmag)*KU_Cz(z)
    g = dt - (1/np.sqrt(u))*(X**3)*KU_Sz(z)
    R = f*r + g*v
    Rmag = np.linalg.norm(R)

    fdot = (np.sqrt(u)/(Rmag*rmag))*(alpha*(X**3)*KU_Sz(z) - X)
    gdot = 1 - ((X**2)/Rmag)*KU_Cz(z)
    V = fdot*r + gdot*v
    Vmag = np.linalg.norm(V)

    return R, V

def Uv_prop(cb,N,tstep,r,v):
    n = 1
    rs = np.zeros((int(N),3))
    vs = np.zeros((int(N),3))
    rs[0] = r
    vs[0] = v

    while n < N:
        
        R,V= rv_from_kepler_U(cb['mu'],tstep,r,v)

        r = R
        v = V

        vs[n] = V
        rs[n] = R
        n = n + 1
        #print(n,'/',N)
        
    return rs,vs