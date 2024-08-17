import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.sin(1/x)

def plotquad(a,b,step):

    xv = np.arange(a,b,step)
    yv = f(xv)
    plt.plot(xv,yv)
    plt.grid()
    plt.show()

    res,err = sp.integrate.quad(f, 0, 1)
    return res,err

a = -1
b = 0
step = 0.0000001

res,err = plotquad(a,b,step)

print('I =',res,' +-',err)