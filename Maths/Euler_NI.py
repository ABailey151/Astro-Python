import numpy as np
import matplotlib.pyplot as plt

x0 = 0
y0 = -1

F = lambda x,y: 3*x + y

h = 0.04

xmax = 0.2
Nmax = int(np.ceil(xmax/0.04)+1)
i = 1

xs = np.zeros([Nmax])
xs[0] = x0
ys = np.zeros([Nmax])
ys[0] = y0

while i < Nmax:
    x = x0 + h
    y = y0 + h*F(x0,y0)

    x0 = x
    y0 = y
    
    xs[i] = x
    ys[i] = y
    i = i + 1


plt.plot(xs,ys)
plt.show()