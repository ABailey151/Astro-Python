from sys import path
path.append('./Tools')
import Planetary_data as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from plot_func import *

G = 6.6743015e-20

def f(t,y,m1,m2):
    x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2 = y

    r = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

    ax1 = G * m2 * ((x2-x1)/(r**3))
    ax2 = G * m1 * ((x1-x2)/(r**3))
    ay1 = G * m2 * ((y2-y1)/(r**3))
    ay2 = G * m1 * ((y1-y2)/(r**3))
    az1 = G * m2 * ((z2-z1)/(r**3))
    az2 = G * m1 * ((z1-z2)/(r**3))

    dydt = np.array([vx1,vy1,vz1,vx2,vy2,vz2,ax1,ay1,az1,ax2,ay2,az2])
    return dydt

R1 = np.array([0, 0, 0])
R2 = np.array([3000, 0, 0])
V1 = np.array([10, -20, 10])
V2 = np.array([0, 40, 0])
m1 = 1e26
m2 = 1e26

y0 = np.concatenate((R1,R2,V1,V2), axis=0)

#Set time period and step length/amount
ttotal = 480
tstep = 1
N = int(np.ceil(ttotal/tstep) + 1)
n = 0

Rs1 = np.zeros([N,3])
Rs2 = np.zeros([N,3])

tb = sp.integrate.ode(f,jac=None).set_integrator('lsoda')
while n < N:

    tb.set_initial_value(y0,0).set_f_params(m1,m2)
    tb.integrate(tstep)

    Rs1[n] = np.array([tb.y[0],tb.y[1],tb.y[2]])
    Rs2[n] = np.array([tb.y[3],tb.y[4],tb.y[5]])
    y0 = tb.y
    n = n + 1

##Plots
fig = plt.figure(figsize=(10,8))
axes = plt.subplot(projection='3d')
p1 = Traj_axes(axes,np.linalg.norm(R2-R1),(Rs1),ttotal,N)
p2 = Traj(axes,Rs2,ttotal,N,col='b')
plt.subplots_adjust(left=0.05,right=0.95,bottom=0.09,top=0.95)
plt.show()