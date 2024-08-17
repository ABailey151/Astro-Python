from sys import path
path.append('./Tools')
#from Planetary_data import *
import Planetary_data as pd
from sv_coe import *
from geo_func import *
from int_props import tb_prop
from plot_func import *
from misc_maths import *
from frame_transform import int_to_fix_rs
import scipy as sp

#Define central body and t0 angle
cb = pd.moon
phi0 = 0 

#Ground target position
lat = -49.9
lon = 0
gt = geo_r_from_lat_lon(cb['radius'],lat,lon)

#Orbit Initial conditions
ra = cb['radius'] + 150
rp = cb['radius'] + 80
a = (ra+rp)/2
emag = (ra-rp)/(ra+rp)
if emag < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(emag**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(emag**2)))
i = np.deg2rad(90)
lon_an = np.deg2rad(0)
pearg = np.deg2rad(270)
ta0 = np.deg2rad(180)

r0,v0 = sv_from_coe(cb['mu'],emag,hmag,i,lon_an,pearg,ta0)
a,h,hmag,i,lon_an,e,emag,pearg,ta0 = coe_from_sv(cb['mu'],r0,v0)

T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))

#Set time period and step length/amount
ttotal = T*20 #Only use T for closed orbits
tstep = 10
N = int(np.ceil(ttotal/tstep))
rs,vs = tb_prop(cb,N,tstep,r0,v0)

#Convert eci coords to ecf
rsecf = int_to_fix_rs(phi0,cb['dphi'],rs,tstep,N)
vsecf = int_to_fix_rs(phi0,cb['dphi'],vs,tstep,N)
hecf = np.zeros_like(rsecf)
hi = 0
while hi < N:
    hecf[hi] = np.cross(rsecf[hi],vsecf[hi])
    hi = hi + 1

#Initialise arrays
ts =  np.linspace(0,ttotal,N) #time
gtr = np.zeros_like(rsecf) #range vectors
gtrmag = np.zeros_like(ts) #range mag array
b = np.zeros_like(ts) #pointing angle

#Calculate range for each ecf position vector 
ti = 0
while ti < N:
    
    gtr[ti] = rsecf[ti] - gt
    gtrmag[ti] = np.linalg.norm(gtr[ti])
    b[ti] = np.rad2deg(np.arccos(np.dot(rsecf[ti],gtr[ti])/(gtrmag[ti]*np.linalg.norm(rsecf[ti])))) - 90

    ti = ti + 1

####################################################################################
#Plots

def make_plots(axs):

    def axes_labels(ax):
        ax.set_xlabel('$km$')
        ax.set_ylabel('$km$')
        ax.set_zlabel('$km$')

    #Plot ECF
    ax1 = axs[0]
    pltecf = Orbit_Plot(ax1,cb['radius'],rsecf,ttotal,N,tr_col='r',cb_col=cb['col'],wireframe=True)
    ax1.plot(gt[0],gt[1],gt[2],marker='X',markersize=10,color='b',label='Ground Target')
    ax1.view_init(elev=0,azim=lon_an-90)
    ax1.set_title((cb['name'] +"-Centered "+cb['name']+"-Fixed Frame"))
    axes_labels(ax1)

    #Plot ECI
    ax2 = axs[1]
    plteci = Orbit_Plot(ax2,cb['radius'],rs,ttotal,N,tr_col='r',cb_col=cb['col'],wireframe=True)
    ax2.plot(gt[0],gt[1],gt[2],marker='X',markersize=10,color='b',label='Ground Target')
    ax2.set_title((cb['name']+"-Centered Intertial Frame"))
    ax2.legend()
    axes_labels(ax2)
Orbfig, Orbaxs = plt.subplots(1, 2, subplot_kw=dict(projection='3d'), figsize=(12,6))
plt.subplots_adjust(left=0.05,right=0.95,bottom=0.06,top=0.930,wspace=0.2)
make_plots(Orbaxs)

#Create the range plot
range_plt = plt.figure(figsize=(8,8))
ax2d1 = range_plt.add_subplot(211)
ax2d1.grid()
ax2d1.plot((ts/60), gtrmag,'r')
plt.title('Spacecraft Range to Ground Target (Pe = '+str(rp-cb['radius'])+'km, \u03C9 ='+str((pearg))+'\u00b0)')
plt.xlabel('Time (min)')
plt.ylabel('Range (km)')

#Angle Plot
ax2d2 = range_plt.add_subplot(212,sharex=ax2d1)
ax2d2.grid()
ax2d2.plot((ts/60), b,'b')
plt.title('Angle from Spacecraft Local Horizon to Ground Target')
plt.xlabel('Time (min)')
plt.ylabel('Pointing Angle (Degrees)')
plt.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.94,hspace=0.25)
plt.show()

#np.savetxt('Pointing.csv',np.concatenate((ts,b),axis=1),delimiter=',')
