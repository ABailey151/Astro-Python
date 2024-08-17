from sys import path
path.append('./Tools')
import Planetary_data as pd
from sv_coe import *
from geo_func import *
from Kep_Anom import *
from plot_func import *
from int_props import *
from frame_transform import int_to_fix_rs

def make_plots(axs,drift=False):

    def axes_labels(ax):
        ax.set_xlabel('$km$')
        ax.set_ylabel('$km$')
        ax.set_zlabel('$km$')

    #Plot ECF
    ax1 = axs[0]
    pltecf = Orbit_Plot(ax1,cb['radius'],rsecf,ttotal,N,tr_col='r',cb_col=cb['col'],wireframe=True)
    #axs[0].view_init(elev=0,azim=lon_an-90)
    ax1.set_title((cb['name']+"-centered "+cb['name']+"-fixed frame"))
    axes_labels(ax1)

    #Plot ECI
    ax2 = axs[1]
    plteci = Orbit_Plot(axs[1],cb['radius'],rs,ttotal,N,tr_col='r',cb_col=cb['col'],wireframe=True)
    ax2.set_title((cb['name']+"-centered intertial frame"))
    ax2.legend()
    axes_labels(ax2)

    if drift == True:

        #Plot Drift
        rsc = rsecf - rsecf[0]
        ax3 = axs[2]
        sc = np.linalg.norm(rsc[-1])
        scc_plt = Traj_axes(ax3,sc,rsc,ttotal,N)
        ax3.set_title('Geostationary Drift')
        axes_labels(ax3)

cb = pd.earth
phi = np.rad2deg(0)

#Ground target position
lat = np.deg2rad(-89.9)
lon = np.deg2rad(0)
gt = geo_r_from_lat_lon(cb['radius'],lat,lon)

ra = cb['radius'] + 3062
rp = cb['radius'] + 300
a = 42166 #(ra+rp)/2
e = 0 #(ra-rp)/(ra+rp)
if e < 1:
    hmag = np.sqrt((cb['mu']*a)*(1-(e**2)))
else:
    hmag = np.sqrt((cb['mu']*(-a))*(1-(e**2)))
i = np.deg2rad(0.01)
lon_an = np.deg2rad(0)
pearg = np.deg2rad(0)
ta0 = np.deg2rad(0)

r0,v0 = sv_from_coe(cb['mu'],e,hmag,i,lon_an,pearg,ta0)
r0mag = np.linalg.norm(r0)
v0mag = np.linalg.norm(v0)
vr0 = (np.dot(r0,v0))/r0mag
#a,h,hmag,i,lon_an,e,emag,pearg,ta0 = coe_from_sv(cb['mu'],r0,v0)

T = ((2*np.pi)/np.sqrt(cb['mu']))*(a**(3/2))

#Set time period and step length/amount
ttotal = 10*24*(60**2) #Only use T for closed orbits
tstep = ttotal/1000
N = int(np.ceil(ttotal/tstep) + 1)
rs,vs = J2_Cowell(cb,N,tstep,r0,v0)

elems = sv_to_elems(cb,rs,vs,N)

#Convert eci coords to ecf
rsecf = int_to_fix_rs(phi,cb['dphi'],rs,tstep,N)

#Initialise time and range magnitude vectors
ts = np.zeros((N,1))
gtr = np.zeros_like(rsecf)
gtrmag = np.zeros_like(ts)

#Calculate range for each ecf position vector 
ti = 0
while ti < N:
    
    gtr[ti] = rsecf[ti] - gt
    gtrmag[ti] = np.linalg.norm(gtr)
    ts[ti] = tstep*ti

    ti = ti + 1

#convert time to minutes
tm = ts/60

elemfig,elemaxs = elem_plot(elems,ts)

fig, axs = plt.subplots(1, 3, subplot_kw=dict(projection='3d'), figsize=(18,6))
plt.subplots_adjust(left=0.01,right=0.99,bottom=0.1,top=0.9,wspace=0.2)
make_plots(axs,drift=True)
plt.show()

#np.savetxt('Range.csv',np.concatenate((ts,gtrmag),axis=1),delimiter=',')