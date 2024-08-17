from sys import path
path.append('./Tools')
import Planetary_data as pd
from sv_coe import *
from geo_func import *
from Kep_Anom import *
from plot_func import *

cb = pd.earth
u = cb['mu']
R_e = cb['radius']

#Initial conditions
r0 = np.array([3231.6687, -5597.41438, -2352.46243]) #np.array([0,42164.1712,0])   
v0 = np.array([0, 0, 10]) #np.array([-np.sqrt((u)/np.linalg.norm(r0)),0,0])

r0mag = np.linalg.norm(r0)
v0mag = np.linalg.norm(v0)
vr0 = (np.dot(r0,v0))/r0mag

#Calculate orbital elements
a,h,hmag,i,lon_an,e,emag,pearg,ta0 = coe_from_sv(u,r0,v0)
rp = ((hmag**2)/u)*(1/(1+emag))
if emag < 1:
    ra = rp*((1+emag)/(1-emag))
print('h (km^2/s) =',hmag,'\na (km) =',a,'\ni (Deg) =',i,'\nRAAN (Deg) =',lon_an,'\ne =',emag,'\nArg of Pe (Deg) =',pearg,'\nTa (Deg) =',ta0)

if emag < 1:
    T = ((2*np.pi)/np.sqrt(u))*(a**(3/2))

#Set time period and step length/amount
ttotal = 0.5*(60**2) #Only use T for closed orbits
tstep = 60
N = int(np.ceil(ttotal/tstep) + 1)
n = 1

#Prop to specified time
rs,vs = Uv_prop(cb,N,tstep,r0,v0)

#Prop full Orbit
if emag < 1 and ttotal < T:
    rls = np.array([rs[(N-1),0],rs[(N-1),1],rs[(N-1),2]])
    rlsmag = np.linalg.norm(rls)
    vls = np.array([vs[(N-1),0],vs[(N-1),1],vs[(N-1),2]])
    vlsmag = np.linalg.norm(vls)
    Nfull = np.ceil((T-ttotal)/tstep) + 1
    rfull,vfull = Uv_prop(cb,Nfull,tstep,rls,vls)

#np.savetxt('rs.csv',np.concatenate((rs,vs),axis=1),delimiter=',')

#Define Plot
fig = plt.figure(figsize=(8,8))
ax = plt.axes(projection="3d")

#Axis and Origin
ax.plot(0,0,0,marker='X',markersize=10)
axu,axv,axw = [[R_e,0,0],[0,R_e,0],[0,0,R_e]]
ax.quiver(0,0,0,axu,axv,axw,length=3,color='#4d4d4d',arrow_length_ratio=0.1)
ax.text(3*R_e,0,0,'X',color='#4d4d4d')
ax.text(0,3*R_e,0,'Y',color='#4d4d4d')
ax.text(0,0,3*R_e,'Z',color='#4d4d4d')

#Partial Orbit
rX = rs[:,0]
rY = rs[:,1]
rZ = rs[:,2]
ax.plot(rX, rY, rZ,color='r')
ax.plot(r0[0],r0[1],r0[2],marker='X',color='r',markersize=10)
ax.plot(rs[(N-1),0],rs[(N-1),1],rs[(N-1),2],marker='X',color='r',markersize=10)
r0x,r0y,r0z = [r0[0],r0[1],r0[2]]
ax.text(r0x,r0y,r0z,'Initial Position',color='r')
th = str(round(ttotal/(60*60),3))
ax.text(rs[(N-1),0],rs[(N-1),1],rs[(N-1),2],'Position after '+th+' Hours',color='r')

#Full Orbit
if emag < 1 and ttotal < T:
    rfX = rfull[:,0]
    rfY = rfull[:,1]
    rfZ = rfull[:,2]
    ax.plot(rfX, rfY, rfZ,color='c')

#Angular Momentum
Hu,Hv,Hw = [h[0]/(5),h[1]/(5),h[2]/(5)]
ax.quiver(0,0,0,Hu,Hv,Hw,color='y',arrow_length_ratio=0.1)
ax.text(Hu,Hv,Hw,'Angular Momentum',color='y')

#Perigee
if emag == 0:
    eu,ev,ew = [r0[0],r0[1],r0[2]]
else:
    eu,ev,ew = [(rp/emag)*e[0],(rp/emag)*e[1],(rp/emag)*e[2]]
ax.quiver(0,0,0,eu,ev,ew,color='g',arrow_length_ratio=0.2)
print_rp = str(round(rp,3))
ax.text(eu,ev,ew,'Perigee: '+str(print_rp)+'km',color='g')
#Apogee
if 0 < emag < 1:
    rau,rav,raw = [-(ra/emag)*e[0],-(ra/emag)*e[1],-(ra/emag)*e[2]]
    ax.quiver(0,0,0,rau,rav,raw,color='g',arrow_length_ratio=0.2*(rp/ra))
    print_ra = str(round(ra,3))
    ax.text(rau,rav,raw,'Apogee: '+str(print_ra)+'km',color='g')

#Earth
_u,_v = np.mgrid[0:2*np.pi:20j,0:np.pi:20j]
_x = R_e * np.cos(_u)*np.sin(_v)
_y = R_e * np.sin(_u)*np.sin(_v)
_z = R_e * np.cos(_v)
ax.plot_wireframe(_x,_y,_z,color=cb['col'])

#M = (u**2)/(h**3) * ((e**2)-1)**(3/2) * ttotal
#F,ta = ta_from_F(e,M)

ax.set_aspect('equal')
plt.show()
