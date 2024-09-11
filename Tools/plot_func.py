import matplotlib.pyplot as plt
import numpy as np
from int_props import tb_prop

#Plot functions
def Orbit_Plot(ax,R_p,rs,ttotal,N,tr_col,cb_col,wireframe=True,legend=True,show_cb=True,line=True):
    '''
    #Define Plot
    plt.style.use('bmh')
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(projection = '3d')
    #ax = plt.axes(projection="3d")
    '''
    #Trajectory

    if line == True:
        ax.plot(rs[:,0], rs[:,1], rs[:,2],color=tr_col)
    else:
        ax.scatter(rs[:,0], rs[:,1], rs[:,2],color=tr_col)
    ax.plot(rs[0,0],rs[0,1],rs[0,2],marker='X',color=tr_col,label='Initial Position',markersize=10)
    th = str(round(ttotal/(60*60),3))
    ax.plot(rs[(-1),0],rs[(-1),1],rs[(-1),2],marker='X',color='g',label='Position after '+th+' Hours',markersize=10)
    #ax.text(rs[0,0],rs[0,1],rs[0,2],'Initial Position',color='r')
    #ax.text(rs[(n-1),0],rs[(n-1),1],rs[(n-1),2],'Position after '+th+' Hours',color='g')

    #Axis and Origin
    ox,oy,oz = [[0,0,0],[0,0,0],[0,0,0]] #Origin point
    axu,axv,axw = [[R_p,0,0],[0,R_p,0],[0,0,R_p]]
    ax.quiver(ox,oy,oz,axu,axv,axw,length=3,color='#4d4d4d',arrow_length_ratio=0.1)
    ax.text(3*R_p,0,0,'X',color='#4d4d4d')
    ax.text(0,3*R_p,0,'Y',color='#4d4d4d')
    ax.text(0,0,3*R_p,'Z',color='#4d4d4d')
    ax.set_title(("Orbit Plot"))

    if legend == True:
        ax.legend()

    if show_cb == True:
        #Primary Body
        _u,_v = np.mgrid[0:2*np.pi:20j,0:np.pi:20j]
        _x = R_p * np.cos(_u)*np.sin(_v)
        _y = R_p * np.sin(_u)*np.sin(_v)
        _z = R_p * np.cos(_v)
        if wireframe == True:
            ax.plot_wireframe(_x,_y,_z,color=cb_col)
        elif wireframe == False:
            ax.plot_surface(_x,_y,_z,color=cb_col)

    ax.set_aspect('equal')

def Traj_axes(ax,scale,rs,ttotal,n):
    
    #Trajectory
    ax.plot(rs[:,0], rs[:,1], rs[:,2],color='r')
    ax.plot(rs[0,0],rs[0,1],rs[0,2],marker='X',color='r',label='Initial Position',markersize=10)
    th = str(round(ttotal/(60*60),3))
    ax.plot(rs[(n-1),0],rs[(n-1),1],rs[(n-1),2],marker='X',color='g',label='Position after '+th+' Hours',markersize=10)
    #ax.text(rs[0,0],rs[0,1],rs[0,2],'Initial Position',color='r')
    #ax.text(rs[(n-1),0],rs[(n-1),1],rs[(n-1),2],'Position after '+th+' Hours',color='g')

    #Axis and Origin
    ox,oy,oz = [[0,0,0],[0,0,0],[0,0,0]] #Origin point
    axu,axv,axw = [[scale,0,0],[0,scale,0],[0,0,scale]]
    ax.quiver(ox,oy,oz,axu,axv,axw,length=1,color='#4d4d4d',arrow_length_ratio=0.1)
    ax.text(scale,0,0,'X',color='#4d4d4d')
    ax.text(0,scale,0,'Y',color='#4d4d4d')
    ax.text(0,0,scale,'Z',color='#4d4d4d')

    ax.set_aspect('equal')

def Traj(ax,rs,ttotal,N,tr_col):

    #Trajectory
    ax.plot(rs[:,0], rs[:,1], rs[:,2],color=tr_col)
    ax.plot(rs[0,0],rs[0,1],rs[0,2],marker='X',color='r',label='Initial Position',markersize=10)
    th = str(round(ttotal/(60*60),3))
    ax.plot(rs[(N-1),0],rs[(N-1),1],rs[(N-1),2],marker='X',color='g',label='Position after '+th+' Hours',markersize=10)

    ax.set_aspect('equal')

def elem_plot(elems,ts,hours=True):
    if hours == True:
        ts = ts/3600
        x_label = 'Time (Hours)'
    else:
        ts = ts/(3600*24)
        x_label = 'Time (Days)'
    #Element Plots
    elemfig, elemaxs = plt.subplots(2, 3, sharex='all', figsize=(16,8))
    elemfig.suptitle('Orbital Elements Delta', fontsize=16)
    plt.subplots_adjust(left=0.04,right=0.96,bottom=0.075,top=0.9,wspace=0.22,hspace=0.275)

    elemaxs[0,0].grid()
    elemaxs[0,0].plot((ts), (elems[:,0]),'b')
    elemaxs[0,0].set_title('Semi-Major Axis (a)')
    elemaxs[0,0].set_xlabel(x_label)
    elemaxs[0,0].set_ylabel('(km)')

    elemaxs[0,1].grid()
    elemaxs[0,1].plot((ts), (elems[:,1]),'b')
    elemaxs[0,1].set_title('Eccentricity (e)')
    elemaxs[0,1].set_xlabel(x_label)

    elemaxs[0,2].grid()
    elemaxs[0,2].plot((ts), np.rad2deg(elems[:,2]),'b')
    elemaxs[0,2].set_title('True Anomaly (Theta)')
    elemaxs[0,2].set_xlabel(x_label)
    elemaxs[0,2].set_ylabel('(Degrees)')

    elemaxs[1,0].grid()
    elemaxs[1,0].plot((ts), np.rad2deg(elems[:,3]),'r')
    elemaxs[1,0].set_title('Right ascension of the Acsending Node (Omega)')
    elemaxs[1,0].set_xlabel(x_label)
    elemaxs[1,0].set_ylabel('(Degrees)')

    elemaxs[1,1].grid()
    elemaxs[1,1].plot((ts), np.rad2deg(elems[:,4]),'r')
    elemaxs[1,1].set_title('Inclination (i)')
    elemaxs[1,1].set_xlabel(x_label)
    elemaxs[1,1].set_ylabel('(Degrees)')

    elemaxs[1,2].grid()
    elemaxs[1,2].plot((ts), np.rad2deg(elems[:,5]),'r')
    elemaxs[1,2].set_title('Argument of Perigee (omega)')
    elemaxs[1,2].set_xlabel(x_label)
    elemaxs[1,2].set_ylabel('(Degrees)')

    return elemfig, elemaxs

def elem_plot_d(elems,ts):
#Element Plots
    elemfig, elemaxs = plt.subplots(2, 3, sharex='all', figsize=(16,8))
    elemfig.suptitle('Orbital Elements Delta', fontsize=16)
    plt.subplots_adjust(left=0.04,right=0.96,bottom=0.075,top=0.9,wspace=0.22,hspace=0.275)

    elemaxs[0,0].grid()
    elemaxs[0,0].plot((ts/60**2), (elems[:,0]-elems[0,0]),'b')
    elemaxs[0,0].set_title('Semi-Major Axis (a)')
    elemaxs[0,0].set_xlabel('Time (Hours)')
    elemaxs[0,0].set_ylabel('(km)')

    elemaxs[0,1].grid()
    elemaxs[0,1].plot((ts/60**2), (elems[:,1]-elems[0,1]),'b')
    elemaxs[0,1].set_title('Eccentricity (e)')
    elemaxs[0,1].set_xlabel('Time (Hours)')

    elemaxs[0,2].grid()
    elemaxs[0,2].plot((ts/60**2), np.rad2deg(elems[:,2]),'b')
    elemaxs[0,2].set_title('True Anomaly (Theta)')
    elemaxs[0,2].set_xlabel('Time (Hours)')
    elemaxs[0,2].set_ylabel('(Degrees)')

    elemaxs[1,0].grid()
    elemaxs[1,0].plot((ts/60**2), np.rad2deg(elems[:,3]-elems[0,3]),'r')
    elemaxs[1,0].set_title('Right ascension of the Acsending Node (Omega)')
    elemaxs[1,0].set_xlabel('Time (Hours)')
    elemaxs[1,0].set_ylabel('(Degrees)')

    elemaxs[1,1].grid()
    elemaxs[1,1].plot((ts/60**2), np.rad2deg(elems[:,4]-elems[0,4]),'r')
    elemaxs[1,1].set_title('Inclination (i)')
    elemaxs[1,1].set_xlabel('Time (Hours)')
    elemaxs[1,1].set_ylabel('(Degrees)')

    elemaxs[1,2].grid()
    elemaxs[1,2].plot((ts/60**2), np.rad2deg(elems[:,5]-elems[0,5]),'r')
    elemaxs[1,2].set_title('Argument of Perigee (omega)')
    elemaxs[1,2].set_xlabel('Time (Hours)')
    elemaxs[1,2].set_ylabel('(Degrees)')

    return elemfig, elemaxs

def Kep_Prop_Plot(ax,cb,ttotal,r0,v0,tr_col='r',step_div=1000,legend=True,show_cb=True):

    tstep = ttotal/step_div
    N = int(np.ceil(ttotal/tstep) + 1)
    rs,vs = tb_prop(cb,N,tstep,r0,v0)

    Orbit_Plot(ax,cb['radius'],rs,ttotal,N,tr_col,cb_col=cb['col'],legend=legend,show_cb=show_cb)

def Kep_traj_Prop_Plot(ax,cb,ttotal,r0,v0,tr_col='r',step_div=1000,legend=True):

    tstep = ttotal/step_div
    N = int(np.ceil(ttotal/tstep) + 1)
    rs,vs = tb_prop(cb,N,tstep,r0,v0)

    Traj(ax,rs,ttotal,N,tr_col)