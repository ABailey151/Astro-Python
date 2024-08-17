from sys import path
path.append('./Tools')
import numpy as np
import Planetary_data as pd

#Import functions where needed
#from sv_coe import *
#from geo_func import *
#from Kep_Anom import *
#from plot_func import *
#from misc_maths import *

def JD(gregdate):
    """
    Epoch format, array([D, M, Y, h, m, s])
    """
        
    UT = gregdate[3] + (gregdate[4]/60) + (gregdate[5]/3600)

    J0 = (367*gregdate[2]) - int((7*(gregdate[2]+int((gregdate[1]+9)/12)))/4) + int((275*gregdate[1])/9) + gregdate[0] + 1721013.5

    return J0,UT

GD1 = np.array([4,10,1957,19,26,24])
GD2 = np.array([12,5,2004,14,45,30])

J01,UT1 = JD(GD1)
JD1 = J01 + (UT1/24)
J02,UT2 = JD(GD2)
JD2 = J02 + (UT2/24)

diff  = JD2 - JD1

print('\nJulian Date 1 =',JD1)
print('Julian Date 2 =',JD2)
print('Difference =',diff,'Days,',diff/365,'Years\n')
