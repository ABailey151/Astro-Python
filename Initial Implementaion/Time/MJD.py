from sys import path
path.append('./Tools')
import numpy as np 
from misc_maths import JulianD

epoch = np.array([25,3,2024,8,0,0])

J0,UT = JulianD(epoch)
JD = J0 + UT/24

MJD = JD - 2400001

print(MJD)