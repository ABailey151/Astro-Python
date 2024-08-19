from sys import path
path.append('./Tools')
from Kep_Anom import *

e = float(input("\nEccentricity = "))
M = float(input("\nMean Anomaly = "))

if e < 1:
    E, ta = ta_from_E(e,M)
    dta = np.rad2deg(ta)
    print("\n-------------------------------------------------------------------")
    print("Eccentric Anomaly =",E,"RAD/s")
    print("True Anomaly =",ta,"RAD /",dta,"Degress")
    print("-------------------------------------------------------------------\n")

elif e ==1:
    ta = para_ta(M)
    dta = np.rad2deg(ta)
    print("\n-------------------------------------------------------------------")
    print("True Anomaly =",ta,"RAD /",dta,"Degress")
    print("-------------------------------------------------------------------\n")

else:
    F, ta = ta_from_F(e,M)
    dta = np.rad2deg(ta)
    print("\n-------------------------------------------------------------------")
    print("Eccentric Anomaly =",F,"RAD/s")
    print("True Anomaly =",ta,"RAD /",dta,"Degress")
    print("-------------------------------------------------------------------\n")
