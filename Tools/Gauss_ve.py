import numpy as np
from spher_harm import j2_rsw

def guass_dt(cb,rmag,hmag,emag,ta,i,pearg):
    p = j2_rsw(cb,rmag,i,pearg,ta)

    u = pearg + ta

    dh_dt = rmag * p[1]

    de_dt = ((hmag/cb['mu']) * np.sin(ta) * p[0]) + ((1/(cb['mu']*hmag)) * ((((hmag**2) + (cb['mu']*rmag)) * np.cos(ta)) + (cb['mu']*emag*rmag)) * p[1])

    dtta_dt = (hmag/(rmag**2)) + ((1/(emag*hmag)) * ((((hmag**2)/cb['mu']) * np.cos(ta) * p[0]) - ((rmag + ((hmag**2)/cb['mu'])) * np.sin(ta) * p[1])))

    dtlonan_dt = (rmag/(hmag*np.sin(i))) * np.sin(u) * p[2]

    di_dt = (rmag/hmag) * np.cos(u) * p[2]

    dtpearg_dt = (-(1/(emag*hmag)) * ((((hmag**2)/cb['mu']) * np.cos(ta) * p[0]) - ((rmag + ((hmag**2)/cb['mu'])) * np.sin(ta) * p[1]))) - (((rmag*np.sin(u))/(hmag*np.tan(i))) * p[2])

    return np.array([dh_dt, de_dt, dtta_dt, dtlonan_dt, di_dt, dtpearg_dt])