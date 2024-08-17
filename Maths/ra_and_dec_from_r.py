from sys import path
path.append('./Tools')
import numpy as np

def ra_and_dec_from_r(r):

    rmag = np.linalg.norm(r)

    l = r[0]/rmag
    m = r[1]/rmag
    n = r[2]/rmag

    dec = np.arcsin(n)
    decdeg = np.rad2deg(dec)

    if m > 0:
        ra = np.arccos((l/np.cos(dec)))
        radeg = np.rad2deg(ra)

    elif m <= 0:
        ra = 2*np.pi - np.arccos((l/np.cos(dec)))
        radeg = np.rad2deg(ra)
    
    return ra, dec

r = np.array([-3000, -6000, -9000])

ra, dec = ra_and_dec_from_r(r)

decdeg = np.rad2deg(dec)
radeg = np.rad2deg(ra)

print("\nDeclination =",dec,"RADs /",decdeg,"°\nRight Ascention =",ra,"RADs /",radeg,"°\n")