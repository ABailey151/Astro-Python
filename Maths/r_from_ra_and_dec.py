from sys import path
path.append('./Tools')
import numpy as np

def r_from_ra_and_dec(ra,dec,rmag):
    ra_rad = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)

    ur = np.array([(np.cos(dec_rad)*np.cos(ra_rad)), (np.cos(dec_rad)*np.sin(ra_rad)), (np.sin(dec_rad))])

    r = rmag * ur

    return r

ra = 300
dec = -60
rmag = 6878.14

r = r_from_ra_and_dec(ra,dec,rmag)

print(r)
