import numpy as np

#Geocentric functions
def geo_r_from_lat_lon(R_p,lat,lon):

    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    
    z = R_p * np.sin(lat)

    x = R_p * np.cos(lon)*np.cos(lat)


    if 0 < lon < np.pi:
        y = np.sqrt((R_p**2)-(x**2)-(z**2))
    elif np.pi < lon < 2*np.pi:
        y = -np.sqrt((R_p**2)-(x**2)-(z**2))
    else:
        y = 0

    geo_r = np.array([x, y, z])

    return geo_r
  
