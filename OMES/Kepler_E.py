import numpy as np

e = float(input('Orbit Eccentricity = '))
M = float(input('Time since Perigee Mean Anomaly = '))

if M < np.pi:
    E = M+(e/2)

else:
    E = M-(e/2)

error = float(input('Final value tolerance = '))

ratio = float(1)

while abs(ratio) > error:
    ratio = (E - e*np.sin(E)-M)/(1-e*np.cos(E))
    E = E - ratio

print('E =',E)