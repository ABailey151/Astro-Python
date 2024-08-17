import numpy as np
import math

a = np.array([1,6,18])
b = np.array([42,-69,98])

adb = np.dot(a,b);

amag = np.linalg.norm(a)
bmag = np.linalg.norm(b)

#Part A: find the agle between a and b

abrads = np.arccos(adb/(amag*bmag))
abtheta = np.rad2deg(abrads)

print('The angle between a and b is', abtheta, 'degrees')

#Part B: Project b onto A

Ba = adb/amag

print('The component of b in the direction of a is', Ba)

#Part C: Project A onto B

Ab = adb/bmag

print('The component of a in the direction of b is', Ab)