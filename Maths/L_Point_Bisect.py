from sys import path
path.append('./Tools')
from misc_maths import *

m1 = 5.974e24
m2 = 7.347e22
pi2 = m2/(m1+m2)

f = lambda S: (1-pi2)*((S + pi2)/((abs(S + pi2))**3)) + pi2*((S + pi2 - 1)/((abs(S + pi2 - 1))**3)) - S
tol = 1e-6

#L1
a1 = 0.7
b1 = 0.9
L1 = bisection(f,a1,b1,tol)
print("\nL1 =",L1)


#L2
a2 = 1
b2 = 1.2
L2 = bisection(f,a2,b2,tol)
print("L2 =",L2)


#L3
a3 = -1.1
b3 = -0.9
L3 = bisection(f,a3,b3,tol)
print("L3 =",L3,"\n")