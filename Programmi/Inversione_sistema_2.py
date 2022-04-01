import numpy as np
import matplotlib.pyplot as plt
from sympy import *




x_0=200
y_0=200
N=20




def f_1(x,y):
    return x**2+y**2

def f_2(x,y):
    return x**2-y**2

funz=[f_1,f_2]

a=Symbol('a')
b=Symbol('b')

simb=[a,b]

'''
z_1=f_1(*simb)

print(z_1)

print(z_1.diff(simb[0]))
print(z_1.diff(simb[1]))

z=funz[1](*simb).diff(simb[1])
f=lambdify((a,b),z, 'numpy')
print(f(x_0,y_0))
'''

F=np.zeros(2)
J=np.zeros((2,2))
var=np.array([x_0, y_0])

x=np.zeros(N)
y=np.zeros(N)

for k in range(0,N):

    for i in range(0,len(F)):

        F[i]=funz[i](*var)

        for j in range(0,len(J[i])):

            z=funz[j](*simb).diff(simb[i])
            f=lambdify((a,b),z, 'numpy')

            J[j][i]=f(*var)

    #print(F)
    #print(J)

    var=var-np.dot(np.linalg.inv(J),F)

    x[k]=var[0]
    y[k]=var[1]

print(var)

iter=np.linspace(0,N-1,N)

plt.plot(iter, x, color="blue", label="x")
plt.plot(iter,y,color="red", label="y")

plt.xlabel("Numero iterazioni")
plt.ylabel("Valori variabili")
plt.legend()
plt.show()

