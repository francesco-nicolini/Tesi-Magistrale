import numpy as np
import matplotlib.pyplot as plt

'''
a=np.array( [(1,2),(1,3)] )
b=np.array([1,0])

print(a,"\n",b)
print(np.dot(a,b))

c=np.linalg.inv(a)

print(a)
print(c)
print(np.dot(a,c))
'''

x_0=1
y_0=2
N=100

def f_1(x,y):
    return x**2+y**2

def der_f(x):
    return 2*x

def f_2(x,y):
    return x**2-y**2

funz=[f_1,f_2]

#print(funzioni[0](1,1))

F=np.zeros(2)
J=np.zeros((2,2))
var=np.array([x_0, y_0])
#print(J)
#print(funz[0](*var))

for k in range(0,N):

    for i in range(0,len(F)):

        F[i]=funz[i](*var)

        for j in range(0,len(J[i])):

            J[j][i]=der_f(var[i])
            if (i==1 and j==1):
                J[j][i]=-J[j][i]

    #print(F)
    #print(J)

    var=var-np.dot(np.linalg.inv(J),F)

print(var)


