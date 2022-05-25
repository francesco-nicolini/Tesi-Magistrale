import numpy as np


N=200
K=10

omega= np.linspace(1, 10, N)

rap= int(N/K)
res= N%K


if(res==0):

    a= omega[rap-1::rap]

else:
    a= omega[res-1+rap::rap]

print(a)
print("lung= ",len(a))
print(a[-1], omega[-1])


