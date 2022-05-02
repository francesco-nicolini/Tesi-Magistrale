import numpy as np

K=10

M=np.ones((K,K))

for i in range(0, K):

    for j in range(0, K):

        if(j!=0 and j!=(K-1)):

            M[i][j]=2

print(M/2)


