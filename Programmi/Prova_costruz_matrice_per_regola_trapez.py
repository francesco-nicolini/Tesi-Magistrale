import numpy as np

K=21

a=np.zeros((K,K))

for i in range(0,K):

    for j in range(0,K):

        if ( (i*j==0 and i+j==K-1) or (i+j==0) or (i+j==2*(K-1)) ):

            a[i][j]=1


# è un elif quindi si possono mettere casi già previsti dall'if precedente dal momento che questi, verificando il primo if, non saranno analizzati dal secondo
        elif ( i*j==0 or i==K-1 or j==K-1 ):

            a[i][j]=2

        else:

            a[i][j]=4
print(a)
