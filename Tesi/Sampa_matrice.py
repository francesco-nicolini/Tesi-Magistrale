import numpy as np

K=4



for i in range(0,K):

    for j in range(0,K):

        if ( (i*j==0 and i+j==K-1) or (i+j==0) or (i+j==2*(K-1)) ):

            testo= "C(m_{0}, m_{1},\\nu)".format(i, j)

            if(j==K-1):
                print(testo + " \\\\")

            else:
                print(testo + " &", end=" ")


# è un elif quindi si possono mettere casi già previsti dall'if precedente dal momento che questi, verificando il primo if, non saranno analizzati dal secondo
        elif ( i*j==0 or i==K-1 or j==K-1 ):

            testo= "2\,C(m_{0}, m_{1},\\nu)".format(i, j)

            if(j==K-1):
                print(testo + " \\\\")

            else:
                print(testo + " &", end=" ")

        else:

            testo= "4\,C(m_{0}, m_{1},\\nu)".format(i, j)

            if(j==K-1):
                print(testo + " \\\\")

            else:
                print(testo + " &", end=" ")

