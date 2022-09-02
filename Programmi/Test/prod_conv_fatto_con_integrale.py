import numpy as np
import math
import matplotlib.pyplot as plt



K= 500

mu= 5
sigma= 2
A= 10





def f_m(m, A, mu, sigma):

    #return (1/np.sqrt(2*math.pi*sigma**2))*np.exp(-(m-mu)**2/(2*sigma**2))
    return A*np.exp(-(m-mu)**2/(2*sigma**2))




masse= np.linspace(0, 30, K)
dM= masse[1] - masse[0]


val_conv= masse

prod_conv= np.zeros(len(val_conv))



for i in range(0, len(val_conv)):

    integrale= 0

    for j in range(0, len(masse)):

        prod= f_m(masse[j], A, mu, sigma)*f_m(val_conv[i]-masse[j], A, mu, sigma)

        if( j==0 or j==(K-1) ):
            integrale+= dM*prod/2

        else:
            integrale+= dM*prod

    prod_conv[i]= integrale


figure= plt.figure()

plt.plot(masse, f_m(masse, A, mu, sigma), linestyle="-", marker="", color="red", label="funzione di partenza")
plt.plot(val_conv, prod_conv, linestyle="-", marker="", color="blue", label="prodotto di convoluzione")


plt.legend()

plt.show()