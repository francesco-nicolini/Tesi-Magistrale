import matplotlib.pyplot as plt
import numpy as np
from sympy import symbols, nsolve, exp


# valore minimo e massimo della massa rappresentati

minimo= 0
massimo= 10



# parametri funzione

a= -1
b= 5
c= 0

bordo= 2.5
sigma= 1




# bordo indica il valore di x a destra del quale si ha una gaussiana e a sinistra una parabola, si indicano i coefficienti della parabola (a, b, c) e la deviazione standard della gaussiana (sigma), la costante moltiplicativa e il valor medio vengono fissate imponendo continuita e continuita della derivata prima

def trova_coefficienti(a, b, c, bordo, sigma):

    mu= symbols("mu")
    D= symbols("D")

    r1= a*bordo**2 + b*bordo + c
    r2= 2*a*bordo + b

    f1= D*exp( -(bordo-mu)**2/(2*sigma**2) ) - r1
    f2= D*exp( -(bordo-mu)**2/(2*sigma**2) )*(-(bordo-mu)/sigma**2) - r2

    return nsolve( (f1, f2), (D, mu), (1, 1) )





def parab_gauss(m, bordo, a, b, c, sigma, D, mu):

    D= float(D)
    mu= float(mu)

    if( m<bordo ):
        return a*m**2+b*m+c

    else:
        return D*np.exp( -(m-mu)**2/(2*sigma**2) )



def parabola (m, a, b, c):

    return a*m**2 + b*m + c



def gaussiana (m, sigma, D, mu):

    D= float(D)
    mu= float(mu)

    return D*np.exp( -(m-mu)**2/(2*sigma**2) )






D, mu= trova_coefficienti(a, b, c, bordo, sigma)


print(D, mu)

masse= np.linspace(minimo, massimo, 1000)



distrib= np.zeros( len(masse) )

for i in range(0, len(distrib)):

    distrib[i]= parab_gauss(masse[i], bordo, a, b, c, sigma, D, mu)


parab= parabola(masse, a, b, c)




gauss= gaussiana(masse, sigma, D, mu)



plt.figure()

plt.title("Distribuzione Parabola e Gaussiana")

plt.plot(masse, parab, linestyle=(0, (4,7)), marker="", color= "green", label="parabola")
plt.plot(masse, gauss, linestyle=(0, (4,7)), marker="", color= "orange", label="gaussiana")
plt.plot(masse, distrib, linestyle="-", marker="", color= "navy", label="distribuzione")

plt.plot([minimo, massimo], [0, 0], marker="", linestyle="-", color="black", linewidth=0.8)

plt.xlabel("masse [U.A.]")
plt.xlabel("distribuzione [U.A.]")

plt.xlim(minimo, massimo)

plt.legend()

plt.show()












