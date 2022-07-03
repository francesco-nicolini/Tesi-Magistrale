import numpy as np
import matplotlib.pyplot as plt
import math


N= 1000

mu= 0
sigma= 1

m_min= -20
m_max= 20


def gaussiana(m, mu, sigma):

    return (1/np.sqrt(2*math.pi*sigma**2))*np.exp(-(m-mu)**2/(2*sigma**2))



masse= np.linspace(m_min, m_max, N)
f_m= gaussiana(masse, mu, sigma)


shift= np.zeros(len(f_m))

for i in range(0, int( len(shift)/2 )):

    shift[i]= f_m[ int( len(shift)/2 ) + i ]

for i in range(int( len(shift)/2 ), len(shift)):

    shift[i]= f_m[ i - int( len(shift)/2 )]




fig= plt.figure()

plt.title("Applicazione Condizioni di Periodicità")

plt.xlabel("m [u.a.]")
plt.ylabel("f(m) [u.a.]")

plt.plot(masse, f_m, linestyle="-", marker="", color="blue", label="Segnale Originale")
plt.plot(masse, shift, linestyle="-", marker="", color="red", label="Segnale con  Condizioni\nPeriodiche")


plt.xlim(masse[0], masse[-1])


plt.legend()
plt.tight_layout()





# CREAZIONE DI UNA FUNZIONE CHE SHIFTA L'ASSE X
def asse_shift( x, K):

    return x + np.sign( int(K/2) - x )*int(K/2)

x= np.linspace(0, N, N)

fig= plt.figure()

plt.title("Funzione che Shifta l'Aasse delle Ascisse")

plt.plot( x, asse_shift(x, len(masse) ), linestyle="-", marker="", color="blue")

plt.xlabel("m [u.a.]")
plt.ylabel("Funzione che Shifta [u.a.]")

plt.xlim(x[0], x[-1])

plt.tight_layout()

plt.show()

















'''
# ESPERIMENTO SECONDARIO: STUDIO FUNZIONE np.fft.fftshift(). NON CENTRA IN REALTA NULLA CON L'IMPOSIZIONE DELLE CONDZIONI AL CONTORNO

masse_shift= np.fft.fftshift(masse)


fig= plt.figure()

plt.title("Confronto tra Asse x Normale e Asse x Su Cui E' Applicato np.fft.fftshift()")

plt.xlabel("m [u.a.]")
plt.ylabel("f(m) [u.a.]")

plt.plot(masse, f_m, linestyle="-", marker="", color="blue", label="Asse x Normale")
plt.plot(masse_shift, f_m, linestyle="-", marker="", color="red", label="Asse x con np.fft.fftshift()")


index_min= np.where(masse_shift==masse[0])
index_max= np.where(masse_shift==masse[-1])


plt.plot( [masse_shift[index_min[0]], masse_shift[index_max[0]]], [f_m[index_min[0]], f_m[index_max[0]] ],  linestyle="-", marker="", color="white", linewidth=1)


plt.xlim(masse[0], masse[-1])

plt.legend()
plt.tight_layout()
plt.show()
'''