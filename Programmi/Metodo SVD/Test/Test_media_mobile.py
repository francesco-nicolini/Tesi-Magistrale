import numpy as np
import matplotlib.pyplot as plt
import scipy.signal

N= 100000

m_min= 0
m_max= 1

lung_sottoin= 5000

wind_size= 500
poly_order= 3




masse= np.linspace(m_min, m_max, N)

F_M= np.ones(N)

random= np.random.rand(N) - 0.5

F_M= F_M + random

'''
print(F_M)
'''

# media mobile 1

F_M_medie= np.zeros( int( len(F_M)/lung_sottoin ) )
masse_medie= np.zeros( int( len(F_M)/lung_sottoin ) )

for i in range(0, int( len(F_M)/lung_sottoin ) ):

    sum= 0

    for j in range( 0, lung_sottoin):

        sum+= F_M[lung_sottoin*i + j]

    media= sum/lung_sottoin

    F_M_medie[i]= media
    masse_medie[i]= masse[lung_sottoin*i + int(j/2)]





# media mobile 2

F_M_medie_1= np.zeros( len(F_M) - lung_sottoin )
masse_medie_1= np.zeros( len(F_M) - lung_sottoin )


for i in range( int(lung_sottoin/2), len(F_M) - int(lung_sottoin/2)):

    sum= 0

    for j in range( -int(lung_sottoin/2), int(lung_sottoin/2)):

         sum+= F_M[i + j]

    media= sum/lung_sottoin

    F_M_medie_1[i-int(lung_sottoin/2)]= media
    masse_medie_1[i-int(lung_sottoin/2)]= masse[i]




# Savitzky Golay

F_M_S_V= scipy.signal.savgol_filter(F_M, wind_size, poly_order)




fig= plt.figure()

plt.plot(masse, F_M, color="blue", linestyle="-", marker="", label="Segnale Originale")
plt.plot(masse_medie, F_M_medie, color="red", linestyle="-", marker="", label="Media Mobile Metodo 1")
plt.plot(masse_medie_1, F_M_medie_1, color="yellow", linestyle="-", marker="", label="Media Mobile Metodo 2")
plt.plot(masse, F_M_S_V, color="black", linestyle="-", marker="", label="Savitzky Golay")




plt.title("Confronto tra il Segnale Originale e la Sua Media Mobile")

plt.xlabel("M [u.a.]")
plt.ylabel("F(M) [u.a.]")

plt.xlim(masse[0], masse[-1])

plt.legend()
plt.tight_layout()

plt.show()