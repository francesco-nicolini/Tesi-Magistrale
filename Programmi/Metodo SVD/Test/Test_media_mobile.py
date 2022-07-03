import numpy as np
import matplotlib.pyplot as plt

N= 10000

m_min= 0
m_max= 1

lung_sottoin= 500




masse= np.linspace(m_min, m_max, N)

F_M= np.ones(N)

random= np.random.rand(N) - 0.5

F_M= F_M + random

'''
print(F_M)
'''


F_M_medie= np.zeros( int( len(F_M)/lung_sottoin ) )
masse_medie= np.zeros( int( len(F_M)/lung_sottoin ) )

for i in range(0, int( len(F_M)/lung_sottoin ) ):

    sum= 0

    for j in range( 0, lung_sottoin):

        sum+= F_M[lung_sottoin*i + j]

    media= sum/lung_sottoin

    F_M_medie[i]= media
    masse_medie[i]= masse[lung_sottoin*i + int(j/2)]





fig= plt.figure()

plt.plot(masse, F_M, color="blue", linestyle="-", marker="", label="Segnale Originale")
plt.plot(masse_medie, F_M_medie, color="red", linestyle="-", marker="", label="Media Mobile")


plt.title("Confronto tra il Segnale Originale e la Sua Media Mobile")

plt.xlabel("M [u.a.]")
plt.ylabel("F(M) [u.a.]")

plt.xlim(masse[0], masse[-1])

plt.legend()
plt.tight_layout()

plt.show()