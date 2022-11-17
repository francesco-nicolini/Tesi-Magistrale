import numpy as np
import matplotlib.pyplot as plt







# COSTANTI DEL PROGRAMMA

N= 2048
df= 0.01
A= 1
mu= 0
sigma= 1/5

n= int(N/2 + 1)
dt= 1/(N*df)







# DEFINIZIONE FUNZIONI

# Gaussiana
def gauss(x, A, mu, sigma):

    return A*np.exp( -(x - mu)**(2)/(2*sigma**(2)) )




# Questa funzione inverte l'ordine di alcuni elementi dell'array in ingresso, nello specifico pone per primi quegli elementi che si trovano oltre la metà. Ciò viene fatto al fine di rendere più bella la rappresentazione grafica degli array ottenuti mediante l'utilizzo di fft, evitando che vengano disegnati dei punti viccini ma non congiunti e punti distanti ma congiunti
def array_graf(a):

    n= len(a)

    if( np.iscomplexobj(a) ):
        b= np.zeros(n, dtype=complex)

    else:
        b= np.zeros(n)

    for i in range(0, int(n/2)):
        b[i] = a[int(n/2) + i]

    for i in range(0, int(n/2)):
        b[int(n/2) + i]= a[i]

    return b







# CREAZIONE DELL'ARRAY DELLE FREQUENZE

f= np.linspace(0, (n-1)*df, n)







# CREAZIONE DELL'ARRAY DEI TEMPI

t= np.zeros(N)

t[0]= 0

t[int(N/2)]= -(N/2)*dt

for i in range(1, int(N/2) ):

    t[i]= dt*i
    t[N-i]= -t[i]







# CREAZIONE DEL SEGALE NEL DOMINIO DEL TEMPO

segn= gauss(t, A, mu, sigma)







# CALCOLO DEL SEGNALE NEL DOMINIO DELLE FREQUENZE

segn_tras= np.fft.rfft(segn)







# GRAFICO DEL SEGNALE NEL DOMINIO DEL TEMPO E NEL DOMINIO DELLA FREQUENZA


fig, ax = plt.subplots(2,1)

plt.suptitle("TRASFORMATA DEL SEGNALE")



ax[0].set_title("Segnale nel dominio del tempo")

t_graf= array_graf(t)
segn_graf= array_graf(segn)

ax[0].plot(t_graf, segn_graf, marker="", linestyle="-", color="blue")
ax[0].plot([t_graf[0], t_graf[-1]], [0, 0], marker="", linestyle="-", color="black", linewidth=0.8)

ax[0].set_xlim(t_graf[0], t_graf[-1])

ax[0].set_xlabel("t [s]")
ax[0].set_ylabel("Segnale [U.A.]")









ax[1].set_title("Segnale nel dominio delle frequenze")

ax[1].plot(f, segn_tras.real, marker="", linestyle="-", color="red")
ax[1].plot([f[0], f[-1]], [0, 0], marker="", linestyle="-", color="black", linewidth=0.8)

ax[1].set_xlim(f[0], f[-1])

ax[1].set_xlabel("f [Hz]")
ax[1].set_ylabel("Segnale [U.A.]")



plt.tight_layout()
plt.show()





















