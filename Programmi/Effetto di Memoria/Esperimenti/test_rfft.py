import numpy as np
import matplotlib.pyplot as plt


def gauss(x, A, mu, sigma):

    return A*np.exp( -(x - mu)**(2)/(2*sigma**(2)) )







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






A= 1
mu= 0
sigma= 1/5

dt=0.1







# per rendere la fft più efficiente possibile conviene scegliere una potenza di 2 come dimenzione n dell'array, inoltre si deve costruire l'aarray delle x in modo che il primo valore contenuto sia lo zeros, il secondo sia pari a meno l'ultimo, il terzo sia pari a meno il penultimo e così via. Così facendo avanza un vlore, ossia l' (n/2 +1)-esimo, che non è pari al negativo di nessun altro valore e che da quello che ho visto viene solitamente posto negativo

t= np.zeros(2048)

t[0]= 0

t[int(len(t)/2)]= -(len(t)/2)*dt

for i in range(1, int(len(t)/2) ):

    t[i]= dt*i
    t[len(t)-i]= -t[i]






segn= gauss(t, A, mu, sigma)

n= len(segn)






plt.figure()

plt.title("Segnale nel dominio del tempo")

t_graf= array_graf(t)
segn_graf= array_graf(segn)

plt.plot(t_graf, segn_graf, linestyle="-", marker="", color="blue")


plt.plot([t_graf[0], t_graf[-1]], [0, 0], linestyle="-", marker="", color="black", linewidth=0.8)

plt.xlim(t_graf[0], t_graf[-1])

plt.xlabel("t [s]")
plt.ylabel("segnale [U.A.]")





tras_rfft= np.fft.rfft(segn)
print("La lunghezza dell'array trasformato è {0}, mentre quella dell'array di partenza è {1}".format(len(tras_rfft), len(segn)) )


df= 1/(dt*n)
f_rfft_min= 0
f_rfft_max= n*df/2

f_rfft= np.linspace(f_rfft_min, f_rfft_max, int(n/2 + 1))





tras_fft= np.fft.fft(segn)
f_fft= np.fft.fftfreq(n, dt)





plt.figure()

plt.title("Segnale nel dominio delle frequenze")

plt.plot(f_rfft, tras_rfft.real, linestyle="-", marker="", color="red", label="rfft")


f_fft_graf=array_graf(f_fft)

print(tras_fft)
tras_fft_graf=array_graf(tras_fft)

plt.plot(f_fft_graf, tras_fft_graf.real, linestyle="--", marker="", color="yellow", label="fft")


plt.plot([f_fft_graf[0], f_fft_graf[-1]], [0, 0], linestyle="-", marker="", color="black", linewidth=0.8)


plt.xlim(f_fft_graf[0], f_fft_graf[-1])

plt.xlabel("f [Hz]")
plt.ylabel("trasformata del segnale [U.A.]")

plt.legend()


plt.show()







