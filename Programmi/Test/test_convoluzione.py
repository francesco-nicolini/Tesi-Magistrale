import numpy as np

K=500

wind_x=K/10
wind_y= 0.001

# definizione della funzione di cui si calcola il prodotto di convoluizone per se stessa

def f_m(m, mu, sigma):

    return (m**2/np.sqrt(2*math.pi*sigma**2))*np.exp(-(m-mu)**2/(2*sigma**2))

asse= np.linspace(0, 30, K)

mu= 10
sigma= 1

a= f_m(asse, mu, sigma)
b= f_m(asse, mu, sigma)

'''
a= a/a.sum()
b= b/b.sum()
'''

'''
asse_a= np.linspace(0, len(a), len(a), endpoint= False)
asse_b= np.linspace(0, len(b), len(b), endpoint= False)
'''

asse_a= asse
asse_b= asse

dM= asse[1] - asse[0]
conv= dM*np.convolve(a, b, mode="full")

asse_conv= np.linspace(0, len(conv), len(conv), endpoint=False)*dM


# Identificazione massimo della convoluzione

max_conv= max(conv)

index= np.argmax(conv)
M_conv_max= asse_conv[index]
print(M_conv_max)
print(max(conv))

print("La somma dei valori del segnale convoluto è: {0}".format(conv.sum()))

plt.plot(asse_a, a, linestyle="-", color="blue", label="Primo segnale")
plt.plot(asse_b, b, linestyle="-", color="orange", label="Secondo segnale")
plt.plot(asse_conv[::1], conv[::1], linestyle="-", color="black", label="Convoluzione")


minimo_x= min(np.amin(asse_a), np.amin(asse_b))
minimo_x= min(minimo_x, np.amin(asse_conv))

massimo_x= max(np.amax(asse_a), np.amax(asse_b))
massimo_x= max(massimo_x, np.amax(asse_conv))

minimo_y= min(np.amin(a), np.amin(b))
minimo_y= min(minimo_y, np.amin(conv))

massimo_y= max(np.amax(a), np.amax(b))
massimo_y= max(massimo_y, np.amax(conv))


asse_x= np.linspace(minimo_x - wind_x, massimo_x + wind_x, 100)

plt.plot(asse_x, np.zeros(len(asse_x)), linestyle="-", color="black")

plt.xlim(minimo_x - wind_x, massimo_x + wind_x)
plt.ylim(minimo_y - wind_y, massimo_y + wind_y)

plt.title("CONVOLUZIONE DI DUE SEGNALI")
plt.xlabel("N")
plt.ylabel("Segnali")

ax= plt.gca()




plt.legend()
plt.tight_layout()

plt.show()


print( "Lunghezza prima funzione= {0}, Lunghezza seconda funzione= {1}, Lunghezza convoluzione= {2}".format(len(a), len(b), len(conv[::2])) )



