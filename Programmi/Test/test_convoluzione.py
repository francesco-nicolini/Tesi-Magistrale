import numpy as np

K=1000

wind_x=K/10
wind_y= 0.001

a= np.linspace(0, 1, K)
b= np.linspace(0, 1, K)**4

a= a/a.sum()
b= b/b.sum()

asse_a= np.linspace(0, len(a), len(a), endpoint= False)
asse_b= np.linspace(0, len(b), len(b), endpoint= False)

conv=np.convolve(a, b, mode="full")
asse_conv= np.linspace(0, len(conv), len(conv), endpoint=False)

print("La somma dei valori del segnale convoluto è: {0}".format(conv.sum()))

plt.plot(asse_a, a, linestyle="-", color="blue", label="Primo segnale")
plt.plot(asse_b, b, linestyle="-", color="orange", label="Secondo segnale")
plt.plot(asse_conv, conv, linestyle="-", color="red", label="Convoluzione")


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


plt.legend()
plt.tight_layout()

plt.show()
