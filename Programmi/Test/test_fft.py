import numpy as np
import matplotlib.pyplot as plt
import math


# definizione segnale
t= np.linspace(-10, 10, 1000)
#sign= np.sin(2*math.pi*1*t)
sign= np.exp(-(t)**2)

# calcolo trasformata
tras= np.fft.fft(sign)

n= sign.size
dt= t[1]-t[0]
freq= np.fft.fftfreq(n, dt)

# calcolo della derivata inversa
inv= np.fft.ifft(tras)

n= tras.size
df= freq[1] - freq[0]
time= np.fft.fftfreq(n, df)
time= np.fft.fftshift(time)


# se scrivo sin(omega*t), il periodo èverifica: omega*T=2*pi, quindi T=2*pi/omega e la frequenza nu che è 1/T vale omega/2*pi
# il comando np.fft.fftfreq restituisce nu


plt.subplot(2,2,1)

plt.title("Segnale nel Dominio del Tempo", fontsize=13)


plt.plot(t, sign, linestyle="-", color="blue")
plt.plot(t, 0*t, linestyle="-", color="black", linewidth=0.8)

plt.xlim(min(t), max(t))

plt.xlabel("t [s]")
plt.ylabel("Segnale")




plt.subplot(2,2,2)


plt.title("Segnale nel Dominio del Frequenze", fontsize=13)


plt.plot(freq, abs(tras), linestyle="-", color="blue")
plt.plot(freq, 0*freq, linestyle="-", color="black", linewidth=0.8)

plt.xlim(min(freq), max(freq))

plt.xlabel("f [Hz]")
plt.ylabel("Segnale")




plt.subplot(2,2,3)

plt.title("Antitrasformata", fontsize=13)


plt.plot(time, inv, marker="", linestyle="-", color="blue")
plt.plot(time, 0*time, linestyle="-", color="black", linewidth=0.8)

plt.xlim(min(time), max(time))

plt.xlabel("t [s]")
plt.ylabel("Segnale")




plt.subplot(2,2,4)


plt.title("Segnale di Partenza e Antitrasformata", fontsize=13)

plt.plot(t, sign, linestyle="-", color="blue", label="segnale originale")
plt.plot(time, inv, linestyle="--", color="orange", label="antitrasformata")


infer= max( min(t), min(time) )
super= min( max(t), max(time) )


plt.plot([infer, super], [0, 0], linestyle="-", color="black", linewidth=0.8)



plt.xlim( infer, super)

plt.xlabel("t [Hz]")
plt.ylabel("Segnali")

plt.legend()





plt.tight_layout()

plt.show()