import numpy as np
import matplotlib.pyplot as plt

# Path del file contenente la soluzione vera
file_name_f_m="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt\\f_m.txt"

K=50

m_min=0
m_max=30


masse=np.linspace(m_min, m_max, K)

iniz=np.zeros(K)

#da 13 a 17 sale e poi da 18 a 22 scende

for i in range(13, 18):
    iniz[i]= 40/5 + iniz[i-1]

for i in range(18, 22):
    iniz[i]= iniz[i-1] - 40/4



masse_graf, f_esatta= np.loadtxt(file_name_f_m, unpack=True)


plt.title("Prova Condizioni Iniziali")
plt.xlabel("massa [M_sole]")
plt.ylabel("f(m)")

plt.plot(masse_graf, f_esatta, color="red", linestyle="-", marker="")

plt.plot(masse, iniz, color="blue", linestyle="", marker=".")


plt.show()