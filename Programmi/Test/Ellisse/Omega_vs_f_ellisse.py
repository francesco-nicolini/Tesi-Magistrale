import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import dblquad



# num contiene il numero di frequenze usate per realizzare il grafico
num=500

# freq_min e fre_max contengono il massimo valore della frequenza e il minimo (l'intervallo delle frequenze è creato in scala logaritmica
freq_min=10**(-8)
freq_max=10**(1)

# percorso del file in cui stampare i risultati
path="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt\\Ellisse"
name_omega= "omega_GW_ellisse.txt"
name_f_m= "f_m_ellisse.txt"

# valore minimo e valore massimo considerati per la massa

m_min=1
m_max=30

# valori dei parametri della funzione f_m (ossia f(m) )

mu=10
sigma=1





# massa del sole in kg
M_s=1.989*10**(30)

# unità astronomica in m
UA=1.496*10**(11)

# costante di gravitazione universale in Nm**2/kg**2
G= 6.7*10**(-11)

# costante di gravitazione universale in U.A.**3/(M_S*s**3)
G= (M_s/UA**3)*G



# VALORE DI ALCUNE GRANDEZZE FISICHE CHE DEFINISCONO IL SISTEMA DEI DUE BUCHI NERI

# Incertezza sul parametro di Hubble
h_70= 0.7

# Valore di omega della materia oscura
ome_DM= 0.25

# Valore di delta_loc
delta_loc= 10**8

# Velocità iniziale tra i due buchi neri in Km/s
v_0= 10

# Costante che moltiplica l'integrale nell'equazione 18 dell'articolo 2109.11376. (1/10)**2 compare poiche nell'equazione compare (nu/10)**2 (nu in Hz), (1/1000)**2 compare poiché nell'equazione compare (dm1/1000)*(dm2/1000) (m1 e m2 aono in masse solari)

cost= 2.39*10**(-13)*h_70*(ome_DM/0.25)**2*(delta_loc/10**8)*(v_0/10)**(-11/7)





# descrive la distribuzione in massa dei buchi neri primordiali
def f_m(m, mu, sigma):

    return (m**2/np.sqrt(2*math.pi*sigma**2))*np.exp(-(m-mu)**2/(2*sigma**2))


# è la funzione che appare nell'integrale di omega escludendo le f(m)
# nu la frequenza e m_1 e m_2 sono le masse dei due BH in masse solari
def I(m_1, m_2):


    z= ( (m_1 + m_2)**(23/21) )/(m_1*m_2)**(5/7)

    return z


# è la funzione da integrare
def integranda(m_1, m_2, mu, sigma):

    z= f_m(m_1, mu, sigma)*f_m(m_2, mu, sigma)*I(m_1, m_2)

    return z

# è la funzione omega_GW
def omega( mu, sigma, m_min, m_max):

    integrale= dblquad(integranda, m_min, m_max, lambda x: m_min, lambda x: m_max, args=[ mu, sigma])[0]

    return cost*integrale

nu=np.logspace(np.log10(freq_min), np.log10(freq_max), base=10, num=num)


omega_GW=np.zeros(len(nu))

integrale= omega( mu=mu, sigma=sigma, m_min=m_min, m_max=m_max)

for i in range(0, len(nu)):

    omega_GW[i]= integrale*nu[i]



# Grafico

fig, ax = plt.subplots()

ax.plot(nu, omega_GW, linestyle="-", color="blue")

plt.title("$\\Omega_{GW}$ in funzione della frequenza")
plt.xlabel("f [Hz]")
plt.ylabel("$\\Omega_{GW}$")
plt.yscale("log")
plt.xscale("log")


#ax.ticklabel_format( axis="y", style="sci", scilimits=(0,0))




# stampa su file

masse=np.linspace(m_min, m_max, 1000)
distrib= f_m(masse, mu, sigma)

file_name_omega= path + "\\" + name_omega
file_name_f_m= path + "\\" + name_f_m

file_omega= open(file_name_omega, "w")
file_f_m= open(file_name_f_m, "w")

np.savetxt(file_omega, np.c_[nu, omega_GW])
file_omega.close()

np.savetxt(file_f_m, np.c_[masse, distrib])
file_f_m.close()





plt.show()
