import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import dblquad



# num contiene il numero di frequenze usate per realizzare il grafico
num=100

# freq_min e fre_max contengono il massimo valore della frequenza e il minimo (l'intervallo delle frequenze è creato in scala logaritmica
freq_min=10**(-8)
freq_max=10**(1)

# percorso del file in cui stampare i risultati
path="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt"
name="Omega_GW.txt"



# massa del sole in kg
M_s=1.989*10**(30)

# unità astronomica in m
UA=1.496*10**(11)

# costante di gravitazione universale in Nm**2/kg**2
G= 6.7*10**(-11)

# costante di gravitazione universale in U.A.**3/(M_S*s**3)
G= (M_s/UA**3)*G



# incertezza sul parametro di Hubble
h_70=0.7

# valore di omega della materia
ome_m=0.3

# valore di omega della materia oscura
ome_DM=0.25

# valore di delta_loc
delta_loc=10**8

# valore del semiasse maggiore in U.A.
a=1

# valore di y definito come e**2 - 1 con e l'eccentricità
y=0.01


# costante che moltiplica la funzione omega nell'equazione 18 dell'articolo 2109.11376. (1/10)**2 compare poiche nell'equazione compare (nu/10)**2 (nu in Hz), (1/1000)**2 compare poiché nell'equazione compare (dm1/1000)*(dm2/1000) (m1 e m2 aono in masse solari)

cost= 9.81*10**(-12)*h_70*(ome_m/0.3)**(-1/2)*(ome_DM/0.25)**2*(delta_loc/10**8)*(a/1)*(y/0.01)*(1/10)**2*(1/1000)**2

'''
print(cost)

nu_0= np.sqrt(a**3/(G*(2)))
print("nu_0=", nu_0)

x_0= 2*math.pi*nu_0*10**(-2)
print("x_0=", x_0)

xi= y - np.arctan(y)
print("xi=", xi)

print(np.exp(2*x_0*xi))

print( "z= ",(1 - y**2 + 4*y**4 + 1.5*(x_0*y**6)/(xi))/( np.exp(2*x_0*xi)*(1 + y**2)**2))
'''

# descrive la distribuzione in massa dei buchi neri primordiali
def f_m(m, mu, sigma):

    return m**2/np.sqrt(2*math.pi*sigma**2)*np.exp(-(m-mu)**2/(2*sigma**2))


# è la funzione che appare nell'integrale di omega escludendo le f(m)
# nu la frequenza e m_1 e m_2 sono le masse dei due BH in masse solari
def I(m_1, m_2, nu):


    nu_0= np.sqrt(a**3/(G*(m_1 + m_2)))
    x_0= 2*math.pi*nu_0*nu

    xi= y - np.arctan(y)

    z= (1 - y**2 + 4*y**4 + 1.5*(x_0*y**6)/(xi))/( np.exp(2*x_0*xi)*(1 + y**2)**2)

    return z


# è la funzione da integrare
def integranda(m_1, m_2, nu, mu, sigma):

    z= f_m(m_1, mu, sigma)*f_m(m_2, mu, sigma)*I(m_1, m_2, nu)

    return z

# è la funzione omega_GW
def omega(nu, mu, sigma, m_min, m_max):

    integrale=dblquad(integranda, m_min, m_max, lambda x: m_min, lambda x: m_max, args=[nu, mu, sigma])[0]

    return cost*nu**2*integrale

nu=np.logspace(np.log10(freq_min), np.log10(freq_max), base=10, num=num)


omega_GW=np.zeros(len(nu))

for i in range(0, len(nu)):

    omega_GW[i]=omega(nu[i], mu=10, sigma=1, m_min=1, m_max=30)


fig, ax = plt.subplots()

ax.plot(nu, omega_GW, linestyle="-", color="blue")

plt.title("$\Omega_{GW}$ in funzione della frequenza")
plt.xlabel("f [Hz]")
plt.ylabel("$\Omega_{GW}$")
plt.yscale("log")
plt.xscale("log")


#ax.ticklabel_format( axis="y", style="sci", scilimits=(0,0))
plt.show()

file_name=path + "\\" + name

file= open(file_name, "w")
np.savetxt(file, np.c_[nu, omega_GW])

file.close()





