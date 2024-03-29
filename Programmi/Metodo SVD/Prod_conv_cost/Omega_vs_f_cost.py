import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad



# num contiene il numero di frequenze usate per realizzare il grafico
num=500

# freq_min e fre_max contengono il massimo valore della frequenza e il minimo (l'intervallo delle frequenze è creato in scala logaritmica
freq_min=10**(-8)
freq_max=10**(1)

# percorso del file in cui stampare i risultati
path="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\Metodo SVD\\Prod_conv_cost\\file_txt"
name_omega= "omega_GW_" + str(num) + "_prod_conv_cost.txt"
name_f_m= "f_m_" + str(num) + "_prod_conv_cost.txt"

# valore minimo e valore massimo considerati per la massa

m_min=1
m_max=20

# valori dei parametri della funzione F_M (ossia F(M))

costante=3





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




# descrive il prodotto di convoluzione della distribuzione in massa dei buchi neri primordiali per se stessa
def F_M(M, costante):

    if (hasattr(M, "__len__")==False):

        return costante

    else:

        array= np.ones(len(M))
        return costante*array

# è la funzione che appare nell'integrale di omega escludendo F(M)
# nu la frequenza e M è la massa della somma dei due BH in masse solari
def I(M, nu):


    nu_0= np.sqrt(a**3/(G*(M)))
    x_0= 2*math.pi*nu_0*nu

    xi= y - np.arctan(y)

    z= (1 - y**2 + 4*y**4 + 1.5*(x_0*y**6)/(xi))/( np.exp(2*x_0*xi)*(1 + y**2)**2)

    return z


# è la funzione da integrare
def integranda(M, nu, costante):

    z= F_M(M, costante)*I(M, nu)

    return z

# è la funzione omega_GW
def omega(nu, costante, m_min, m_max):

    integrale=quad(integranda, m_min, m_max, args=(nu, costante))[0]

    return cost*nu**2*integrale

nu=np.logspace(np.log10(freq_min), np.log10(freq_max), base=10, num=num)


omega_GW= np.zeros(len(nu))

for i in range(0, len(nu)):

    omega_GW[i]=omega(nu[i], costante, m_min, m_max)







# Grafico di Omega

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
distrib= F_M(masse, costante)

print(len(masse), len(distrib))
file_name_omega= path + "\\" + name_omega
file_name_f_m= path + "\\" + name_f_m

file_omega= open(file_name_omega, "w")
file_f_m= open(file_name_f_m, "w")

np.savetxt(file_omega, np.c_[nu, omega_GW])
file_omega.close()

np.savetxt(file_f_m, np.c_[masse, distrib])
file_f_m.close()



# Grafico di F_M

fig, ax = plt.subplots()


ax.plot(masse, F_M(masse, costante), linestyle="-", color="blue")

plt.title("$\\Omega_{GW}$ in funzione della frequenza")
plt.xlabel("f [Hz]")
plt.ylabel("$\\Omega_{GW}$")

bordi= 10

plt.xlim( m_min-bordi, m_max+bordi)
plt.ylim( -1, costante+2)

asse= np.linspace(m_min-bordi, m_max+bordi, 1000)
ax.plot( asse, 0*asse, linestyle="-", color="black")
ax.plot( [m_min, m_min], [0, costante], color="blue", linestyle="--")
ax.plot( [m_max, m_max], [0, costante], color="blue", linestyle="--")



plt.show()

