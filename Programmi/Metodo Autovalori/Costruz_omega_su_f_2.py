import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import dblquad



# num contiene il numero di frequenze usate per realizzare il grafico
num=50

# freq_min e fre_max contengono il minimo e il massimo valore dell'intervallo della frequenza usato per visualizzare tutta la funzione omega (l'intervallo delle frequenze è creato in scala logaritmica
freq_tutta_min=10**(-8)
freq_tutta_max=10**(1)


# alpha_min e alpha_max contengono il valore minimo e il valore massimo di alpha (definito come l'inverso del quadrato della frequenza)
alpha_min= 10**(0)
alpha_max= 10**(2)






# percorso del file in cui stampare i risultati
path="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\Metodo Autovalori\\file_txt"
name_omega_tutta_minore= "omega_GW_tutta_minore.txt"
name_omega_tutta_maggiore= "omega_GW_tutta_maggiore.txt"
name_omega= "omega_GW.txt"
name_f_m= "f_m.txt"

# valore minimo e valore massimo considerati per la massa
m_min= alpha_min
m_max= alpha_max

# valori dei parametri della funzione f_m ( ossia f(m) )

mu=10
sigma=1



# GRANDEZZE FISICHE

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






# DEFINIZIONE DELLE FUNZIONI


# descrive la distribuzione in massa dei buchi neri primordiali
def f_m(m, mu, sigma):

    return (m**2/np.sqrt(2*math.pi*sigma**2))*np.exp(-(m-mu)**2/(2*sigma**2))


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





# CREAZIONE DEI DUE ARREY DI OEMGA

# creazione estremi delle frequenze considerate nella finestra

# freq_min e fre_max contengono il minimo e il massimo valore dell'intervallo della frequenza usato per la costruzione della funzione omega nell'intervallo di interesse (l'intervallo delle frequenze è creato in scala logaritmica

freq_min= 1/(np.sqrt(alpha_max))
freq_max= 1/(np.sqrt(alpha_min))

alpha= np.linspace(alpha_min, alpha_max, endpoint=False, num=num)

nu= 1/np.sqrt(alpha)

nu= nu[::-1]

print("nu= ",nu)
print("\nalpha= ", alpha)


riscalamento= 3

nu_tutta_minore= np.logspace(np.log10(freq_tutta_min), np.log10(nu[0]), base=10, num= int(num/riscalamento))

omega_GW_tutta_minore= np.zeros(len(nu_tutta_minore))

for i in range(0, len(nu_tutta_minore)):

    omega_GW_tutta_minore[i]= omega(nu_tutta_minore[i], mu=mu, sigma=sigma, m_min=m_min, m_max=m_max)



nu_tutta_maggiore= np.logspace(np.log10(nu[-1]), np.log10(freq_tutta_max), base=10, num= int(num/riscalamento))

omega_GW_tutta_maggiore= np.zeros(len(nu_tutta_maggiore))

for i in range(0, len(nu_tutta_maggiore)):

    omega_GW_tutta_maggiore[i]= omega(nu_tutta_maggiore[i], mu=mu, sigma=sigma, m_min=m_min, m_max=m_max)





omega_GW= np.zeros(len(nu))


print("La frequenza varia all'interno di tale range: [ {:.2} - {:.2} ]".format(nu[0], nu[-1]))


for i in range(0, len(nu)):

    omega_GW[i]= omega(nu[i], mu=mu, sigma=sigma, m_min=m_min, m_max=m_max)

'''
# CREAZIONE DEI DUE ARREY DI OEMGA

nu_tutta= np.logspace(np.log10(freq_tutta_min), np.log10(freq_tutta_max), base=10, num= int(num/10))

omega_GW_tutta= np.zeros(len(nu_tutta))

for i in range(0, len(nu_tutta)):

    omega_GW_tutta[i]= omega(nu_tutta[i], mu=mu, sigma=sigma, m_min=m_min, m_max=m_max)




# MASCHERE PER I GRAFICI


mask_minore= ( nu_tutta < freq_min )

nu_tutta_minore= nu_tutta[mask_minore]
omega_GW_tutta_minore= omega_GW_tutta[mask_minore]

nu_tutta_minore= np.concatenate([nu_tutta_minore, [freq_min]])

val_sup= omega(freq_min, mu=mu, sigma=sigma, m_min=m_min, m_max=m_max)
omega_GW_tutta_minore= np.concatenate([omega_GW_tutta_minore, [val_sup]])



mask_maggiore= ( nu_tutta > freq_max )

nu_tutta_maggiore= nu_tutta[mask_maggiore]
omega_GW_tutta_maggiore= omega_GW_tutta[mask_maggiore]

nu_tutta_maggiore= np.concatenate([[freq_max], nu_tutta_maggiore])

val_inf= omega(freq_max, mu=mu, sigma=sigma, m_min=m_min, m_max=m_max)
omega_GW_tutta_maggiore= np.concatenate([[val_inf], omega_GW_tutta_maggiore])
'''

# GRAFICO

fig, ax = plt.subplots()

color_tut= "cornflowerblue"
color_wind= "midnightblue"


ax.plot(nu_tutta_minore, omega_GW_tutta_minore, linestyle="-", color=color_tut)
ax.plot(nu_tutta_maggiore, omega_GW_tutta_maggiore, linestyle="-", color=color_tut, label="Tutto il Range")
ax.plot(nu, omega_GW, linestyle="-", color=color_wind, label="Finestra di Interesse")



inf_omega= min(omega_GW_tutta_minore)

ax.plot([nu[0], nu[0]], [inf_omega, omega_GW[0]], linestyle="--", marker="", color=color_wind, alpha=1)
plt.text( nu[0]/1.2, inf_omega*10, "f= {:.2} Hz".format(freq_min), ha="right", va="top", color=color_wind, fontsize=9)

ax.plot([nu[-1], nu[-1]], [inf_omega, omega_GW[-1]], linestyle="--", marker="", color=color_wind, alpha=1)
plt.text( nu[-1]*1.2, inf_omega*10, "f= {:.2} Hz".format(freq_max), ha="left", va="top", color=color_wind, fontsize=9)




plt.title("$\\Omega_{GW}$ in funzione della frequenza")
plt.xlabel("f [Hz]")
plt.ylabel("$\\Omega_{GW}$")
plt.yscale("log")
plt.xscale("log")
plt.xlim(freq_tutta_min, freq_tutta_max)
plt.ylim(inf_omega, 10**(-15))

plt.legend()




# stampa su file

masse=np.linspace(m_min, m_max, 1000)
distrib= f_m(masse, mu, sigma)

file_name_omega_tutta_minore= path + "\\" + name_omega_tutta_minore
file_name_omega_tutta_maggiore= path + "\\" + name_omega_tutta_maggiore
file_name_omega= path + "\\" + name_omega
file_name_f_m= path + "\\" + name_f_m

file_omega_tutta_minore= open(file_name_omega_tutta_minore, "w")
file_omega_tutta_maggiore= open(file_name_omega_tutta_maggiore, "w")
file_omega= open(file_name_omega, "w")
file_f_m= open(file_name_f_m, "w")


np.savetxt(file_omega_tutta_minore, np.c_[nu_tutta_minore, omega_GW_tutta_minore])
file_omega_tutta_minore.close()

np.savetxt(file_omega_tutta_maggiore, np.c_[nu_tutta_maggiore, omega_GW_tutta_maggiore])
file_omega_tutta_maggiore.close()

np.savetxt(file_omega, np.c_[nu, omega_GW])
file_omega.close()

np.savetxt(file_f_m, np.c_[masse, distrib])
file_f_m.close()


plt.show()