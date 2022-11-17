import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



# valori dei parametri della funzione f_m (ossia f(m) )

mu=5
sigma=1

# path e nome dei file in cui è salvato lo spettro in potenza e quello in cui sono salvati i valori della frequenza e delle due masse
path= "C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\Effetto di Memoria\\Calcolo_diretto_dello_spettro_in_potenza\\Risultati"

file_name_spettro= "DP_in_funzione_f_m_1_m_2"

file_name_f_m_1_m_2= "valori_di_f_m_1_m_2"







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

# valore dell'eccentricità
e=1.005

# H_0 è la costante di Hubble


# rho_c è la densità critica


# v_0 è la velocità asintotica

# costante che compare nell'integrale

cost= 25.4*10**(-8)*h_70**4*(ome_DM/0.25)**2*(delta_loc/10**8)*(e**2-1)/(v_0/c)**3
cost= cost/(rho_c*np.sqrt(ome_m))







# DEFINIZIONE FUNZIONI

# descrive la distribuzione in massa dei buchi neri primordiali
def f_m(m, mu, sigma):

    return (1/np.sqrt(2*math.pi*sigma**2))*np.exp(-(m-mu)**2/(2*sigma**2))







# LETTURA DA FILE DELL'ARRAY CONTENENTE LE MASSE m_1, LE MASSE m_2, LE FREQUENZE f E LO SPETTRO DI POTENZA. QUEST'ULTIMO E' UN ARRAY CON 3 INDICI, IL PRIMO RELATIVO AD f, IL SECONDO AD m_1 E IL TERZO A m_2


file_name_f_m_1_m_2= file_name_f_m_1_m_2 + "____e_" + str(e) + "_a_" + str(a)
file_name_f_m_1_m_2= path + "\\" + file_name_f_m_1_m_2

f, m_1, m_2= np.loadtxt(file_name_f_m_1_m_2, unpack=True)



file_name_spettro= file_name_spettro + "____e_" + str(e) + "_a_" + str(a)
file_name_spettro= path + "\\" + file_name_spettro

df_spettro = pd.read_csv(file_name_spettro, header=None, index_col=[0,1])
spettro= np.array(df_spettro.agg(list, 1).groupby(level=0).agg(list).tolist())







# CALCOLO DELL'INTEGRALE NEL REDSHIFT


# I contiene il risultato di tale integrazione e dipende da tre indici, il primo per la frequenza, il secondo per m_1 e il terzo per m_2

I=np.zeros([ len(f_oss), len(m_1), len(m_2)])

for i in range(0, len(f_oss)):

    mask= f > f_oss[i]

    f_masked= f[mask]
    spettro_masked= spettro[mask,:,:]

    for j in range(len(m_1)):

        for k in range(len(m_2)):

            integranda_f= spettro_masked[:,j,k]/(f_masked)**(2.5)

            I[i,j,k]= f_oss[i]**(3.5)*np.trapz( integranda, f_masked)







# CALCOLO DELL'INTEGRALE NELLE MASSE



# INTEGRALE IN m_2

f_m_1= f_m(m_1, mu, sigma)
f_m_2= f_m(m_2, mu, sigma)

itegrale_m_2= np.zeros([ len(f_oss), len(m_1)])

for i in range(0, len(f_oss)):

    for j in range(0, len(m_1))

    integranda_m_2= f_m_2*( m_2 + m_1[j] )**2*I[i,j,:]/(m_2)**2

    integrale_m_2[i,j] = np.trapz( integranda_m_2, m_2)



# INTEGRALE IN m_1

Omega_GW= np.zeros(len(f_oss))

for i in range(0, len(f_oss)):

    integranda_m_1= f_m_1*integrale_m_2[i,:]/m_1**2

    Omega_GW[i]= np.trapz( integranda_m_1, m_1)







# GRAFICO DI OMEGA_GW

plt.figure()

plt.title("$\Omega_{GW}$ in funzione della frequenza osservata")

plt.plot(f_oss, Omega_GW, linestyle="-", marker="", color="blue")
plt.plot( [ min(f_oss), max(f_oss)], [ 0, 0], linestyle="-", marker="", color="black", linewidth=0.8)

plt.xlabel("$f^{oss}$ [Hz]")
plt.ylabel("$\Omega_{GW}$")

plt.xlim( min(f_oss), max(f_oss) )

plt.yscale("log")

plt.show()


















