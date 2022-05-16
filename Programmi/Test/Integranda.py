import numpy as np
import matplotlib.pyplot as plt
import math
import keyboard



# OPZIONI VARIE

# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)
K=50


# Estremo superiore e estremo inferiore delle masse considerate nello studio della fuzione integranda

m_min=1
m_max=30

dm= (m_max - m_min)/(K-1)

# Se option è pari a "read" omega_GW è letto da file altrimenti viene creato tramite una funzione
option="read"

# Path del file se option è uguale a "read"
file_name_omega="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt\\omega_GW.txt"


# Estremo superiore e estremo inferiore delle frequenze considerate (da selezionare solo se option è diverso da "read")
freq_min=10**(-8)
freq_max=10**(1)


# camp è il numero di frequenze (+1) tra un campionamento e l'altro quando si realizza il grafico
camp= 2


# se la variabile omega è pari a True viene realizzato il grafico della funzione integranda riscalata per Omega_GW (letto da file, quindi applicabile solo se option="read"), se invece assume un qualsiasi altro valore viene realizzato il grafico dell'integranda non riscalata
omega= 1










# COSTANTI FISICHE

# Massa del sole in kg
M_s=1.989*10**(30)

# Unità astronomica in m
UA=1.496*10**(11)

# Costante di gravitazione universale in Nm**2/kg**2
G= 6.7*10**(-11)

# Costante di gravitazione universale in U.A.**3/(M_S*s**2)
G= (M_s/UA**3)*G







# VALORE DI ALCUNE GRANDEZZE FISICHE CHE DEFINISCONO IL SISTEMA DEI DUE BUCHI NERI

# Incertezza sul parametro di Hubble
h_70=0.7

# Valore di omega della materia
ome_M=0.3

# Valore di omega della materia oscura
ome_DM=0.25

# Valore di delta_loc
delta_loc=10**8

# Valore del semiasse maggiore in U.A.
a=1

# Valore di y definito come e**2 - 1 con e l'eccentricità
y=0.01

xi= y - np.arctan(y)

# Costante che moltiplica l'integrale nell'equazione 18 dell'articolo 2109.11376. (1/10)**2 compare poiche nell'equazione compare (nu/10)**2 (nu in Hz), (1/1000)**2 compare poiché nell'equazione compare (dm1/1000)*(dm2/1000) (m1 e m2 aono in masse solari)

cost= 9.81*10**(-12)*h_70*(ome_M/0.3)**(-1/2)*(ome_DM/0.25)**2*(delta_loc/10**8)*(a/1)*(y/0.01)*(1/10)**2*(1/1000)**2






# CREAZIONE DELL'ARRAY CONTENENTE LE FREQUENZE E DELL'ARRAY CONTENENTE I VALORI DI OMEGA_GW

# freq contiene la lista delle frequenze
# omega_GW contiene la lista dei valori dello spettro in potenza del fondo


if (option=="read"):

    freq, omega_GW= np.loadtxt(file_name_omega, unpack=True)

    print("Il vettore contenente omega_GW ha una dimensione pari a {0}\n".format(len(omega_GW)))


else:

    def funz_omeg(nu):

        n=len(nu)

        return 10**(-22)*np.ones(n)


    freq=np.logspace(np.log10(freq_min), np.log10(freq_max), K)
    omega_GW=funz_omeg(freq)

'''
# Grafico di omega_GW

fig, ax = plt.subplots()

ax.plot(freq, omega_GW, linestyle="-", color="blue")

plt.title("$\\Omega_{GW}$ in funzione della frequenza")
plt.xlabel("f [Hz]")
plt.ylabel("$\\Omega_{GW}$")
plt.yscale("log")
plt.xscale("log")
'''







# La funzione usata per costruire la matrice, M è la massa totale del sistema, mentre nu è una frequenza

def integ(M, nu):

    nu_0= np.sqrt(a**3/(G*(M)))
    x_0= 2*math.pi*nu_0*nu

    z= (1 - y**2 + 4*y**4 + 1.5*(x_0*y**6)/(xi))/( np.exp(2*x_0*xi)*(1 + y**2)**2)

    return cost*(nu**2)*z





# Grafico della funzione integranda al variare della massa per diverse frequenze

# l'if che segue serve a far si che venfa sempre realizzato il grafico corrispondente all'ultima frequenza contenuta nel vettore freq

resto= len(freq)%camp

if( resto==0 ):
    freq_graf= freq[camp-1::camp]

    if( omega==True ):
        omega_graf= omega_GW[camp-1::camp]

else:
    freq_graf= freq[resto-1::camp]

    if( omega==True ):
        omega_graf= omega_GW[resto-1::camp]

'''
print("\n\nLa massima frequenza studiata è pari a {:.2e} Hz".format(max(freq_graf)))
'''
masse= np.linspace(1, 100000, 100000)

fig= plt.figure()

plt.xlabel("M [M_sole]", fontsize=12)
plt.xlim(min(masse)-1, max(masse))

for i in range(0, len(freq_graf)):

    nu= freq_graf[i]

    integranda= integ(masse, nu)

    color= "C{0}".format(i)

    if( 1<nu<10 ):
        label= "{:.2} Hz".format(nu)

    else:
        label= "{:.2e} Hz".format(nu)

    if( omega==True ):

        Omega= omega_graf[i]
        print("frequenza={:.2e}  omega= {:.2e}  ultimo valore integranda= {:.2e}  rapporto= {:.2e}".format(nu, Omega, integranda[-1], integranda[-1]/Omega))

        plt.title("Funzione Integranda Riscalata per $\\Omega_{GW}$ al Variare della Massa e per Diverse Frequenze", fontsize=15)

        plt.ylabel("Funzione Integranda Riscalata", fontsize=12)

        plt.plot(masse, integranda/Omega, color=color, linestyle="-", marker="", label=label)


    else:

        print("frequenza={:.2e}  ultimo valore integranda= {:.2e}".format(nu, integranda[-1]))

        plt.title("Funzione Integranda al Variare della Massa per Diverse Frequenze", fontsize=15)

        plt.ylabel("Funzione Integranda", fontsize=12)

        plt.plot(masse, integranda, color=color, linestyle="-", marker="", label=label)


plt.legend( title="frequenze", loc="upper left", bbox_to_anchor=(1, 1.0125))
plt.show()
























