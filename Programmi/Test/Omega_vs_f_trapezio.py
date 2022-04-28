import numpy as np
import matplotlib.pyplot as plt
import math



# OPZIONI VARIE

# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)
K=1000

# Estremo superiore e estremo inferiore delle masse considerate nell'integrale e scarto tra due masse consecutive
m_min=1
m_max=30

dm= (m_max - m_min)/(K-1)

# Se option è pari a "read" omega_GW di Omega_vs_f.py viene letto da file e disegnato
option="read"

# Path del file se option è uguale a "read"
file_name_omega="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt\\omega_GW.txt"


# Estremo superiore e estremo inferiore delle frequenze considerate (da selezionare solo se option è diverso da "read")
freq_min=10**(-8)
freq_max=10**(1)

# Parametri della distribuzione logaritmica di massa

mu= 10
sigma= 1






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

    freq_1, omega_GW_1= np.loadtxt(file_name_omega, unpack=True)






freq= np.logspace(np.log10(freq_min), np.log10(freq_max), K)










# DEFINIZIONE DELLE FUNZIONI CHE COMPAIONO NELLA RELAZIONE DI RIFERIMENTO UNA VOLTA CHE E' STATO DISCRETIZZATO L'INTEGRALE


# Lista delle masse indagate nella discretizzazione dell'integrale (sono in masse solari)

masse=np.linspace(m_min, m_max, K)


# descrive la distribuzione in massa dei buchi neri primordiali
def f_m(m, mu, sigma):

    return (m**2/np.sqrt(2*math.pi*sigma**2))*np.exp(-(m-mu)**2/(2*sigma**2))



# La funzione usata per costruire le matrici (si tiene conto anche del dm**2 in essa)

def integ(m_1, m_2, nu):

    nu_0= np.sqrt(a**3/(G*(m_1 + m_2)))
    x_0= 2*math.pi*nu_0*nu

    z= dm**2*(1 - y**2 + 4*y**4 + 1.5*(x_0*y**6)/(xi))/( np.exp(2*x_0*xi)*(1 + y**2)**2)

    return cost*(nu**2)*z



# La funzione che crea le matrici (si è usata la regola del trapezio, vedi file Prova_costruz_matrice_per_regola_trapez.py)

def fun_mat(nu, K):

    a=np.zeros((K,K))

    for i in range(0,K):

        for j in range(0,K):

            if ( (i+j==0) or (i+j==2*(K-1)) or (i*j==0 and i+j==K-1) ):

                a[i][j]=integ(masse[i], masse[j], nu)/4


            # è un elif quindi si possono mettere casi già previsti dall'if
            #precedente dal momento che questi, verificando quest'ultimo, non saranno
            #analizzati dal secondo
            elif ( i*j==0 or i==K-1 or j==K-1 ):

                a[i][j]=integ(masse[i], masse[j], nu)/2

            else:

                a[i][j]=integ(masse[i], masse[j], nu)

    return a



# La lista matrix è riempita con delle matrici delle forme quadrariche

matrix=[]

for i in range(0, K):


    matrix.append(fun_mat(freq[i], K))



# Tramite le matrici in matrix sono create le funzioni che vanno a comporre il sistema e che sono inserite nella lista funz

funz=[]

for i in range(0,K):

    def f(*var,i=i):

        vec_dx=np.dot(var,matrix[i])
        return np.dot(vec_dx,var)-omega_GW[i]/cost

    funz.append(f)




var= f_m(masse, mu, sigma)



omega_GW= np.zeros(K)

for i in range(0, K):

    omega_GW[i]= funz[i](*var)













# Grafico di omega_GW

fig, ax = plt.subplots()

ax.plot(freq, omega_GW, linestyle="-", color= "blue", label= "Regola Trapezio")

if( option=="read" ):
    ax.plot(freq_1, omega_GW_1, linestyle="--", color="yellow", label= "Funzione Python")

plt.title("$\\Omega_{GW}$ in funzione della frequenza")
plt.xlabel("f [Hz]")
plt.ylabel("$\\Omega_{GW}$")
plt.yscale("log")
plt.xscale("log")
plt.legend()

plt.tight_layout()
plt.show()