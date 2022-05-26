import numpy as np
from numpy.linalg import eig
import matplotlib.pyplot as plt
import math




# OPZIONI VARIE


# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)
K=5

# Estremo superiore e estremo inferiore delle masse considerate nell'integrale e scarto tra due masse consecutive
m_min=1
m_max=30

dm= (m_max - m_min)/(K-1)

# Se option è pari a "read" omega_GW è letto da file altrimenti viene creato tramite una funzione
option="read"

# Path del file se option è uguale a "read"
file_name_omega="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt\\omega_GW.txt"

# Se option è pari a "read" ponendo disegna uguale a True verra realizzato il grafico della soluzione esatta nello stesso piano in cui vengono rappresentate le soluzioni trovate minimizzando
disegna=True

# Path del file se disegna è uguale a True
file_name_f_m="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt\\f_m.txt"

# Estremo superiore e estremo inferiore delle frequenze considerate (da selezionare solo se option è diverso da "read")
freq_min=10**(-8)
freq_max=10**(1)






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

    if (disegna==True):

        #  masse_graf contiene le masse per effettuare il grafico di f_esatta
        # f_esatta la forma della distribuzione letta da file (è la soluzione esatta)

        masse_graf, f_esatta= np.loadtxt(file_name_f_m, unpack=True)

else:

    def funz_omeg(nu):

        n=len(nu)

        return 10**(-22)*np.ones(n)


    freq=np.logspace(np.log10(freq_min), np.log10(freq_max), K)
    omega_GW=funz_omeg(freq)


# Campionamento del vettore omega_GW e del vettore freq in modo che abbia una dimensione pari a K


lung= len(omega_GW)


rap= int( lung/K )
res= lung%K


if(res==0):

    omega_GW= omega_GW[rap-1::rap]
    freq= freq[rap-1::rap]

else:
    omega_GW= omega_GW[res-1+rap::rap]
    freq= freq[res-1+rap::rap]


print("lunghezza Omega_GW= {0}, numero variabili= {1}".format(len(omega_GW), K))







# DEFINIZIONE DELLE FUNZIONI CHE COMPAIONO NELLA RELAZIONE DI RIFERIMENTO UNA VOLTA CHE E' STATO DISCRETIZZATO L'INTEGRALE


# Lista dei valori della massa totale utilizzati per effettuare la discretizzazione dell'integrale (sono in masse solari)

masse= np.linspace(m_min, m_max, K)

dM= masse[1] - masse[0]



# La funzione usata per costruire la matrice, M è la massa totale del sistema, mentre nu è una frequenza

def integ(M, nu):

    nu_0= np.sqrt(a**3/(G*(M)))
    x_0= 2*math.pi*nu_0*nu

    z= (1 - y**2 + 4*y**4 + 1.5*(x_0*y**6)/(xi))/( np.exp(2*x_0*xi)*(1 + y**2)**2)

    return cost*(nu**2)*z



# La funzione che crea la matrice del sistema (si è usata la regola del trapezio, vedi file Costruzuzione_matric_regola_trapez_1D.py)
# freq è il vettore contenente tutte le frequenze considerate.
# M è il vettore contenente i valori della massa totale utilizzati nella discretizzazione (sia freq che M devono avere dimenzione K)

def fun_mat(freq, M):

    global K

    a=np.zeros((K,K))

    for i in range(0,K):

        nu= freq[i]

        for j in range(0,K):

            massa= M[j]


            if ( j==0 or j==(K-1) ):
                a[i][j]= integ(massa, nu)/2

            else:
                a[i][j]= integ(massa, nu)



    dM= (M[-1] - M[0])/(K-1)

    return dM*a


# RICERCA AUTOVALORI E AUTOVETTORI DELLA MATRICE

matrix= fun_mat(freq, masse)/omega_GW

lamb, v= eig(matrix)



# creazione di una matrice 2*K in cui ogni riga contiene un'autovalore e il suo corrispondente autovettore. Tale matrice è realizzata così da avere gli autovettori legati ai loro corrispondenti autovettori e quindi non avere problemi quando poi questi vengono ordinati in senso decrescente

autovett= []

for i in range(0, len(lamb)):

    rig=[]
    autovett.append(rig)



for i in range(0, len(lamb)):

    autovett[i].append(lamb[i].real)
    autovett[i].append(v[i].real)



autovett= np.array(autovett)
autovett= autovett[(-autovett)[:,0].argsort()]



for i in range(0, len(lamb)):

    print(autovett[i])








