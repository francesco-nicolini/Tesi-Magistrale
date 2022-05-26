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

matrix= fun_mat(freq, masse)

lamb, v= eig(matrix)



# creazione di una matrice 2*K in cui ogni riga contiene un'autovalore e il suo corrispondente autovettore. Tale matrice è realizzata così da avere gli autovettori legati ai loro corrispondenti autovettori e quindi non avere problemi quando poi questi vengono ordinati in senso decrescente

autovett= []

for i in range(0, len(lamb)):

    rig=[]
    autovett.append(rig)



for i in range(0, len(lamb)):

    autovett[i].append(lamb[i].real)
    autovett[i].append(v[i].real)



autovett= np.array(autovett, dtype=object)
autovett= autovett[(-abs(autovett))[:,0].argsort()]



for i in range(0, len(lamb)):

    print(autovett[i])



# RICERCA DEI COEFFICIENTI DELLA COMBINAZIONE LINEARE DI AUTOSTATI DELLA MATRICE CHE MEGLIO APPROSSIMA omega_GW

num_autost= 0

# itero la ricerca di tali coefficienti aumentando ogni volta il numero di autostai considerati. La ricerca termina quando la differenza tra omega_GW e tale combinazione (mem_err) risulta più grande dell'iterazione precedente.

mem_err= 10**(10)


while(1):

    num_autost+= 1

    rap= int( K/num_autost )
    res= K%num_autost


    if(res==0):
        ini= rap-1

    else:
        ini= res-1+rap


    omega_per_appros= omega_GW[ini::rap]


    # Creazione della matrice

    mat= np.zeros((num_autost, num_autost))

    for j in range(0, num_autost):

        autostato= autovett[j][1]
        autostato= autostato[ini::rap]

        for i in range(0, num_autost):

            mat[i][j]= autostato[i]


    # Indivuduazione coefficienti e calcolo errore

    coeff= np.linalg.solve(mat, omega_per_appros)


    omega_appros= np.zeros( len(autovett[0][1]) )

    for i in range(0, num_autost):

        omega_appros= coeff[i]*autovett[i][1]


    err= 0

    for i in range(0,K):

        err+= np.abs(omega_appros[i]-omega_GW[i])

    if ( err> mem_err):
        break

    mem_err= err



print("L'approssimazione migliore si ottiene con {0} autovettori".format(num_autost))




# CONFRONTO

dif_max= 0
i_max= 0

for i in range(0, len(omega_GW)):

    diff= abs( (omega_appros[i] - omega_GW[i])/omega_GW[i]  )

    if( diff> dif_max ):
        dif_max= diff
        i_max= i


print("Il valore massimo della differenza relativa ( modulo di omega ottenuto con convoluzione meno omega \"esatto\" diviso quest'ultimo) è pari a {:.2e} e corrisponede alla frequenza {:.2e}".format(dif_max, freq[i_max]))






# GRAFICO PER IL CONFRONTO DELL'APPROSSIMAZIONE MIGLIORE CON omega_GW

fig, ax = plt.subplots(2)

ax[0].plot(freq, omega_GW, linestyle=(0, (1, 1)), color="blue", label="Soluzione Esatta")
ax[0].plot(freq, omega_appros, linestyle="--", color="red", label="Approssimazione Migliore")


ax[0].set_title("$\\Omega_{GW}$ in funzione della frequenza", fontsize=14)
ax[0].set_xlabel("f [Hz]", fontsize=10)
ax[0].set_ylabel("$\\Omega_{GW}$", fontsize=10)
ax[0].set_xlim(min(freq), max(freq))
ax[0].set_xscale("log")
ax[0].set_yscale("log")

ax[0].legend()


scarto= (omega_appros-omega_GW)/omega_GW

ax[1].plot(freq, abs(scarto), linestyle="-", color="blue")

ax[1].plot( freq[i_max], dif_max, marker=".", color="blue")
ax[1].text( freq[i_max], dif_max, "{:.2e}".format(dif_max), horizontalalignment="right")

ax[1].set_title("Modulo della Differenza Relativa tra le $\\Omega_{GW}$ ", fontsize=14)
ax[1].set_xlabel("f [Hz]", fontsize=10)
ax[1].set_ylabel("$|\\Delta(\\Omega_{GW})/\\Omega_{GW}|$", fontsize=10)
ax[1].set_xlim(min(freq), max(freq))
#ax[1].set_ylim(10**(-13), 10**(-5))
ax[1].set_xscale("log")
ax[1].set_yscale("log")




plt.tight_layout()
plt.show()





















