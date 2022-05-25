import numpy as np
import matplotlib.pyplot as plt
import math
import keyboard



# OPZIONI VARIE

# N contiene il numero delle iterazioni
N=5000

# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)
K=500

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

            #print("massa=  {0}, nu={1}\n".format( massa, nu))
            #print("tipo massa=  {0}, tipo nu={1}\n".format( type(massa), type(nu)))


            if ( j==0 or j==(K-1) ):

                #print("integ= {0}".format(integ(M[0], freq[0])[int(K/2)]/2))
                #print("tipo integ= {0}".format(type(integ(M[0], freq[0])/2)))

                a[i][j]= integ(massa, nu)/2

            else:

                a[i][j]= integ(massa, nu)



    dM= (M[-1] - M[0])/(K-1)

    return dM*a




# COSTRUZIONE DEL PRODOTTO DI CONVOLUZIONE


# definizione di f(m)

def f_m(m, mu, sigma):

    return (m**2/np.sqrt(2*math.pi*sigma**2))*np.exp(-(m-mu)**2/(2*sigma**2))



# convoluzione
'''
masse= np.linspace(0, 30, K)

mu= 10
sigma= 1

a= f_m(masse, mu, sigma)
b= a

conv= dM*np.convolve(a, b, mode="full")
asse_conv= np.linspace(0, len(conv), len(conv), endpoint=False)*dM

figure=plt.figure()

plt.plot(asse_conv[::2], conv[::2])
'''

mu= 10
sigma= 1


dM= masse[1] - masse[0]


val_conv= masse

conv= np.zeros(len(val_conv))



for i in range(0, len(val_conv)):

    integrale= 0

    for j in range(0, len(masse)):

        prod= f_m(masse[j], mu, sigma)*f_m(val_conv[i]-masse[j], mu, sigma)

        if( j==0 or j==(K-1) ):
            integrale+= dM*prod/2

        else:
            integrale+= dM*prod

    conv[i]= integrale

# grafico prodotto di convoluzione

plt.title("Prodotto di convoluzione", fontsize=14)

plt.xlabel("M [M_solare]", fontsize=10)
plt.ylabel("f(m)*f(m)", fontsize=10)

plt.plot(val_conv, conv, linestyle="-", marker="", color="blue")



# RISOLUZIONE DEL SISTEMA

masse= np.linspace(m_min, m_max, K)



matrix= fun_mat(freq, masse)

omega_conv= np.dot(matrix, conv)


# CONFRONTO

dif_max= 0
i_max= 0

for i in range(0, len(omega_GW)):

    diff= abs( (omega_conv[i] - omega_GW[i])/omega_GW[i]  )

    if( diff> dif_max ):
        dif_max= diff
        i_max= i


print("Il valore massimo della differenza relativa ( modulo di omega ottenuto con convoluzione meno omega \"esatto\" diviso quest'ultimo) è pari a {:.2e} e corrisponede alla frequenza {:.2e}".format(dif_max, freq[i_max]))



# GRAFICI

fig, ax = plt.subplots(2)

ax[0].plot(freq, omega_GW, linestyle=(0, (1, 1)), color="blue", label="Soluzione Esatta")
ax[0].plot(freq, omega_conv, linestyle="--", color="red", label="Calcolo con Prodotto\ndi Convoluzione")


ax[0].set_title("$\\Omega_{GW}$ in funzione della frequenza", fontsize=14)
ax[0].set_xlabel("f [Hz]", fontsize=10)
ax[0].set_ylabel("$\\Omega_{GW}$", fontsize=10)
ax[0].set_xlim(min(freq), max(freq))
ax[0].set_xscale("log")
ax[0].set_yscale("log")

ax[0].legend()


scarto= (omega_conv-omega_GW)/omega_GW

ax[1].plot(freq, abs(scarto), linestyle="-", color="blue")

ax[1].plot( freq[i_max], dif_max, marker=".", color="blue")
ax[1].text( freq[i_max], dif_max, "{:.2e}".format(dif_max), horizontalalignment="right")

ax[1].set_title("Differenza Relativa tra le $\\Omega_{GW}$ ", fontsize=14)
ax[1].set_xlabel("f [Hz]", fontsize=10)
ax[1].set_ylabel("$|\\Delta\\Omega_{GW}/\\Omega_{GW}|$", fontsize=10)
ax[1].set_xlim(min(freq), max(freq))
ax[1].set_ylim(10**(-13), 10**(-5))
ax[1].set_xscale("log")
ax[1].set_yscale("log")




plt.tight_layout()
plt.show()



























