import numpy as np
import matplotlib.pyplot as plt
import math
import keyboard



# OPZIONI VARIE

# N contiene il numero delle iterazioni
N=5000

# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)
K=50


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

    print("Il vettore contenente omega_GW ha una dimensione pari a {0}\n".format(len(omega_GW)))

    if (disegna==True):

        #  masse_graf contiene le masse per effettuare il grafico di f_esatta
        # f_esatta la forma della distribuzione letta da file (è la soluzione esatta)

        masse_graf, f_esatta= np.loadtxt(file_name_f_m, unpack=True)

        print("f_esatta è lungo {0}".format(len(f_esatta)))

else:

    def funz_omeg(nu):

        n=len(nu)

        return 10**(-22)*np.ones(n)


    freq=np.logspace(np.log10(freq_min), np.log10(freq_max), K)
    omega_GW=funz_omeg(freq)


# Grafico di omega_GW

fig, ax = plt.subplots()

ax.plot(freq, omega_GW, linestyle="-", color="blue")

plt.title("$\\Omega_{GW}$ in funzione della frequenza")
plt.xlabel("f [Hz]")
plt.ylabel("$\\Omega_{GW}$")
plt.yscale("log")
plt.xscale("log")











# DEFINIZIONE DELLE FUNZIONI CHE COMPAIONO NELLA RELAZIONE DI RIFERIMENTO UNA VOLTA CHE E' STATO DISCRETIZZATO L'INTEGRALE


# Lista dei valori della massa totale utilizzati per effettuare la discretizzazione dell'integrale (sono in masse solari)

masse=np.linspace(m_min, m_max, K)



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

    a=np.ones((K,K))

    for i in range(0,K):

        nu= freq[i]

        for j in range(0,K):

            if ( j==0 or j==(K-1) ):

                a[i][j]= integ(M[j], nu)/2

            else:

                a[i][j]= integ(M[j], nu)


    dM= (M[-1] - M[0])/(K-1)

    return dM*a








# RISOLUZIONE DEL SISTEMA


matrix= fun_mat(freq, masse)

F_M= np.linalg.solve( matrix, omega_GW)

if ( np.allclose(np.dot(matrix, F_M), omega_GW) ):

    print("Soluzione individuata")




# GRAFICO CON IL CONFORNTO DELLA CONVOLUZIONE DELLA FUNZIONE ESATTA


if ( disegna==True ):

    asse_F_M= np.linspace(0, len(F_M), len(F_M), endpoint=False)

    conv= 2*np.convolve( f_esatta, f_esatta, mode="full")
    asse_conv= np.linspace(0, len(conv), len(conv), endpoint=False)


    fig=plt.figure()

    plt.plot(asse_F_M, F_M, linestyle="-", color="blue", label="Soluzione Sistema")
    plt.plot(asse_conv, conv, linestyle="-", color="orange", label="Convoluzione Soluzione Corretta")


    minimo_x= min(np.amin(asse_F_M), np.amin(asse_conv))
    massimo_x= max(np.amax(asse_F_M), np.amax(asse_conv))

    minimo_y= min(np.amin(F_M), np.amin(conv))
    massimo_y= max(np.amax(F_M), np.amax(conv))


    wind_x= 100
    wind_y= 1


    asse_x= np.linspace(minimo_x - wind_x, massimo_x + wind_x, 100)
    plt.plot(asse_x, np.zeros(len(asse_x)), linestyle="-", color="black", linewidth=0.8)


    plt.xlim(minimo_x - wind_x, massimo_x + wind_x)
    plt.ylim(minimo_y - wind_y, massimo_y + wind_y)

    plt.title("CONFRONTO DELLA SOLUZIONE DEL SISTEMA E DELLA CONVOLUZIONE DELLA SOLUZIONE CORRETTA")
    plt.xlabel("N")
    plt.ylabel("f(m)*f(m)")


    plt.legend()
    plt.tight_layout()







# DETERMINAZIONE DI f(m) E SUA RAPPRESENTAZIONE

# calcolo della trasformata
tras= np.fft.fft(F_M)

n= F_M.size
dm= masse[1] - masse[0]
k_masse= np.fft.fftfreq( n, d=dm)


# calcolo di f_m
f_m= np.fft.ifft( np.sqrt(tras/2) )

n= tras.size
dk= k_masse[1] - k_masse[0]
masse_f_m= np.fft.fftfreq( n, d=dk)
masse_f_m= np.fft.fftshift(masse_f_m)

fig= plt.figure()

plt.plot(masse_f_m, abs(f_m), linestyle="-", color="blue", label="Soluzione Individuata")


if ( disegna==True ):

    plt.plot(masse_graf, f_esatta, linestyle="-", color="orange", label="Soluzione Esatta")


plt.title("FUNZIONE LOGARITMICA DI MASSA")
plt.xlabel("Massa [M_sole]")
plt.ylabel("f(m)")

plt.legend()

plt.tight_layout()
plt.show()





