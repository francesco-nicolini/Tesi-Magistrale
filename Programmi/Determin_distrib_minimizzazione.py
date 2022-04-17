import numpy as np
import matplotlib.pyplot as plt
import math



# OPZIONI VARIE

# N contiene il numero delle iterazioni
N=200

# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)
K=50

# b è l'ampiezza del passo
b=10**(6)

# Estremo superiore e estremo inferiore delle masse considerate nell'integrale e scarto tra due masse consecutive
m_min=1
m_max=30

dm= (m_max - m_min)/(K-1)

# Se option è pari a "read" omega_GW è letto da file altrimenti viene creato tramite una funzione
option="read"

# Path del file se option è uguale a "read"
file_name="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt\\omega_GW.txt"

# Estremo superiore e estremo inferiore delle frequenze considerate (da selezionare solo se option è diverso da "read"
freq_min=10**(-8)
freq_max=10**(0)






# VALORI INIZIALI DELLA DISTRIBUZIONE DA CERCARE

# Contiene il valore di partenza delle variabile
val_iniz=np.ones(K)


if (K%2==0):

    for i in range(1,int(K/2)):

        val_iniz[i]=val_iniz[i-1]+2/K
        val_iniz[K-1-i]=val_iniz[i]

if (K%2==1):

    for i in range(1, int((K-1)/2)):

        val_iniz[i]=val_iniz[i-1]+2/K
        val_iniz[K-1-i]=val_iniz[i]

    val_iniz[int((K-1)/2)]= val_iniz[int((K-1)/2)-1]  + 2/K







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







# CREAZIONE DEL FILE CONTENENTE LE FREQUENZE E DELL'ARRAY CONTENENTE I VALORI DI OMEGA_GW

# freq contiene la lista delle frequenze
# omega_GW contiene la lista dei valori dello spettro in potenza del fondo

if (option=="read"):

    freq, omega_GW=np.loadtxt(file_name, unpack=True)

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

plt.show()





# DEFINIZIONE DELLE FUNZIONI CHE COMPAIONO NELLA RELAZIONE DI RIFERIMENTO UNA VOLTA CHE E' STATO DISCRETIZZATO L'INTEGRALE


# Lista delle masse indagate nella discretizzazione dell'integrale (sono in masse solari)

masse=np.linspace(m_min, m_max, K)



# La funzione usata per costruire le matrici (si tiene conto anche del dm**2 in essa)

def integ(m_1, m_2, nu):

    nu_0= np.sqrt(a**3/(G*(m_1 + m_2)))
    x_0= 2*math.pi*nu_0*nu

    z= dm**2*(1 - y**2 + 4*y**4 + 1.5*(x_0*y**6)/(xi))/( np.exp(2*x_0*xi)*(1 + y**2)**2)

    return (nu**2)*z



# La funzione che crea le matrici (si è usata la regola del trapezio, vedi file Prova_costruz_matrice_per_regola_trapez.py)

def fun_mat(nu, K):

    a=np.zeros((K,K))

    for i in range(0,K):

        for j in range(0,K):

            if ( (i+j==0) or (i+j==2*(K-1)) or (i*j==0 and i+j==K-1) ):

                a[i][j]=integ(masse[i], masse[j], nu)/4


            # è un elif quindi si possono mettere casi già previsti dall'if
            #precedente dal momento che questi, verificando il primo if, non saranno
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



# Funzioni per il calcolo delle derivate parziali, queste sono inserite nella lista chiamata deriv (l è l'indice della variabile rispetto alla quale si vuole calcolare la derivata parziale)

deriv=[]

for i in range(0, K):

    def f_deriv(*var, l, i=i):

        sum= 2*var[l]*matrix[i][l][l]

        for j in range (0, K):

            if (j!=l):

                sum=sum + 2*var[j]*matrix[i][j][l]

        return sum

    deriv.append(f_deriv)






# CREAZIONE DELLA FUNZIONE DA MINIMIZZARE E CALCOLO DELLE SUE DERIVATE PARZIALI

# Questa è la funzione di cui si derve cercare il minimo
def funz_da_minim(*var):

    somma=0

    for i in range(0, len(funz)):
        somma=somma + (funz[i](*var))**2

    return somma


# Funzione per il calcolo delle derivate parziali della funzione funz_da_minim (index è l'indice della variabile rispetto alla quale si vuole calcolare la derivata parziale)

def gradiente(*var, index):

    somma=0

    for i in range(0, len(funz)):

        somma= somma + funz[i](*var)*deriv[i](*var, l=index)

    return 2*somma







# ESECUZIONE DELL'ALGORITMO DI MINIMIZZAZIONE

print("Funzione calcolata con i valori iniziali:", funz_da_minim(*val_iniz))


# Vettori introdotti per poter fare i grafici finali

graf=[]

for i in range(0,K):

    a=np.zeros(N+1)
    graf.append(a)


val=val_iniz



# array che contiene i valori del nuovo punto (esempio in due dimenzioni: se mettessi val[0]=val[0] - b*derivate(*val) andrebbe bene, ma val[1]=val[1] - b*derivate(*val) darebbe problemi poichè val a destra conterrebbe il vecchio valore di val[1] e il nuovo valore di val[0] e non i vecchi valori di entrambi

provv= np.ones(len(val))



# Ciclo con le iterazioni

for i in range(0,N):

    for j in range(0, len(val)):
        graf[j][i]=val[j]


    for j in range(0, len(val)):
        provv[j]= val[j] - b*gradiente(*val, index=j)

    val= provv



# Inserimento delle soluzioni trovate nelle liste per la realizzazione dei grafici

for i in range(0, len(val)):
    graf[i][N]=val[i]



# Stampa dei risultati finali

print("\nSoluzioni:\n",val,"\n")

print("Funzione calcolate con i valori delle variabili individuati:")
print("f =", funz_da_minim(*val))







# REALIZZAZIONE DEI GRAFICI

# Grafico degli andamenti delle soluzioni all'aumentare del numero di iterazioni

if (K<30):

    plt.figure()

    plt.title("Valore delle Variabili al Variare del Numero di Iterazioni")
    plt.xlabel("Numero iterazioni")
    plt.ylabel("Valori variabili")

    iter=np.linspace(0,N,N+1)

    for i in range(0,K):
        plt.plot(iter, graf[i], color="C{0}".format(i), label="f(m_{0})".format(i), linestyle="-", marker="", markersize=3)

    plt.legend()
    plt.show()



# Grafico della soluzione ottenuta in funzione della massa

plt.figure()

plt.title("f(m) in funzione di m")
plt.xlabel("massa [$M_sole]")
plt.ylabel("f(m)")


plt.plot(masse, val, color="blue", linestyle="", marker="o")

plt.tight_layout()
plt.show()
