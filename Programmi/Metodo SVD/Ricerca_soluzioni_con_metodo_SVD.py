import numpy as np
from numpy.linalg import svd
import matplotlib.pyplot as plt
import math
import keyboard



# OPZIONI VARIE


# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)
K=500

# num_sing contiene il numero di valori singolari che si vogliono considerare per la risoluzione del problema
num_sing= 29

# Estremo superiore e estremo inferiore delle masse considerate nell'integrale e scarto tra due masse consecutive
m_min=1
m_max=30


# Se option è pari a "read" omega_GW è letto da file altrimenti viene creato tramite una funzione
option="read"


# Path del file se option è uguale a "read", num è il numero di valori della frequenza considerati nel file che si vuole aprire
num= 500

file_name_omega="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\Metodo SVD\\file_txt\\omega_GW_" + str(num) + ".txt"


# Se option è pari a "read" ponendo disegna uguale a True verra realizzato il grafico della soluzione esatta nello stesso piano in cui vengono rappresentate le soluzioni trovate minimizzando
disegna=True


# Path del file se disegna è uguale a True
file_name_f_m="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\Metodo SVD\\file_txt\\f_m_" + str(num) + ".txt"


# Estremo superiore e estremo inferiore delle frequenze considerate (da selezionare solo se option è diverso da "read")
freq_min=10**(-8)
freq_max=10**(1)



# Estremi della finestra in massa che si considera una volta ottenuta la funzione F(M) (i valori fuori si escludono poiche sono caratterizzati dalla presenza di artefatti)
mask_min= 15.15
mask_max= 25.15



# Numero di zeri aggiunti a destra e a sinistra dell'array contenente f(M) dopo che è stato selezionato nella sola finestra di cui si è parlato sopra
num_zeri= 200



# Se opzione_smooth è pari a media_mobile, la funzione F(M) viene smussata tramite il metodo della media mobile con finestre indipendenti. Se invece è pari a media_mobile_1, si utilizza il metodo della media mobile con finestre dipendenti.  Se invece è pari a S_V, allora lo smussamento avviene mediante un filtro di Savitzky Golay
opzione_smooth= "media_mobile_1"



# Se opzione_smooth è pari a media_mobile o a media_mobile_1, allora è necessario indicare lung_sottoin, ossia la lunghezza dei sottointervalli di cui si calcola la media
lung_sottoin= 20



# Se opzione_smooth è pari a S_V, allora è necessario indicare wind_size, ossia la dimensione della finestra per la determinazione del polinomio locale, e poly_order, ossia l'ordine del polinomio
wind_size= 10
poly_order= 3









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

        print("\nf_esatta è lungo {0}".format(len(f_esatta)))

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






# Ricampionamento del vettore freq e omega in modo tale che abbia le giuste dimensioni (ossia K)

lung= len(freq)

rap= int(lung/K)
res= lung%K


if(res==0):

    freq= freq[rap-1::rap]
    omega_GW= omega_GW[rap-1::rap]

else:
    freq= freq[res-1+rap::rap]
    omega_GW= omega_GW[res-1+rap::rap]


print("Dopo il ricampionamento il file contenente omega_GW ha dimensione {0}, mentre il numero delle variabili è pari a {0}".format( len(omega_GW), K))




# Divisione di omega_GW per la costante. Tale operazione viene effettuata semplicemnete per incrementare l'ordine di grandezza dei numeri con cui si ha a che fare all'interno del programma

omega_GW= omega_GW/cost










# DEFINIZIONE DELLE FUNZIONI CHE COMPAIONO NELLA RELAZIONE DI RIFERIMENTO UNA VOLTA CHE E' STATO DISCRETIZZATO L'INTEGRALE


# Lista dei valori della massa totale utilizzati per effettuare la discretizzazione dell'integrale (sono in masse solari)

masse=np.linspace(m_min, m_max, K)



# La funzione usata per costruire la matrice, M è la massa totale del sistema, mentre nu è una frequenza

def integ(M, nu):

    nu_0= np.sqrt(a**3/(G*(M)))
    x_0= 2*math.pi*nu_0*nu

    z= (1 - y**2 + 4*y**4 + 1.5*(x_0*y**6)/(xi))/( np.exp(2*x_0*xi)*(1 + y**2)**2)

    return (nu**2)*z



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










# FATTORIZZAZIONE SVD

matrix= fun_mat(freq, masse)


# La funzione svd effettua la fattorizzazione SVD, S e D sono le matrici unitarie, mentre v è un vettore contenente, in ordine decrescente, i valori singolari

S, v, D= svd(matrix)


print("\n\nLista dei valori singolari:\n")

for i in range(0, len(v)):
    print("{:}  {:.3}".format( i, v[i]))










# OTTENIEMTO DELLA SOLUZIONE

v_i= 1/v


for i in range(num_sing, len(v_i)):
    v_i[i]=0


V_i= np.diag(v_i)



F_M= np.dot(D.transpose(), V_i)
F_M= np.dot(F_M, S.transpose())
F_M= np.dot(F_M, omega_GW)










# CALCOLO PRODOTTO DI CONVOLUZIONE PARTENDO DAL RISULTATO ESATTO

# definizione di f(m)

def f_m(m, mu, sigma):

    return (m**2/np.sqrt(2*math.pi*sigma**2))*np.exp(-(m-mu)**2/(2*sigma**2))


mu= 10
sigma= 1


dM= masse[1] - masse[0]


val_conv= masse

conv= np.zeros(len(val_conv))



for i in range(0, len(val_conv)):

    integrale= 0

    for j in range(0, len(masse)):

        prod= f_m(masse[j], mu, sigma)*f_m(val_conv[i]-masse[j], mu, sigma)

        if( j==0 or j==len(masse)-1 ):
            integrale+= dM*prod/2

        else:
            integrale+= dM*prod

    conv[i]= integrale









# SELEZIONE FINESTRA IN CUI LA SOLUZIONE TROVATA COINCIDE CON LA SOLUZIONE CORRETTA E AGGIUNTA DI ZERI


mask= ( masse>mask_min ) & ( masse<mask_max )

masse= masse[mask]
F_M= F_M[mask]


dm= masse[1] - masse[0]

masse_sotto= np.linspace( masse[0] - dm - (num_zeri-1)*dm, masse[0] - dm, num=num_zeri)
masse= np.concatenate((masse_sotto, masse))

masse_sopra= np.linspace(masse[-1] + dm, masse[-1] + dm + (num_zeri-1)*dm, num=num_zeri)
masse= np.concatenate((masse, masse_sopra))


array_zeri= np.zeros(num_zeri)

F_M= np.concatenate((array_zeri, F_M))
F_M= np.concatenate((F_M, array_zeri))












# SMUSSAMENTO TRAMITE LA TECNICA DELLA MEDIA MOBILE (DIVIDO IN FINESTRE INDIPENDENTI, QUINDI SI RIDUCE IL NUMERO DEI PUNTI ALLA FINE)
# divido in sottointervalli e per ognuna calcolo la media, le finestra successiva non condivide alcun punto con quella precedente, quindi la dimensione finale è pari alla dimensione iniziale diviso la lunghezza della finestra

if (opzione_smooth=="media_mobile"):

    masse_medie= np.zeros(int( len(F_M)/lung_sottoin ))
    F_M_medie= np.zeros(int( len(F_M)/lung_sottoin ))

    for i in range(0, int( len(F_M)/lung_sottoin ) ):

        sum= 0

        for j in range( 0, lung_sottoin):

            sum+= F_M[lung_sottoin*i + j]

        media= sum/lung_sottoin

        F_M_medie[i]= media
        masse_medie[i]= masse[lung_sottoin*i + int(j/2)]


    masse= masse_medie
    F_M= F_M_medie









# SMUSSAMENTO TRAMITE LA TECNICA DELLA MEDIA MOBILE (DIVIDO IN FINESTRE CHE CONDIVIDONO DIVERSI PUNTU)
# creo una finestra e calcolo media, la finestra successiva si crea partendo dalla precedente facendo scorrere di un valore a destra, quindi vengono condivisi un numero di punti pari alla lunghezza della finestra meno uno. La dimensione finale è pari a alla dimensione iniziale meno la lunghezza della finestra



if (opzione_smooth=="media_mobile_1"):

    masse_medie_1= np.zeros( len(F_M) - lung_sottoin )
    F_M_medie_1= np.zeros( len(F_M) - lung_sottoin )


    for i in range( int(lung_sottoin/2), len(F_M) - int(lung_sottoin/2)):

        sum= 0

        for j in range( -int(lung_sottoin/2), int(lung_sottoin/2)):

            sum+= F_M[i + j]

        media= sum/lung_sottoin

        F_M_medie_1[i-int(lung_sottoin/2)]= media
        masse_medie_1[i-int(lung_sottoin/2)]= masse[i]


    masse= masse_medie_1
    F_M= F_M_medie_1










# SMUSSAMENTO TRAMITE FILTRO DI SAVITZKY GOLAY

if (opzione_smooth=="S_V"):

    F_M= scipy.signal.savgol_filter(F_M, wind_size, poly_order)










# GRAFICO DELLA SOLUZIONE INDIVIDUATA E CONFRONTO CON SOLUZIONE ESATTA

fig, ax= plt.subplots()

ax.set_title("Soluzione Individuata Considerando {0} Valori Singolari".format(num_sing))

ax.plot([masse[0],masse[-1]], [0,0], linestyle="-", color="black", marker="", linewidth=0.75, alpha=1)

ax.plot(masse, F_M, linestyle="-", color="blue", marker="", label="Soluzione Ottenuta")
ax.plot(val_conv, conv, linestyle="-", color="red", marker="", label="Soluzione Esatta: $\\mu$= {0}, $\\sigma$= {1}".format( mu, sigma))


ax.set_xlabel("M [M_sun]")
ax.set_ylabel("F(M)")
ax.set_xlim(max(masse[0], val_conv[0]), min(masse[-1], val_conv[-1]))


ax.legend()
plt.tight_layout()











# calcolo della trasformata del prodotto di convoluzione e sua rappresentazione grafica
tras= np.fft.fft(F_M)

n= F_M.size
dm= masse[1] - masse[0]
k_masse= np.fft.fftfreq( n, d=dm)

asse_x= np.fft.fftshift(k_masse)

asse_x= k_masse

fig= plt.figure()

plt.subplot(2,2,1)

plt.title("Trasformata del Prodotto di Convoluzione (Componenete Reale)")

plt.plot(asse_x, tras.real)
plt.plot(asse_x, 0*asse_x, linestyle="-", color="black", linewidth=0.8)

plt.xlim(min(asse_x), max(asse_x))

plt.xlabel("Frequenza")
plt.ylabel("Parte Reale")


plt.subplot(2,2,3)

plt.title("Trasformata del Prodotto di Convoluzione (Componenete Immaginaria)")

plt.plot(asse_x, tras.imag)
plt.plot(asse_x, 0*asse_x, linestyle="-", color="black", linewidth=0.8)

plt.xlim(min(asse_x), max(asse_x))

plt.xlabel("Frequenza")
plt.ylabel("Parte Immaginaria")

plt.subplot(2,2,2)

plt.title("Trasformata del Prodotto di Convoluzione (Modulo)")

plt.plot(asse_x, abs(tras))
plt.plot(asse_x, 0*asse_x, linestyle="-", color="black", linewidth=0.8)

plt.xlim(min(asse_x), max(asse_x))

plt.xlabel("Frequenza")
plt.ylabel("Modulo")

plt.subplot(2,2,4)

plt.title("Trasformata del Prodotto di Convoluzione (Fase)")

plt.plot(asse_x, np.angle(tras))
plt.plot(asse_x, 0*asse_x, linestyle="-", color="black", linewidth=0.8)

plt.xlim(min(asse_x), max(asse_x))

plt.xlabel("Frequenza")
plt.ylabel("Fase")

plt.tight_layout()


# calcolo di f_m

uni= np.ones(len(tras))

'''
for i in range(0, int(len(uni)/2)):

    var= np.random.rand()

    estraz= int(var + 0.5) - (1 - int(var + 0.5))

    uni[i]= estraz
    uni[len(uni)-1-i]= estraz

#tras=tras*uni
'''

f_m= np.fft.ifft( np.sqrt(tras) )

n= tras.size
dk= k_masse[1] - k_masse[0]
masse_f_m= np.fft.fftfreq( n, d=dk)
masse_f_m= np.fft.fftshift(masse_f_m)

fig= plt.figure()

plt.subplot(2,1,1)

plt.plot(masse_f_m, f_m, linestyle="-", color="blue", label="Soluzione Individuata")

if ( disegna==True ):

    def f_m_funzione(m, mu, sigma):

        return (m**2/np.sqrt(2*math.pi*sigma**2))*np.exp(-(m-mu)**2/(2*sigma**2))

    mu= 10
    sigma= 1

    plt.plot(masse_f_m, f_m_funzione(masse_f_m, mu, sigma), linestyle="-", color="orange", label="Soluzione Esatta")


plt.title("FUNZIONE LOGARITMICA DI MASSA")
plt.xlabel("Massa [M_sole]")
plt.ylabel("f(m)")

plt.legend()


plt.subplot(2,1,2)

plt.title("Fase")

plt.plot(masse_f_m, np.angle(f_m), linestyle="", marker=".", color="blue", label="Soluzione Individuata")



plt.tight_layout()
plt.show()
