import numpy as np
from numpy.linalg import svd
import matplotlib.pyplot as plt
import math
from scipy.optimize import minimize
from scipy.optimize import Bounds
from scipy.optimize import dual_annealing


# OPZIONI VARIE


# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)
K=1000

# num_sing contiene il numero di valori singolari che si vogliono considerare per la risoluzione del problema
num_sing= 25

# Estremo superiore e estremo inferiore delle masse considerate nell'integrale e scarto tra due masse consecutive
m_min=1
m_max=25


# Se option è pari a "read" omega_GW è letto da file altrimenti viene creato tramite una funzione
option="read"


# Path del file se option è uguale a "read", num è il numero di valori della frequenza considerati nel file che si vuole aprire
num= 1000

file_name_omega="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\Metodo SVD\\Gaussiana_doppia\\file_txt\\omega_GW_" + str(num) + ".txt"


# Se option è pari a "read" ponendo disegna uguale a True verra realizzato il grafico della soluzione esatta nello stesso piano in cui vengono rappresentate le soluzioni trovate minimizzando
disegna=True


# Path del file se disegna è uguale a True
file_name_f_m="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\Metodo SVD\\Gaussiana_doppia\\file_txt\\f_m_" + str(num) + ".txt"



# Estremo superiore e estremo inferiore delle frequenze considerate (da selezionare solo se option è diverso da "read")
freq_min=10**(-8)
freq_max=10**(1)



# Estremi della finestra in massa che si considera una volta ottenuta la funzione F(M) (i valori fuori si escludono poiche sono caratterizzati dalla presenza di artefatti)
mask_min= 0
mask_max= 60



# Se scelta_valori_singolari è posta pari ad "auto", allora il numero dei valori singolari considerato per la determinazione della soluzione viene scelto imponendo che il rappporto tra la somma di questi ultimi e la somma di tutti sia pari a 1-epsilon con epsilon scelto. Se scelta_valori_singolari assume un qualsiasi altro valore allora il numero di valori singolari considerato è quello fornito.
scelta_valori_singolari="auto"
epsilon= 1e-15



# Numero di zeri aggiunti a destra e a sinistra dell'array contenente f(M) dopo che è stato selezionato nella sola finestra di cui si è parlato sopra
num_zeri= 100



# Se opzione_smooth è pari a "media_mobile_1", la funzione F(M) viene smussata tramite il metodo della media mobile con finestre indipendenti. Se invece è pari a "media_mobile_1", si utilizza il metodo della media mobile con finestre dipendenti.  Se invece è pari a "S_V", allora lo smussamento avviene mediante un filtro di Savitzky Golay
opzione_smooth= "media_mobile_1"



# Se opzione_smooth è pari a media_mobile, allora è necessario indicare lung_sottoin, ossia la lunghezza dei sottointervalli di cui si calcola la media
lung_sottoin= 10



# Se opzione_smooth è pari a S_V, allora è necessario indicare wind_size, ossia la dimensione della finestra per la determinazione del polinomio locale, e poly_order, ossia l'ordine del polinomio
wind_size= 10
poly_order= 3



# Se elimina_negativi è pari a "parabola", allora, negli intervalli di massa in cui F(M) assume valori negativi, si sostituisce a questa una parabola definita positiva della forma a(x-b)**2. Tale operazione si può effettuare in quanto almeno un punto subito a destra o subito a sinistra dell'intervallo è nullo per costruzione, quindi b e pari al corrispondente valore in massa di quest'ultomo, mentre a è dato dal rapporto tra il valore di F(M) nell'altro punto appena fuori dall'intervallo e il quadrato tra la differenza del valore sull'ascisse corrispondente e il parametro b
elimina_negativi="parabola"



# Se elimina_negativi="parabola", si trasformano in parabola tutti i punti degli intervalli in cui F(M) è minore di zero. Se gli indici degli estremi di questi intervalli sono i e j, allora il valore dell'indice  dell'estremo con il valore di F(M) più basso non viene modificato, mentre l'altro può essere spostato di shift_da_bordo (verso sinistra se questo indice è i e verso destra se invece è j)
shift_da_bordo= 10




# L'opzione funzione definisce il tipo di funzione da minimizzare nella determinazione di f_m: se pari a "semplice" è pari semplicemente alla somma dei quadrati delle differenze tra f_m*f_m e F(M); se è pari a "derivata_prima" si aggiunge un termine pari alla somma dei quadrati delle differenze tra elementi successivi di f_m (derivata prima discretizzata) per cost_prima; se è pari a "derivata_seconda" si aggiunge la somma dei quadrati della derivata seconda discretizzata moltiplicata per cost_seconda
funzione= "derivata_seconda"
cost_prima= 0.00000001
cost_seconda= 0.01




# se si pone valori_iniziali pari a gaussiana, allora il programma individua la gaussiana che meglio approssima F_M, trova quindi i parametri (ampiezza, media e deviazione standard) della gaussiana il cui prodotto di convoluzione con se stessa restituisce l'altra gaussiana e la utilizza come condizione iniziale. Con qualunque altro valore rende la scelta delle consizioni iniziali personalizzabile
valori_iniziali="gaussiana"




# f_m_val_iniz contiene la lista dei valori iniziali per f_m richiesti per procedere con la minimizazione
'''
f_m_val_iniz=[1]*( )
'''







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

import os

print(os.getcwd())





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


    dM= (M[1] - M[0])

    return dM*a










# FATTORIZZAZIONE SVD

matrix= fun_mat(freq, masse)


# La funzione svd effettua la fattorizzazione SVD, S e D sono le matrici unitarie, mentre v è un vettore contenente, in ordine decrescente, i valori singolari

S, v, D= svd(matrix)


print("\n\nLista dei valori singolari:\n")

for i in range(0, len(v)):
    print("{:}  {:.3}".format( i, v[i]))










# GRAFICO DEI VALORI SINGOLARI IN FUNZIONE DEL LORO INDICE
plt.figure()

asse_x= np.linspace(1, K+1, K, endpoint=False)

plt.title("Valori Singolari in Funzione del Loro Indice")
plt.plot( asse_x, v, marker="o", color="blue", linestyle="")
plt.xscale("log")
plt.yscale("log")

plt.xlabel("i")
plt.ylabel("valori singolari")









# DETERMINAZIONE DEL NUMERO DI VALORI SINGOLARI

if(scelta_valori_singolari=="auto"):

    somma= 0
    totale= v.sum()

    for i in range(0, len(v)):

        somma+= v[i]

        if ( (somma/totale)>=(1-epsilon) ):

            num_sing= i
            break










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

def f_m_funzione(m, A_0, mu_0, sigma_0, A_1, mu_1, sigma_1):

    if ( m<0 ):
        return 0

    else:
        return A_0*np.exp(-(m-mu_0)**2/(2*sigma_0**2)) + A_1*np.exp(-(m-mu_1)**2/(2*sigma_1**2))


A_0= 10
mu_0= 5
sigma_0= 1

A_1= 10
mu_1= 10
sigma_1= 1


dM= masse[1] - masse[0]

array= np.zeros(len(masse))

for i in range(0, len(masse)):

    array[i]=  f_m_funzione(masse[i], A_0, mu_0, sigma_0, A_1, mu_1, sigma_1)

conv= dM*np.convolve( array, array, mode="full")

val_conv= np.linspace(0, len(conv), len(conv), endpoint=False)*dM + 2*masse[0]









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

    masse_medie= []
    F_M_medie= []

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










# SOSTITUZIONE PARTI NEGATIVE DI F_M CON PARABOLE

if(elimina_negativi=="parabola"):


    F_M_par= np.zeros(len(F_M))

    for i in range(0, len(F_M)):

        F_M_par[i]= F_M[i]



    # ricerca indici degli estremi degli intervalli in cui F_M è negativa

    indici= np.where( F_M_par<0 )

    indici= indici[0]


    inizio=[]
    fine=[]

    if ( (indici[1]-indici[0])>1 ):

        # il primo indice definisce un intervallo con un solo punto
        inizio.append(indici[0])
        fine.append(indici[0])

    else:

        inizio.append(indici[0])


    for i in range( 1, len(indici)-1):


        if ( ((indici[i] - indici[i-1])>1) and ((indici[i+1] - indici[i])>1) ):

            # l'entirvallo è costituito da un solo punto

            inizio.append(indici[i])
            fine.append(indici[i])

        elif ( ((indici[i] - indici[i-1])>1) and ((indici[i+1] - indici[i]) == 1) ):

            # questo indice segna l'inizio di un intervallo

            inizio.append(indici[i])

        elif ( ((indici[i+1] - indici[i])>1) and ((indici[i] - indici[i-1])==1) ):

            # questo indice segna la fine di un intervallo

            fine.append(indici[i])




    if ( (indici[-1]-indici[-2])>1 ):

        # l'ultimo indice definisce un intervallo con un solo punto
        inizio.append(indici[-1])
        fine.append(indici[-1])

    else:

        fine.append(indici[-1])


    print(len(inizio), len(fine))

    # stampa degli estremi degli intervalli di massa in cui F(M) è minore di zero
    print("\nGli intervalli in massa in cui F(M) è minore di zero sono:")

    for i in range(0, len(fine)):

        print("[{:.3},{:.3}]".format(masse[inizio[i]], masse[fine[i]]))



    # sostituzione degli intervalli con parabole

    # ind_min e ind_max sono gli indici degli estremi dell'intervallo in cui la funzione è minor di zero

    def coefficienti (ind_min, ind_max):

        x_0= masse[ind_min-1]
        x_1= masse[ind_max+1]
        y_0= F_M_par[ind_min-1]
        y_1= F_M_par[ind_max+1]

        if ( y_0<y_1 ):
            b= x_0
            a= y_1/(x_1 - b)**2

        else:
            b= x_1
            a= y_0/(x_0 - b)**2

        print("\na= {:.4}, b={:.4}".format(a, b))

        return a, b


    def parabola(x, a, b):

        return a*(x - b)**2


    for k in range(0, len(fine)):



        if ( F_M_par[inizio[k] -1 ]<F_M_par[fine[k] + 1] ):
            ind_sx= inizio[k]
            ind_dx= fine[k] + shift_da_bordo

        else:
            ind_sx= inizio[k] - shift_da_bordo
            ind_dx= fine[k]


        a, b= coefficienti(ind_sx, ind_dx)

        for i in range (ind_sx, ind_dx):

            F_M_par[i]= parabola(masse[i], a, b)

    indici_rest= np.where( F_M_par < 0 )

    for i in indici_rest:

        F_M_par[i]= 0

    F_M= F_M_par

    print("\n\n\n")







# GRAFICO DELLA SOLUZIONE INDIVIDUATA E CONFRONTO CON SOLUZIONE ESATTA

fig, ax= plt.subplots()

ax.set_title("Soluzione Individuata Considerando {0} Valori Singolari".format(num_sing))

ax.plot([masse[0],masse[-1]], [0,0], linestyle="-", color="black", marker="", linewidth=0.75, alpha=1)

ax.plot(masse, F_M, linestyle="-", color="blue", marker="", label="Soluzione Ottenuta")
ax.plot(val_conv, conv, linestyle="-", color="red", marker="", label="Soluzione Esatta:\n$A_0$={0}, $\\mu_0$= {1}, $\\sigma_0$= {2}\n$A_1$={3}, $\\mu_1$= {4}, $\\sigma_1$= {5}".format( A_0, mu_0, sigma_0, A_1, mu_1, sigma_1))


ax.set_xlabel("M [M_sun]")
ax.set_ylabel("F(M)")
ax.set_xlim( max(masse[0], val_conv[0]), min(masse[-1], val_conv[-1]))


ax.legend()
plt.tight_layout()


plt.show()







# RICERCA DI f(m) MEDIANTE UN ALGORITMO DI MINIMIZZAZIONE


# dal momento che si vuole effettuare la minimizzazione considerando esclusivamente i valori positivi della massa, seleziono esclusivamente i valori di F_M corrispondenti

mask= masse>0

masse= masse[mask]
F_M= F_M[mask]

if( len(F_M)%2!=0 ):

    masse= masse[:-1]
    F_M= F_M[:-1]



# L'array da minimizzare corrispondente ad f_m deve essere definito in un intervallo in massa che ha come estremi 0 e la meta del valore massimo dell'intervallo in cui è definito F_M, inoltre questa seconda lista di masse si deve costruire con il doppio della sensibilità (la differenza tra valori successivi della massa è pari a metà rispetto alla lista di F_M). In questo modo, la lista del prodotto di convoluizone avrà il doppio della dimensione della lista di F_M, tuttavia mediando elementi successivi si avrà uno stesso numero di componenti, inoltre le due liste saranno definite nello stesso intervallo in massa

masse_f_m= np.linspace(masse[0]/2, masse[-1]/2, len(masse))
dm_f= masse_f_m[1] - masse_f_m[0]


# Scelta dei valori iniziali

if (valori_iniziali=="gaussiana"):

    def param_gauss(masse, funz):


        ind_max= np.argmax(funz)
        massimo= funz[ind_max]

        meta= massimo/2

        for i in range(ind_max, 0, -1):

            if ( funz[i]<=meta ):
                break

        ind_sx= i


        for i in range(ind_max, len(funz)):

            if ( funz[i]<=meta ):
                break

        ind_dx= i

        sigma= (masse[ind_dx] - masse[ind_sx])/2

        print("La semi larghezza a metà altezza è pari a {:.3}".format(sigma))

        sigma= sigma/np.sqrt(2)

        A= np.sqrt( massimo/(np.sqrt(np.pi)*sigma) )

        mu= masse[ind_max]/2

        return A, mu, sigma



    def gauss(x, A, mu, sigma):

        return A*np.exp(-( (x-mu)/(np.sqrt(2)*sigma) )**2)




    A_gauss, mu_gauss, sigma_gauss= param_gauss(masse, F_M)

    f_m_val_iniz= gauss(masse_f_m, A_gauss, mu_gauss, sigma_gauss)


else:
    '''
    def rettangoli(m, mu_0, mu_1, sigma, ampiezza):

        if ( mu_0-sigma<m<mu_0+sigma ):
            return ampiezza

        elif( mu_1-sigma<m<mu_1+sigma ):
            return ampiezza

        else:
            return 0

    f_m_val_iniz= np.zeros(len(masse_f_m))

    for i in range(0, len(f_m_val_iniz)):

        f_m_val_iniz[i]= rettangoli(masse_f_m[i], mu_0=5, mu_1=10, sigma=1, ampiezza=10)

    f_m_val_iniz=[0.1]*len(masse_f_m)
    '''

    f_m_val_iniz= np.zeros(len(masse_f_m))

    for i in range(0, len(f_m_val_iniz)):

        f_m_val_iniz[i]= f_m_funzione(masse_f_m[i], A_0, mu_0, sigma_0, A_1, mu_1, sigma_1)






# E' la funzione da minimizzare, f_m contiene le variabili, dm_f contiene la differenza tra due elementi successivi nella lista delle masse corrispondenti a f_m, mentre F_M è la soluzione trovata precedentemente mediante la fattorizzazione SVD che deve essere confrontata con il prodoto di convoluzione di f_m per se stessa

if (funzione=="semplice"):

    def funz_minim(f_m, dm_f, F_M):

        prod= dm_f*np.convolve(f_m, f_m, mode="full")

        prod= np.append(prod, [0])

        prod_media= np.zeros(len(f_m))

        for i in range( 0, len(prod_media)):

            prod_media[i]= ( prod[2*i] + prod[2*i+1] )/2

        somma= 0

        for i in range(0, len(F_M)):

            somma+= ( prod_media[i] - F_M[i] )**2

        return somma



if (funzione=="derivata_prima"):

    def funz_minim(f_m, dm_f, F_M):

        prod= dm_f*np.convolve(f_m, f_m, mode="full")

        prod= np.append(prod, [0])

        prod_media= np.zeros(len(f_m))

        for i in range( 0, len(prod_media)):

            prod_media[i]= ( prod[2*i] + prod[2*i+1] )/2

        somma= 0

        for i in range(0, len(F_M)):

            somma+= ( prod_media[i] - F_M[i] )**2


        deriv= 0

        for i in range(1, len(f_m)):

            deriv+= (f_m[i] - f_m[i-1])**2

        return somma + cost_prima*deriv



if (funzione=="derivata_seconda"):

    def funz_minim(f_m, dm_f, F_M):

        prod= dm_f*np.convolve(f_m, f_m, mode="full")

        prod= np.append(prod, [0])

        prod_media= np.zeros(len(f_m))

        for i in range( 0, len(prod_media)):

            prod_media[i]= ( prod[2*i] + prod[2*i+1] )/2

        somma= 0

        for i in range(0, len(F_M)):

            somma+= ( prod_media[i] - F_M[i] )**2


        deriv= 0

        for i in range(1, len(f_m) -1):

            deriv+= (f_m[i+1] - 2*f_m[i] + f_m[i-1])**2

        return somma + cost_seconda*deriv










# Minimizzazione

# imposizione dei limiti per le diverse varuabili (devono essere tutte maggiori di zero)

minimi=[0]*len(f_m_val_iniz)
massimi=[np.inf]*len(f_m_val_iniz)


bounds= Bounds(minimi, massimi)

risultati= minimize(funz_minim, f_m_val_iniz, args=(dm_f, F_M), method="TNC", bounds= bounds, options={'disp': True})


'''
bounds=[]

for i in range(0, len(f_m_val_iniz)):

    a= [0, 10e6]

    bounds.append(a)


risultati= dual_annealing(funz_minim, args=(dm_f, F_M),  bounds= bounds)

'''

f_m_risult= risultati.x










# SMUSSAMENTO TRAMITE LA TECNICA DELLA MEDIA MOBILE (DIVIDO IN FINESTRE CHE CONDIVIDONO DIVERSI PUNTI)

# creo una finestra e calcolo media, la finestra successiva si crea partendo dalla precedente facendo scorrere di un valore a destra, quindi vengono condivisi un numero di punti pari alla lunghezza della finestra meno uno. La dimensione finale è pari a alla dimensione iniziale meno la lunghezza della finestra



if (opzione_smooth=="media_mobile_1"):

    masse_f_m_medie= np.zeros( len(F_M) - lung_sottoin )
    f_m_medie= np.zeros( len(F_M) - lung_sottoin )


    for i in range( int(lung_sottoin/2), len(f_m_risult) - int(lung_sottoin/2)):

        sum= 0

        for j in range( -int(lung_sottoin/2), int(lung_sottoin/2)):

            sum+= f_m_risult[i + j]

        media= sum/lung_sottoin

        masse_f_m_medie[i-int(lung_sottoin/2)]= masse_f_m[i]
        f_m_medie[i-int(lung_sottoin/2)]= media

    '''
    masse= masse_f_m_medie
    f_m_risult= f_m_medie
    '''









# GRAFICO DEI RISULTATI OTTENUTI

fig= plt.figure()


plt.plot(masse_f_m, f_m_risult, linestyle="-", marker="", color="blue", label="Soluzione Individuata")

if (opzione_smooth=="media_mobile_1"):

    stringa="Soluzione individuata con\nmedia mobile con intervalli\ndi lunghezza {0}".format(lung_sottoin)
    plt.plot(masse_f_m_medie, f_m_medie, linestyle="-", marker="", color="red", label=stringa)


if ( disegna==True ):

    distrib=np.zeros(len(masse_f_m))

    for i in range(0, len(masse_f_m)):

        distrib[i]=  f_m_funzione(masse_f_m[i], A_0, mu_0, sigma_0, A_1, mu_1, sigma_1)

    plt.plot(masse_f_m, distrib, linestyle="-", color="orange", label="Soluzione Esatta")




plt.title("FUNZIONE LOGARITMICA DI MASSA")
plt.xlabel("Massa [M_sole]")
plt.ylabel("f(m)")

plt.legend()









# CONFRONTO TRA IL PRODOTTO DI CONVOLUZIONE DI f_m_risult CON SE STESSO E F_M

dM= masse_f_m[1] - masse_f_m[0]
conv= dM*np.convolve(f_m_risult, f_m_risult, mode="full")

masse_conv= np.linspace(0, len(conv), len(conv), endpoint=False)*dM + 2*masse_f_m[0]


fig= plt.figure()

plt.plot(masse_conv, conv, linestyle="-", color="blue", label="f(m)_m*f(m)_m")
plt.plot(masse, F_M, linestyle="-", color="orange", marker="", label="F(M)")

plt.title("PRODOTTO DI CONVOLUZIONE DELLA FUNZIONE LOGARITMICA DI MASSA")
plt.xlabel("Massa [M_sole]")
plt.ylabel("f(m)*f(m)")

plt.xlim( max(masse[0], masse_conv[0]), min(masse[-1], masse_conv[-1]))











plt.legend()

plt.tight_layout()
plt.show()