import numpy as np
import matplotlib.pyplot as plt
import math
import keyboard



# OPZIONI VARIE


# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)
K=50

# Estremo superiore e estremo inferiore delle masse considerate nell'integrale e scarto tra due masse consecutive
m_min=1
m_max=10**2


# Se option è pari a "read" omega_GW è letto da file altrimenti viene creato tramite una funzione
option="read"

# Path dei file se option è uguale a "read"
path_omega="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\Metodo Autovalori\\file_txt"

name_omega_tutta_minore= "omega_GW_tutta_minore.txt"
name_omega_tutta_maggiore= "omega_GW_tutta_maggiore.txt"
name_omega= "omega_GW.txt"


# Se option è pari a "read" ponendo disegna uguale a True verra realizzato il grafico della soluzione esatta nello stesso piano in cui vengono rappresentate le soluzioni trovate minimizzando
disegna=True

# Path del file se disegna è uguale a True
path_f="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\Metodo Autovalori\\file_txt"

name_f_m= "f_m.txt"


# Estremo superiore e estremo inferiore delle frequenze considerate (da selezionare solo se option è diverso da "read")
freq_min= 0.18
freq_max= 1.8













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

# poiche si seleziona solo una finestra in frequenze, gli array il cui nome termina con minore e maggiore sono utilizzati per visualizzare il grafico di omega_GW nella sua interezza

if (option=="read"):

    file_name_omega_tutta_minore= path_omega + "//" + name_omega_tutta_minore
    file_name_omega_tutta_maggiore= path_omega + "//" + name_omega_tutta_maggiore
    file_name_omega= path_omega + "//" + name_omega

    freq_tutta_minore, omega_GW_tutta_minore= np.loadtxt(file_name_omega_tutta_minore, unpack= True)
    freq_tutta_maggiore, omega_GW_tutta_maggiore= np.loadtxt(file_name_omega_tutta_maggiore, unpack=True)
    freq, omega_GW= np.loadtxt(file_name_omega, unpack=True)

    print("Il vettore contenente omega_GW ha una dimensione pari a {0}\n".format(len(omega_GW)))

    if (disegna==True):

        #  masse_graf contiene le masse per effettuare il grafico di f_esatta
        # f_esatta la forma della distribuzione letta da file (è la soluzione esatta)

        file_name_f_m= path_f + "//" + name_f_m

        masse_graf, f_esatta= np.loadtxt(file_name_f_m, unpack=True)

        '''
        print("f_esatta è lungo {0}".format(len(f_esatta)))
        '''
else:

    def funz_omeg(nu):

        n=len(nu)

        return 10**(-22)*np.ones(n)


    freq=np.logspace(np.log10(freq_min), np.log10(freq_max), K)
    omega_GW=funz_omeg(freq)


# Grafico di omega_GW

fig, ax = plt.subplots()

color_tut= "lightsteelblue"
color_wind= "midnightblue"


ax.plot(freq_tutta_minore, omega_GW_tutta_minore, linestyle="-", color=color_tut)
ax.plot(freq_tutta_maggiore, omega_GW_tutta_maggiore, linestyle="-", color=color_tut, label="Tutto il Range")
ax.plot(freq, omega_GW, linestyle="-", color=color_wind, label="Finestra di Interesse")



inf_omega= min(omega_GW_tutta_minore)

ax.plot([freq[0], freq[0]], [inf_omega, omega_GW[0]], linestyle="--", marker="", color=color_wind, alpha=1)
plt.text( freq[0]/1.2, inf_omega*10, "f= {:.2} Hz".format(freq[0]), ha="right", va="top", color=color_wind, fontsize=9)

ax.plot([freq[-1], freq[-1]], [inf_omega, omega_GW[-1]], linestyle="--", marker="", color=color_wind, alpha=1)
plt.text( freq[-1]*1.2, inf_omega*10, "f= {:.2} Hz".format(freq[-1]), ha="left", va="top", color=color_wind, fontsize=9)




plt.title("$\\Omega_{GW}$ in funzione della frequenza")
plt.xlabel("f [Hz]")
plt.ylabel("$\\Omega_{GW}$")
plt.yscale("log")
plt.xscale("log")
plt.xlim(freq_tutta_min, freq_tutta_max)
plt.ylim(inf_omega, 10**(-15))

plt.legend()

plt.show()


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










# CREAZIONE DELLE LISTE CONTENENTI LE MASSE, I VALORI DI ALPHA OMEGA DIVISO PER LA FREQUENZA LA QUADRATO

# Lista dei valori della massa totale utilizzati per effettuare la discretizzazione dell'integrale (sono in masse solari)

masse= np.logspace( np.log10(m_min), np.log10(m_max), K)

alpha= 1/freq**2

alpha= alpha[::-1]

rap_m= masse[1]/masse[0]
rap_alpha= alpha[1]/alpha[0]


print("\nIl rapporto tra elementi successivi della seguenza delle masse è pari a {:.3}, mentre il rapporto tra l'inverso dei quadrati di elementi successivi nella seguenza delle frequenze (alpha) è pari a {:.3}".format( rap_m, rap_alpha))


omega_risc= omega_GW/freq**2
omega_risc= omega_risc[::-1]









# DEFINIZIONE DELLA FUNZIONE CHE COMPARE NELLA RELAZIONE DI RIFERIMENTO UNA VOLTA CHE E' STATO DISCRETIZZATO L'INTEGRALE



# La funzione usata per costruire la matrice, M è la massa totale del sistema, mentre alpha è l'inerso al quadrato della frequenza

def integ(M, alpha):

    nu_0= np.sqrt(a**3/(G*(M)))
    x_0= 2*math.pi*nu_0*1/np.sqrt(alpha)

    z= (1 - y**2 + 4*y**4 + 1.5*(x_0*y**6)/(xi))/( np.exp(2*x_0*xi)*(1 + y**2)**2)

    return cost*z


# CREAZIONE DELLA MATRICE E TEST DI SIMMETRIA


# Creazione della matrice connessa alla tecnica di discretizzazione usata. Nel caso della regola del trapezio è una matrice diagonale le cui componenti sono pari a 1/2, tuttavia al fine di rendere la matrice finale simmetrica, si crea una matrice diagonale D le cui componenti sono pari a 1/sqrt(2). In questo modo il problema diviene (D*matrix*D)*(D*F)= (D*omega_riscalato), dove matrix è la matrice che si ottiene discretizzando la funzione integranda che dipemde sia dalla frequenza che della massa e che risulta simmetrica, mentre F è il vettore incognito

a= np.sqrt(1/2)*np.ones(K)

D= np.diag(a)






# Creazione della matrice della funzione integranda nota e test di simmetria

print("\n\n\n\nCostruzione matrice")

print("\n\nMatrice di cui sono determinate tutte le componenti")

matrix= np.zeros( (K,K) )


for i in range(0,K):

    for j in range(0,K):

        matrix[i][j]= integ( masse[i], alpha[j])





# test di simmetria

trasposta= matrix.transpose()


if ( (matrix==trasposta).all() ):
    print("\nLa matrice è simmetrica")

else:
    print("\nLa matrice non è simmetrica")


massimo= 0

for i in range(0, K):

    for j in range(0, K):

        diff= abs( (matrix[i][j] - trasposta[i][j])/matrix[i][j] )

        if ( diff>massimo ):
            massimo= diff

print("Il massimo della differenza relativa in valore assoluto tra la matrice e la sua trasposta è {:.2e}".format(massimo))






print("\n\nMatrice che costruisco simmetrica:")

matrix= np.zeros( (K,K) )

for i in range(0,K):

    matrix[i][i]= integ( masse[i], alpha[i])

    for j in range(i+1,K):

        matrix[i][j]= integ( masse[i], alpha[j])

        matrix[j][i]= matrix[i][j]


trasp= matrix.transpose()



if ( (matrix==trasp).all() ):
    print("\nLa matrice è simmetrica")

else:
    print("\nLa matrice non è simmetrica")





# COSTRUZIONE DELLA MATRICE FINALE E DELL'ARRAY DEI TERMINI NOTI FINALE

mat_fin= np.dot( D, matrix)
mat_fin= np.dot( mat_fin, D)


omega_fin= np.dot( D, omega_risc)


# N.B. l'arrai che si trova non è F ma D*F




