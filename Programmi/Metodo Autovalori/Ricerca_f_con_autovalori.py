import numpy as np
from numpy.linalg import eigh
import matplotlib.pyplot as plt
import math
import keyboard



# OPZIONI VARIE


# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)
K=50

# num_autoval contiene il numero di autostati considerati per approssimare la funzione omega/f**2
num_autoval= 50

# Estremo superiore e estremo inferiore delle masse considerate nell'integrale e scarto tra due masse consecutive
m_min=10**(0)
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


# Il programma indica il numero e in caso di necessita la lista dei moduli di un prodotto scalare di un autovettore per un autovettore con un valore superiore a quello contenuto nrlla variabile soglia
soglia= 10**(-10)









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
plt.xlim(min(freq_tutta_minore), max(freq_tutta_maggiore))
plt.ylim(inf_omega, 10**(-15))

plt.legend()




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










# CREAZIONE DELLE LISTE CONTENENTI LE MASSE, I VALORI DI ALPHA E OMEGA DIVISO PER LA FREQUENZA LA QUADRATO

# Lista dei valori della massa totale utilizzati per effettuare la discretizzazione dell'integrale (sono in masse solari)

masse= np.linspace( m_min, m_max, endpoint=False, num=K)

alpha= 1/freq**2

alpha= alpha[::-1]

rap_m= masse[1]/masse[0]
rap_alpha= alpha[1]/alpha[0]


print("\nIl rapporto tra elementi successivi della seguenza delle masse è pari a {:.3}, mentre il rapporto tra l'inverso dei quadrati di elementi successivi nella seguenza delle frequenze (alpha) è pari a {:.3}".format( rap_m, rap_alpha))


omega_risc= omega_GW/freq**2
omega_risc= omega_risc[::-1]


print("alpha:", alpha)

print("\nmasse:", masse)







# DEFINIZIONE DELLA FUNZIONE CHE COMPARE NELLA RELAZIONE DI RIFERIMENTO UNA VOLTA CHE E' STATO DISCRETIZZATO L'INTEGRALE



# La funzione usata per costruire la matrice, M è la massa totale del sistema, mentre alpha è l'inerso al quadrato della frequenza

def integ(M, alpha):

    nu_0= np.sqrt(a**3/(G*(M)))
    x_0= 2*math.pi*nu_0*1/np.sqrt(alpha)

    z= (1 - y**2 + 4*y**4 + 1.5*(x_0*y**6)/(xi))/( np.exp(2*x_0*xi)*(1 + y**2)**2)

    return z








# CREAZIONE DELLA MATRICE E TEST DI SIMMETRIA


# Creazione della matrice connessa alla tecnica di discretizzazione usata. Nel caso della regola del trapezio è una matrice diagonale le cui componenti sono pari a 1/2, tuttavia al fine di rendere la matrice finale simmetrica, si crea una matrice diagonale D le cui componenti sono pari a 1/sqrt(2). In questo modo il problema diviene (D*matrix*D)*(D*F)= (D*omega_riscalato), dove matrix è la matrice che si ottiene discretizzando la funzione integranda che dipemde sia dalla frequenza che della massa e che risulta simmetrica, mentre F è il vettore incognito

a= np.sqrt(1/2)*np.ones(K)

D= np.diag(a)






# Creazione della matrice della funzione integranda nota e test di simmetria

print("\n\n\n\nCostruzione matrice")

print("\n\nMatrice di cui sono determinate tutte le componenti")

matrix= np.zeros( (K,K) )


m= masse[0]
a= alpha[0]

#print( integ(m, a))


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


omega_fin= np.dot( D, omega_risc)/cost


# N.B. l'arrai che si trova non è F ma D*F





# DETERMINAZIONE DEGLI AUTOVALORI E DEGLI AUTOVETTORI

lamb, v= eigh(mat_fin)



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

    autovett[i][1]= autovett[i][1]/np.sqrt((np.dot(autovett[i][1], autovett[i][1])))


print("\n\n\nLista degli autovalori individuati ordinati in senso decrescente del loro modulo:\n")


for i in range(0, len(lamb)):

    print(autovett[i][0])






# RICERCA DEI COEFFICIENTI DELLA COMBINAZIONE LINEARE DI AUTOSTATI DELLA MATRICE CHE MEGLIO APPROSSIMA OMEGA_GW E REALIZZAZIONE DEL GRAFICO DELLA FUNZIONE STIMATA




coeffic= np.zeros( num_autoval)

for i in range(0, num_autoval):

    coeffic[i]= np.dot(autovett[i][1], omega_fin)/np.dot(autovett[i][1], autovett[i][1])



omega_stima= np.zeros(len(omega_fin))

for i in range(0, num_autoval):

    omega_stima+= coeffic[i]*autovett[i][1]



fig, ax = plt.subplots(2)

ax[0].plot(alpha, omega_fin, linestyle=(0, (1, 1)), color="blue", label="Soluzione Esatta")
ax[0].plot(alpha, omega_stima, linestyle="--", color="red", label="Approssimazione trovata\ncon {0} autostati".format(num_autoval))


ax[0].set_title("$\\Omega_{GW}$ in funzione della frequenza", fontsize=14)
ax[0].set_xlabel("$\\alpha$ [1/Hz**2]", fontsize=10)
ax[0].set_ylabel("$\\Omega_{GW}$", fontsize=10)
ax[0].set_xlim(min(alpha), max(alpha))
#ax[0].set_ylim(10**(-100), 10**(-8))
ax[0].set_xscale("log")
ax[0].set_yscale("log")

ax[0].legend()


scarto= (omega_fin-omega_stima)/omega_fin

ax[1].plot(alpha, abs(scarto), linestyle="-", color="blue")

'''
ax[1].plot( freq[i_max], dif_max, marker=".", color="blue")
ax[1].text( freq[i_max], dif_max, "{:.2e}".format(dif_max), horizontalalignment="right")
'''

ax[1].set_title("Modulo della Differenza Relativa tra le $\\Omega_{GW}$ ", fontsize=14)
ax[1].set_xlabel("$\\alpha$ [1/Hz**2]", fontsize=10)
ax[1].set_ylabel("$|\\Delta(\\Omega_{GW})/\\Omega_{GW}|$", fontsize=10)
ax[1].set_xlim(min(alpha), max(alpha))
#ax[1].set_ylim(10**(-13), 10**(-5))
ax[1].set_xscale("log")
ax[1].set_yscale("log")




plt.tight_layout()









'''
print("\n\nComponenti di omega stimata:\n")

for i in range(0, len(omega_stima)):
    print(omega_stima[i])
'''

print("\n\nCoefficienti:\n")

for i in range(0, len(coeffic)):
    print(coeffic[i])

'''
print("\n\nautostati:\n")

for i in range(0, len(coeffic)):

    print("\nAutovettore numero {0}:\n".format(i))

    for j in range(0, len(autovett[i][1])):

        print(autovett[i][1][j])
'''





# TEST SUGLI AUTOVETTORI INDIVIDUATI



massimo=0

for i in range(0, len(autovett)):

    coeff= np.dot(autovett[i][1], autovett[i][1])

    if ( massimo<coeff ):
        massimo=coeff
        i_massimo= i

print("\n\nIl modulo massimo è {:.2} e si ha per l'autostato {:}".format(massimo, i_massimo))



massimo= 0

for i in range(0, len(autovett)):

    for j in range(i+1, len(autovett)):

        coeff_ort= np.dot(autovett[i][1], autovett[j][1])

        if ( massimo<abs(coeff_ort) ):
            massimo= coeff_ort
            i_massimo= i
            j_massimo= j


print("\n\n\nMassimo del valore assoluto del prodotto scalare di un'autovettore per un diverso autovettore: {:.2}\nTale risultato si ottiene per la coppia ({:},{:})".format(massimo, i_massimo, j_massimo ) )




sopra_sogl=[]
i_sup_sogl=[]
j_sup_sogl=[]


for i in range(0, len(autovett)):

    for j in range(i+1, len(autovett)):

        coeff_ort= np.dot(autovett[i][1], autovett[j][1])

        if ( abs(coeff_ort)>soglia ):

            sopra_sogl.append(coeff_ort)
            i_sup_sogl.append(i)
            j_sup_sogl.append(j)


print("\n\n\nCi sono {:} valori assoluti di prodotti scalari di autovettori diversi tra loro che sono superiori alla soglia di {:.2} :\n".format( len(sopra_sogl), soglia))

'''
print("\n\n\nLista dei valori assoluti di prodotti scalari di autovettori diversi tra loro che sono superiori alla soglia di {:.2} e corrispettive coppie:\n".format(soglia))


for i in range(0, len(sopra_sogl)):

    print("{:.2} ({:},{:})".format(sopra_sogl[i], i_sup_sogl[i], j_sup_sogl[i]))
'''




# CALCOLO F(M)

F_M= np.zeros(K)

for i in range(0, num_autoval):

    F_M+= (coeffic[i]/autovett[i][0])*autovett[i][1]


fig, ax= plt.subplots()

plt.plot(masse, F_M)
plt.xlabel("M [M_sun]")
plt.ylabel("F(M)")





plt.show()





