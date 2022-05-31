import numpy as np
import matplotlib.pyplot as plt
import math
import keyboard



# OPZIONI VARIE

# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle componenti della matrice)
K=50

# Estremo superiore e estremo inferiore delle masse considerate nello studio della fuzione integranda

m_min=1
m_max=10**(2)

# Estremo superiore e estremo inferiore delle frequenze considerate (da selezionare solo se option è diverso da "read")
freq_min=0.18
freq_max=1.8











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






# COSTRUZIONE DEGLI ARRAY CONTENENTI LE MASSE E LE FREQUENZE

m_min= np.log10(m_min)
m_max= np.log10(m_max)

masse= np.logspace(m_min, m_max, K)



alpha_min= 1/(freq_max)**2
alpha_max= 1/(freq_min)**2

print("alpha varia nel range [ {:.4} - {:.4} ]".format(alpha_min, alpha_max))

alpha_min= np.log10(alpha_min)
alpha_max= np.log10(alpha_max)

alpha= np.logspace(alpha_min, alpha_max, K)






# COSTRUZIONE DELLA MATRICE E TEST DI SEMMETRICITA'


# La funzione usata per costruire la matrice, M è la massa totale del sistema, mentre alpha è definito come l'inverso del quadrato della frequenza

def integ(M, alpha):

    nu_0= np.sqrt(a**3/(G*(M)))
    x_0= 2*math.pi*nu_0*(1/np.sqrt(alpha))

    z= (1 - y**2 + 4*y**4 + 1.5*(x_0*y**6)/(xi))/( np.exp(2*x_0*xi)*(1 + y**2)**2)

    return cost*z



print("\n\n\nMatrice di cui sono determinate tutte le componenti")


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






print("\n\n\n\nMatrice che costruisco simmetrica:")

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




















