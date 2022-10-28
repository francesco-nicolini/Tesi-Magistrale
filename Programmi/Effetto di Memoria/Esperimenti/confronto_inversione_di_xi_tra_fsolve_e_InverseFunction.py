import numpy as np
from scipy.optimize import fsolve
from pynverse import inversefunc



# PARAMETRI DA FISSARE


# e è l'eccentricità
e= 1.005

# a è il semi asse maggiore in U.A.
a= 1

# m_1 è la massa del primo ogetto in masse solari
m_1= 1

# m_2 è la massa del secondo ogetto in masse solari
m_2= 1




# COSTANTI VARIE


# massa del sole in kg
M_s=1.989*10**(30)

# unità astronomica in m
UA=1.496*10**(11)

# costante di gravitazione universale in Nm**2/kg**2
G= 6.7*10**(-11)

# costante di gravitazione universale in U.A.**3/(M_S*s**2)
G= (M_s/UA**3)*G

# M è la massa totale del sistema in masse solari
M=m_1+m_2

# mu è la massa ridotta del sistema in masse solari
mu= m_1*m_2/M


# alpha è il prodotto tra la costante di gravitazione aniversale G, la massa totale M, e la massa ridotta mu
alpha= G*M*mu

# A è una costante che compare nelle componenti del quadrupolo e nelle unita di misura scelta si misura in secondi
A=np.sqrt(mu*a**(3)/alpha)




# VALORE DI u DI CUI CALCOLARE IL CORRISPONDENTE xi

u= -1e1


# METODO BASATO SU fsolve

print("Metodo basato su fsolve:")

# FUNZIONE DA INVERTIRE

def u_xi(xi, u):

    return u - A*(e*np.sinh(xi) - xi)



# INVERSIONE DELL'EQUAZIONE


if(u < 0):
    initial_guess = -np.log(-2*u/A)

elif (u==0):
    initial_guess = 0

else:
    initial_guess= np.log(2*u/A)


xi_soluz_fs= fsolve(u_xi, initial_guess, args=(u))


print("\nLa soluzione ottenuta con la funzione fsolve è {:.3}".format(xi_soluz_fs[0]))


# serve per verificare che si è ottenuto un risultato corretto
print("\nLa differenza tra il valore del tempo proprio u={:.2} e il valore della funzione\ncon la soluzione trovata è pari a {:.2}".format( u, u_xi(xi_soluz_fs[0], u) ) )





#  METODO BASATO SU InverseFunction

print("\n\n\nMetodo basato su InverseFunction:")

# FUNZIONE DA INVERTIRE

def funz(xi):

    return A*(e*np.sinh(xi) - xi)


# INVERSIONE DELL'EQUAZIONE


xi_soluz_if= inversefunc(funz, y_values=u)
print("\nLa soluzione ottenuta con la funzione InverseFunction è {:.3}".format(xi_soluz_if))


print("\nLa differenza tra il valore del tempo proprio u= {:.2} e il valore della funzione\ncon la soluzione trovata è pari a {:.2}".format( u, u - funz(xi_soluz_if) ) )

