import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve



# VARIABILI PROGRAMMA


window= 1e20
N= 1e4
path=
file_name=




# PARAMETRI FISICI


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

# A è una costante che compare nelle componenti del quadrupolo
A=np.sqrt(mu*a**(3)/alpha)




# FUNZIONE DA INVERTIRE

def u_xi(xi, u):

    return A*(e*np.sinh(xi) - xi) - u




# INVERSIONE DELL'EQUAZIONE


def invert(u):

    if(u < 0):
        initial_guess = -np.log(-2*u/A)

    elif (u==0):
        initial_guess = 0

    else:
        initial_guess= np.log(2*u/A)


    xi_soluz= fsolve(u_xi, initial_guess, args=(u))

    return xi_soluz[0]





# DETERMINAZIONE DELLE xi

U= np.linspace(-window, window, int(N))
xi= np.zeros(len(U))


for i in range(0, len(U)):
    xi[i]= invert(U[i])





# GRAFICO DI xi IN FUNZIONE DI u

plt.figure()

plt.title("$\\xi$ in funzione di U")

plt.plot(U, xi, linestyle="-", marker="", color="blue")

plt.xlabel("U [s]")
plt.ylabel("$\\xi$")

plt.show()





# STAMPA SU FILE

file_name= path + "\\" + file_name

file= open(file_name, "w")

np.savetxt(file, np.c_[U, xi])
file.close()

