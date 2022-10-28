import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad



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

# costante di gravitazione universale in U.A.**3/(M_S*s**3)
G= (M_s/UA**3)*G

# M è la massa totale del sistema in masse solari
M=m_1+m_2

# mu è la massa ridotta del sistema in masse solari
mu= m_1*m_2/M


# alpha è il prodotto tra la costante di gravitazione aniversale G, la massa totale M, e la massa ridotta mu
alpha= G*M*mu

# A è una costante che compare nelle componenti del quadrupolo e nelle unita di misura scelta si misura in secondi
A=np.sqrt(mu*a**(3)/alpha)




# DEFINIZIONE DELLE FUNZIONI DA UTILIZZARE



# DERIVATE DI xi

# xi_1 e la derivata prima di xi
def xi_1(xi):

    return 1/(A*(e*np.cosh(xi) - 1))


# xi_2 e la derivata seconda di xi
def xi_2(xi):

    return -(e*np.sinh(xi))/(A**(2)*(e*np.cosh(xi) - 1)**(3))


# xi_3 e la derivata terza di xi
def xi_3(xi):

    return e*(2*e*np.sinh(xi)**2 + np.cosh(xi) - e)/(A**(3)*(e*np.cosh(xi) - 1)**(5))



# DERIVATE CUBICHE DELLE COMPONENTI DEL QUADRUPOLO

# M_11_3 è la derivata terza della componente 11 del quadrupolo
def M_11_3(xi):
    return (   m*a**2*(  ( (3-e**2)*np.sinh(2*xi) - 4*e*np.sinh(xi) )*xi_3(xi)
                     + 6*( (3-e**2)*np.cosh(2*xi) - 2*e*np.cosh(xi) )*xi_1(xi)*xi_2(xi)
                     + 4*( (3-e**2)*np.sinh(2*xi) - e*np.sinh(xi) )*( xi_1(xi) )**3  )    )


# M_12_3 è la derivata terza della componente 12 del quadrupolo
def M_12_3(xi):
    return (   -3*m*a**2*np.sqrt( (e**2)-1 )*(  ( np.cosh(2*xi) - e*np.cosh(xi) )*xi_3(xi)
                                            + 3*( 2*np.sinh(2*xi) - e*np.sinh(xi) )*xi_1(xi)*xi_2(xi)
                                              + ( 4*np.cosh(2*xi) - e*np.cosh(xi) )*( xi_1(xi) )**3  )    )


# M_22_3 è la derivata terza della componente 22 del quadrupolo
def M_22_3(xi):
    return (   m*a**2*(  ( (2*e**2-3)*np.sinh(2*xi) + 2*e*np.sinh(xi) )*xi_3(xi)
                     + 6*( (2*e**2-3)*np.cosh(2*xi) + e*np.cosh(xi) )*xi_1(xi)*xi_2(xi)
                     + 2*( 2*(2*e**2-3)*np.sinh(2*xi) + e*np.sinh(xi) )*( xi_1(xi) )**3  )    )


# M_33_3 è la derivata terza della componente 33 del quadrupolo
def M_33_3(xi):
    return (   -m*a**2*(  ( (e**2)*np.sinh(2*xi) + 2*e*np.sinh(xi) )*xi_3(xi)
                      + 6*( (e**2)*np.cosh(2*xi) + e*np.cosh(xi) )*xi_1(xi)*xi_2(xi)
                      + 2*( 2*(e**2)*np.sinh(2*xi) + e*np.sinh(xi) )*( xi_1(xi) )**3  )    )





# LETURA DA FILE DI U E xi

# l'integrale da calcolare dipende da xi solo per quanto riguarda il suo estremo superiore. Quest, corrispondenti ai valori del tempo ritardato che si vuole esplorare, sono letti da un file esterno

U, xi= np.




# CALCOLO DEGLI INTEGRALI

# gli integrali da calolare sono combinazioni di sei integrali di riferimento, quello di (M_11^(3))**2, quello di (M_12^(3))**2, quello di M_11^(3)*M_12^(3), quello di M_12^(3)*M_22^(3), quello di (M_22^(3))**2 e infine quello di (M_33^(3))**2



# CALCOLO DELL'INTEGRALE DI (M_11^(3))**2

# definizione funzione integranda

def integr_11_quad(xi):

    return ( M_11_3(xi) )**2/( xi_1(xi) )






# CALCOLO DELL'INTEGRALE DI (M_12^(3))**2

# definizione funzione integranda

def integr_12_quad(xi):

    return ( M_12_3(xi) )**2/( xi_1(xi) )






# CALCOLO DELL'INTEGRALE DI M_11^(3)*M_12^(3)

# definizione funzione integranda

def integr_11_12(xi):

    return ( M_11_3(xi)*M_12_3(xi) )/( xi_1(xi) )






# CALCOLO DELL'INTEGRALE DI M_12^(3)*M_22^(3)

# definizione funzione integranda

def integr_12_22(xi):

    return ( M_12_3(xi)*M_22_3(xi) )/( xi_1(xi) )






# CALCOLO DELL'INTEGRALE DI (M_22^(3))**2

# definizione funzione integranda

def integr_22_quad(xi):

    return ( M_22_3(xi) )**2/( xi_1(xi) )






# CALCOLO DELL'INTEGRALE DI (M_33^(3))**2

# definizione funzione integranda

def integr_33_quad(xi):

    return ( M_33_3(xi) )**2/( xi_1(xi) )