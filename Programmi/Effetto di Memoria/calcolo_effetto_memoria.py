import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


# VARIABII PROGRAMMA

# nella variabile componenti indicare quali componenti si vogliono calcolare. Le possibilita sono "H_11", "H_11", "H_11", "H_11". Se si vogliono indicare più componenti queste vanno separate da una virgola, ad esempio se vi vuole conoscere H_11 e H_22 basta scrivere componenti="H_11, H_22".
componenti="H_11, H_12, H_22, H_33"

# percorso del file da cui leggere i valori di xi e dei corrisponendenti U
file_U_xi=""

# path e nome dei file in cui salvare le componenti di H_ij calcolate
path_H=""
file_name_H_11="H_11"
file_name_H_12="H_12"
file_name_H_22="H_22"
file_name_H_33="H_33"





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


# velocità della luce nel vuoto in m/s
c= 299792458

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

U, xi= np.loadtxt(file_U_xi, unpack=True)





# CALCOLO DEGLI INTEGRALI

# gli integrali da calolare sono combinazioni di sei integrali di riferimento, quello di (M_11^(3))**2, quello di (M_12^(3))**2, quello di M_11^(3)*M_12^(3), quello di M_12^(3)*M_22^(3), quello di (M_22^(3))**2 e infine quello di (M_33^(3))**2



# CALCOLO DELL'INTEGRALE DI (M_11^(3))**2

if ("H_11" in componenti):

    # definizione funzione integranda

    def integranda_11_11(xi):

        return ( M_11_3(xi) )**2/( xi_1(xi) )


    # calcolo dell'integrale

    integrale_11_11= np.zeros(len(U))

    for i in range(0, len(U)):
        integrale_11_11[i]= quad(integranda_11_11, -np.inf, xi[i] )








# CALCOLO DELL'INTEGRALE DI (M_12^(3))**2

if ( ("H_11" in componenti) or ("H_22" in componenti) ):

    # definizione funzione integranda

    def integranda_12_12(xi):

        return ( -2*G/(7*c**(5)) )*( M_12_3(xi) )**2/( xi_1(xi) )


    # calcolo dell'integrale

    integrale_12_12= np.zeros(len(U))

    for i in range(0, len(U)):
        integrale_12_12[i]= quad(integranda_12_12, -np.inf, xi[i] )








# CALCOLO DELL'INTEGRALE DI M_11^(3)*M_12^(3)

if ("H_12" in componenti):

    # definizione funzione integranda

    def integranda_11_12(xi):

        return ( -2*G/(7*c**(5)) )*( M_11_3(xi)*M_12_3(xi) )/( xi_1(xi) )


    # calcolo dell'integrale

    integrale_11_12= np.zeros(len(U))

    for i in range(0, len(U)):
        integrale_11_12[i]= quad(integranda_11_12, -np.inf, xi[i] )








# CALCOLO DELL'INTEGRALE DI M_12^(3)*M_22^(3)

if ("H_12" in componenti):

    # definizione funzione integranda

    def integranda_12_22(xi):

        return ( -2*G/(7*c**(5)) )*( M_12_3(xi)*M_22_3(xi) )/( xi_1(xi) )


    # calcolo dell'integrale

    integrale_12_22= np.zeros(len(U))

    for i in range(0, len(U)):
        integrale_12_22[i]= quad(integranda_12_22, -np.inf, xi[i] )








# CALCOLO DELL'INTEGRALE DI (M_22^(3))**2

if ("H_22" in componenti):

    # definizione funzione integranda

    def integranda_22_22(xi):

        return ( -2*G/(7*c**(5)) )*( M_22_3(xi) )**2/( xi_1(xi) )


    # calcolo dell'integrale

    integrale_22_22= np.zeros(len(U))

    for i in range(0, len(U)):
        integrale_22_22[i]= quad(integranda_22_22, -np.inf, xi[i] )








# CALCOLO DELL'INTEGRALE DI (M_33^(3))**2

if ("H_33" in componenti):

    # definizione funzione integranda

    def integranda_33_33(xi):

        return ( -2*G/(7*c**(5)) )*( M_33_3(xi) )**2/( xi_1(xi) )


    # calcolo dell'integrale

    integrale_33_33= np.zeros(len(U))

    for i in range(0, len(U)):
        integrale_33_33[i]= quad(integranda_33_33, -np.inf, xi[i] )








# DETERMINAZIONE DELLE COMPONENTI DI H

if ("H_11" in componenti):
    H_11= integrale_11_11 + integrale_12_12

if ("H_12" in componenti):
    H_12= integrale_11_12 + integrale_12_22

if ("H_22" in componenti):
    H_22= integrale_12_12 + integrale_22_22

if ("H_33" in componenti):
    H_33= integrale_33_33








# GRAFICO DELLE COMPOENTI DI H_ij IN FUNZIONE DI U

plt.figure()

plt.title("Componenti di $H_{ij}$ in funzione del termpo ritardato U")

if ("H_11" in componenti):
    plt.plot(U, H_11, linestyle="-", marker="", color="C0", label="$H_{11}$")

if ("H_12" in componenti):
    plt.plot(U, H_12, linestyle="-", marker="", color="C1", label="$H_{12}$")

if ("H_22" in componenti):
    plt.plot(U, H_22, linestyle="-", marker="", color="C2", label="$H_{22}$")

if ("H_33" in componenti):
    plt.plot(U, H_33, linestyle="-", marker="", color="C3", label="$H_{33}$")

plt.xlabel("U [s]")
plt.ylabel("$H_{ij}$")

plt.xlim(U[0], U[-1])

plt.legend()
plt.show()








# SALVATAGGIO DELLE COMPONENENTI DI H_ij SU FILE

if ("H_11" in componenti):
    file_name_H_11= path_H + "\\" + file_name_H_11
    file_H_11= open(file_name_H_11, "w")
    np.savetxt(file_H_11, np.c_[U, H_11])
    file_H_11.close()

if ("H_12" in componenti):
    file_name_H_12= path_H + "\\" + file_name_H_12
    file_H_12= open(file_name_H_12, "w")
    np.savetxt(file_H_12, np.c_[U, H_12])
    file_H_12.close()

if ("H_22" in componenti):
    file_name_H_22= path_H + "\\" + file_name_H_22
    file_H_22= open(file_name_H_22, "w")
    np.savetxt(file_H_22, np.c_[U, H_22])
    file_H_22.close()

if ("H_33" in componenti):
    file_name_H_33= path_H + "\\" + file_name_H_33
    file_H_33= open(file_name_H_33, "w")
    np.savetxt(file_H_33, np.c_[U, H_33])
    file_H_33.close()


