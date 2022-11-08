import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import hankel1


# VARIABII PROGRAMMA


# percorso del file da cui leggere i valori di xi e dei corrisponendenti U
file_U_xi=""

# path e nome dei file in cui salvare lo spettro in potenza
path_P=""
file_name_H_11="H_11"







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

# velocità della luce nel vuoto in m/s
c= 299792458

# velocità della luce nel vuoto in UA/s
c= c/UA

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

# cost è una costante che compare nella formula per lo spettro in potenza
cost= -( 4*G**(2) )/( 1260*np.pi**(2)*c**(10) )

# nu_0 è una costante che compare nelle trasformate di Fourier delle componenti del quadrupolo

nu_0 = np.sqrt( a**(3)/(G*M) )







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
    return (   mu*a**2*(  ( (3-e**2)*np.sinh(2*xi) - 4*e*np.sinh(xi) )*xi_3(xi)
                     + 6*( (3-e**2)*np.cosh(2*xi) - 2*e*np.cosh(xi) )*xi_1(xi)*xi_2(xi)
                     + 4*( (3-e**2)*np.sinh(2*xi) - e*np.sinh(xi) )*( xi_1(xi) )**3  )    )



# M_12_3 è la derivata terza della componente 12 del quadrupolo
def M_12_3(xi):
    return (   -3*mu*a**2*np.sqrt( (e**2)-1 )*(  ( np.cosh(2*xi) - e*np.cosh(xi) )*xi_3(xi)
                                            + 3*( 2*np.sinh(2*xi) - e*np.sinh(xi) )*xi_1(xi)*xi_2(xi)
                                              + ( 4*np.cosh(2*xi) - e*np.cosh(xi) )*( xi_1(xi) )**3  )    )



# M_22_3 è la derivata terza della componente 22 del quadrupolo
def M_22_3(xi):
    return (   mu*a**2*(  ( (2*e**2-3)*np.sinh(2*xi) + 2*e*np.sinh(xi) )*xi_3(xi)
                     + 6*( (2*e**2-3)*np.cosh(2*xi) + e*np.cosh(xi) )*xi_1(xi)*xi_2(xi)
                     + 2*( 2*(2*e**2-3)*np.sinh(2*xi) + e*np.sinh(xi) )*( xi_1(xi) )**3  )    )



# M_33_3 è la derivata terza della componente 33 del quadrupolo
def M_33_3(xi):
    return (   -mu*a**2*(  ( (e**2)*np.sinh(2*xi) + 2*e*np.sinh(xi) )*xi_3(xi)
                      + 6*( (e**2)*np.cosh(2*xi) + e*np.cosh(xi) )*xi_1(xi)*xi_2(xi)
                      + 2*( 2*(e**2)*np.sinh(2*xi) + e*np.sinh(xi) )*( xi_1(xi) )**3  )    )




# TRASFORMATE DI FOURIER DELLE COMPONENTI DEL QUADRUPOLO

# hankel1_primo è funzione di Hankel primata. nu è l'argoemnto, mentre eta è l'ordine. La funzione hankel1 appartiene alla libreria scipy e fornisce i valori della funzione di Hankel del primo tipo.
def hankel1_primo(nu, eta):

    return 0.5*( hankel1(eta-1, nu) - hankel1(eta+1, nu) )



# M_11_fourier è la trasformata di Fourier della componente 11 del quadrupolo
def M_11_fourier(f):

    nu= f*nu_0
    return  (a**(2)*mu*np.pi)/(4*f)*( 16*e*hankel1_primo(1j*nu*e, 1j*nu) + (e**(2) - 3)*hankel1_primo(1j*nu*e*0.5, 1j*nu) )



# M_12_fourier è la trasformata di Fourier della componente 12 del quadrupolo
def M_12_fourier(f):

    nu= f*nu_0
    return  -(3*a**(2)*mu*np.pi)/(4*f*e)*np.sqrt(e**(2) - 1)*( 4*e*hankel1(1j*nu, 1j*nu*e) - hankel1(1j*nu, 1j*nu*e*0.5) )



# M_22_fourier è la trasformata di Fourier della componente 22 del quadrupolo
def M_22_fourier(f):

    nu= f*nu_0
    return  -(a**(2)*mu*np.pi)/(4*f)*( 8*e*hankel1_primo(1j*nu*e, 1j*nu) - (3 - 2*e**(2))*hankel1_primo(1j*nu*e*0.5, 1j*nu) )



# M_33_fourier è la trasformata di Fourier della componente 33 del quadrupolo
def M_33_fourier(f):

    nu= f*nu_0
    return  (a**(2)*mu*np.pi)/(4*f)*( 8*e*hankel1_primo(1j*nu*e, 1j*nu) + e**(2)*hankel1_primo(1j*nu*e*0.5, 1j*nu) )







# DEFINIZIONE DELL'ARRAY DELLE FREQUENZE CHE SI VOGLIONO STUDIARE

freq=







# CALCOLO DELLE TRASFORMATE DI FOURIER DEI QUADRATI DELLE DERIVATE CUBICHE DELLE COMPONENTI DEL QUADRUPOLO

# costruzione dell'array dei tempi corrispondente a freq

n= freq.size
df= freq[1] - freq[0]
time= np.fft.fftfreq(n, df)
time= np.fft.fftshift(time)

xi= inver_t_xi(time)

# trasformata del quadrato della derivata cubica della componenti 11 del quadrupolo

antitras_11= ( M_11_3(xi) )**2

# FORSE E' MEGLIO RFFT, FFT DA UN ARRAY LUNGO COME L'ARRAY DI INPUT, PER RFFT INVECE E' DIVERSO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# uso ifft per come è definita la trasformata discreta di forier in numpy

tras_11= np.fft.ifft(antitras_11)

