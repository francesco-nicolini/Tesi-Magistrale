import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import hankel1







# VARIABII PROGRAMMA


# percorso del file da cui leggere i valori di xi e dei corrisponendenti U
file_U_xi=""

# path e nome dei file in cui salvare lo spettro in potenza
path_DP=
file_name_DP=

# N è la dimenzione dell'array dei tempi. E' conveniente scegliere per questa variabile una potenza di 2, così da rendere più efficiente il calcolo della fft
N= 2048

# n è la dimenzione dell'array delle frequenze. Poichè si utilizza la funzione np.fft.rfft(), l'array in uscita da questa ha dimenzione N/2 + 1 e, per N pari, la sua prima componente corrisponde alla frequenza nulla, mentre l'ultima alla frequenza di Nyquist, pari a metà della frequenza di campionamento. Vengono inoltre calcolati solo i valori corrisponndenti a frequenze positive, in quanto, essendo presente in ingresso un'array reale, la trasformata risulta essere hermitiana.
n= N/2 + 1

# freq_min è la frequenza minima considerata
freq_min= 0

# freq_max è la frequenza massima considerata
freq_max= 10**(3)







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
A= np.sqrt(mu*a**(3)/alpha)

# cost è una costante che compare nella formula per lo spettro in potenza
cost= -( 2*G**(2) )/( 315*np.pi**(2)*c**(10) )

# nu_0 è una costante che compare nelle trasformate di Fourier delle componenti del quadrupolo
nu_0 = np.sqrt( a**(3)/(G*M) )







# DEFINIZIONE DELLE FUNZIONI DA UTILIZZARE


# FUNZIONI PER PASSARE DAL TEMPO A xi

# FUNZIONE DA INVERTIRE

def t_xi(xi, t):

    return A*(e*np.sinh(xi) - xi) - t




# INVERSIONE DELL'EQUAZIONE


def inver_t_xi(t):

    if(t < 0):
        initial_guess = -np.log(-2*t/A)

    elif (t==0):
        initial_guess = 0

    else:
        initial_guess= np.log(2*t/A)


    xi_soluz= fsolve(t_xi, initial_guess, args=(t))

    return xi_soluz[0]














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

freq= np.linspace(freq_min, freq_max, n)

df= freq[1] - freq[0]







# COSTRUZIONE DELL'ARRAY DEI TEMPI CORRISPONEDENTI ALL'ARRAY freq

# per rendere la fft più efficiente possibile conviene scegliere una potenza di 2 come dimenzione N dell'array, inoltre si deve costruire l'array delle ascisse in modo che il primo valore contenuto sia lo zero, il secondo sia pari a meno l'ultimo, il terzo sia pari a meno il penultimo e così via. Così facendo avanza un valore, ossia l' (N/2 +1)-esimo, che non è pari al negativo di nessun altro che viene solitamente posto negativo

dt= 1/(N*df)

time= np.zeros(N)

time[0]= 0

time[int(N/2)]= -(N/2)*dt

for i in range(1, int(N/2)):

    time[i]= dt*i
    time[N-i]= -t[i]







# CALCOLO DEGLI xi CORRISPONEDENTI AI TEMPI E RAPPRESENTAZIONE GRAFICA


xi= np.zeros(len(time))

for i in range(0, len(time)):

xi[i]= inver_t_xi(time[i])


# grafico di xi in funzione di t

plt.figure()

plt.title("$\\xi$ in funzione del tempo")

plt.plot(time, xi, linestyle="-", marker="", color="blue")
plt.plot([time[0], time[-1]], [0, 0], linestyle="-", marker="", color="black", linewidth= 0.8)

plt.xlabel("t [s]")
plt.ylabel("$\\xi$")

plt.show()







# CALCOLO DELLE ANTITRASORMATE


antitras_11= ( M_11_3(xi) )**2
antitras_12= ( M_12_3(xi) )**2
antitras_22= ( M_22_3(xi) )**2
antitras_33= ( M_33_3(xi) )**2







# CALCOLO DELLE TRASFORMATE DI FOURIER DEI QUADRATI DELLE DERIVATE CUBICHE DELLE COMPONENTI DEL QUADRUPOLO


# trasformata del quadrato della derivata cubica della componenti 11 del quadrupolo
tras_11= np.fft.rfft(antitras_11)

# trasformata del quadrato della derivata cubica della componenti 12 del quadrupolo
tras_12= np.fft.rfft(antitras_12)

# trasformata del quadrato della derivata cubica della componenti 22 del quadrupolo
tras_22= np.fft.rfft(antitras_22)

# trasformata del quadrato della derivata cubica della componenti 33 del quadrupolo
tras_33= np.fft.rfft(antitras_33)







# CALCOLO DELLA CORREZIONE ALLO SPETTRO IN FREQUENZA

DP= np.zeros(n, dtype=complex)

for i in range(0, n):

    DP[i]= tras_11[i]*M_11_fourier(-freq[i]) + tras_22[i]*M_22_fourier(-freq[i]) + tras_33[i]*M_33_fourier(-freq[i]) + tras_12[i]*( M_11_fourier(-freq[i]) + M_22_fourier(-freq[i]) )

    DP[i]= cost*(-1j*freq[i]**(3))*DP[i]


# converzione in Joule diviso Hertz

DP= (M_s*UA**2)*DP







# GRAFICO DELLA CORREZIONE ALLO SPETTRO IN FREQUNZA

fig=plt.figure()

plt.title("Correzione allo spettro di potenza")


plt.plot(freq, DP, marker="", linestyle="-", color="blue")
plt.plot( [freq[0], freq[-1]], [0, 0], marker="", linestyle="-", color="black", linewidth=0.8)


plt.xlim(freq[0], freq[-1])

plt.xlabel("f [Hz]")
plt.ylabel("$\\Delta$P [J/Hz]")







# SALVATAGGIO SU FILE

file_name= path_DP + "\\" + file_name_DP

file= open(file_name, "w")

np.savetxt(file, np.c_[freq, DP])
file.close()















