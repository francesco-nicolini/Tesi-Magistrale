import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import hankel1
from scipy.optimize import fsolve
import pandas as pd




# VARIABII PROGRAMMA


# path e nome dei file in cui salvare lo spettro in potenza e del file in cui salvare la frequenze e i valori delle due masse
path= "C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\Effetto di Memoria\\Calcolo_diretto_dello_spettro_in_potenza\\Risultati"

file_name_DP= "DP_in_funzione_f_m_1_m_2"

file_name_f= "valori_di_f"
file_name_m_1= "valori_di_m_1"
file_name_m_2= "valori_di_m_2"


# N è la dimenzione dell'array dei tempi. E' conveniente scegliere per questa variabile una potenza di 2, così da rendere più efficiente il calcolo della fft
N= 2048

# n è la dimenzione dell'array delle frequenze. Poichè si utilizza la funzione np.fft.rfft(), l'array in uscita da questa ha dimenzione N/2 + 1 e, per N pari, la sua prima componente corrisponde alla frequenza nulla, mentre l'ultima alla frequenza di Nyquist, pari a metà della frequenza di campionamento. Vengono inoltre calcolati solo i valori corrisponndenti a frequenze positive, in quanto, essendo presente in ingresso un'array reale, la trasformata risulta essere hermitiana.
n= int( N/2 + 1 )

# freq_min è la frequenza minima considerata
freq_min= 0

# freq_max è la frequenza massima considerata
freq_max= 20







# PARAMETRI FISICI


# e è l'eccentricità
e= 1.05

# distanza minima di avvicnamento in U.A.
r_min= 0.05

# a è il semi asse maggiore in U.A.
a= r_min/(e-1)

# m_1_min e m_1_max sono il valore minimo e il valore massimo per la massa del primo ogetto in masse solari,n_m_1 è la dimenzione dell'array delle masse del primo ogetto
m_1_min= 1
m_1_max= 10
n_m_1= 100

# m_2_min e m_2_max sono il valore minimo e il valore massimo per la massa del secondo ogetto in masse solari,,n_m_2 è la dimenzione dell'array delle masse del secondo ogetto
m_2_min= 1
m_2_max= 10
n_m_2= 100







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

# cost è una costante che compare nella formula per lo spettro in potenza
cost= -( 4*G**(2) )/( 315*c**(10) )

'''
# nu_0 è una costante che compare nelle trasformate di Fourier delle componenti del quadrupolo
nu_0 = np.sqrt( a**(3)/(G*M) )
'''





# DEFINIZIONE DELLE FUNZIONI DA UTILIZZARE


# FUNZIONI PER PASSARE DAL TEMPO A xi

# FUNZIONE DA INVERTIRE

def t_xi(xi, t, A):

    return A*(e*np.sinh(xi) - xi) - t




# INVERSIONE DELL'EQUAZIONE


def inver_t_xi(t, A):

    if(t < 0):
        initial_guess = -np.log(-2*t/A)

    elif (t==0):
        initial_guess = 0

    else:
        initial_guess= np.log(2*t/A)


    xi_soluz= fsolve(t_xi, initial_guess, args=(t, A))

    return xi_soluz[0]







# DERIVATE DI xi

# xi_1 e la derivata prima di xi
def xi_1(xi, A):

    return 1/(A*(e*np.cosh(xi) - 1))


# xi_2 e la derivata seconda di xi
def xi_2(xi, A):

    return -(e*np.sinh(xi))/(A**(2)*(e*np.cosh(xi) - 1)**(3))


# xi_3 e la derivata terza di xi
def xi_3(xi, A):

    return e*(2*e*np.sinh(xi)**2 + np.cosh(xi) - e)/(A**(3)*(e*np.cosh(xi) - 1)**(5))




# DERIVATE CUBICHE DELLE COMPONENTI DEL QUADRUPOLO

# M_11_3 è la derivata terza della componente 11 del quadrupolo
def M_11_3(xi, mu, A):
    return (   mu*a**2*(  ( (3-e**2)*np.sinh(2*xi) - 4*e*np.sinh(xi) )*xi_3(xi, A)
                     + 6*( (3-e**2)*np.cosh(2*xi) - 2*e*np.cosh(xi) )*xi_1(xi, A)*xi_2(xi, A)
                     + 4*( (3-e**2)*np.sinh(2*xi) - e*np.sinh(xi) )*( xi_1(xi, A) )**3  )    )



# M_12_3 è la derivata terza della componente 12 del quadrupolo
def M_12_3(xi, mu, A):
    return (   -3*mu*a**2*np.sqrt( (e**2)-1 )*(  ( np.cosh(2*xi) - e*np.cosh(xi) )*xi_3(xi, A)
                                            + 3*( 2*np.sinh(2*xi) - e*np.sinh(xi) )*xi_1(xi, A)*xi_2(xi, A)
                                              + ( 4*np.cosh(2*xi) - e*np.cosh(xi) )*( xi_1(xi, A) )**3  )    )



# M_22_3 è la derivata terza della componente 22 del quadrupolo
def M_22_3(xi, mu, A):
    return (   mu*a**2*(  ( (2*e**2-3)*np.sinh(2*xi) + 2*e*np.sinh(xi) )*xi_3(xi, A)
                     + 6*( (2*e**2-3)*np.cosh(2*xi) + e*np.cosh(xi) )*xi_1(xi, A)*xi_2(xi, A)
                     + 2*( 2*(2*e**2-3)*np.sinh(2*xi) + e*np.sinh(xi) )*( xi_1(xi, A) )**3  )    )



# M_33_3 è la derivata terza della componente 33 del quadrupolo
def M_33_3(xi, mu, A):
    return (   -mu*a**2*(  ( (e**2)*np.sinh(2*xi) + 2*e*np.sinh(xi) )*xi_3(xi, A)
                      + 6*( (e**2)*np.cosh(2*xi) + e*np.cosh(xi) )*xi_1(xi, A)*xi_2(xi, A)
                      + 2*( 2*(e**2)*np.sinh(2*xi) + e*np.sinh(xi) )*( xi_1(xi, A) )**3  )    )





# Questa funzione inverte l'ordine di alcuni elementi dell'array in ingresso, nello specifico pone per primi quegli elementi che si trovano oltre la metà. Ciò viene fatto al fine di rendere più bella la rappresentazione grafica degli array ottenuti mediante l'utilizzo di fft, evitando che vengano disegnati dei punti viccini ma non congiunti e punti distanti ma congiunti
def array_graf(a):

    n= len(a)

    if( np.iscomplexobj(a) ):
        b= np.zeros(n, dtype=complex)

    else:
        b= np.zeros(n)

    for i in range(0, int(n/2)):
        b[i] = a[int(n/2) + i]

    for i in range(0, int(n/2)):
        b[int(n/2) + i]= a[i]

    return b







# COSTRUZIONE ARRAY MASSE

m_1= np.linspace(m_1_min, m_1_max, n_m_1)
m_2= np.linspace(m_2_min, m_2_max, n_m_2)







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
    time[N-i]= -time[i]









'''
# grafico di xi in funzione di t

plt.figure()

plt.title("$\\xi$ in funzione del tempo")

plt.plot(time, xi, linestyle="-", marker="", color="blue")
plt.plot([time[0], time[-1]], [0, 0], linestyle="-", marker="", color="black", linewidth= 0.8)

plt.xlabel("t [s]")
plt.ylabel("$\\xi$")

plt.show()
'''






# CALCOLO DELLA CORREZIONE ALLO SPETTRO IN FREQUENZA



DP= np.zeros( [ len(freq), len(m_1), len(m_2) ], dtype=complex)




for i in range(0, len(m_1)):

    if(i<10):
        print(" {:.0f} % ".format(100*i/len(m_1)) )

    else:
        print("{:.0f} % ".format(100*i/len(m_1)) )

    for j in range(0, len(m_2)):

        # M è la massa totale del sistema in masse solari
        M=m_1[i]+m_2[j]

        # mu è la massa ridotta del sistema in masse solari
        mu= m_1[i]*m_2[j]/M


        # alpha è il prodotto tra la costante di gravitazione aniversale G, la massa totale M, e la massa ridotta mu
        alpha= G*M*mu

        # A è una costante che compare nelle componenti del quadrupolo e nelle unita di misura scelta si misura in secondi
        A= np.sqrt(mu*a**(3)/alpha)


        # CALCOLO DEGLI xi CORRISPONEDENTI AI TEMPI E RAPPRESENTAZIONE GRAFICA


        xi= np.zeros(len(time))

        for l in range(0, len(time)):

            xi[l]= inver_t_xi(time[l], A)

        xi= array_graf(xi)



        # CALCOLO DELLE TRASFORMATE DI FOURIER DEI QUADRATI DELLE DERIVATE CUBICHE DELLE COMPONENTI DEL QUADRUPOLO


        tras_11_cubo = np.fft.rfft( M_11_3(xi, mu, A)*np.hanning(len(xi)) )
        #tras_12_cubo = np.fft.rfft( M_12_3(xi, mu, A)*np.hanning(len(xi)) )
        tras_22_cubo = np.fft.rfft( M_22_3(xi, mu, A)*np.hanning(len(xi)) )
        tras_33_cubo = np.fft.rfft( M_33_3(xi, mu, A)*np.hanning(len(xi)) )






        # CALCOLO DELLE ANTITRASORMATE


        antitras_11_quadro= ( M_11_3(xi, mu, A) )**2*np.hanning(len(xi))
        antitras_12_quadro= ( M_12_3(xi, mu, A) )**2*np.hanning(len(xi))
        antitras_22_quadro= ( M_22_3(xi, mu, A) )**2*np.hanning(len(xi))
        antitras_33_quadro= ( M_33_3(xi, mu, A) )**2*np.hanning(len(xi))





        # CALCOLO DELLE TRASFORMATE DI FOURIER DEI QUADRATI DELLE DERIVATE CUBICHE DELLE COMPONENTI DEL QUADRUPOLO


        # trasformata del quadrato della derivata cubica della componenti 11 del quadrupolo
        tras_11= np.fft.rfft(antitras_11_quadro)

        # trasformata del quadrato della derivata cubica della componenti 12 del quadrupolo
        tras_12= np.fft.rfft(antitras_12_quadro)

        # trasformata del quadrato della derivata cubica della componenti 22 del quadrupolo
        tras_22= np.fft.rfft(antitras_22_quadro)

        # trasformata del quadrato della derivata cubica della componenti 33 del quadrupolo
        tras_33= np.fft.rfft(antitras_33_quadro)



        for k in range(0, len(freq)):


            DP[k,i,j]= tras_11[k]*np.conjugate( tras_11_cubo[k] )  + tras_22[k]*np.conjugate( tras_22_cubo[k] ) + tras_33[k]*np.conjugate( tras_33_cubo[k] ) + 3*tras_12[k]*( np.conjugate( tras_11_cubo[k] ) + np.conjugate( tras_22_cubo[k] ) )

            DP[k,i,j]= cost*DP[k,i,j]







DP= abs(DP)


# SALVATAGGIO SU FILE

string_e= f'{e:.3f}'
string_e= string_e.replace(".", "_")

string_a= f'{a:.3f}'
string_a= string_a.replace(".", "_")



file_name_DP= file_name_DP + "____e_" + string_e + "_a_" + string_a
file_name_DP= path + "\\" + file_name_DP + ".txt"

df = pd.DataFrame( DP.reshape((-1, DP.shape[-1])), index= pd.MultiIndex.from_product( [range(DP.shape[0]), range(DP.shape[1])] ) ).to_csv(file_name_DP)



file_name_f= file_name_f + "____e_" + string_e + "_a_" + string_a
file_name_f= path + "\\" + file_name_f + ".txt"

file_f= open(file_name_f, "w")

np.savetxt(file_f, freq)
file_f.close()



file_name_m_1= file_name_m_1 + "____e_" + string_e + "_a_" + string_a
file_name_m_1= path + "\\" + file_name_m_1 + ".txt"

file_m_1= open(file_name_m_1, "w")

np.savetxt(file_m_1, m_1)
file_m_1.close()



file_name_m_2= file_name_m_2 + "____e_" + string_e + "_a_" + string_a
file_name_m_2= path + "\\" + file_name_m_2 + ".txt"

file_m_2= open(file_name_m_2, "w")

np.savetxt(file_m_2, m_2)
file_m_2.close()





