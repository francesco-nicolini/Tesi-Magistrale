import numpy as np
import matplotlib.pyplot as plt
import math
import keyboard



# OPZIONI VARIE

# N contiene il numero delle iterazioni
N=5000

# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)
K=50


# Estremo superiore e estremo inferiore delle masse considerate nell'integrale e scarto tra due masse consecutive
m_min=1
m_max=30

dm= (m_max - m_min)/(K-1)

# Se option è pari a "read" omega_GW è letto da file altrimenti viene creato tramite una funzione
option="read"

# Path del file se option è uguale a "read"
file_name_omega="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt\\omega_GW.txt"

# Se option è pari a "read" ponendo disegna uguale a True verra realizzato il grafico della soluzione esatta nello stesso piano in cui vengono rappresentate le soluzioni trovate minimizzando
disegna=True

# Path del file se disegna è uguale a True
file_name_f_m="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt\\f_m.txt"


# Se metodo è uguale a "semplice" il minimo viene ricercato usando esclusivamente il gradiente.
# Se metodo è pari a "barriera" il minimo vie ne cercato con il gradiente modificato da una funzione data dal prodotto tra rho e la somma di tutti gli esponenziali di argomento meno rho per il valore di una delle variabile su cui si sta minimizzando e (rho è un coefficiente che deve essere fornito e che dovreebe essere grande)
# Se metodo è pari a "barriera + semplice" effettua prima la minimizzazione con il metodo "barriera" e quindi utilizza i valori trovati come valori iniziali per il metodo semplice.
# Se metodo è pari a "annealing" la ricerca del minimo viene effettuata mediante la funzione scipy.optimize.dual_annealing
# Se metodo è "SPA" (Semplice con Passo Aggiornato) si procede come nel caso "semplice", con la differenza però che si conforonta il nuovo valore della funzione con il avlore al paasso precedente e se la prima risulta maggiore il valore delle soluzioni non si aggiorna e il passo viene diviso per agg, altrimenti si ha l'aggiornamento e il passo viene moltiplicato per agg. La ricerca termina quando il passo diviene inferiore a soglia. Se b superiore a massimo_b, b smette di essere aggiornato fino a che non inizia a diminuire. dim_loop è il numero di passi precedenti a quello corrente che vengono considerati al fine di capire se si è raggiunto un loop, come descritto nel commento a memory nella funzione min_agg()

# b è l'ampiezza del passo
b= 10**(-5)

metodo= "SPA"

rho= 100

agg_dis= 2
agg_sal=1.1

soglia= 10**(-20)

massimo_b= 10**(300)

dim_loop= 100


# delle prime num_dis variabili viene realizzato il grafico del loro andamento al variare dell'iterazione

num_dis= 15

# Estremo superiore e estremo inferiore delle frequenze considerate (da selezionare solo se option è diverso da "read")
freq_min=10**(-8)
freq_max=10**(1)






# VALORI INIZIALI DELLA DISTRIBUZIONE DA CERCARE

# Contiene il valore di partenza delle variabile


val_iniz=10*np.ones(K)

'''
for i in range(12, 17):
    val_iniz[i]= 40/5 + val_iniz[i-1]

for i in range(17, 21):
    val_iniz[i]= val_iniz[i-1] - 40/4



for i in range(12, 19):

    val_iniz[i]= 40

'''



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
h_70= 0.7

# Valore di omega della materia oscura
ome_DM= 0.25

# Valore di delta_loc
delta_loc= 10**8

# Velocità iniziale tra i due buchi neri in Km/s
v_0= 10

# Costante che moltiplica l'integrale nell'equazione 18 dell'articolo 2109.11376. (1/10)**2 compare poiche nell'equazione compare (nu/10)**2 (nu in Hz), (1/1000)**2 compare poiché nell'equazione compare (dm1/1000)*(dm2/1000) (m1 e m2 aono in masse solari)

cost= 2.39*10**(-13)*h_70*(ome_DM/0.25)**2*(delta_loc/10**8)*(v_0/10)**(-11/7)







# CREAZIONE DELL'ARRAY CONTENENTE LE FREQUENZE E DELL'ARRAY CONTENENTE I VALORI DI OMEGA_GW

# freq contiene la lista delle frequenze
# omega_GW contiene la lista dei valori dello spettro in potenza del fondo


if (option=="read"):

    freq, omega_GW= np.loadtxt(file_name_omega, unpack=True)

    if (disegna==True):

        #  masse_graf contiene le masse per effettuare il grafico di f_esatta
        # f_esatta la forma della distribuzione letta da file (è la soluzione esatta)

        masse_graf, f_esatta= np.loadtxt(file_name_f_m, unpack=True)

else:

    def funz_omeg(nu):

        n=len(nu)

        return 10**(-22)*np.ones(n)


    freq=np.logspace(np.log10(freq_min), np.log10(freq_max), K)
    omega_GW=funz_omeg(freq)


# Grafico di omega_GW

fig, ax = plt.subplots()

ax.plot(freq, omega_GW, linestyle="-", color="blue")

plt.title("$\\Omega_{GW}$ in funzione della frequenza")
plt.xlabel("f [Hz]")
plt.ylabel("$\\Omega_{GW}$")
plt.yscale("log")
plt.xscale("log")

plt.show()



print("La dimensione del vettore contenente i valori di Omega_GW è pari a {0}".format(len(omega_GW)))



# DEFINIZIONE DELLE FUNZIONI CHE COMPAIONO NELLA RELAZIONE DI RIFERIMENTO UNA VOLTA CHE E' STATO DISCRETIZZATO L'INTEGRALE


# Lista delle masse indagate nella discretizzazione dell'integrale (sono in masse solari)

masse=np.linspace(m_min, m_max, K)



# La funzione usata per costruire le matrici (si tiene conto anche del dm**2 in essa)

def integ(m_1, m_2, nu):

    z= dm**2*( ( m_1 + m_2 )**(23/21) )/( m_1*m_2 )**(5/7)

    return (nu**(2/3))*z



# La funzione che crea le matrici (si è usata la regola del trapezio, vedi file Prova_costruz_matrice_per_regola_trapez.py)


def fun_mat(nu, K):

    a=np.zeros((K,K))

    for i in range(0,K):

        for j in range(0,K):

            if ( (i+j==0) or (i+j==2*(K-1)) or (i*j==0 and i+j==K-1) ):

                a[i][j]=integ(masse[i], masse[j], nu)/4


            # è un elif quindi si possono mettere casi già previsti dall'if
            #precedente dal momento che questi, verificando quest'ultimo, non saranno
            #analizzati dal secondo
            elif ( i*j==0 or i==K-1 or j==K-1 ):

                a[i][j]=integ(masse[i], masse[j], nu)/2

            else:

                a[i][j]=integ(masse[i], masse[j], nu)

    return a



# La lista matrix è riempita con delle matrici delle forme quadrariche

matrix=[]

for i in range(0, K):


    matrix.append(fun_mat(freq[i], K))



# Tramite le matrici in matrix sono create le funzioni che vanno a comporre il sistema e che sono inserite nella lista funz

funz=[]

for i in range(0,K):

    def f(*var,i=i):

        vec_dx=np.dot(var,matrix[i])
        return np.dot(vec_dx,var)-omega_GW[i]/cost

    funz.append(f)



# Funzioni per il calcolo delle derivate parziali, queste sono inserite nella lista chiamata deriv (l è l'indice della variabile rispetto alla quale si vuole calcolare la derivata parziale)

deriv=[]

for i in range(0, K):

    def f_deriv(*var, l, i=i):

        sum= 2*var[l]*matrix[i][l][l]

        for j in range (0, K):

            if (j!=l):

                sum=sum + 2*var[j]*matrix[i][j][l]

        return sum

    deriv.append(f_deriv)







# CREAZIONE DELLA FUNZIONE DA MINIMIZZARE E CALCOLO DELLE SUE DERIVATE PARZIALI

# Questa è la funzione di cui si derve cercare il minimo

def funz_da_minim(*var):

    somma=0

    for i in range(0, len(funz)):
        somma=somma + (funz[i](*var))**2

    return somma



# Funzione per il calcolo delle derivate parziali della funzione funz_da_minim (index è l'indice della variabile rispetto alla quale si vuole calcolare la derivata parziale)

def gradiente(*var, index):

    somma=0

    for i in range(0, len(funz)):

        somma= somma + funz[i](*var)*deriv[i](*var, l=index)

    return 2*somma







# VARIE FUNZIONI PER EFFETTUARE LA MINIMIZZAZIONE

# funzione che si basa esclusivamente sul gradiente. funz è la funzione da minimizzare, gradinete il gradiente di quest'ultima, b è lacostante che moltiplica il gradiente quando si calcola il passo successivo dell'iterazione, iniz_val sono i valori iniziali,  N è il numero delle iterazioni da effettuare, K è il numero delle variabili (sia N che K sono variabili globali quindi non devono essere fornite alla funzione). La funzione restituisce i valori che minimizzano e un lista contenente i valori trovat ad ogni iterazione in modo da poter effettuare dei grafici con l'andamento delle soluzioni al variare del numero di iterazioni. Restituisce prima le soluzioni e poi la lista.

def min_sempl(funz, gradiente, b, *iniz_val):

    global K, N


    print("Funzione calcolata con i valori iniziali:", funz(*iniz_val))


    # Vettori introdotti per poter fare i grafici finali

    graf=[]

    for i in range(0,K):

        a=np.zeros(N+1)
        graf.append(a)


    val=np.array(iniz_val)



    # array che contiene i valori del nuovo punto (esempio in due dimenzioni: se mettessi val[0]=val[0] - b*derivate(*val) andrebbe bene, ma val[1]=val[1] - b*derivate(*val) darebbe problemi poichè val a destra conterrebbe il vecchio valore di val[1] e il nuovo valore di val[0] e non i vecchi valori di entrambi


    provv= np.ones(len(val))




    # Ciclo con le iterazioni

    for i in range(0,N):

        for j in range(0, len(val)):
            graf[j][i]=val[j]


        for j in range(0, len(val)):

            grad= gradiente(*val, index=j)
            provv[j]= val[j] - b*grad

            #print(grad, val[j])


        #input("---------------------")

        for j in range(0, len(val)):

            val[j]=provv[j]





    # Inserimento delle soluzioni trovate nelle liste per la realizzazione dei grafici

    for i in range(0, len(val)):
        graf[i][N]=val[i]



    # Stampa dei risultati finali

    print("\nSoluzioni:\n",val,"\n")

    print("Funzione calcolate con i valori delle variabili individuati mediante il metodo \"semplice\":")

    print("f =", funz(*val))
    print("\nRapporto tra la funzione calcolata dopo aver eseguito l'algoritmo e la funzione calcolata con i valori iniziali:")
    print("rapporto=",funz(*val)/funz(*iniz_val))


    return val, graf





# Funzione per eviatare che vengano ottenuti valori negativi. Sitratta di un esponenziale.

def barriera(rho, *val):

    global K

    sum=0

    for i in range(0, K):

        sum= sum + np.exp(-rho*val[i])

    return rho*sum





# Funzione che si basa sul gradiente modificato da una funzione esponenziale in moda da evitare che si ottengano soluzioni con valori negativi. funz è la funzione da minimizzare, gradinete il gradiente di quest'ultima, barriera è la funzione che modifica il gradiente, rho è il parametro che compare in quest'ultima, b è lacostante che moltiplica il gradiente quando si calcola il passo successivo dell'iterazione, iniz_val sono i valori iniziali,  N è il numero delle iterazioni da effettuare, K è il numero delle variabili (sia N che K sono variabili globali quindi non devono essere fornite alla funzione). La funzione restituisce i valori che minimizzano e un lista contenente i valori trovat ad ogni iterazione in modo da poter effettuare dei grafici con l'andamento delle soluzioni al variare del numero di iterazioni. Restituisce prima le soluzioni e poi la lista.

def min_barriera(funz, gradiente, bar, rho, b, *iniz_val):

    global K, N


    print("Funzione calcolata con i valori iniziali:", funz(*iniz_val))


    # Vettori introdotti per poter fare i grafici finali

    graf=[]

    for i in range(0,K):

        a=np.zeros(N+1)
        graf.append(a)


    val=iniz_val



    # array che contiene i valori del nuovo punto (esempio in due dimenzioni: se mettessi val[0]=val[0] - b*derivate(*val) andrebbe bene, ma val[1]=val[1] - b*derivate(*val) darebbe problemi poichè val a destra conterrebbe il vecchio valore di val[1] e il nuovo valore di val[0] e non i vecchi valori di entrambi


    provv= np.ones(len(val))


    # Ciclo con le iterazioni

    for i in range(0,N):

        for j in range(0, len(val)):
            graf[j][i]=val[j]


        for j in range(0, len(val)):
            provv[j]= val[j] - b*( gradiente(*val, index=j) - bar(*val, rho) )

        val= provv



    # Inserimento delle soluzioni trovate nelle liste per la realizzazione dei grafici

    for i in range(0, len(val)):
        graf[i][N]=val[i]



    # Stampa dei risultati finali

    print("\nSoluzioni:\n",val,"\n")

    print("Funzione calcolate con i valori delle variabili individuati con il metodo \"barriera\":")
    print("f =", funz(*val))

    print("\nRapporto tra la funzione calcolata dopo aver eseguito l'algoritmo e la funzione calcolata con i valori iniziali:")
    print("rapporto=",funz(*val)/funz(*iniz_val))


    return val, graf




# Funzione che si basa sul gradiente come per min_sempl() ma in cui il passo è aggiornato ogni volta (vedi scelta modo all'inzio del programma). funz è la funzione da minimizzare, gradinete il gradiente di quest'ultima, b è il valore iniziale per lo scalre che moltiplica il gradiente quando si calcola il passo successivo dell'iterazione, iniz_val sono i valori iniziali, K è il numero delle variabili, soglia è il valore che deve raggiungere lo scalare affinche la ricerca termini, agg è il valore usato per aggiornare il passo, massimo_b è il valore che il passo non deve superare per evitare overflow e dim_loop è la dimensione dell'array che serve a capire se si ha un loop ( K, soglia, agg, massimo_b e dim_loop sono variabili globali quindi non devono essere fornite alla funzione). La funzione restituisce i valori che minimizzano e un lista contenente i valori trovat ad ogni iterazione in modo da poter effettuare dei grafici con l'andamento delle soluzioni al variare del numero di iterazioni. Restituisce prima le soluzioni e poi la lista.


def min_agg(funzione, gradiente, b, *iniz_val):

    global K, soglia, agg, massimo_b, dim_loop, funz


    print("\n\nPremere t per far terminare il programma e osservare i risultati ottenuti senza che la condizione sul passo sia stata eseguita.\n\n")

    print("Premere l per passare alla modalità lettura.\n\n")

    print("Premere r per ritornare alla modalità rapida.\n\n")

    print("Premere s per stampare ad ogni iterazione alcune informazioni ottenute.\n\n")

    print("Premere n per non stampare.\n\n")

    print("Premere g per visualizzare il grafico della soluzione trovata senza interrompere la ricerca.\n\n")


    print("Funzione calcolata con i valori iniziali:", funzione(*iniz_val))


    # Vettori introdotti per poter fare i grafici finali

    graf=[]

    for i in range(0,K):

        a=[]
        graf.append(a)


    val=np.array(iniz_val)


    for j in range(0, len(val)):
        graf[j].append(val[j])


    # costruzione dell'array che controlla che non si appia un loop: questo è dato da una sequenza di 0 e 1 alternati

    # il comando sotto crea una lista lunga dim_loop e le cui componenti sono tutte pari a 0

    loop= [0]*dim_loop

    for i in range(1, dim_loop):

        loop[i]= 1 - loop[i-1]


    # costruzione dell'array che tiene conto, per ognuna delle ultime dim_loop iterazioi, se la funzione calcolata nei nuovi valori è maggiore (in ale caso la componente corrispondente è 0) o minore (in tale caso la componente corrispondente è 1) di quella calcolata con i valori trovati al passo precedente. Se l'algoritmo arriva al punto in cui entra in un loop, ossia ad un iterazione il passo viene aumentato, a quello successivo diminuito, a quello successivo ancora aumentato e cosi via, allora le componenti di questo array sarano uguali a ad 1 e 0 in modo alternato e si avra pertanto un'uguaglianza con l'array loop

    memory=[2]*dim_loop


    # array che contiene i valori del nuovo punto (esempio in due dimenzioni: se mettessi val[0]=val[0] - b*derivate(*val) andrebbe bene, ma val[1]=val[1] - b*derivate(*val) darebbe problemi poichè val a destra conterrebbe il vecchio valore di val[1] e il nuovo valore di val[0] e non i vecchi valori di entrambi


    provv= np.ones(len(val))


    # count tiene conto del numero di iterazioni che sono state effettuate

    count= 0


    # se flag_l è apri ad 1 il programma passa alla modalità di lettura in cui ogni volta che un'iterazione è completatata il programma si blocca fino a che non si preme un pulsante (diverso da t o r)

    flag_l= 0

    # se flag_s è pari ad 1 il programma stampa alcune informazioni ad ogni iterazione

    flag_s= 1



    # se si realizza il grafico della soluzione attualmente individuata, viene anche disegnato il grafico della differenza tra Omega_GW (omega_ora) calcolato con la soluzione attuale e la forma esatta e i valori di Omega_GW determinati a partire dai valori iniziali (omega_iniz) a cui vengono sempre sottratti i valori esatti.

    omega_iniz= np.zeros(K)
    omega_ora= np.zeros(K)


    for i in range(0, K):
        omega_iniz[i]= funz[i](*val)


    # Ciclo con le iterazioni

    while(b > soglia):


        for j in range(0, len(val)):

            provv[j]= val[j] - b*gradiente(*val, index=j)



        if ( funzione(*val) < funzione(*provv)):

            b= b/agg_dis

            #list.pop(index) elimina la componente della lista di indice pari a index
            memory.pop(0)
            memory.append(0)


        else:

            if( b<massimo_b):
                b= b*agg_sal


            for j in range(0, len(val)):
                val[j]= provv[j]


            memory.pop(0)
            memory.append(1)


            # Inserimento delle soluzioni trovate nelle liste per la realizzazione dei grafici
            for j in range(0, len(val)):
                graf[j].append(val[j])



        count= count + 1





        if ( memory==loop ):
            if ( flag_s==1 ):
                print("E' stato individuato un loop")

            b=b/(agg**2)


        if( keyboard.is_pressed("s") ):
            print("\nModalità stampa\n")
            flag_s=1

        if( keyboard.is_pressed("n") ):
            print("\nModalità senza stampa\n")
            flag_s=0


        # comando testato in test_semplice_agg.py, se si preme il pulsante lt flag viene posto pari ad 1 e il programma passa alla modalità di lettura

        if( keyboard.is_pressed("l") ):
            print("\nModalità lettura\n")
            flag_l= 1


        # se si preme il pulsante r flag viene posto pari ad 0 e il programma ritorna alla modalità rapida

        if( keyboard.is_pressed("r") ):
            print("\nModalità rapida\n")
            flag_l= 0

        # se si preme g viene realizzato il grafico della soluzione trovata all'iterazione corrente senza interrompere la ricerca.

        if( keyboard.is_pressed("g") ):


            title="Soluzione all'Iterazione Numero {0}".format(count)
            plt.suptitle(title, fontsize=13)



            # grafico della soluzione


            plt.subplot(2,1,1)

            plt.title("Valori di f(m) Individuati", fontsize=12)


            plt.plot(masse, val, color="orange", linestyle="", marker=".", markersize=7, label="soluzione individuata")
            plt.plot(masse, val_iniz, color="lightcoral", linestyle="", marker=".", markersize= 7, label="valori iniziali")

            if (disegna==True):

                plt.plot(masse_graf, f_esatta, color="midnightblue", linestyle="-", marker="", label="soluzione corretta")


            plt.xlabel("massa [M_sole]", fontsize=10)
            plt.ylabel("f(m)", fontsize=10)

            plt.legend()

            # grafico di Omega_GW

            plt.subplot(2,1,2)

            for i in range(0, K):
                omega_ora[i]= funz[i](*val)

            plt.title("Confronto tra $\\Omega_{GW}$ Calcolato con i Valori Individuati e $\\Omega_{GW}$ Esatto", fontsize=12)



            plt.plot(freq, omega_iniz, color="red", linestyle="-", label="con valori iniziali")
            plt.plot(freq, omega_ora, color="blue", linestyle="-", label="con valori attuali")


            plt.xlabel("f [Hz]", fontsize=10)
            plt.ylabel("$\\Delta\\Omega_{GW}$", fontsize=10)

            plt.yscale("log")
            plt.xscale("log")


            plt.tight_layout()
            plt.legend()

            plt.show()




        if ( flag_s==1 ):

            print("iterazione= {:}, b= {:.1e}, rapporto= {:.10e}".format(count, b, funzione(*val)/funzione(*iniz_val)))
            '''
            print(memory)
            '''


        if ( flag_l==1 and flag_s==1 ):
            input("-------------------------------------------")



        # termina il loop se viene premuto il tasto t
        if( keyboard.is_pressed("t") ):
            print("\n\nCiclo interrotto\n\n")
            break


    # Stampa dei risultati finali

    print("\nSono state efftuate {0} iterazioni.".format(count))

    print("\nSoluzioni:\n",val,"\n")

    print("Funzione calcolate con i valori delle variabili individuati mediante il metodo \"Semplice con passo aggoirnato\":")
    print("f =", funzione(*val))

    print("\nRapporto tra la funzione calcolata dopo aver eseguito l'algoritmo e la funzione calcolata con i valori iniziali:")
    print("rapporto=",funzione(*val)/funzione(*iniz_val))

    print("\nComponente in modulo piu grande del gradiente calcolato con le soluzioni finali: ")

    mass_grad= 0

    for i in range(0, len(val)):

        a= abs( gradiente(*val, index=i) )

        if ( a>mass_grad ):
            mass_grad= a

    print(mass_grad)



    return val, graf






# FUNZIONI PER LA REALIZZAZIONE DI GRAFICI

# grafici realizza sia i grafici degli andamenti delle soluzioni trovate al avriare del numero di iterazioni (solo se K<30, altrimenti fa il grafico solo delle prime 10), sia il grafico della f(m).
# graf è la lista di array che contengono i valori delle varie variabili all'aumentare del numero di iterazioni, val è un array contenente le solzioni trovate, tipologia è una stringa che si aggiunge al titolo per chiarire con quale processo i valori che minimizzano sono stati individuati (se tiplogoa è pari a "nulla", non si aggiunge nulla al titolo)
# delle prime num_dis variabile viene realizzato un grafico al variare del numero dell'iterazione

def grafici(graf, val, tipologia ):

    global K, num_dis
    global masse

    global disegna, masse_graf, f_esatta


    dimenz= len(graf[0])

    # Grafico degli andamenti delle soluzioni all'aumentare del numero di iterazioni

    if (K<=num_dis):
        num_variab= K
        title= "Valore delle Variabili al Variare del Numero di Iterazioni"

    else:
        num_variab= num_dis
        title= "Valore delle Prime Dieci Variabili al Variare del Numero di Iterazioni"

    plt.figure()

    if (tipologia!="nulla"):
        title= title + " (" + tipologia + ")"


    plt.title(title)
    plt.xlabel("Numero iterazioni")
    plt.ylabel("Valori delle variabili")

    iter=np.linspace(0, dimenz, num= dimenz, endpoint=False )

    for i in range(0,num_variab):
        plt.plot(iter, graf[i], color="C{0}".format(i), label="f(m_{0})".format(i), linestyle="-", marker="", markersize=3)



    plt.legend()
    plt.tight_layout()





    # Grafico della soluzione ottenuta in funzione della massa

    plt.figure()

    if (tipologia=="nulla"):
        title= "f(m) in funzione di m"

    else:
        title= "f(m) in funzione di m"
        title= title + " (" + tipologia + ")"

    plt.title(title)
    plt.xlabel("massa [M_sole]")
    plt.ylabel("f(m)")


    plt.plot(masse, val, color="orange", linestyle="", marker=".", markersize=7, label="soluzione individuata")
    plt.plot(masse, val_iniz, color="lightcoral", linestyle="", marker=".", markersize= 7, label="valori iniziali")

    if (disegna==True):

        plt.plot(masse_graf, f_esatta, color="midnightblue", linestyle="-", marker="", label="soluzione corretta")

    plt.tight_layout()
    plt.legend()

    return








# ESECUZIONE DELL'ALGORITMO DI MINIMIZZAZIONE E REALIZZAZIONE DEI GRAFICI

if (metodo=="semplice"):

    val, graf= min_sempl(funz_da_minim, gradiente, b, *val_iniz)
    grafici(graf, val, tipologia=metodo)



elif (metodo=="barriera"):

    val, graf= min_barriera(funz_da_minim, gradiente, barriera, rho, b, *val_iniz)
    grafici(graf, val, tipologia=metodo)



elif (metodo=="barriera + semplice"):


    val_1, graf_1= min_barriera(funz_da_minim, gradiente, barriera, rho, b, *val_iniz)
    tipologia= "barriera"
    grafici(graf, val, tipologia= tipologia)

    val_2, graf_2= min_sempl(funz_da_minim, gradiente, b, *val_1)
    grafici(graf, val, tipologia= metodo)



elif (metodo=="SPA"):

    val, graf=min_agg(funz_da_minim, gradiente, b, *val_iniz)
    tipologia= "semplice con passo aggiornato"
    grafici(graf, val, tipologia= metodo)



else:

    print("E' stato attribuito alla variabile metodo un valore che non rientra nei casi elencati all'inizio, si procede per tanto con l'effettuare la minimizzazione mediante il metodo \"semplice\". ")

    val, graf= min_sempl(funz_da_minim, gradiente, b, val_iniz)

    tipologia="semplice"
    grafici(graf, val, tipologia= tipologia)



plt.show()




