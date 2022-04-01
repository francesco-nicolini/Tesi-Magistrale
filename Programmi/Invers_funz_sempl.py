import numpy as np
import matplotlib.pyplot as plt
from sympy import *



# N contiene il numero delle iterazioni

N=30



# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)

K=20



# estremo superiore e estremo inferiore delle frequenze considerate nell'integrale

freq_min=1
freq_max=100



# estremo superiore e estremo inferiore delle masse considerate nell'integrale e scarto tra due masse consecutive

m_min=1
m_max=10

dm= (m_max - m_min)/(K-1)



# contiene il valore di partenza delle variabile

var=np.linspace(0, 10, K)



# lista delle masse indagate nella discretizzazione dell'integrale (sono in masse solari)

freq=np.logspace(np.log10(freq_min), np.log10(freq_max), K)



# lista delle masse indagate nella discretizzazione dell'integrale (sono in masse solari)

masse=np.linspace(m_min, m_max, K)



# funzione per il calcolo di omega in funzione della frequenza secondo un certo modello

def funz_omeg(nu):

    n=len(nu)

    return 1*np.ones(n)



# la funzione usata per costruire le matrici

def integ(m_1, m_2, nu):

    z=(m_1*m_2) + nu

    return z



# la funzione che crea le matrici

def fun_mat(nu, K):

    a=np.zeros((K,K))

    for i in range(0, K):
        for j in range(0,K):

            a[i][j]=integ(masse[i], masse[j], nu)

    return a



# omeg contiene i valori di omega per le frequenze analizzate

omeg=funz_omeg(freq)



# La lista matrix è riempita con delle matrici delle forme quadrariche

matrix=[]

for i in range(0, K):

    matrix.append(fun_mat(freq[i], K))



# Tramite le matrici in matrix sono create le funzioni che vanno a comporre il sistema e che sono inserite nella lista funz

funz=[]

for i in range(0,K):

    def f(*var,i=i):
        vec_dx=np.dot(var,matrix[i])
        return np.dot(vec_dx,var)-omeg[i]


    funz.append(f)



# Tali variabile sono introdotte per poter svolgere le derivate delle funzioni

simb=[]

for k in range(0,K):
    simb.append(Symbol('x_{0}'.format(k)))



# Inizializzazione di vari oggetti (F e J contengono rispettivament il valore di ogni funzione del sistema e lo jacobiano calcolati per le soluzioni all'n-esima interazione. var contiene le soluzione all'n-esima interazione

F=np.zeros(K)
J=np.zeros((K,K))



# Vettori introdotti per poter fare i grafici finali

graf=[]

for i in range(0,K):

    a=np.zeros(N+1)
    graf.append(a)



# determinazione dello Jacobiano in forma simbolica

J_simb=[]

for i in range(0,K):

    a=[]

    for j in range(0,K):

        z=funz[i](*simb).diff(simb[j])
        fu=lambdify(simb,z, 'numpy')

        a.append(fu)

    J_simb.append(a)



# k è l'indice dell'iterazione
for k in range(0,N):

    print(k/N*100, "%")

# riempie i vettori per i grafici
    for i in range(0,K):
        graf[i][k]=var[i]

# applicazione metodo di Newton
    for i in range(0,len(F)):

        F[i]=funz[i](*var)

        for j in range(0,len(J[i])):

            J[j][i]=J_simb[j][i](*var)

    delta=np.linalg.solve(J,F)
    var=var - delta

# inserisce negli array per i grafici i risultati finali ottenuti

for i in range(0,K):
    graf[i][N]=var[i]



# Stampa dei risultati finali

print("\nSoluzioni:\n",var,"\n")


if (K<30):

    print("Funzioni calcolate con i valori delle variabili individuati:")
    for i in range(0,K):
        print("f_{0} =".format(i), funz[i](*var))



# Grafico degli andamenti delle soluzioni all'aumentare del numero di iterazioni

if (K<30):

    plt.figure()

    plt.title("Valore delle Variabili al Variare del Numero di Iterazioni")
    plt.xlabel("Numero iterazioni")
    plt.ylabel("Valori variabili")

    iter=np.linspace(0,N,N+1)
    for i in range(0,K):
        plt.plot(iter, graf[i], color="C{0}".format(i), label=simb[i], linestyle="-", marker="", markersize=3)

    plt.legend()
    plt.show()



# Grafico della soluzione ottenuta in funzione della massa

plt.figure()

plt.title("f(m) in funzione di m")
plt.xlabel("massa [$M_\odot\$]")
plt.ylabel("f(m)")


plt.plot(masse, var, color="blue", linestyle="-", marker="")

plt.show()




