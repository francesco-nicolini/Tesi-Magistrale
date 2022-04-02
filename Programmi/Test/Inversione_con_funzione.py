import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import *

# N contiene il numero di iterazioni
N=20


# val contiene i valori di partenza dei parametri
val=np.array([1,100])

# matrix è una lista contenente le matrci corrispondenti alle varie equazioni
matrix=[]

a=np.array( [(1,3),(2,1)] )
matrix.append(a)
a=np.array( [(1,-1),(-1,100)] )
matrix.append(a)

K=len(matrix[0][0])

# cost è una lista contenente i termini noti delle equazioni
cost=[1,0]

# creazione della lista di simboli
simb=[]

for i in range(0,K):

    a=Symbol("x_{0}".format(i))
    simb.append(a)

# f_m è la distribuzione in massa
def f_m(m, *par):
    '''
    return 1/(sqrt(par[1]**2))*exp(-(m-par[0])**2/(2*(par[1])**2))

    2*math.pi*
    '''
    return cos(par[0]*m)+par[1]

# m contiene i valori delle masse considerate nella discretizzazione dell'integrale
m=np.linspace(0,1,K)


# vet è una funzione che restuisce un vettore le cui componenti corrispondono a f_m calcolato per
# ogni valore di m.
def vet(*a):

    varia=[]

    for i in range(0,K):

        varia.append(f_m(m[i],*a))

    return varia

# funz è una lista contenente le funzioni che costituiscono il sistema
# tramite le matrici in matrix sono create tali funzioni
funz=[]

for i in range(0,K):

    def f(*par,i=i):

        varia=[]

        for j in range(0,K):

            varia.append(f_m(m[j],*par))

        vec_dx=np.dot(varia,matrix[i])
        return np.dot(vec_dx,varia)-cost[i]

    funz.append(f)


# J_symb è una matrice contenente la forma simbolica dello jacobiano
J_simb=[]

for i in range(0,K):

    a=[]

    for j in range(0,K):

        z=funz[i](*simb).diff(simb[j])
        fu=lambdify(simb,z, 'numpy')

        a.append(fu)

    J_simb.append(a)


# Inizializzazione di vari oggetti (F e J contengono rispettivament il valore di ogni funzione
# del sistema e lo jacobiano calcolati per le soluzioni all'n-esima interazione
# var contiene le soluzione all'n-esima interazione
F=np.zeros(K)

J=np.zeros((K,K))

# Vettori introdotti per poter fare i grafici finali
graf=[]

for i in range(0,K):

    a=np.zeros(N+1)
    graf.append(a)

print(val)

# Iterazioni per la determinazione delle soluzioni
for k in range(0,N):

# riempie i vettori per i grafici
    for i in range(0,K):
        graf[i][k]=val[i]

# applicazione metodo di Newton
    for i in range(0,K):

        F[i]=funz[i](*val)

        for j in range(0,len(J[i])):

            J[j][i]=J_simb[j][i](*val)

    delta=np.linalg.solve(J,F)
    val= val - delta

for i in range(0,K):
    graf[i][N]=val[i]

# Stampa dei risultati finali
print("Soluzioni:\n",val,"\n")

print("Funzioni calcolate con i valori delle variabili individuati:")
for i in range(0,K):
    print("f_{0} =".format(i), funz[i](*val))


# Grafico
plt.title("Valore delle Variabili al Variare del Numero di Iterazioni")
plt.xlabel("Numero iterazioni")
plt.ylabel("Valori variabili")

iter=np.linspace(0,N,N+1)
for i in range(0,K):
    plt.plot(iter, graf[i], color="C{0}".format(i), label=simb[i], linestyle="-", marker="", markersize=3)

plt.legend()
plt.show()