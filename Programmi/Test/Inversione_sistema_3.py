import numpy as np
import matplotlib.pyplot as plt
from sympy import *



# contiene il valore di partenza delle variabile
var=np.array([18,1,12,1])
#var=np.array([30,10])
# N contiene il numero delle iterazioni
N=30



# Inizializazione liste, cost conterra i coefficienti che nelle varie funzioni non moltiplicano
# delle variabili, funz conterrà le funzioni del sistema
funz=[]
matrix=[]
cost=[2,0,3,2]
#cost=[30,20]

# La lista Matrix è riempita con delle matrici delle forme quadrariche


a=np.array( [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)] )
matrix.append(a)
a=np.array( [(1,0,0,1),(0,0,0,0),(0,0,-1,0),(1,0,0,1)] )
matrix.append(a)
a=np.array([(1,1,1,1),(1,1,1,1),(1,1,1,1),(1,1,1,1)])
matrix.append(a)
a=np.array([(1,1,1,2),(1,0,1,1),(1,1,0,1),(2,1,1,1)])
matrix.append(a)
'''

a=np.array( [(1,2),(2,1)] )
matrix.append(a)
a=np.array( [(1,0),(0,-15)] )
matrix.append(a)
'''
# K contiene il numero di valori della funzione che si vogliono calcolare (ossia il numero delle variabili)
K=len(matrix)

# Tramite le matrici in matrix sono create le funzioni che vanno a comporre il sistema
for i in range(0,K):

    def f(*var,i=i):
        vec_dx=np.dot(var,matrix[i])
        return np.dot(vec_dx,var)-cost[i]


    funz.append(f)


# Tali variabile sono introdotte per poter svolgere le derivate delle funzioni
simb=[]

for k in range(0,K):
    simb.append(Symbol('x_{0}'.format(k)))

'''
print(simb,"\n")


# Stampa le funzioni per verificare che siano quelle desiderate
for i in range(0,K):

    print("f_{0}=".format(i),funz[i](*simb))

    for j in range(0,K):

        z=funz[i](*simb).diff(simb[j])
        print("derivata rispetto ad",simb[j],"=",z)

    print("\n")
'''

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

for i in range(0,K):
    graf[i][N]=var[i]

# Stampa dei risultati finali
print("Soluzioni:\n",var,"\n")

print("Funzioni calcolate con i valori delle variabili individuati:")
for i in range(0,K):
    print("f_{0} =".format(i), funz[i](*var))

# Grafico
plt.title("Valore delle Variabili al Variare del Numero di Iterazioni")
plt.xlabel("Numero iterazioni")
plt.ylabel("Valori variabili")

iter=np.linspace(0,N,N+1)
for i in range(0,K):
    plt.plot(iter, graf[i], color="C{0}".format(i), label=simb[i], linestyle="-", marker="", markersize=3)

plt.legend()
plt.show()
