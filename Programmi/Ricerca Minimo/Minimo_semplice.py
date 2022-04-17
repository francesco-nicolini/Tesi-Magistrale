import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# numero di iterazioni
N=100000

# ampiezza del passo
b=10**(-1)


# DEFINIZIONE DELLE VARIABILI

# funzione con un minimo locale positivo e un minimo globale negativo

def funz(*x):

    return 1.2*(x[0]+0.2)**(4)+x[1]**(4)-2*(x[0]+0.2)**(3)+(x[0]+0.2)-3

# derivata parziale rispetto la prima variabile

def deriv_x(*x):

    return 1.2*4*(x[0]+0.2)**(3)-2*3*(x[0]+0.2)**(2) + 1

# derivata parziale rispetto la seconda variabile

def deriv_y(*x):

    return 4*x[1]**3

# lista contenente le derivate

derivate=[]
derivate.append(deriv_x)
derivate.append(deriv_y)





# INDIVIDUAZIONE DEL MINIMO

# valori iniziali
val_iniz=np.array([-1,-1])

# val contiene i valori ad ogni iterazione
val=val_iniz



# lista che contiene i valori trovati ad ogni iterazione e usata per disegnare i punti corrispondenti nei grafici
graf=[]

for i in range(0, len(val_iniz)):
    a=np.zeros(N+1)
    graf.append(a)



# array che contiene i valori del nuovo punto (se mettessi val[0]=val[0] - b*derivate(*val) andrebbe bene, ma val[1]=val[1] - b*derivate(*val) darebbe problemi poichè val a destra conterrebbe il vecchio valore di val[1] e il nuovo valore di val[0]

provv= np.ones(len(val))



# iterazione per l'individuazione del punto di minimo

for i in range(0,N):

    for j in range(0, len(val)):
        graf[j][i]=val[j]


    for j in range(0, len(val)):
        provv[j]= val[j] - b*derivate[j](*val)

    val= provv



for i in range(0, len(val)):
    graf[i][N]=val[i]

'''
for i in range(0, len(graf[0])):
    print("({0},{1})".format(graf[0][i], graf[1][i]))
'''

print("Il punto di minimo trovato è ({0},{1})".format(val[0], val[1]))




# REALIZZAZIONE DEL GRAFICO

x=np.linspace(-2,2,1000)
y=np.linspace(-2,2,1000)


fig=plt.figure()

ax=fig.add_subplot(111, projection="3d")

X, Y=np.meshgrid(x, y)

var=[]
var.append(X)
var.append(Y)

ax.plot_surface(var[0], var[1], funz(*var))
ax.plot(graf[0], graf[1], funz(*graf), color="red", linewidth=6)

ax.set_zlim(-3,-1)

plt.show()