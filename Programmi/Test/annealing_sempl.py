import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import dual_annealing
from mpl_toolkits.mplot3d import Axes3D

# numero di iterazioni
N=100

# ampiezza del passo
b=10**(-1)


# DEFINIZIONE DELLE VARIABILI

# funzione con un minimo locale positivo e un minimo globale negativo

def funz(x):

    return 1.2*(x[0]+0.2)**(4)+x[1]**(4)-2*(x[0]+0.2)**(3)+(x[0]+0.2)-3


x=np.array([1, 2])
print(funz(x))

# derivata parziale rispetto la prima variabile

def deriv_x(*x):

    return 1.2*4*(x[0]+0.2)**(3)-2*3*(x[0]+0.2)**(2) + 1



res= dual_annealing(funz, ([-1, 1], [-1, 1]))


risultati=res.x



# REALIZZAZIONE DEL GRAFICO

x=np.linspace(-2,2,1000)
y=np.linspace(-2,2,1000)


fig=plt.figure()

ax=fig.add_subplot(111, projection="3d")

X, Y=np.meshgrid(x, y)

var=[]
var.append(X)
var.append(Y)

ax.plot_surface(var[0], var[1], funz(var))
ax.plot(risultati[0], risultati[1], funz(risultati), color="red", marker="o")

ax.set_zlim(-3,-1)

plt.show()
