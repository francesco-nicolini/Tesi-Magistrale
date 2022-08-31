import numpy as np
import matplotlib.pyplot as plt

# x_max e y_max indicano l'ascissa e l'ordinata del massimo
x_max= 3
y_max= 10

# x_inf e x_sup indicano l'estremo superiore e l'estremo inferiore del supporto
x_inf= 0
x_sup= 20

# minimo e massimo sono gli estremi dell'intervallo delle x rappresentato nel grafico
minimo=-30
massimo= 30


def triangolo(m, x_max, y_max, x_inf, x_sup):


    if ( x_inf<x_max<x_sup):

        if( m<x_inf or m>x_sup ):
            return 0

        elif ( x_inf<=m<=x_max ):

            incl= y_max/(x_max - x_inf)
            q= -incl*x_inf

            return incl*m + q

        elif ( x_max<m<=x_sup ):

            incl= y_max/(x_max - x_sup)
            q= -incl*x_sup

            return incl*m + q


    else:
        print("Errore nella scelta dei parametri")
        return

m= np.linspace(minimo, massimo, 1000)

tri= np.ones(len(m))

for i in range(0, len(m)):

    tri[i]= triangolo(m[i], x_max, y_max, x_inf, x_sup)


plt.figure()

plt.title("Distribuzione a Triangolo")

plt.plot([minimo, massimo], [0, 0], linestyle="-", marker="", color="black", linewidth=0.8)
plt.plot(m, tri, linestyle="-", marker="", color="blue")

plt.xlabel("masse [u.a.]")
plt.ylabel("f(m) [u.a.]")

plt.xlim(minimo, massimo)

plt.show()
