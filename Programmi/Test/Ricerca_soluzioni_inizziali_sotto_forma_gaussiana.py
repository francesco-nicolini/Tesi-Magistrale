import numpy
import matplotlib.pyplot as plt

a= -1
b= 20
c= -90


def f(x, a, b, c):

    par= a*x**2+b*x+c

    if ( par<0 ):
        return 0

    else:
        return par




def param_gauss(masse, funz):


    ind_max= np.argmax(funz)
    massimo= funz[ind_max]

    meta= massimo/2

    for i in range(ind_max, 0, -1):

        if ( funz[i]<=meta ):
            break

    ind_sx= i


    for i in range(ind_max, len(funz)):

        if ( funz[i]<=meta ):
            break

    ind_dx= i

    sigma= (masse[ind_dx] - masse[ind_sx])/2

    print("La semi larghezza a metà altezza è pari a {:.3}".format(sigma))

    sigma= sigma/np.sqrt(2)

    A= np.sqrt( massimo/(np.sqrt(np.pi)*sigma) )

    mu= masse[ind_max]/2

    return A, mu, sigma, ind_sx, ind_dx



def gauss(x, A, mu, sigma):

    return A*np.exp(-( (x-mu)/(np.sqrt(2)*sigma) )**2)





masse= np.linspace(0, 30, 1000)

funz= np.zeros(len(masse))

for i in range(0, len(masse)):

    funz[i]= f(masse[i], a, b, c)



A, mu, sigma, ind_sx, ind_dx= param_gauss(masse, funz)


gaussiana= gauss(masse, np.sqrt(np.pi)*sigma*A**2, 2*mu, np.sqrt(2)*sigma)





plt.figure()


plt.title("TEST APPROSSIMAZIONE CON GAUSSIANA")

plt.plot(masse, funz, linestyle="-", marker="", color="blue", label="Funzione di partenza")
plt.plot(masse, gaussiana, linestyle="-", marker="", color="red", label="Approssimazione")

plt.plot([masse[ind_sx], masse[ind_dx]], [funz[ind_sx], funz[ind_sx]], linestyle="--", marker="", color="green")

plt.plot([masse[0], masse[-1]], [0, 0], linestyle="-", marker="", color="black", linewidth=0.8)

plt.xlabel("masse [U.A.]")
plt.ylabel("F(M)")

plt.xlim(masse[0], masse[-1])

plt.legend()

plt.show()







