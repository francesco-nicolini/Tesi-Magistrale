import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import dblquad


def f_m(m,mu,sigm):

    return 1/np.sqrt(2*math.pi*sigma**2)*np.exp(-(m-mu)**2/(2*sigma**2))

def integranda(x,y,nu):

    return x*y**np.log(nu)

def fun_conf(nu):

    return 2/(np.log(nu)+1)*2**(np.log(nu)+1)

nu=np.linspace(1,10,100)

omega=np.zeros(len(nu))

for i in range(0,len(nu)):

    omega[i]=dblquad(integranda,0,2, lambda x: 0, lambda x: 2, args=[nu[i]])[0]


plt.plot(nu, omega, marker='o', color="blue", markersize=2.5 ,linestyle="")
plt.plot(nu, fun_conf(nu), color="red", linestyle="-")

plt.show()