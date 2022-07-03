import numpy as np


num_zeri= 5

masse= np.linspace(1, 10, 10)

print(masse)

dm= masse[1] - masse[0]

masse_sotto= np.linspace( masse[0] - dm - (num_zeri-1)*dm, masse[0] - dm, num=num_zeri)

masse= np.concatenate((masse_sotto, masse))

print(masse)

masse_sopra= np.linspace(masse[-1] + dm, masse[-1] + dm + (num_zeri-1)*dm, num=num_zeri)

masse= np.concatenate((masse, masse_sopra))

print(masse)


a= np.array([11, 13, 22, 3, 4, 21, 5, 4, 7, 4])
zeri= np.zeros(num_zeri)

a= np.concatenate((zeri, a))
a= np.concatenate((a, zeri))

print(a)

print(len(masse), len(a))


