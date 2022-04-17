import numpy as np
import pylab


path="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt"
name="prova.txt"

a=np.array([1,2,3,4])
b=np.array([11,12,13,14])


file_name=path + "\\" + name
file= open(file_name, "w")
np.savetxt(file, np.c_[a, b])

file.close()

x, y= np.loadtxt(file_name, unpack=True)
print(x)
print(y)