import numpy as np

path="C:\\Users\\39366\\Dropbox\\PC\\Desktop"

file_name= "prova"

a= np.array([1,2,3,4,5,6,7,8,9,10])
b= np.array([10,20,30,40,50,60,70,80,90,100])
c= np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.01])


file_name= path + "\\" + file_name

file= open(file_name, "w")

np.savetxt(file, np.c_[a, b, c])
file.close()


A, B, C= np.loadtxt(file_name, unpack=True)

print(A,"\n")
print(B,"\n")
print(C)
