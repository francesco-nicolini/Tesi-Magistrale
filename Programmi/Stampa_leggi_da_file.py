import numpy as np
import pylab

'''
path="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt"
name="prova_1.txt"

a=np.array([1,2,3,4])
b=np.array([11,12,13,14])


file_name=path + "\\" + name
file= open(file_name, "w")
np.savetxt(file, np.c_[a, b])

file.close()

x, y= np.loadtxt(file_name, unpack=True)
print(x)
print(y)
'''


path="C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\file_txt"
name="prova_2.txt"

a=np.array([1.1234,2.523,3.35234,4])
b=np.array([11,12,13,14,15])


file_name=path + "\\" + name
file= open(file_name, "w")


if (len(a) < len(b)):


    abcd=[]
    print(abcd)

    for i in range(0, len(a)):


        car=str(a[i])[0]
        for j in range(1, len(str(a[i]))):

            car= car + str(a[i])[j]

        abcd.append(car)


    for i in range(len(a), len(b)):

        abcd.append(" ")

    a=abcd


else:

    abcd=np.chararray(len(a), unicode=True)

    for i in range(0, len(b)):

        abcd[i]=str(b[i])

    for i in range(len(b), len(a)):

        abcd[i]="\t"

    b=abcd

'''
print(a[-1]+"ccccc")
for i in range(0, len(a)):

    print(a[i], b[i])


print("\t" + "a")

'''
print(a)

np.savetxt(file, np.c_[a, b], fmt="%s")

file.close()



'''
x, y= np.loadtxt(file_name, unpack=True)
print(x)
print(y)
'''