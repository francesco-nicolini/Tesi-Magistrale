import numpy as np
from numpy.linalg import svd




# Questa funzione stampa la matrice in modo che tutte i numeri dopo i segni siano allineati lungo la stessa colonna. a è la matrice da stampare, mentre cifr_sign è il numero di cifre significative delle componenti che si vuole visualizzare.


def print_matrix(a, cifr_sign):

    formato="{:." + str(cifr_sign) + "}"

    segni=[]

    for j in range(0, len(a[0])):

        flag=0

        for i in range(0, len(a)):

            if ( a[i][j]<0 ):
                segni.append("-")
                flag=1
                break

        if ( flag==0 ):
            segni.append("+")





    lung_max=[]

    for j in range(0, len(a)):

        massimo= 0

        for i in range(0, len(a)):

            if( a[i][j]==int(a[i][j]) ):
                b= len("{0}".format(a[i][j]))

            else:
                b= len(formato.format(a[i][j]))


            if( segni[j]=="-" and a[i][j]>0 ):
                b+=1


            if( massimo<b ):
                massimo= b

        lung_max.append(massimo)





    for i in range(0, len(a)):

        if( a[i][0]==int(a[i][0]) ):

            if( segni[0]=="-" and a[i][0]>0 ):
                print(" {0}".format(a[i][0]), end=" ")

            else:
                print("{0}".format(a[i][0]), end=" ")

        else:
            if( segni[0]=="-" and a[i][0]>0 ):
                print(" "+formato.format(a[i][0]), end=" ")

            else:
                print(formato.format(a[i][0]), end=" ")


        for j in range(1, len(a)):


            if( a[i][j]==int(a[i][j]) ):

                if( a[i][j-1]>0 and segni[j-1]=="-" ):
                    b= " {0}".format(a[i][j-1])

                else:
                    b= "{0}".format(a[i][j-1])

                diff= lung_max[j-1]-len(b)

                space=" "

                for k in range(0, diff):

                    space+=" "

                if( segni[j]=="-" and a[i][j]>0 ):
                    space+=" "

                print(space + "{0}".format(a[i][j]), end=" ")


            else:

                if( a[i][j-1]>0 and segni[j-1]=="-" ):
                    b= " "+formato.format(a[i][j-1])

                else:
                    b= formato.format(a[i][j-1])

                diff= lung_max[j-1]-len(b)

                space=" "

                for k in range(0, diff):

                    space+=" "

                if( segni[j]=="-" and a[i][j]>0 ):
                    space+=" "

                print(space + formato.format(a[i][j]), end=" ")


        print("")


cifr_sign= 3

a= np.array([[1,2,3],[4,5,6],[7,8,100]])



print_matrix(a, cifr_sign)
print("\n\n")



C, s, V= svd(a)


print("Valori Singolari:")

for i in range(0, len(s)):
    print("{:.2}".format(s[i]))

print("\n\n")


print("Matrice C:")
print_matrix(C, cifr_sign)
print("\n\n")



S=np.diag(s)

print("Matrice S:")
print_matrix(S, cifr_sign)
print("\n\n")





print("Matrice V:")
print_matrix(V, cifr_sign)
print("\n\n")






a_1= np.dot(C,S)
a_1= np.dot(a_1,V)


print("Matrice a_1:")
print_matrix(a_1, cifr_sign)
print("\n\n")

