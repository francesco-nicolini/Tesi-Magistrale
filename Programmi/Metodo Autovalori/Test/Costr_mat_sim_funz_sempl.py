import numpy as np
import matplotlib.pyplot as plt


a_min= 10
a_max= 10**3

b_min= 10**(-8)
b_max= 10**(-6)


N= 400



a_min= np.log10(a_min)
a_max= np.log10(a_max)

b_min= np.log10(b_min)
b_max= np.log10(b_max)





a= np.logspace(a_min, a_max, N)
b= np.logspace(b_min, b_max, N)


print("Matrice di cui calcolo tutte le componenti:")


matrix= np.zeros( (N,N) )


for i in range(0, N):

    for j in range(0, N):

        matrix[i][j]= a[i]*b[j]


#print(matrix)


trasposta= matrix.transpose()



if ( (matrix==trasposta).all() ):
    print("\nLa matrice è simmetrica")

else:
    print("\nLa matrice non è simmetrica")


massimo= 0

for i in range(0, N):

    for j in range(0, N):

        diff= abs( (matrix[i][j] - trasposta[i][j])/matrix[i][j] )

        if ( diff>massimo ):
            massimo= diff

print("Il massimo della differenza relativa in valore assoluto tra la matrice e la sua trasposta è {:.2e}".format(massimo))


print("\n\n\nMatrice che costruisco simmetrica:")
matrix= np.zeros( (N,N) )

for i in range(0,N):

    matrix[i][i]= a[i]*b[i]

    for j in range(i+1,N):

        matrix[i][j]= a[i]*b[j]

        matrix[j][i]= matrix[i][j]


trasp= matrix.transpose()



if ( (matrix==trasp).all() ):
    print("\nLa matrice è simmetrica")

else:
    print("\nLa matrice non è simmetrica")




