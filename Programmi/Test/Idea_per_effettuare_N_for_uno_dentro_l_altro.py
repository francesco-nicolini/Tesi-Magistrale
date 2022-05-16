import numpy as np

K= 3
a= 2034

print( [int(x) for x in str(a)])

print(10*[0])

numero= 0


while(numero<2**(K)):

    lista_1= [ int(x) for x in str(bin(numero))[2:] ]

    lung= len(lista_1)

    if( K>lung ):

        for i in range(0, len(lista_1)):

            if(not lista_1[i]):
                lista_1[i]= -1

        lista_2=(K-lung)*[-1]

        lista=lista_2 + lista_1

    else:

        lista= lista_1

        for i in range(0, len(lista)):

            if(not lista[i]):
                lista[i]= -1

    print(lista)
    numero= numero + 1