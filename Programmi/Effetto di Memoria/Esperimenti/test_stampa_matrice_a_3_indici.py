import numpy as np
import pandas as pd



path= "C:\\Users\\39366\\Dropbox\\PC\\Documents\\GitHub\\Tesi-Magistrale\\Programmi\\Effetto di Memoria\\Esperimenti\\file_txt"
name= "stampa_matrice_3_indici"


a= np.zeros( [2,2,2], dtype=complex )

a[0,0,0]= 1+1j
a=abs(a)
print(a)

file_name= path + "\\" + name

df = pd.DataFrame( a.reshape((-1, a.shape[-1])), index= pd.MultiIndex.from_product( [range(a.shape[0]), range(a.shape[1])] ) ).to_csv(file_name)





df = pd.read_csv(file_name, header=None, index_col=[0,1])
A= np.array(df.agg(list, 1).groupby(level=0).agg(list).tolist())

print("\n\n\n")

print(A)
print(type(A))