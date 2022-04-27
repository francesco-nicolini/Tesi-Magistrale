import numpy as np

import keyboard

a=0

while(1):


    if (keyboard.is_pressed("s")):
        print("\n\nletto\n\n")
        break

    a=a+1


print(a)


num=10

b=np.zeros(num)

for i in range(1, num):

    b[i]=abs(b[i-1] - 1)

print(b)

c=0

if(c):
    comand=print
    stringa="ciao"

else:
    comand=input
    stringa="-------"

comand(stringa)


