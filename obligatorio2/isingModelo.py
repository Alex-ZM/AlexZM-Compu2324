import random
import numpy as np
import os
import time


# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 120    # Dimensión de la cuadrícula
T = 2.828    # Temperatura T = [0,5]
t = 20000000  # Tiempo
skip = 10000

# CREACIÓN DE LA MATRIZ DE ESPINES s
s = np.random.choice([-1,+1], size=(N,N)).astype(np.int8)


def monteCarlo(N):
    i = random.randint(0,N-1)  # Se elije una partícula aleatoria
    j = random.randint(0,N-1)  #
    if i==N-1:      #
        up = i-1    #
        down = 0    #
    elif i==0:      #
        up = N-1    # Periodicidad vertical
        down = 1    #
    else:           #
        up = i-1    #
        down = i+1  #
    if j==N-1:       # 
        left = j-1   # 
        right = 0    # 
    elif j==0:       # 
        left = N-1   # Periodicidad horizontal
        right = 1    # 
    else:            # 
        left = j-1   # 
        right = j+1  # 
    deltaE = 2*s[i,j]*(s[up,j]+s[down,j]+s[i,left]+s[i,right])
    pE = np.exp(-deltaE/T)
    if 1<pE:     # 
        p = 1    # Evaluación de p
    else:        # 
        p = pE   # 
    n = random.random()  # Número aleatorio entre 0 y 1
    return n,p,i,j


wd = os.path.dirname(__file__)  # Directorio de trabajo
rd = "ising_data.dat"           # Directorio relativo
fichero = open(os.path.join(wd,rd), "w")  


ini = time.time()
for w in range(t):

    n,p,i,j = monteCarlo(N)

    if n<p:                 # Cambio del espín (i,j)
        s[i,j] = -s[i,j]    # 
    
    # Se guarda el estado de la red en este instante
    if w%skip==0:
        fichero.write("\n")
        for i in range(N):
            fichero.write(' '.join(map(str, s[i])) + "\n")
        

fin = time.time()
print("\n|| Red "+str(N)+"x"+str(N)+"\n|| T = "+str(T)+"\n|| "+str(t)+" iteraciones (~"+f"{(t/N**2):.0f}"+" pmc)")
print("----> Tiempo de ejecución: "+f"{(fin-ini):.2f}"+" s\n")

fichero.close()