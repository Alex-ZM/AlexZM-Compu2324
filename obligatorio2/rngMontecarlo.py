import random
import numpy as np
import os
import time
from numba import jit


# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 10    # Dimensión de la cuadrícula
T = 0.5    # Temperatura T = [0,5]
t = 50000  # Tiempo
skip = 1

# CREACIÓN DE LA MATRIZ DE ESPINES s
s = np.random.choice([-1,+1], size=(N,N)).astype(np.int8)


#@jit(nopython=True, fastmath=True, cache=True)
def bordes(N):
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
    return i,j,up,down,left,right


#@jit(nopython=True, fastmath=True, cache=True)
def funcDeltaE(i,j,up,down,left,right):
    deltaE = 2*s[i,j]*(s[up,j]+s[down,j]+s[i,left]+s[i,right])
    return deltaE


#@jit(nopython=True, fastmath=True, cache=True)
def flipCheck(deltaE):
    pE = np.exp(-deltaE/T)
    if 1<pE:     # 
        p = 1    # Evaluación de p
    else:        # 
        p = pE   # 
    n = random.random()  # Número aleatorio entre 0 y 1
    return n,p


wd = os.path.dirname(__file__)  # Directorio de trabajo
rd = "ising_data.dat"           # Directorio relativo
fichero = open(os.path.join(wd,rd), "w")  


ini = time.time()
for w in range(0,t):

    i,j,up,down,left,right = bordes(N)

    dE = funcDeltaE(i,j,up,down,left,right)

    n,p =flipCheck(dE)

    if n<p:                 # Cambio del espín (i,j)
        s[i,j] = -s[i,j]    # 
    
    # Se guarda el estado de la red en este instante
    if w%skip==0:
        fichero.write("\n")
        for i in range(0,N):
            fichero.write(' '.join(map(str, s[i])) + "\n")
        

fin = time.time()
print("|| N = "+str(N)+"\n|| T = "+str(T)+" s\n|| "+str(t)+" iteraciones")
print("Tiempo de ejecución: "+str(fin-ini))

fichero.close()
