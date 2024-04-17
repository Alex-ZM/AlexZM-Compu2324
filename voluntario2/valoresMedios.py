import random
import numpy as np
import os
import time


# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 16    # Dimensión de la cuadrícula
T = 2    # Temperatura T = [0,5]
t = 10**6*N**2  # Tiempo

magnProm = np.array([])
E = np.array([])

# CREACIÓN DE LA MATRIZ DE ESPINES s
s = np.random.choice([-1,+1], size=(N,N)) 

ini = time.time()
for _ in range(0,t):

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
        p = 1   # Evaluación de p
    else:        # 
        p = pE    # 

    n = random.random()  # Número aleatorio entre 0 y 1

    if n<p:                 # Cambio del espín (i,j)
        s[i,j] = -s[i,j]    # 
    

    if t%(N*N*100)==0:  # Cada 100pmc calcula promedios:

        magnProm.append(np.sum(s)/N**2)

        Ei = 0
        for i in range(0,N):
            for j in range(0,N):
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

                    Ei = Ei-0.5*(s[i,j]+(s[up,j]+s[down,j]+s[i,left]+s[i,right]))
        E.append(Ei)

        
        
        
fin = time.time()
print("|| N = "+str(N)+"\n|| T = "+str(T)+" s\n|| "+str(t)+" iteraciones")
print("Tiempo de ejecución: "+str(fin-ini))
