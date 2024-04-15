import random
import numpy as np
import time


# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 5  # Dimensión de la cuadrícula
T = 5  # Temperatura T = [0,5]
t = 10  # Tiempo

# CREACIÓN DE LA MATRIZ DE ESPINES s
s = np.random.choice([-1,+1], size=(N,N))
print(s)


for _ in range(0,t):


    for n in range(0,N):

        if n==N:        #
            up = n-1    #
            down = 0    #
        elif n==0:      #
            up = N      # Periodicidad vertical
            down = n+1  #
        else:           #
            up = n-1    #
            down = n+1  #


        for m in range(0,N):

            if m==N:         # 
                left = m-1   # 
                right = 0    # 
            elif m==0:       # 
                left = N     # Periodicidad horizontal
                right = m+1  # 
            else:            # 
                left = m-1   # 
                right = m+1  # 

            deltaE = 2*s[n,m]*(s[up,m]+s[down,m]+s[n,left]+s[n,right])

            if 1<np.exp(-deltaE/T):     # 
                p = np.exp(-deltaE/T)   # Evaluación de p
            else:                       # 
                p = 1                   # 

            a = random.random()  # Número aleatorio entre 0 y 1

            if a<p:                 # Cambio del espín (n,m)
                s[n,m] = -s[n,m]    # 



# print(str(s[0,0]))