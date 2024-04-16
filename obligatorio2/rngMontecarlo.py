import random
import numpy as np


# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 30    # Dimensión de la cuadrícula
T = 1.0    # Temperatura T = [0,5]
t = 1000  # Tiempo

# CREACIÓN DE LA MATRIZ DE ESPINES s
s = np.random.choice([-1,+1], size=(N,N))
print(s)

fichero = open("ising_data.dat", "w")

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

    if 1<np.exp(-deltaE/T):     # 
        p = np.exp(-deltaE/T)   # Evaluación de p
    else:                       # 
        p = 1                   # 

    n = random.random()  # Número aleatorio entre 0 y 1

    if n<p:                 # Cambio del espín (i,j)
        s[i,j] = -s[i,j]    # 
    
    # Se guarda el estado de la red en este instante
    for _ in range(0,N):
        for i in range(0,N):
            fichero.write(' '.join(map(str, s[i])) + "\n")
        fichero.write("\n")

fichero.close()
