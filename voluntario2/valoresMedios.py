import random
import numpy as np
import os
import time


# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 4           # Dimensión de la cuadrícula
T = 3.0         # Temperatura T = [0,5]
t = 10**6*N**2  # Tiempo

magnProm = []  # Magnetización promedio de cada macroestado
E = []  # Energía de cada macroestado
ECuad = []

# CREACIÓN DE LA MATRIZ DE ESPINES s
s = np.ones((N,N),dtype=int) 

cienPMC = 100*N**2

ini = time.time()
for v in range(1,10002):

    for _ in range(0,cienPMC):
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

        if n<p:                 # Cambio del espín (i,j)
            s[i,j] = -s[i,j]    # 
        

    # Cada 100pmc calcula promedios:
    magnProm.append(np.sum(s)/N**2)

    E_i = 0
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

                E_i = E_i-0.5*(s[i,j]+(s[up,j]+s[down,j]+s[i,left]+s[i,right]))
    E.append(E_i)


for i in range(len(E)):
    ECuad.append(E[i]**2)

wd = os.path.dirname(__file__)          # Directorio de trabajo
rd = "valoresMedios_data.dat"           # Directorio relativo

# Guardamos los resultados en el fichero
fichero = open(os.path.join(wd,rd), "a")  
fichero.write("|| Magnetización promedio total: " + str(np.mean(magnProm)))
fichero.write("\n|| La energía media es: " + str(np.mean(E)/(2*N)))
fichero.write("\n|| El calor específico es: " + str((np.mean(ECuad)-np.mean(E)**2)/(T*N**2)))
fin = time.time()
fichero.write("\n|| Red "+str(N)+"x"+str(N)+"\n|| T = "+str(T)+" K\n|| "+str(t)+" iteraciones")
fichero.write("\n|| Tiempo de ejecución: "+str(fin-ini)+"\n\n")
fichero.close()
