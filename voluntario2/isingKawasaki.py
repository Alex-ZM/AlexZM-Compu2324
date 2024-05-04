    ##################################################################
    ### SIMULACIÓN DEL MODELO DE ISING CON LA DINÁMICA DE KAWASAKI ###
    ##################################################################

import random
import numpy as np
import os
import time

###################################################################################################################

# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 12    # Dimensión de la cuadrícula
T = 2.828    # Temperatura T = [0,5]
t = 20000  # Tiempo
skip = 10

# CREACIÓN DE LA MATRIZ DE ESPINES s
s = np.random.choice([-1,+1], size=(N,N)).astype(np.int8)

###################################################################################################################

def condContorno(i,j):
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
    return up,down,left,right

## CÁLCULO DE LA ENERGÍA DEL SISTEMA EN EL INSTANTE INICIAL
#E = 0
#for i in range(N):
#    for j in range(N):
#        up,down,left,right = condContorno(i,j)
#        E = E-0.5*(s[i,j]+(s[up,j]+s[down,j]+s[i,left]+s[i,right]))


###################################################################################################################

wd = os.path.dirname(__file__)  # Directorio de trabajo
rd = "ising_data.dat"           # Directorio relativo
fichero = open(os.path.join(wd,rd), "w")  

ini = time.time()
for w in range(t):

    i = random.randint(0,N-1)   # Se elijen dos partículas aleatorias vecinas
    j = random.randint(0,N-1)   # p1=(i,j) ; p2=(u,v)
    u = i + random.choice(-1,1) #
    v = j + random.choice(-1,1) #
    
    up1,down1,left1,right1 = condContorno(i,j)
    up2,down2,left2,right2 = condContorno(u,v)

    deltaE = 2*s[i,j]*(s[up1,j]+s[down1,j]+s[i,left1]+s[i,right1])###########  RECALCULAR NUEVA ENERGÍA Y CAMBIAR ESTA
    pE = np.exp(-deltaE/T)
    if 1<pE:     # 
        p = 1    # Evaluación de p
    else:        # 
        p = pE   # 
    n = random.random()  # Número aleatorio entre 0 y 1

    # PERMUTACIÓN DE ESPINES s1,s2 = s2,s1
    if n<p:  
        s[i,j], s[u,v] = s[u,v], s[i,j]     
    
    # Se guarda el estado de la red en este instante
    if w%skip==0:
        fichero.write("\n")
        for i in range(N):
            fichero.write(' '.join(map(str, s[i])) + "\n")


fin = time.time()
print("\n|| Red "+str(N)+"x"+str(N)+"\n|| T = "+str(T)+"\n|| "+str(t)+" iteraciones (~"+f"{(t/N**2):.0f}"+" pmc)")
print("----> Tiempo de ejecución: "+f"{(fin-ini):.2f}"+" s\n")

fichero.close()
