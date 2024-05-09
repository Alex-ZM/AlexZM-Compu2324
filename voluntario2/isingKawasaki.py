
    ####################################################################
    ###  SIMULACIÓN DEL MODELO DE ISING CON LA DINÁMICA DE KAWASAKI  ###
    ####################################################################

import random
import numpy as np
import os
import time

###################################################################################################################

# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 32    # Dimensión de la cuadrícula
T =1    # Temperatura T = [0,5]
t = 200000  # Tiempo
skip = 200

# CREACIÓN DE LA MATRIZ DE ESPINES s
s = np.random.choice([-1,+1], size=(N+2,N)).astype(np.int8)
for j in range(N):
    s[0,j] = 1
    s[N+1,j] = -1


###################################################################################################################

def condContorno(i,j):
    if j==N-1:       # 
        left = j-1   # 
        right = 0    # 
    elif j==0:       # 
        left = N-1   # Periodicidad horizontal
        right = 1    # 
    else:            # 
        left = j-1   # 
        right = j+1  # 
    return left,right

## CÁLCULO DE LA ENERGÍA DEL SISTEMA EN EL INSTANTE INICIAL
#E = 0
#for i in range(N):
#    for j in range(N):
#        up,down,left,right = condContorno(i,j)
#        E = E-0.5*(s[i,j]+(s[up,j]+s[down,j]+s[i,left]+s[i,right]))


###################################################################################################################

wd = os.path.dirname(__file__)  # Directorio de trabajo
rd = "isingKawasaki_data.dat"           # Directorio relativo
fichero = open(os.path.join(wd,rd), "w")  

ini = time.time()
for w in range(t):

    i = random.randint(1,N-1)   # Se elijen dos partículas aleatorias vecinas
    j = random.randint(0,N-1)   # p1=(i,j) ; p2=(u,v)
    flip = np.random.choice([-1,1])
    if flip == -1: 
        u = i + 1
        v = j
    else: 
        u = i
        if j != N-1:
            v = j + 1 
        else:
            v = 0

    # Si las dos partículas tienen igual espín, no hacer nada
    # Si las dos partículas tienen diferente espín, calcular energía y probabilidad
    if s[i,j] != s[u,v]:  
        up1,up2 = i-1,i-1
        down1,down2 = i+1,i+1
        left1,right1 = condContorno(i,j)
        left2,right2 = condContorno(u,v)

        if flip == -1:
            deltaE = 2*s[i,j]*(s[i,right1]+s[i,left1]+s[up1,j]-s[u,right2]-s[u,left2]-s[down2,v]) # Pareja vertical
        else:
            deltaE = 2*s[i,j]*(s[up1,j]+s[down1,j]+s[i,left1]-s[up2,v]-s[down2,v]-s[u,right2]) # Pareja horizontal

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
        for i in range(1,N+1):
            fichero.write(' '.join(map(str, s[i])) + "\n")


fin = time.time()
print("\n|| Red "+str(N)+"x"+str(N)+"\n|| T = "+str(T)+"\n|| "+str(t)+" iteraciones (~"+f"{(t/N**2):.0f}"+" pmc)")
print("----> Tiempo de ejecución: "+f"{(fin-ini):.2f}"+" s\n")

fichero.close()
