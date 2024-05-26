
    ####################################################################
    ###  SIMULACIÓN DEL MODELO DE ISING CON LA DINÁMICA DE KAWASAKI  ###
    ####################################################################

import random
import numpy as np
import os
import time
from numba import jit

###################################################################################################################

# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 64    # Dimensión de la cuadrícula
T = 1    # Temperatura T = [0,5]
pmc = 400  # Número de pasos Montecarlo
t = pmc*N**2
skip = 40*pmc

# SELECCIÓN DE LAS CONDICIONES INICIALES
# Opciones: "magnNula", "magnAleatoria"
condIni = "magnNula"  

# CREACIÓN DE LA MATRIZ DE ESPINES s
if condIni == "magnAleatoria":
    s = np.random.choice([-1,+1], size=(N+2,N)).astype(np.int8)
    for j in range(N):
        s[0,j] = -1
        s[N+1,j] = 1
elif condIni == "magnNula":
    s = np.ones((N+2,N)).astype(np.int8)
    for _ in range(int(N**2/2)):
        i = np.random.randint(1,N)
        j = np.random.randint(0,N-1)
        while s[i,j] == -1:
            i = np.random.randint(1,N)
            j = np.random.randint(0,N-1)
        s[i,j] = -1
    for j in range(N):
        s[0,j] = -1
        s[N+1,j] = 1
###################################################################################################################

#@jit(nopython=True)
def condContorno(j,N):
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

# CÁLCULO DE LA MAGNETIZACIÓN DEL SISTEMA
magn = 0
for i in range(1,N+1):
    for j in range(N):
        magn += s[i,j]

def magnArriba(s,N):
    magnUp = 0
    for i in range(1,int((N+1)/2)):
        for j in range(N):
            magnUp += s[i,j]
    return magnUp/(N**2/2)

def magnAbajo(s,N):
    magnDown = 0
    for i in range(int((N+1)/2),N+1):
        for j in range(N):
            magnDown += s[i,j]
    return magnDown/(N**2/2)

## CÁLCULO DE LA ENERGÍA DEL SISTEMA EN EL INSTANTE INICIAL
#E = 0
#for i in range(N):
#    for j in range(N):
#        left,right = condContorno(j,N)
#        E = E-0.5*(s[i,j]+(s[up,j]+s[down,j]+s[i,left]+s[i,right]))

#@jit(nopython=True)
def isingKawasaki(flip,s,i,j,u,v,N,T):
    if s[i,j] != s[u,v]:  
        up1,up2 = i-1,i-1
        down1,down2 = i+1,i+1
        left1,right1 = condContorno(j,N)
        left2,right2 = condContorno(v,N)

        if flip == -1:
            deltaE = 2*s[i,j]*(s[i,right1]+s[i,left1]+s[up1,j]-s[u,right2]-s[u,left2]-s[down2,v]) # Pareja vertical
        else:
            deltaE = 2*s[i,j]*(s[up1,j]+s[down1,j]+s[i,left1]-s[up2,v]-s[down2,v]-s[u,right2]) # Pareja horizontal

        pE = np.exp(-deltaE/T)
        if 1<pE:     # 
            p = 1    # Evaluación de p
        else:        # 
            p = pE   # 
        n = np.random.rand()  # Número aleatorio entre 0 y 1

        # PERMUTACIÓN DE ESPINES s1,s2 = s2,s1
        if n<p:  
            s[i,j], s[u,v] = s[u,v], s[i,j]     


###################################################################################################################

wd = os.path.dirname(__file__)  # Directorio de trabajo
sd = "isingKawasaki_data.dat"          
magnd = "isingKawasaki_magn.dat"       
fichero = open(os.path.join(wd,sd), "w")  
ficheroMagn = open(os.path.join(wd,magnd), "w")  

ini = time.time()
for w in range(t):

    # Se elijen dos partículas aleatorias vecinas [ p1=(i,j) ; p2=(u,v) ]
    i = np.random.randint(1,N-1)
    j = np.random.randint(0,N-1)   
    flip = np.random.choice([-1,1])
    # Vecinas verticales
    if flip == -1:  
        u = i + 1
        v = j
    # Vecinas horizontales
    else:           
        u = i
        if j == N-1:
            v = 0
        else:
            v = j + 1   
    
    isingKawasaki(flip,s,i,j,u,v,N,T)

    # Se guarda el estado de la red en este instante
    if w%skip==0:
        fichero.write("\n")
        for i in range(1,N+1):
            fichero.write(' '.join(map(str, s[i])) + "\n")

        for i in range(1,N+1):
            ficheroMagn.write(str(magnArriba(s,N)) + "," + str(magnAbajo(s,N)) + "\n")


fin = time.time()
print("\n|| Red "+str(N)+"x"+str(N)+"\n|| T = "+str(T)+"\n|| "+str(t)+" iteraciones (~"+f"{(t/N**2):.0f}"+" pmc)")
print("----> Tiempo de ejecución: "+f"{(fin-ini):.2f}"+" s\n")

fichero.close()
