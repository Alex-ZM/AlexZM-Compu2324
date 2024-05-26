
    ####################################################################
    ###  SIMULACIÓN DEL MODELO DE ISING CON LA DINÁMICA DE KAWASAKI  ###
    ####################################################################

import numpy as np
import os
import time
from numba import jit

###################################################################################################################

# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 32    # Dimensión de la cuadrícula
T = 1    # Temperatura T = [0,5]
pmc = 5000  # Número de pasos Monte Carlo
t = pmc*N**2
skip = 40*pmc

# SELECCIÓN DE LAS CONDICIONES INICIALES
# Opciones: "magnNula", "magnAleatoria"
condIni = "magnNula"  

###################################################################################################################
# CREACIÓN DE LA MATRIZ DE ESPINES s
if condIni == "magnAleatoria":
    s = np.random.choice([-1,+1], size=(N,N)).astype(np.int8)
    for j in range(N):
        s[0,j] = -1
        s[N-1,j] = 1
elif condIni == "magnNula":
    s = np.ones((N,N)).astype(np.int8)
    for _ in range(int(N**2/2)-N):
        i = np.random.randint(1,N-1)
        j = np.random.randint(0,N-1)
        while s[i,j] == -1:
            i = np.random.randint(1,N-1)
            j = np.random.randint(0,N-1)
        s[i,j] = -1
    for j in range(N):
        s[0,j] = -1
        s[N-1,j] = 1
###################################################################################################################

# CONDICIONES DE CONTORNO - BORDES VERTICALES CON ESPÍN FIJO
@jit(nopython=True)
def condContornoV(i,N):
    if i==N-1:      #
        up = i-1    #
        down = i    #
    elif i==0:      # 
        up = i      # Condiciones de contorno verticales
        down = i+1  # 
    else:           # 
        up = i-1    # 
        down = i+1  # 
    return up,down

# CONDICIONES DE CONTORNO - PERIODICIDAD HORIZONTAL (CILINDRO)
@jit(nopython=True)
def condContornoH(j,N):
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
#magn = 0
#for i in range(N):
#    for j in range(N):
#        magn += s[i,j]

# MAGNETIZACIÓN PROMEDIO DE LA MITAD SUPERIOR DEL SISTEMA
def magnArriba(s,N):
    magnUp = 0
    for i in range(int((N-1)/2)):
        for j in range(N):
            magnUp += s[i,j]
    return magnUp/(N**2/2)

# MAGNETIZACIÓN PROMEDIO DE LA MITAD INFERIOR DEL SISTEMA
def magnAbajo(s,N):
    magnDown = 0
    for i in range(int((N-1)/2),N):
        for j in range(N):
            magnDown += s[i,j]
    return magnDown/(N**2/2)


# ALGORITMO DE METROPOLIS - ISING KAWASAKI
@jit(nopython=True)
def isingKawasaki(flip,s,i,j,u,v,N,T):
    if s[i,j] != s[u,v]:  
        up1,down1 = condContornoV(i,N)
        up2,down2 = condContornoV(u,N)
        left1,right1 = condContornoH(j,N)
        left2,right2 = condContornoH(v,N)

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

        # Permutación de espines s1,s2 = s2,s1
        if n<p:  
            s[i,j], s[u,v] = s[u,v], s[i,j]    


## CÁLCULO DE LA ENERGÍA DEL SISTEMA EN EL INSTANTE INICIAL
#E = 0
#for i in range(N):
#    for j in range(N):
#        left,right = condContornoH(j,N)
#        E = E-0.5*(s[i,j]+(s[up,j]+s[down,j]+s[i,left]+s[i,right]))

# ENERGÍA MEDIA POR PARTÍCULA EN EL INSTANTE t
@jit(nopython=True)
def energiaMediaPorParticula(s,N):
    E = 0
    for i in range(N):
        for j in range(N):
            up,down = condContornoV(i,N)
            left,right = condContornoH(j,N)
            E += s[i,j]*(s[up,j]+s[down,j]+s[i,left]+s[i,right])
    return -0.5*E/N**2


###################################################################################################################

wd = os.path.dirname(__file__)  # Directorio de trabajo
dirDatos = "isingKawasaki_data.dat"          
dirMag = "isingKawasaki_magn.dat"      
dirE = "isingKawasaki_energiapp" 
fichero = open(os.path.join(wd,dirDatos), "w")  
ficheroMagn = open(os.path.join(wd,dirMag), "w")  
ficheroE = open(os.path.join(wd,dirE), "w")  

ini = time.time()
for w in range(t):

    # Se elijen dos partículas aleatorias vecinas [ p1=(i,j) ; p2=(u,v) ]
    i = np.random.randint(0,N-2)
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
        for i in range(N):
            fichero.write(' '.join(map(str, s[i])) + "\n")

        for i in range(N):
            ficheroMagn.write(str(magnArriba(s,N)) + "," + str(magnAbajo(s,N)) + "\n")

        ficheroE.write(str(energiaMediaPorParticula(s,N)) + "\n")

fin = time.time()
print("\n|| Red "+str(N)+"x"+str(N)+"\n|| T = "+str(T)+"\n|| "+str(t)+" iteraciones (~"+f"{(t/N**2):.0f}"+" pmc)")
print("----> Tiempo de ejecución: "+f"{(fin-ini):.2f}"+" s\n")

fichero.close()
