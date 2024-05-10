
    #####################################################################################
    ### SIMULACIÓN CON DINÁMICA MOLECULAR DE UN GAS CON UN POTENCIAL DE LENNARD-JONES ###
    #####################################################################################

import numpy as np
import random
import time
import os
from numba import jit

###################################################################################################################

# Definimos algunas constantes
h = 0.002    
nIter = 3000     
nParticulas = 20
tEjecIni = time.time()
L = 10
margen = 0.5
skip = 10
hMedios = h/2
lMedios = L/2

# Definimos (y reescalamos) los parámetros iniciales de las partículas
m = np.full(nParticulas,1).astype(np.int8)
r = np.zeros((nParticulas,2))
v = np.zeros((nParticulas,2))
a = np.zeros((nParticulas,2))
w = np.zeros((nParticulas,2))
evR = np.zeros((nParticulas,2))
evV = np.zeros((nParticulas,2))
evA = np.zeros((nParticulas,2))
T = np.zeros((nParticulas,1))
V = np.zeros((nParticulas,1))

###################################################################################################################

# CONDICIONES INICIALES - POSICIONES EQUIESPACIADAS EN CUADRÍCULA
particulasPorFila = int(np.floor(np.sqrt(nParticulas)))
separacionInicialH = L/particulasPorFila
if nParticulas > particulasPorFila**2:
    separacionInicialV = L/(particulasPorFila+1)
else: separacionInicialV = separacionInicialH

for f in range(particulasPorFila):
    for i in range(particulasPorFila):  
        r[f*particulasPorFila+i,1] = i*separacionInicialH  # Equiespaciado horizontal
particulasUltimaFila = nParticulas-particulasPorFila**2
if particulasUltimaFila != 0:
    for i in range(particulasUltimaFila):
        r[particulasPorFila**2+i,1] = i*separacionInicialH  # Equiespaciado horizontal (última fila)

for c in range(particulasPorFila):
    for i in range(particulasPorFila):
        r[c*particulasPorFila+i,0] = c*separacionInicialV  # Equiespaciado vertical
if particulasUltimaFila != 0:
    for i in range(particulasUltimaFila):
        r[particulasPorFila**2+i,0] = particulasPorFila*separacionInicialV  # Equiespaciado vertical (última fila)

# CONDICIONES INICIALES - VELOCIDADES ALEATORIAS
for i in range(nParticulas):
    v[i,0] = 2*random.random()-1
    v[i,1] = 1-np.abs(v[i,0])

# CONDICIONES INICIALES - POSICIONES ALEATORIAS
for i in range(nParticulas):
    r[i,0] += (2*random.random()-1)*margen  # Desplaza aleatoriamente en la componente X (+-margen)
    r[i,1] += (2*random.random()-1)*margen  # Desplaza aleatoriamente en la componente Y (+-margen)

# CONDICIONES DE CONTORNO PERIÓDICO - RECOLOCACIÓN DE PARTÍCULAS 
@jit(nopython=True,fastmath=True)
def bordes(v,p,L):
    evX = v[p,0]
    evY = v[p,1]
    if evX > L:
        evX = evX%L
    if evX < 0:
        evX = evX%L
    if evY > L:
        evY = evY%L
    if evY < 0:
        evY = evY%L
    return np.array([evX,evY])

# CONDICIONES DE CONTORNO PERIÓDICO - DISTANCIA MÍNIMA ENTRE DOS PARTÍCULAS 
def distanciaToroide(vector,p,j,L,lMedios):
    distX = vector[p,0]-vector[j,0]
    distY = vector[p,1]-vector[j,1]
    if np.abs(distX) > lMedios:
        distX = -(L-np.abs(distX))*(distX/np.abs(distX))
    if np.abs(distY) > lMedios:
        distY = -(L-np.abs(distY))*(distY/np.abs(distY))
    return np.array([distX,distY])

# VERLET - ACELERACIÓN (t)
@jit(nopython=True,fastmath=True)
def aceleracion(r,nParticulas,L,lMedios):
    for p in range(nParticulas):
        aux1 = np.array([0.0,0.0])
        auxa = r
        for j in range(nParticulas):
            if p != j:
                R = distanciaToroide(r,p,j,L,lMedios)
                normaR = np.linalg.norm(R)
                aux1 = aux1 + (48/normaR**13 - 24/normaR**7)*R/normaR
        auxa[p] = aux1
    return auxa

# VERLET - W
def uvedoble(w,v,a):
    for p in range(nParticulas):
        w[p] = v[p]+hMedios*a[p]
    return w

# VERLET - POSICION (t+h)
def evPosicion(evR,r,w):
    for p in range(nParticulas):
        evR[p] = r[p]+h*w[p]
    return evR

# VERLET - VELOCIDAD (t+h)
def evVelocidad(w,evA):
    for p in range(nParticulas):
        evV[p] = w[p]+hMedios*evA[p]
    return evV


###################################################################################################################

# Ahora solo queda programar el bucle y guardar los resultados de cada iteración en el
# formato correcto y dentro de un fichero, para poder representarlos luego.
wd = os.path.dirname(__file__)  # Directorio de trabajo
datosPath = os.path.join(wd,"posParticulas.dat")  
TPath = os.path.join(wd,"energiaCinetica.dat")
VPath = os.path.join(wd,"energiaPotencial.dat")
EtotPath = os.path.join(wd,"energiaTotal.dat")
ficheroPlot = open(datosPath, "w")
ficheroT = open(TPath, "w")
ficheroV = open(VPath, "w")
ficheroEtot = open(EtotPath, "w")

# BUCLE - ALGORITMO DE VERLET
for t in range(nIter):

    if t%skip==0:  # Guarda la posición de las partículas cada "skip" iteraciones
        print(t)
        for p in range(nParticulas):
            ficheroPlot.write(str(r[p][0]) + ", " + str(r[p][1]) + "\n")  # Posiciones de las partículas -> Fichero
            #ficheroT.write(str(T[p])+"\n")  # Energía cinética -> Fichero
            #ficheroV.write(str(V[p])+"\n")  # Energía potencial -> Fichero
        ficheroPlot.write("\n") 

    for j in range(nParticulas):
            r[j] = bordes(r,j,L) 

    a = aceleracion(r,nParticulas,L,lMedios)
    w = uvedoble(w,v,a)
    evR = evPosicion(evR,r,w)
    for j in range(nParticulas):
            evR[j] = bordes(evR,j,L) 
    evA = aceleracion(evR,nParticulas,L,lMedios)
    evV = evVelocidad(w,evA)


    # EVOLUCIÓN TEMPORAL
    for p in range(nParticulas):
        r[p] = evR[p]
        v[p] = evV[p]

    # CÁLCULO DE ENERGÍAS
    if t%skip == 0:
        sumaT = 0
        sumaV = 0
        for p in range(nParticulas):

            # ENERGÍA CINÉTICA
            T[p] = 0.5*np.linalg.norm(v[p])**2
            sumaT += T[p]
            
            # ENERGÍA POTENCIAL
            V[p] = 0
            for j in range(nParticulas):
                if p != j:
                    R = np.linalg.norm(distanciaToroide(r,p,j,L,lMedios))
                    V[p] += (R**(-12)-R**(-6))
            V[p] = 2*V[p]
            sumaV += V[p]

        ficheroT.write(str(sumaT)+"\n")
        ficheroV.write(str(sumaV)+"\n")
        ficheroEtot.write(str(sumaT+sumaV)+"\n")
        
###################################################################################################################

# Por último, escribimos algunos datos de interés al final del fichero
ficheroPlot.write("# Se han realizado "+str(nIter)+" iteraciones con h = "+str(h)+" y skip "+str(skip)+"\n")
tEjecFin = time.time()
ficheroPlot.write("# Tiempo de ejecución: "+str(tEjecFin-tEjecIni))
