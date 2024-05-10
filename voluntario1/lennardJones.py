
    #####################################################################################
    ### SIMULACIÓN CON DINÁMICA MOLECULAR DE UN GAS CON UN POTENCIAL DE LENNARD-JONES ###
    #####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import random
import time
import os
from numba import jit

###################################################################################################################

# Definimos algunas constantes
h = 0.002    
nIteraciones = 50000     
nParticulas = 20
L = 10
lMedios = L/2
margen = 0.5
skip = 10

eCinetica = []
ePotencial = []
eTotal = []

# CONDICIONES DE CONTORNO PERIÓDICO - DISTANCIA MÍNIMA ENTRE DOS PARTÍCULAS 
def distanciaToroideGlobal(vector,t,p,j,L,lMedios):
    distX = vector[t,p,0]-vector[t,j,0]
    distY = vector[t,p,1]-vector[t,j,1]
    if np.abs(distX) > lMedios:
        distX = -(L-np.abs(distX))*(distX/np.abs(distX))
    if np.abs(distY) > lMedios:
        distY = -(L-np.abs(distY))*(distY/np.abs(distY))
    return np.array([distX,distY])


@jit(nopython=True,fastmath=True)
def verlet(h,nIteraciones,nParticulas,skip,L,margen):
    # PARÁMETROS INICIALES DEL SISTEMA
    r = np.zeros((nParticulas,2))
    v = np.zeros((nParticulas,2))
    a = np.zeros((nParticulas,2))
    w = np.zeros((nParticulas,2))
    evR = np.zeros((nParticulas,2))
    evV = np.zeros((nParticulas,2))
    evA = np.zeros((nParticulas,2))
    posiciones = np.zeros((int(nIteraciones/skip),nParticulas,2))
    velocidades = np.zeros((int(nIteraciones/skip),nParticulas,2))

    hMedios = h/2
    lMedios = L/2


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


    for t in range(nIteraciones):

        if t%skip == 0:
            for p in range(nParticulas):
                posiciones[int(t/skip),p] = r[p]
                velocidades[int(t/skip),p] = v[p]
    

        for j in range(nParticulas):
            r[j] = bordes(r,j,L)

        for p in range(nParticulas):
            aux1 = np.array([0.0,0.0])
            for j in range(nParticulas):
                if p != j:
                    R = distanciaToroide(r,p,j,L,lMedios)
                    normaR = np.linalg.norm(R)
                    aux1 = aux1 + (48/normaR**13 - 24/normaR**7)*R/normaR
            a[p] = aux1

        for p in range(nParticulas):
            w[p] = v[p]+hMedios*a[p]

        for p in range(nParticulas):
            evR[p] = r[p]+h*w[p]

        for j in range(nParticulas):
            evR[j] = bordes(evR,j,L) 

        for p in range(nParticulas):
            aux2 = np.array([0.0,0.0])
            for j in range(nParticulas):
                if p != j:
                    R = distanciaToroide(evR,p,j,L,lMedios)
                    normaR = np.linalg.norm(R)
                    aux2 = aux2 + (48/normaR**13 - 24/normaR**7)*R/normaR
            evA[p] = aux2

        for p in range(nParticulas):
            evV[p] = w[p]+hMedios*evA[p]

        # EVOLUCIÓN TEMPORAL
        for p in range(nParticulas):
            r[p] = evR[p]
            v[p] = evV[p]

    return posiciones,velocidades

###################################################################################################################

wd = os.path.dirname(__file__)  # Directorio de trabajo
datosPath = os.path.join(wd,"posParticulas.dat")  
ficheroPlot = open(datosPath, "w")

# CÁLCULO DE POSICIONES, VELOCIDADES Y PERÍODOS MEDIANTE LA FUNCIÓN "verlet()"
tEjecIni = time.time()
r,v = verlet(h,nIteraciones,nParticulas,skip,L,margen)
tEjecFin = time.time()


# ESCRITURA DE DATOS EN EL FICHERO
tFicherosIni = time.time()
for t in range(int(nIteraciones/skip-1)):

    # CÁLCULO ENERGÍAS
    sumaT = 0
    sumaV = 0
    for p in range(nParticulas):  

        sumaT += 0.5*np.linalg.norm(v[t,p])**2

        ePotencialAux = 0
        for j in range(nParticulas):
            if p != j:
                R = np.linalg.norm(distanciaToroideGlobal(r,t,p,j,L,lMedios))
                ePotencialAux += (R**(-12)-R**(-6))
        sumaV += 2*ePotencialAux

    # ESCRITURA EN FICHERO
    eCinetica.append(sumaT)
    ePotencial.append(sumaV)
    eTotal.append(sumaT+sumaV)
    for p in range(nParticulas):
        ficheroPlot.write(str(r[t,p,0]) + ", " + str(r[t,p,1]) + "\n")
    ficheroPlot.write("\n") 

plt.plot(eCinetica, label= "T")
plt.plot(ePotencial, label="V")
plt.plot(eTotal, label="E = T+V")
plt.legend()
plt.show()

# Por último, escribimos algunos datos de interés al final del fichero
ficheroPlot.write("# Se han realizado "+str(nIteraciones)+" iteraciones con h = "+str(h)+" y skip "+str(skip)+"\n")
tEjecFin = time.time()
ficheroPlot.write("# Tiempo de ejecución: "+str(tEjecFin-tEjecIni))

tFicherosFin = time.time()
