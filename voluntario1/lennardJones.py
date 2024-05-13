
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
skip = 10
nIteraciones = 70000     # 70000
nParticulas = 20         # 20
L = 10
lMedios = L/2
margen = 0.05

reposo = False
soloDesplHoriz = False
moduloVelocidad = 1

eCinetica = []
ePotencial = []
eTotal = []
normaV = np.zeros((int(nIteraciones/skip),nParticulas),dtype=np.float32)

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
def verlet(h,nIteraciones,nParticulas,skip,L,margen,reposo,moduloVelocidad,soloDesplHoriz):

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
    fuerzaParedes = np.zeros((int(nIteraciones/skip)))
    fuera = np.zeros((nParticulas,2))

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
    if reposo == False:
        if soloDesplHoriz == True:
            for p in range(nParticulas):
                v[p,0] = moduloVelocidad*random.random()
                v[p,1] = 0
        else:
            for p in range(nParticulas):
                v[p,0] = 2*np.random.rand()-1
                v[p,1] = np.random.choice(np.array([-1,1]))*np.sqrt(1-v[p,0]**2)

    # CONDICIONES INICIALES - POSICIONES ALEATORIAS
    for i in range(nParticulas):
        r[i,0] += (2*random.random()-1)*margen  # Desplaza aleatoriamente en la componente X (+-margen)
        r[i,1] += (2*random.random()-1)*margen  # Desplaza aleatoriamente en la componente Y (+-margen)


    # CONDICIONES DE CONTORNO PERIÓDICO - RECOLOCACIÓN DE PARTÍCULAS 
    def bordes(v,p,L):
        evX = v[p,0]
        evY = v[p,1]
        fueraX = 0
        fueraY = 0
        if evX > L:
            evX = evX%L
            fueraX = 1
        if evX < 0:
            evX = evX%L
            fueraX = -1
        if evY > L:
            evY = evY%L
            fueraY = 1
        if evY < 0:
            evY = evY%L
            fueraY = -1
        return np.array([evX,evY]),np.array([fueraX,fueraY])

    # CONDICIONES DE CONTORNO PERIÓDICO - DISTANCIA MÍNIMA ENTRE DOS PARTÍCULAS 
    def distanciaToroide(vector,p,j,L,lMedios):
        distX = vector[p,0]-vector[j,0]
        distY = vector[p,1]-vector[j,1]
        if np.abs(distX) > lMedios:
            distX = -(L-np.abs(distX))*(distX/np.abs(distX))
        if np.abs(distY) > lMedios:
            distY = -(L-np.abs(distY))*(distY/np.abs(distY))
        return np.array([distX,distY])

    for p in range(nParticulas):
                posiciones[0,p] = r[p]
                velocidades[0,p] = v[p]


    # BUCLE - ALGORITMO DE VERLET
    for t in range(1,nIteraciones):

        # GUARDAR EN LOS VECTORES RESULTADO
        if t%skip == 0:
            for p in range(nParticulas):
                posiciones[int(t/skip),p] = r[p]
                velocidades[int(t/skip),p] = v[p]
                fuerzaParedes[int(t/skip)] += (fuera[p,0]*v[p,0]+fuera[p,1]*v[p,1])/L


        # ALGORITMO DE VERLET + COMPROBACIÓN DE CONDICIONES DE CONTORNO

        fuera = np.zeros((nParticulas,2))

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
            evR[j],fuera[j] = bordes(evR,j,L)

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

    return posiciones,velocidades,fuerzaParedes

###################################################################################################################

wd = os.path.dirname(__file__)  # Directorio de trabajo
datosPath = os.path.join(wd,"posParticulas.dat")
ficheroPlot = open(datosPath, "w")

# CÁLCULO DE POSICIONES, VELOCIDADES Y PERÍODOS MEDIANTE LA FUNCIÓN "verlet()"
tEjecIni = time.time()
r,v,fuerzaParedes = verlet(h,nIteraciones,nParticulas,skip,L,margen,reposo,moduloVelocidad,soloDesplHoriz)
tEjecFin = time.time()


# ESCRITURA DE DATOS EN EL FICHERO
temperatura = 0
tFicherosIni = time.time()
for t in range(int(nIteraciones/skip-1)):

    # CÁLCULO ENERGÍAS
    sumaVelocidades = 0
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

        normaV[t,p] = np.linalg.norm(v[t,p])


        temperatura += v[t,p,0]**2 + v[t,p,1]**2
    

    # ESCRITURA EN FICHERO
    eCinetica.append(sumaT)
    ePotencial.append(sumaV)
    eTotal.append(sumaT+sumaV)
    for p in range(nParticulas):
        ficheroPlot.write(str(r[t,p,0]) + ", " + str(r[t,p,1]) + "\n")
    ficheroPlot.write("\n") 

temperatura = 2*temperatura/int(nIteraciones/skip)

# CÁLCULO DE LAS VELOCIDADES PARA EL HISTOGRAMA
promedioVelocidades = np.zeros(nParticulas)
for t in range(int(nIteraciones/(2*skip)),int(nIteraciones/skip-1)):
    for p in range(nParticulas):
        promedioVelocidades[p] += np.linalg.norm(v[t,p])
for p in range(nParticulas):
    promedioVelocidades[p] = promedioVelocidades[p]/(nIteraciones/(2*skip))


# CÁLCULO DE LA PRESIÓN                                      !!!!!!!!!!!!! REVISAR !!!!!!!!!!!!!!!!
presion = np.zeros(int(nIteraciones*h/skip))                #!!!!!!!!!!!!! REVISAR !!!!!!!!!!!!!!!!
for i in range(int(nIteraciones*h/skip)):                   #!!!!!!!!!!!!! REVISAR !!!!!!!!!!!!!!!!
    for t in range(i*int(1/h),int(1/h)*(i+1)):
        presion[i] += fuerzaParedes[t]


# REPRESENTACIÓN GRÁFICA DE LOS RESULTADOS
ax1 = plt.subplot(2,2,1)
plt.hist(normaV[0])
plt.title("Velocidades iniciales",fontsize=11)
plt.xlabel("|v|",fontsize=9)
plt.ylabel("Frecuencia",fontsize=9)

ax2 = plt.subplot(2,2,2)
plt.hist(promedioVelocidades,bins=15)
plt.title("Promedio de velocidades: t=t_f/2 a t=t_f",fontsize=11)
plt.xlabel("|v|",fontsize=9)
plt.ylabel("Frecuencia",fontsize=9)

ax4 = plt.subplot(2,1,2)
plt.plot(eCinetica, label= "T")
plt.plot(ePotencial, label="V")
plt.plot(eTotal, label="E = T+V")
plt.title("Energías en función del tiempo")
plt.xlabel("tiempo")
plt.ylabel("Energía")
plt.legend()

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.5)
plt.show()

plt.plot(presion)
plt.title("Presión en función del tiempo")
plt.xlabel("t")
plt.ylabel("Presión")
plt.show()

# Por último, escribimos algunos datos de interés al final del fichero
ficheroPlot.write("# Se han realizado "+str(nIteraciones)+" iteraciones con h = "+str(h)+", "+str(nParticulas)+" partículas y skip "+str(skip)+"\n")
tEjecFin = time.time()
ficheroPlot.write("# Tiempo de ejecucion: "+str(tEjecFin-tEjecIni)+"\n")
ficheroPlot.write("# Temperatura: "+str(temperatura)+"\n")

tFicherosFin = time.time()
