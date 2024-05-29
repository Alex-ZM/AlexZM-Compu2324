
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
h = 0.002                # 0.002
skip = 10                # 10
nIteraciones = 100000     # 70000
nParticulas = 16
L = 4
lMedios = L/2
margen = 0.05

#-------------------------------------------------------#
# Configuración de las condiciones iniciales (CAMBIAR)
reposo = True
soloDesplHoriz = False
moduloVelocidad = 1
redHexagonal24 = False

# Configuración para los apartados 6 y 7
acelerarRapido = False  # Apartado 6
acelerarLento = True    # Apartado 7
#-------------------------------------------------------#

# Configuración para observar el comportamiento de la red hexagonal
if redHexagonal24:
    h = 0.0005
    nIteraciones = 20000
    nParticulas = 24 
    L = 4
    lMedios = L/2
    skip = 8  

eCinetica = []
ePotencial = []
eTotal = []
normaV = np.zeros((int(nIteraciones/skip),nParticulas),dtype=np.float32)

###################################################################################################################

# CONDICIONES DE CONTORNO PERIÓDICO - DISTANCIA MÍNIMA ENTRE DOS PARTÍCULAS 
def distanciaGlobal(vector,t1,t2,p,j,L,lMedios):
    distX = vector[t2,p,0]-vector[t1,j,0]
    distY = vector[t2,p,1]-vector[t1,j,1]
    if np.abs(distX) > lMedios:
        distX = -(L-np.abs(distX))*(distX/np.abs(distX))
    if np.abs(distY) > lMedios:
        distY = -(L-np.abs(distY))*(distY/np.abs(distY))
    return np.array([distX,distY])


# BUCLE DE VERLET - POTENCIAL LENNARD-JONES
@jit(nopython=True,fastmath=True)
def verlet(h,nIteraciones,nParticulas,skip,L,margen,reposo,moduloVelocidad,soloDesplHoriz,redHexagonal24,acelerarRapido,acelerarLento):

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
    fuerzaParedes = np.zeros(nIteraciones)
    fuera = np.zeros((nParticulas,2))

    hMedios = h/2
    lMedios = L/2


    # CONDICIONES INICIALES - POSICIONES EQUIESPACIADAS EN CUADRÍCULA
    if redHexagonal24 == False:
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
    
    # CONDICIONES INICIALES - DISTRIBUCIÓN HEXAGONAL HOMOGÉNEA
    else:
        huecosPorFila = 6
        espaciado = L/huecosPorFila
        y = 0
        for f in range(huecosPorFila):
            if f%2==0:
                r[4*f] = np.array([0.0,y])
                r[4*f+1] = np.array([2*espaciado,y])
                r[4*f+2] = np.array([3*espaciado,y])
                r[4*f+3] = np.array([5*espaciado,y])
            else:
                r[4*f] = np.array([0.5*espaciado,y])
                r[4*f+1] = np.array([1.5*espaciado,y])
                r[4*f+2] = np.array([3.5*espaciado,y])
                r[4*f+3] = np.array([4.5*espaciado,y])
            y += espaciado

                        
    # CONDICIONES INICIALES - VELOCIDADES ALEATORIAS
    if reposo == False:
        if soloDesplHoriz == True:
            for p in range(nParticulas):
                v[p,0] = moduloVelocidad*random.random()
                v[p,1] = 0
        else:
            for p in range(nParticulas):
                v[p,0] = moduloVelocidad*(2*np.random.rand()-1)
                v[p,1] = np.random.choice(np.array([-1,1]))*np.sqrt(moduloVelocidad**2-v[p,0]**2)

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
    def distancia(vector,p,j,L,lMedios):
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
        for p in range(nParticulas):
            fuerzaParedes[int(t)] += (fuera[p,0]*v[p,0]+fuera[p,1]*v[p,1])/L

        # ACELERAR LAS PARTÍCULAS CADA CIERTO TIEMPO
        if acelerarRapido:
            if t*h==20 or t*h==30 or t*h==35 or t*h==45:
                for p in range(nParticulas):
                    v[p] *= 1.5
        elif acelerarLento:
            if (t*h)%60==0:
                for p in range(nParticulas):
                    v[p] *= 1.1

        # ALGORITMO DE VERLET + COMPROBACIÓN DE CONDICIONES DE CONTORNO

        fuera = np.zeros((nParticulas,2))

        for p in range(nParticulas):
            aux1 = np.array([0.0,0.0])
            for j in range(nParticulas):
                if p != j:
                    R = distancia(r,p,j,L,lMedios)
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
                    R = distancia(evR,p,j,L,lMedios)
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
r,v,fuerzaParedes = verlet(h,nIteraciones,nParticulas,skip,L,margen,reposo,moduloVelocidad,soloDesplHoriz,redHexagonal24,acelerarRapido,acelerarLento)
tEjecFin = time.time()


# ESCRITURA DE DATOS EN EL FICHERO
temperaturaProm = 0
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
                R = np.linalg.norm(distanciaGlobal(r,t,t,p,j,L,lMedios))
                ePotencialAux += (R**(-12)-R**(-6))
        sumaV += 2*ePotencialAux

        normaV[t,p] = np.linalg.norm(v[t,p])


        temperaturaProm += v[t,p,0]**2 + v[t,p,1]**2
    

    # ESCRITURA EN FICHERO
    eCinetica.append(sumaT)
    ePotencial.append(sumaV)
    eTotal.append(sumaT+sumaV)
    for p in range(nParticulas):
        ficheroPlot.write(str(r[t,p,0]) + ", " + str(r[t,p,1]) + "\n")
    ficheroPlot.write("\n") 

temperaturaProm = 0.5*temperaturaProm/(nIteraciones/skip)

# CÁLCULO DEL PROMEDIO DE VELOCIDADES
promedioVelocidades = np.zeros(int(nIteraciones/(2*skip)))
promedioVelocidadesX = np.zeros(int(nIteraciones/(2*skip)))
promedioVelocidadesY = np.zeros(int(nIteraciones/(2*skip)))
for t in range(0,int(nIteraciones/(2*skip))-nParticulas,nParticulas):
    for p in range(nParticulas):
        promedioVelocidades[t+p] = np.linalg.norm(v[t+int(nIteraciones/(2*skip)),p])
        promedioVelocidadesX[t+p] = v[t+int(nIteraciones/(2*skip)),p,0]
        promedioVelocidadesY[t+p] = v[t+int(nIteraciones/(2*skip)),p,1]

# CÁLCULO DE LA PRESIÓN                                      
presionPromedio = 0                                    
presion = np.zeros(int(nIteraciones*h))                
for i in range(int(nIteraciones*h)):                   
    for t in range(i*int(1/h),int(1/h)*(i+1)):
        presion[i] += fuerzaParedes[t]
    presionPromedio += presion[i]
presionPromedio = presionPromedio/(nIteraciones*h)

# CÁLCULO DE LA TEMPERATURA
temperatura = np.zeros(int(nIteraciones*h))
for i in range(int(nIteraciones*h)):                   
    for t in range(i*int(1/(skip*h)),int(1/(skip*h))*(i+1)):
        for p in range(nParticulas):
            temperatura[i] += v[t,p,0]**2 + v[t,p,1]**2
    temperatura[i] = 0.5*temperatura[i]*skip*h

# CÁLCULO DE LA FLUCTUACIÓN RESPECTO A LA POSICIÓN INICIAL DE UNA PARTÍCULA p
desvPos = np.zeros(int(nIteraciones*h))
p = np.random.randint(0,nParticulas)
for i in range(int(nIteraciones*h)):
    for t in range(i*int(1/(skip*h)),int(1/(skip*h))*(i+1)):
        desvPos[i] += np.linalg.norm(distanciaGlobal(r,0,t,p,p,L,lMedios))**2
    desvPos[i] = desvPos[i]*skip*h

# CÁLCULO DE LA SEPARACIÓN CUADRÁTICA MEDIA ENTRE UN PAR DE ÁTOMOS
desvCuad = np.zeros(int(nIteraciones*h))
p1 = np.random.randint(0,nParticulas)
p2 = np.random.randint(0,nParticulas)
while p2==p1:
    p2 = np.random.randint(0,nParticulas)
for i in range(int(nIteraciones*h)):
    for t in range(i*int(1/(skip*h)),int(1/(skip*h))*(i+1)):
        desvCuad[i] += np.linalg.norm(distanciaGlobal(r,t,t,p1,p2,L,lMedios))**2
    desvCuad[i] = desvCuad[i]*skip*h


# REPRESENTACIÓN GRÁFICA - VELOCIDADES Y ENERGÍAS
plt.rcParams["figure.figsize"] = (10,8.5)

ax1 = plt.subplot(3,2,1)
plt.hist(normaV[0])
plt.title("Velocidades iniciales",fontsize=11)
plt.xlabel("|v|",fontsize=9)
plt.ylabel("Frecuencia",fontsize=9)

ax2 = plt.subplot(3,2,2)
plt.hist(promedioVelocidades,bins=100)
plt.title("Promedio de velocidades: t=t_f/2 a t=t_f",fontsize=11)
plt.xlabel("|v|",fontsize=9)
plt.ylabel("Frecuencia",fontsize=9)

ax3 = plt.subplot(3,2,3)
plt.hist(promedioVelocidadesX,bins=100)
plt.title("Promedio de velocidades (eje X): t=t_f/2 a t=t_f",fontsize=11)
plt.xlabel("|v_x|",fontsize=9)
plt.ylabel("Frecuencia",fontsize=9)

ax4 = plt.subplot(3,2,4)
plt.hist(promedioVelocidadesY,bins=100)
plt.title("Promedio de velocidades (eje Y): t=t_f/2 a t=t_f",fontsize=11)
plt.xlabel("|v_y|",fontsize=9)
plt.ylabel("Frecuencia",fontsize=9)

ax5 = plt.subplot(3,1,3)
plt.plot(eCinetica, label= "T")
plt.plot(ePotencial, label="V")
plt.plot(eTotal, label="E = T+V")
plt.title("Energías en función del tiempo")
plt.xlabel("tiempo")
plt.ylabel("Energía")
plt.legend()

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.5)
plt.show()

# REPRESENTACIÓN GRÁFICA - PRESIÓN
plt.rcParams["figure.figsize"] = (5,4)
plt.plot(presion, label="Presión")
plt.title("Presión en función del tiempo")
plt.xlabel("t")
plt.ylabel("Presión")
plt.axhline(y=presionPromedio, color='r', linestyle='-', label="Presión promedio: "+f"{(presionPromedio):.3f}")
plt.legend()
plt.show()

# REPRESENTACIÓN GRÁFICA - DESVIACIÓN PROMEDIO, TEMPERATURA Y SEP. CUAD. MEDIA
plt.rcParams["figure.figsize"] = (5,8.5)

ax6 = plt.subplot(3,1,1)
plt.plot(desvPos, label="Desviación promedio")
plt.title("Desviación promedio de la posición de una partícula")
plt.xlabel("t")
plt.ylabel("Desviación promedio")
plt.legend()

ax7 = plt.subplot(3,1,2)
plt.plot(temperatura, label="Temperatura")
plt.title("Temperatura en función del tiempo")
plt.xlabel("t")
plt.ylabel("Temperatura")
plt.legend()

ax8 = plt.subplot(3,1,3)
plt.plot(temperatura, label="Desv. Cuad. Media")
plt.title("Desv. cuad. media en función del tiempo")
plt.xlabel("t")
plt.ylabel("Desviación cuadrática media")
plt.legend()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.5)
plt.show()


# Por último, escribimos algunos datos de interés al final del fichero
ficheroPlot.write("# Se han realizado "+str(nIteraciones)+" iteraciones con h = "+str(h)+", "+str(nParticulas)+" particulas y skip "+str(skip)+"\n")
tEjecFin = time.time()
ficheroPlot.write("# Tiempo de ejecucion: "+str(tEjecFin-tEjecIni)+"\n")
ficheroPlot.write("# Temperatura: "+str(temperaturaProm)+"\n")
ficheroPlot.write("# Presion: "+f"{(presionPromedio):.3f}"+"\n")

tFicherosFin = time.time()
