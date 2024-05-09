
    ############################################################
    ###  ALGORITMO DE VERLET - SIMULACIÓN DEL SISTEMA SOLAR  ###
    ############################################################

import numpy as np
import random
import time
import os
from numba import jit

###################################################################################################################

# Definimos algunas constantes
masaSolar = 1.98855*10**30 
UA = 1.496*10**11  
G = 6.67*10**(-11)  
h = 0.0001          
nIteraciones = 1000000    
skip = 1000         
guardarVelocidades = False  # <--- Elije si guardar también las velocidades (CAMBIAR)
nPlanetas = 10

hMedios = h/2

def reescalarV(v):  # Función para reescalar t
    return v*np.sqrt(UA/(G*masaSolar))

m = np.array(
    [[1],                       
    [330.2*10**21/masaSolar],  
    [4868.5*10**21/masaSolar], 
    [5973.6*10**21/masaSolar], 
    [641.85*10**21/masaSolar], 
    [1.899*10**27/masaSolar],  
    [0.568*10**27/masaSolar],  
    [0.087*10**27/masaSolar],  
    [0.102*10**27/masaSolar], 
    [12.5*10**21/masaSolar]])  

r = np.array(
    [[0,10**-15],               
    [57.9*10**9/UA,10**-15],   
    [108.2*10**9/UA,10**-15],  
    [1,10**-15],               
    [227.9*10**9/UA,10**-15],  
    [778.6*10**9/UA,10**-15],  
    [1433.5*10**9/UA,10**-15], 
    [2872.5*10**9/UA,10**-15], 
    [4495.1*10**9/UA,10**-15], 
    [5870*10**9/UA,10**-15]])   

v = np.array(
    [[0,0],               
    [0,reescalarV(47890)], 
    [0,reescalarV(35030)], 
    [0,reescalarV(29790)], 
    [0,reescalarV(24130)], 
    [0,reescalarV(13100)], 
    [0,reescalarV(9700)], 
    [0,reescalarV(6800)], 
    [0,reescalarV(5400)], 
    [0,reescalarV(4700)]])

###################################################################################################################

# DEFINICIÓN DE FUNCIÓN - ALGORITMO DE VERLET
@jit(nopython=True,fastmath=True)
def verlet(m,r,v,nIteraciones,nPlanetas,skip):

    T = np.ones(10)
    posiciones = np.zeros((int(nIteraciones/skip),nPlanetas,2))
    velocidades = np.zeros((int(nIteraciones/skip),nPlanetas,2))
    a = np.zeros((nPlanetas,2))
    w = np.zeros((nPlanetas,2))
    evR = np.zeros((nPlanetas,2))
    evV = np.zeros((nPlanetas,2))
    evA = np.zeros((nPlanetas,2))

    hMedios = h/2
    
    for t in range(nIteraciones):

        # GUARDA LOS RESULTADOS EN LOS VECTORES "posiciones" Y "velocidades"
        if t%skip == 0:
            for p in range(nPlanetas):
                posiciones[int(t/skip),p] = r[p]
                velocidades[int(t/skip),p] = r[p]

        # CÁLCULOS DEL ALGORITMO DE VERLET
        for p in range(nPlanetas):
            aux1 = np.array([0.0,0.0])
            for j in range(nPlanetas):
                if p != j:
                    R = np.subtract(r[p],r[j])
                    aux1 = aux1 - (m[j]*R)/(np.linalg.norm(R))**3
            a[p] = aux1

        for p in range(nPlanetas):
            w[p] = v[p] + hMedios*a[p]

        for p in range(nPlanetas):
            evR[p] = r[p] + h*w[p]
        
        for p in range(nPlanetas):
            aux2 = np.array([0.0,0.0])
            for j in range(nPlanetas):
                if p != j:
                    R = np.subtract(evR[p],evR[j])
                    aux2 = aux2 - (m[j]*R)/(np.linalg.norm(R))**3
            evA[p] = aux2

        for p in range(nPlanetas):
            evV[p] = w[p] + hMedios*evA[p]

        # EVOLUCIÓN TEMPORAL
        for p in range(nPlanetas):
            r[p] = evR[p]
            v[p] = evV[p]

        for i in range(1,nPlanetas):
            if T[i]==1 and r[i][1]<0:
                T[i] = 2*t  # Guarda el período de cada planeta

    return posiciones,velocidades,T

###################################################################################################################

wd = os.path.dirname(__file__)  # Directorio de trabajo
datosPath = os.path.join(wd,"posPlanetas.dat")  
ficheroPlot = open(datosPath, "w")

# CÁLCULO DE POSICIONES, VELOCIDADES Y PERÍODOS MEDIANTE LA FUNCIÓN "verlet()"
tEjecIni = time.time()
r,v,T = verlet(m,r,v,nIteraciones,nPlanetas,skip)
tEjecFin = time.time()

# ESCRITURA DE DATOS EN EL FICHERO
tFicherosIni = time.time()
for t in range(int(nIteraciones/skip-1)):

    for p in range(nPlanetas):
        ficheroPlot.write(str(r[t,p,0]) + ", " + str(r[t,p,1]) + "\n")  # Posiciones de las partículas -> Fichero
    ficheroPlot.write("\n") 

    #for p in range(nPlanetas):  # EVOLUCIÓN TEMPORAL + CÁLCULO ENERGÍAS

        #T[p] = 0.5*np.linalg.norm(v[p])**2
        #V[p] = 0
    
        #for j in range(nPlanetas):
        #    if p != j:
        #        R = np.linalg.norm(distanciaToroide(r,p,j))
        #        V[p] += (R**(-12)-R**(-6))
        #V[p] = 4*V[p]
tFicherosFin = time.time()
        
###################################################################################################################

# Por último, escribimos algunos datos de interés al final del fichero
ficheroPlot.write("# Se han realizado "+str(nIteraciones)+" iteraciones con h = "+str(h)+" y skip "+str(skip)+"\n"+"#"+"\n")
ficheroPlot.write("# Tiempo de ejecucion - ALGORITMO DE VERLET .............. "+str(tEjecFin-tEjecIni)+"\n")
ficheroPlot.write("# Tiempo de ejecucion - ESCRITURA DE DATOS EN FICHEROS ... "+str(tFicherosFin-tFicherosIni)+"\n"+"#"+"\n")
ficheroPlot.write("# T(1) = "+str(T[1]/T[3]*365.256)+" dias terrestres (vs. 87.969)\n")
ficheroPlot.write("# T(2) = "+str(T[2]/T[3]*365.256)+" dias terrestres (vs. 224.699)\n")
ficheroPlot.write("# T(3) = "+str(T[3]/T[3]*365.256)+" dias terrestres (vs. 365.256)\n")
ficheroPlot.write("# T(4) = "+str(T[4]/T[3]*365.256)+" dias terrestres (vs. 686.979)\n")
ficheroPlot.write("# T(5) = "+str(T[5]/T[3]*365.256)+" dias terrestres (vs. 4332.589)\n")
ficheroPlot.write("# T(6) = "+str(T[6]/T[3]*365.256)+" dias terrestres (vs. 10759.23)\n"+"#"+"\n")

ficheroPlot.close()
