# Simulación del sistema solar
import numpy as np
import time
import os
from numba import njit


# Definimos algunas constantes (Sistema Internacional - reescalado)
masaSolar = 1.98855*10**30  # Masa del Sol
UA = 1.496*10**11  # Distancia Tierra-Sol
G = 6.67*10**(-11)  # Cte de Gravitación Universal
h = 0.0001          # <---------- Paso temporal, inverso a la precisión (CAMBIAR)
nIter = 10000     # <---------- Número de iteraciones (CAMBIAR)
skip = 50          # <---------- Cada cuántas iteraciones guarda datos en los ficheros (CAMBIAR)
guardarVelocidades = False  # <--- Elije si guardar también las velocidades (CAMBIAR)
nPlanetas = 8
tEjecIni = time.time()

hMedios = h/2

def reescalarV(v):  # Función para reescalar t
    return v*np.sqrt(UA/(G*masaSolar))


# Definimos (y reescalamos) los parámetros iniciales de los planetas
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
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
    [12.5*10**21/masaSolar],
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
    [5870*10**9/UA,10**-15],
    [(5870*10**9+(5870-4491)*10**9)/UA,0],
    [(5870*10**9+2*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+3*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+4*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+5*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+6*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+7*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+8*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+9*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+10*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+11*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+12*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+13*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+14*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+15*(5870-4491)*10**9)/UA,0],
    [(5870*10**9+16*(5870-4491)*10**9)/UA,0]])   

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
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)],
    [0,reescalarV(4700)]])

T = np.array(
    [[1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1],
    [1]])

vFichero = []
rFichero = []

################################################################################################
#@njit
def sistSolar(nIter,nPlanetas,skip):
    t=0
    a = np.zeros((nPlanetas,2))
    w = np.zeros((nPlanetas,2))
    evR = np.zeros((nPlanetas,2))
    evA = np.zeros((nPlanetas,2))
    for n in range(nIter):
        if n%skip==0:  # Guardamos la posición de los planetas cada "skip" iteraciones

            if guardarVelocidades:  # Si lo elegimos, guardamos las velocidades en un fichero
                for i in range(nPlanetas):
                    vFichero.append(v[i])  # Guarda las velocidades de cada frame en un vector
                    
            for i in range(nPlanetas):
                rFichero.append(r[i])  # Guarda las posiciones de cada frame en un vector


        for pl in range(nPlanetas):
            aux1 = np.array([0.0,0.0])  # a
            for j in range(nPlanetas):
                if pl != j:
                    R = np.subtract(r[pl],r[j])
                    aux1 = aux1 - (m[j]*R)/(np.linalg.norm(R))**3
                a[pl] = aux1

            for j in range(nPlanetas):  # w
                w[pl] = v[pl]+hMedios*a[pl]

            for j in range(nPlanetas):  # evR
                evR[pl] = r[pl]+h*w[pl]

            aux2 = np.array([0.0,0.0])  # evA
            for j in range(nPlanetas):  
                if pl != j:
                    R = np.subtract(evR[pl],evR[j])
                    aux2 = aux2 - (m[j]*R)/(np.linalg.norm(R))**3
                evA[pl] = aux2

            for j in range(nPlanetas):
                r[pl] = evR[pl]
                v[pl] = w[pl]+hMedios*evA[pl]
                

        for pl in range(1,nPlanetas):
            if T[pl]==1 and r[pl,1]<0:
                T[pl] = 2*t  # Guarda el período de cada planeta
        t = t + h  
                
################################################################################################

sistSolar(nIter,nPlanetas,skip)

# Ahora solo queda programar el bucle y guardar los resultados de cada iteración en el
# formato correcto y dentro de un fichero, para poder representarlos luego.
wd = os.path.dirname(__file__)      # Directorio de trabajo
planetasPath = os.path.join(wd,"planets_dataN.dat") # Nombre del fichero de datos
velocidadesPath = os.path.join(wd,"velocidadesN") # Nombre del fichero de datos
ficheroPlot = open(planetasPath, "w")
if guardarVelocidades:
    ficheroVelocidades = open(velocidadesPath, "w")

for i in range(len(vFichero)):
    if guardarVelocidades:  # Si lo elegimos, guardamos las velocidades en un fichero
            for i in range(nPlanetas):
                ficheroVelocidades.write(str(vFichero[i]) + "\n")  # Calcula e introduce las velocidades de los planetas en el fichero
                
    for i in range(nPlanetas):
        ficheroPlot.write(str(r[i] + "\n"))  # Calcula e introduce las posiciones de los planetas en el fichero
    ficheroPlot.write("\n")  # Para separar los grupos de datos por instante temporal

# Por último, escribimos algunos datos de interés al final del fichero
ficheroPlot.write("# Se han realizado "+str(nIter)+" iteraciones con h = "+str(h)+" y skip "+str(skip)+"\n")
#ficheroPlot.write("# T(1) = "+str(T[1]/T[3]*365.256)+" días terrestres (vs. 87.969)\n")
#ficheroPlot.write("# T(2) = "+str(T[2]/T[3]*365.256)+" días terrestres (vs. 224.699)\n")
#ficheroPlot.write("# T(3) = "+str(T[3]/T[3]*365.256)+" días terrestres (vs. 365.256)\n")
#ficheroPlot.write("# T(4) = "+str(T[4]/T[3]*365.256)+" días terrestres (vs. 686.979)\n")
#ficheroPlot.write("# T(5) = "+str(T[5]/T[3]*365.256)+" días terrestres (vs. 4332.589)\n")
#ficheroPlot.write("# T(6) = "+str(T[6]/T[3]*365.256)+" días terrestres (vs. 10759.23)\n")

tEjecFin = time.time()
ficheroPlot.write("# Tiempo de ejecución: "+str(tEjecFin-tEjecIni))

ficheroPlot.close()
if guardarVelocidades:
    ficheroVelocidades.close()
