# Simulación del sistema solar
import numpy as np
import random
import time
import os


# Definimos algunas constantes (Sistema Internacional - reescalado)
h = 0.0001    
nIter = 1000000     
nParticulas = 5
tEjecIni = time.time()
L = 10
margen = 0.5
skip = 20
t = 0
hMedios = h/2

# Definimos (y reescalamos) los parámetros iniciales de las partículas
m = np.full(nParticulas,1).astype(np.int8)
r = np.zeros((nParticulas,2))
v = np.zeros((nParticulas,2))
a = np.zeros((nParticulas,2))
w = np.zeros((nParticulas,2))
evR = np.zeros((nParticulas,2))
evA = np.zeros((nParticulas,2))

###################################################################################################################

# POSICIONAMIENTO EQUIESPACIADO DE LOS PLANETAS
particulasPorFila = int(np.floor(np.sqrt(nParticulas)))
separacionInicial = L/particulasPorFila
for f in range(particulasPorFila):
    for i in range(particulasPorFila):  
        r[f*particulasPorFila+i,1] = i*separacionInicial  # Equiespaciado horizontal
particulasUltimaFila = nParticulas-particulasPorFila**2
if particulasUltimaFila != 0:
    for i in range(particulasUltimaFila):
        r[particulasPorFila**2+i,1] = i*separacionInicial  # Equiespaciado horizontal (última fila)

for c in range(particulasPorFila):
    for i in range(particulasPorFila):
        r[c*particulasPorFila+i,0] = c*separacionInicial  # Equiespaciado vertical
if particulasUltimaFila != 0:
    for i in range(particulasUltimaFila):
        r[particulasPorFila**2+i,0] = particulasPorFila*separacionInicial  # Equiespaciado vertical (última fila)

# LIGEROS DESPLAZAMIENTOS ALEATORIOS
for i in range(nParticulas):
    r[i,0] += (2*random.random()-1)*margen  # Desplaza aleatoriamente en la componente X (+-margen)
    r[i,1] += (2*random.random()-1)*margen  # Desplaza aleatoriamente en la componente Y (+-margen)

###################################################################################################################

# Ahora solo queda programar el bucle y guardar los resultados de cada iteración en el
# formato correcto y dentro de un fichero, para poder representarlos luego.
wd = os.path.dirname(__file__)                   # Directorio de trabajo
datosPath = os.path.join(wd,"posParticulas.dat")  # Nombre del fichero de datos
ficheroPlot = open(datosPath, "w")

# BUCLE - ALGORITMO DE VERLET
for t in range(nIter):

    if t%skip==0:  # Guardamos la posición de los planetas cada "skip" iteraciones    
        for i in range(nParticulas):
            ficheroPlot.write(str(r[i][0]) + ", " + str(r[i][1]) + "\n")  # Introduce las posiciones de los planetas en el fichero
        ficheroPlot.write("\n")  # Para separar los grupos de datos por instante temporal

    for part in range(nParticulas):
        aux1 = np.array([0.0,0.0])  # a
        for j in range(nParticulas):
            if part != j:
                R = np.subtract(r[part],r[j])
                aux1 = aux1 - (m[j]*R)/(np.linalg.norm(R))**3
            a[part] = aux1

        for j in range(nParticulas):  # w
            w[part] = v[part]+hMedios*a[part]

        for j in range(nParticulas):  # evR
            evR[part] = r[part]+h*w[part]

        aux2 = np.array([0.0,0.0])  # evA
        for j in range(nParticulas):  
            if part != j:
                R = np.subtract(evR[part],evR[j])
                aux2 = aux2 - (m[j]*R)/(np.linalg.norm(R))**3
            evA[part] = aux2

        for j in range(nParticulas):
            r[part] = evR[part]
            v[part] = w[part]+hMedios*evA[part]

###################################################################################################################

# Por último, escribimos algunos datos de interés al final del fichero
ficheroPlot.write("# Se han realizado "+str(nIter)+" iteraciones con h = "+str(h)+" y skip "+str(skip)+"\n")
tEjecFin = time.time()
ficheroPlot.write("# Tiempo de ejecución: "+str(tEjecFin-tEjecIni))
