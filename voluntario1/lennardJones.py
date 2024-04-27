# Simulación del sistema solar
import numpy as np
import random
import time
import os


# Definimos algunas constantes (Sistema Internacional - reescalado)
h = 0.0001          # <---------- Paso temporal, inverso a la precisión (CAMBIAR)
nIter = 1000000     # <---------- Número de iteraciones (CAMBIAR)
t = 0
nParticulas = 28
tEjecIni = time.time()
L = 10
margen = 0.5
hMedios = h/2

# Definimos (y reescalamos) los parámetros iniciales de los planetas
m = np.full(nParticulas,1).astype(np.int8)

r = np.zeros((nParticulas,2))

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
    r[i,0] += (2*random.random() - 1)*margen  # Desplaza aleatoriamente en la componente X (+-margen)
    r[i,1] += (2*random.random() - 1)*margen  # Desplaza aleatoriamente en la componente Y (+-margen)


def a(i):  # Valor de la aceleración del planeta i en el instante actual
    aFinal = np.array([0.0,0.0])
    for j in range(0, nParticulas):
        if i != j:
            R = np.subtract(r[i], r[j])
            aFinal = aFinal - (m[j]*R)/(np.linalg.norm(R))**3
    return aFinal


def w(i):  #Valor de w del planeta i
    return (v[i] + (hMedios)*a(i))


def evR(i):  # Evolución temporal de la posición del planeta i
    return (r[i] + h*w(i))


def evA(i):  # Evolución temporal de la aceleración del planeta i
    aFinal = np.array([0.0,0.0])
    for j in range(0,nPlanetas):
        if i != j:
            R = np.subtract(evR(i), evR(j))
            aFinal = aFinal - (m[j]*R)/(np.linalg.norm(R))**3
    return aFinal


def evV(i):  # Evolución temporal de la velocidad del planeta i
    return (w(i) + (hMedios)*evA(i))


# Ahora solo queda programar el bucle y guardar los resultados de cada iteración en el
# formato correcto y dentro de un fichero, para poder representarlos luego.
wd = os.path.dirname(__file__)      # Directorio de trabajo
planetasPath = os.path.join(wd,"planets_data.dat") # Nombre del fichero de datos
posicionesPath = os.path.join(wd,"posiciones_data") # Nombre del fichero de datos
velocidadesPath = os.path.join(wd,"velocidades") # Nombre del fichero de datos
ficheroPlot = open(planetasPath, "w")
if guardarVelocidades:
    ficheroPosiciones = open(posicionesPath, "w")
    ficheroVelocidades = open(velocidadesPath, "w")
for j in range(nIter):

    if j%skip==0:  # Guardamos la posición de los planetas cada "skip" iteraciones

        if guardarVelocidades:  # Si lo elegimos, guardamos las velocidades en un fichero
            for i in range(0, nPlanetas):
                ficheroVelocidades.write(str(v[i][0]) + " " + str(v[i][1]) + "\n")  # Calcula e introduce las velocidades de los planetas en el fichero
                ficheroPosiciones.write(str(r[i][0]) + " " + str(r[i][1]) + "\n")  # Calcula e introduce las posiciones de los planetas en el fichero
                
        for i in range(0, nPlanetas):
            ficheroPlot.write(str(r[i][0]) + ", " + str(r[i][1]) + "\n")  # Calcula e introduce las posiciones de los planetas en el fichero
        ficheroPlot.write("\n")  # Para separar los grupos de datos por instante temporal

    for i in range(1,nPlanetas):
        if T[i]==1 and r[i][1]<0:
            T[i] = 2*t  # Guarda el período de cada planeta

        r[i] = evR(i)  # 
        v[i] = evV(i)  # Avance temporal: t = t+h
    t = t + h          #

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

if guardarVelocidades:
    ficheroPosiciones.close()
    ficheroVelocidades.close()
