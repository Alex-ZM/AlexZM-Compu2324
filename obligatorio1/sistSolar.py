# Simulación del sistema solar
import numpy as np
import time


# Definimos algunas constantes (Sistema Internacional - reescalado)
masaSolar = 1.98855*10**30  # Masa del Sol
UA = 1.496*10**11  # Distancia Tierra-Sol
G = 6.67*10**(-11)  # Cte de Gravitación Universal
h = 0.0001          # <---------- Paso temporal, inverso a la precisión (CAMBIAR)
nIter = 10000     # <---------- Número de iteraciones (CAMBIAR)
skip = 50          # <---------- Cada cuántas iteraciones guarda datos en los ficheros (CAMBIAR)
guardarVelocidades = False  # <--- Elije si guardar también las velocidades (CAMBIAR)
t = 0
nPlanetas = 4
tEjecIni = time.time()


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
    [1]])


def a(i):  # Valor de la aceleración del planeta i en el instante actual
    aFinal = np.array([0.0,0.0])
    for j in range(0, nPlanetas):
        if i != j:
            R = np.subtract(r[i], r[j])
            aFinal = aFinal - (m[j] * R)/(np.linalg.norm(R))**3
    return aFinal


def w(i):  #Valor de w del planeta i
    return (v[i] + (h/2)*a(i))


def evR(i):  # Evolución temporal de la posición del planeta i
    return (r[i] + h*w(i))


def evA(i):  # Evolución temporal de la aceleración del planeta i
    aFinal = np.array([0.0,0.0])
    for j in range(0,nPlanetas):
        if i != j:
            R = np.subtract(evR(i), evR(j))
            aFinal = aFinal - (m[j] * R)/(np.linalg.norm(R))**3
    return aFinal


def evV(i):  # Evolución temporal de la velocidad del planeta i
    return (w(i) + (h/2)*evA(i))


# Ahora solo queda programar el bucle y guardar los resultados de cada iteración en el
# formato correcto y dentro de un fichero, para poder representarlos luego.
ficheroPlot = open("planets_data.dat", "w")
ficheroPosiciones = open("posiciones.dat", "w")
ficheroVelocidades = open("velocidades.dat", "w")
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
ficheroPlot.write("# Se han realizado "+str(nIter)+" iteraciones con h = "+str(h)+"\n")
ficheroPlot.write("# T(1) = "+str(T[1]/T[3]*365.256)+" días terrestres (vs. 87.969)\n")
ficheroPlot.write("# T(2) = "+str(T[2]/T[3]*365.256)+" días terrestres (vs. 224.699)\n")
ficheroPlot.write("# T(3) = "+str(T[3]/T[3]*365.256)+" días terrestres (vs. 365.256)\n")
#ficheroPlot.write("# T(4) = "+str(T[4]/T[3]*365.256)+" días terrestres (vs. 686.979)\n")
#ficheroPlot.write("# T(5) = "+str(T[5]/T[3]*365.256)+" días terrestres (vs. 4332.589)\n")
#ficheroPlot.write("# T(6) = "+str(T[6]/T[3]*365.256)+" días terrestres (vs. 10759.23)\n")

tEjecFin = time.time()
ficheroPlot.write("Tiempo de ejecución: "+str(tEjecFin-tEjecIni))

ficheroPosiciones.close()
ficheroVelocidades.close()
