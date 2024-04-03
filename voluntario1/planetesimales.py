#
# Formación de los planetas a partir de planetesimales
#
import numpy as np
import random

# Definimos algunas constantes
masaSolar = 1.98855*10**30      # Masa del Sol (kg)
radioSolar = 6.96*10**8         # Radio del sol (m)
radioSistemaSolar = 4.5*10**12  # Radio del sistema solar (m)            (CAMBIAR)
UA = 1.496*10**11               # Distancia Tierra-Sol (m)
G = 6.67*10**(-11)              # Cte de Gravitación Universal
numeroCuerpos = 2500            # Número de planetesimales inicial       (CAMBIAR)
h = 0.05                         # Paso temporal, inverso a la precisión  (CAMBIAR)
nIter = 1000
velocidadMaxima = 50000         # Velocidad maxima de los planetesimales (CAMBIAR)

# Creamos los vectores con los datos iniciales de los cuerpos
r = np.full(numeroCuerpos, 75*10**7, dtype=np.float64)
m = np.full(numeroCuerpos, 6.3*10**22)
d = np.array([[radioSistemaSolar/np.sqrt(2),radioSistemaSolar/np.sqrt(2)]])
v = np.array([[0, 0]])
a = np.array([[0, 0]])
for i in range(1, numeroCuerpos):
    d=np.append(d, [[random.random()*radioSistemaSolar,random.random()*radioSistemaSolar]], axis=0)
    v=np.append(v, [[random.random()*velocidadMaxima,random.random()*velocidadMaxima]], axis=0)
    a=np.append(a, [[0,0]], axis=0)


def acel(i):  # Valor de la aceleración del planetesimal en el instante t
    return (-masaSolar * d[i])/(np.linalg.norm(d[i]))**3


def w(i):  # Valor de w del planetesimal i
    return (v[i] + (h/2)*acel(i))


def evD(i):  # Evolución temporal de la posición del planeta i
    return (r[i] + h*w(i))


def evA(i):  # Valor de la aceleración del planetesimal en el instante t
    return (-masaSolar * evD(i))/(np.linalg.norm(evD(i)))**3


def evV(i):  # Evolución temporal de la velocidad del planeta i
    return (w(i) + (h/2)*evA(i))


def colision(i,j):  # Crea un nuevo planetesimal suma de i y j
    v[i] = (m[i]*v[i]+m[j]*v[j])/(m[i]+m[j])  # Velocidad resultante (conserva el momento)
    r[i] = r[i]+r[j]  # (r[i]**3+r[j]**3)**(0.333)   # Radio del planetesimal suma (conserva el volumen)
    m[i] = (m[i]+m[j])  # Masa del planetesimal suma
    np.delete(r, j)          # 
    np.delete(m, j)          # 
    np.delete(d, j, axis=0)  # Elimina el planetesimal j (residuo)
    np.delete(v, j, axis=0)  #
    np.delete(a, j, axis=0)  #


def checkColision():  # Comprueba si ha habido colision entre i y j
    for i in range(numeroCuerpos):
        for j in range(numeroCuerpos):
            dist = np.linalg.norm(np.subtract(d[i],d[j]))  # Distancia entre i y j
            if dist<(r[i]+r[j]):
                colision(i,j)  # Simula la colisión de los planetesimales si chocan


# Ahora solo queda programar el bucle y guardar los resultados de cada iteración en el
# formato correcto y dentro de un fichero, para poder representarlos luego.
ficheroPosiciones = open("sSolar_data.dat", "w")
for _ in range(nIter):

    for i in range(len(r)):

        ficheroPosiciones.write(str(d[i][0]) + ", " + str(d[i][1]) + "\n")  # Calcula entroduce las posiciones de los planetas en el fichero
    ficheroPosiciones.write("\n")  # Para separar los grupos de datos por instante temporal

    for i in range(len(r)):
        d[i] = evD(i)  #
        a[i] = evA(i)  # Avance temporal: t = t+h
        v[i] = evV(i)  #

    checkColision()  # Comprobamos si hay alguna colision

ficheroPosiciones.close