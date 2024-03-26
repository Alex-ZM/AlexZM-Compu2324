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
h = 0.5                         # Paso temporal, inverso a la precisión  (CAMBIAR)
velocidadMaxima = 50000         # Velocidad maxima de los planetesimales (CAMBIAR)

# Creamos los vectores con los datos iniciales de los cuerpos
r = np.full(numeroCuerpos, 75*10**7)
m = np.full(numeroCuerpos, 6.3*10**22)
d = np.array([[radioSistemaSolar/np.sqrt(2),radioSistemaSolar/np.sqrt(2)]])
v = np.array([[0, 0]])
for i in range(1, numeroCuerpos):
    d=np.append(d, [[random.random()*radioSistemaSolar,random.random()*radioSistemaSolar]], axis=0)
    v=np.append(v, [[random.random()*velocidadMaxima,random.random()*velocidadMaxima]], axis=0)


def a(i):  # Valor de la aceleración del planetesimal en el instante t
    return (-masaSolar * d[i])/(np.linalg.norm(d[i]))**3


def w(i):  # Valor de w del planetesimal i
    return (v[i] + (h/2)*a(i))


def evD(i):  # Evolución temporal de la posición del planeta i
    return (r[i] + h*w(i))


def evA(i):  # Valor de la aceleración del planetesimal en el instante t
    return (-masaSolar * evD[i])/(np.linalg.norm(evD[i]))**3


def evV(i):  # Evolución temporal de la velocidad del planeta i
    return (w(i) + (h/2)*evA(i))


def checkColision():  # Comprueba si ha habido colision entre los plnts i y j
    for i,j in numeroCuerpos:
        if # Comprobar si chocan i y j
