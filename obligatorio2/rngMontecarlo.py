import random
import numpy as np
import time


# DEFINICIÓN DE CONSTANTES Y PARÁMETROS
N = 5  # Dimensión de la cuadrícula
T = 5  # Temperatura T = [0,5]


def rngVectorParejas(N):  # Genera un vector de N parejas de números aleatorios entre 0 y N
    vectorParejas = np.random.choice([-1,+1], size=(N, N))
    print(vectorParejas)  # ---- TEMPORAL, BORRAR ESTA LÍNEA MÁS TARDE ----
    return vectorParejas


def rngX():  # Genera un float entre 0 y 1
    x = random.random()
    print(x)  # ---- TEMPORAL, BORRAR ESTA LÍNEA MÁS TARDE ----
    return x


# Generamos los números aleatorios
rngVectorParejas(N)
rngX()
