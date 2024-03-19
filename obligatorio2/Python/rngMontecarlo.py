import random
import numpy as np


def rngVectorParejas(N):  # Genera un vector de N parejas de números aleatorios entre 0 y N
    vectorParejas = np.random.randint(N, size=(N, N))
    print(vectorParejas)  # ---- TEMPORAL, BORRAR ESTA LÍNEA MÁS TARDE ----
    return vectorParejas


def rngX():  # Genera un float entre 0 y 1
    x = random.random()
    print(x)  # ---- TEMPORAL, BORRAR ESTA LÍNEA MÁS TARDE ----
    return x


# Input del tamaño de la cuadrícula
N = int(input("\nIntroduce N (tamaño de la cuadrícula): "))

# Generamos los números aleatorios
rngVectorParejas(N)
rngX()
